"""
ENAQT N-Site Chain Scaling Analysis
=====================================
Refactored from enaqt_nsite_chain.py (625 lines).

Functions:
- run_scaling:     N-site scaling with saturation model + bootstrap CIs.
- run_fmo_benchmark: FMO-7 validation with explicit caveats.
- time_trajectory: Time-domain dynamics via solve_ivp (dead-code bug fixed).
- fit_models:      Internal helper for linear + power-saturation fitting.

Critical fixes vs original:
1. Uses L_base + L_deph decomposition from core for efficiency.
2. Adds saturation model E(N) = a * N^b * exp(-c*N).
3. Adds bootstrap confidence intervals to all fits.
4. Adds diagnostic statistics (Durbin-Watson, Cook's D).
5. Fixes dead code bug in time_trajectory (deta_correct overwritten).
6. Reports both enhancement ratio AND absolute yield.
7. FMO section includes explicit caveats about Lindblad limitations.

References:
- Rebentrost et al., New J. Phys. 11, 033003 (2009)
- Adolphs & Renger 2006 (FMO Hamiltonian)
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Tuple

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from scipy.stats import t as t_dist

# ---------------------------------------------------------------------------
# Imports from enaqt.core (merged separately)
# ---------------------------------------------------------------------------
from enaqt.core import (
    liouvillian_parts,
    analytical_yield,
    hamiltonian_nsite,
    hamiltonian_funnel,
    hamiltonian_flat,
    hamiltonian_disordered,
    hamiltonian_fmo7,
    optimal_dephasing,
    get_out_dir,
    save_json,
)

# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------
_DELTA: float = 1.0
_KAPPA: float = 0.1
_GAMMA: float = 0.01
_BIAS: float = 5.0
_DISORDER_SIGMA: float = 2.0


# ============================================================================
#  INTERNAL HELPERS
# ============================================================================

def _bootstrap_ci(
    x: np.ndarray,
    y: np.ndarray,
    fit_func,
    p0,
    n_bootstrap: int = 1000,
    confidence: float = 0.95,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute bootstrap confidence intervals for a curve fit.

    Parameters
    ----------
    x, y : data arrays.
    fit_func : callable(x, *params) -> y_pred.
    p0 : initial guess for fit parameters.
    n_bootstrap : number of bootstrap resamples.
    confidence : confidence level (default 0.95 => 95% CI).

    Returns
    -------
    p_opt : best-fit parameters.
    ci_lower : lower CI bound for predictions on a fine grid.
    ci_upper : upper CI bound for predictions on a fine grid.
    """
    # Best fit on full data
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            p_opt, _ = curve_fit(fit_func, x, y, p0=p0, maxfev=10000)
        except RuntimeError:
            p_opt = np.array(p0)

    # Bootstrap resampling
    rng = np.random.default_rng(42)
    n = len(x)
    x_fine = np.linspace(x.min(), x.max(), 200)
    boot_preds = np.zeros((n_bootstrap, len(x_fine)))

    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        xb, yb = x[idx], y[idx]
        try:
            pb, _ = curve_fit(fit_func, xb, yb, p0=p0, maxfev=10000)
            boot_preds[b] = fit_func(x_fine, *pb)
        except (RuntimeError, ValueError):
            boot_preds[b] = fit_func(x_fine, *p_opt)

    alpha = 1 - confidence
    ci_lower = np.percentile(boot_preds, 100 * alpha / 2, axis=0)
    ci_upper = np.percentile(boot_preds, 100 * (1 - alpha / 2), axis=0)
    return p_opt, ci_lower, ci_upper


def _durbin_watson(residuals: np.ndarray) -> float:
    """Durbin-Watson statistic for autocorrelation detection.

    DW ≈ 2  → no autocorrelation.
    DW < 1.5 → possible positive autocorrelation.
    DW > 2.5 → possible negative autocorrelation.
    """
    diff = np.diff(residuals)
    ss_res = np.sum(residuals ** 2)
    if ss_res == 0:
        return 2.0
    return float(np.sum(diff ** 2) / ss_res)


def _cooks_distance(
    x: np.ndarray, y: np.ndarray, fit_func, p_opt: np.ndarray
) -> np.ndarray:
    """Cook's distance for identifying influential points.

    D_i > 4 / n  is a common rule-of-thumb threshold.
    """
    n = len(x)
    y_hat_full = fit_func(x, *p_opt)
    mse_full = np.mean((y - y_hat_full) ** 2)
    if mse_full == 0:
        return np.zeros(n)

    cooks = np.zeros(n)
    for i in range(n):
        x_loo = np.delete(x, i)
        y_loo = np.delete(y, i)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                p_loo, _ = curve_fit(fit_func, x_loo, y_loo, p0=p_opt, maxfev=10000)
            except (RuntimeError, ValueError):
                p_loo = p_opt
        y_hat_loo = fit_func(x, *p_loo)
        cooks[i] = float(np.sum((y_hat_full - y_hat_loo) ** 2) / (len(p_opt) * mse_full))
    return cooks


def _r_squared(y: np.ndarray, y_hat: np.ndarray) -> float:
    """Coefficient of determination R^2."""
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    if ss_tot == 0:
        return 1.0
    return float(1 - ss_res / ss_tot)


def _aic(n: int, rss: float, k: int) -> float:
    """Akaike Information Criterion.

    AIC = n * ln(RSS / n) + 2k
    Lower AIC = better model.
    """
    if rss <= 0 or n <= k:
        return float("inf")
    return float(n * np.log(rss / n) + 2 * k)


# ============================================================================
#  fit_models
# ============================================================================

def fit_models(
    n_values: np.ndarray,
    enhancements: np.ndarray,
    n_bootstrap: int = 1000,
) -> dict:
    """Fit linear model + power-saturation model with bootstrap CIs.

    Models:
    1. Linear:  E(N) = a * N + b
    2. Saturation: E(N) = a * N**b * exp(-c * N)

    Parameters
    ----------
    n_values : array of chain lengths N.
    enhancements : array of ENAQT enhancement ratios.
    n_bootstrap : number of bootstrap resamples for CIs.

    Returns
    -------
    dict with:
        - linear_params, linear_r2, linear_ci_lower, linear_ci_upper,
          linear_aic, linear_durbin_watson, linear_cooks_d
        - sat_params, sat_r2, sat_ci_lower, sat_ci_upper,
          sat_aic, sat_durbin_watson, sat_cooks_d
        - preferred_model ('linear' or 'saturation')
        - x_fine : grid for CI curves
    """
    x = np.asarray(n_values, dtype=float)
    y = np.asarray(enhancements, dtype=float)
    x_fine = np.linspace(x.min(), x.max(), 200)

    # ---------- Linear fit ----------
    linear_func = lambda xv, a, b: a * xv + b
    p0_linear = [np.polyfit(x, y, 1)[0], np.polyfit(x, y, 1)[1]]

    p_linear, ci_l_lin, ci_u_lin = _bootstrap_ci(
        x, y, linear_func, p0_linear, n_bootstrap=n_bootstrap
    )
    y_hat_linear = linear_func(x, *p_linear)
    r2_linear = _r_squared(y, y_hat_linear)
    dw_linear = _durbin_watson(y - y_hat_linear)
    cooks_linear = _cooks_distance(x, y, linear_func, p_linear)
    aic_linear = _aic(len(x), np.sum((y - y_hat_linear) ** 2), 2)

    # ---------- Saturation fit: E(N) = a * N^b * exp(-c * N) ----------
    def sat_func(xv, a, b, c):
        return a * (xv ** b) * np.exp(-c * xv)

    # Initial guess: a ≈ max(y), b ≈ 1, c ≈ small
    p0_sat = [max(y.max(), 1.0), 1.0, 0.01]
    bounds_sat = ([0.0, 0.0, 0.0], [np.inf, 5.0, 10.0])

    p_sat = p0_sat.copy()
    ci_l_sat = np.full_like(x_fine, np.nan)
    ci_u_sat = np.full_like(x_fine, np.nan)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p_sat, _ = curve_fit(sat_func, x, y, p0=p0_sat, bounds=bounds_sat, maxfev=20000)
        _, ci_l_sat, ci_u_sat = _bootstrap_ci(
            x, y, sat_func, p0_sat, n_bootstrap=n_bootstrap
        )
    except (RuntimeError, ValueError):
        pass  # keep initial guess

    y_hat_sat = sat_func(x, *p_sat)
    r2_sat = _r_squared(y, y_hat_sat)
    dw_sat = _durbin_watson(y - y_hat_sat)
    cooks_sat = _cooks_distance(x, y, sat_func, np.array(p_sat))
    aic_sat = _aic(len(x), np.sum((y - y_hat_sat) ** 2), 3)

    preferred = "saturation" if aic_sat < aic_linear else "linear"

    return {
        # Linear
        "linear_params": {"a": float(p_linear[0]), "b": float(p_linear[1])},
        "linear_r2": r2_linear,
        "linear_ci_lower": ci_l_lin.tolist(),
        "linear_ci_upper": ci_u_lin.tolist(),
        "linear_aic": aic_linear,
        "linear_durbin_watson": dw_linear,
        "linear_cooks_d": cooks_linear.tolist(),
        # Saturation
        "sat_params": {"a": float(p_sat[0]), "b": float(p_sat[1]), "c": float(p_sat[2])},
        "sat_r2": r2_sat,
        "sat_ci_lower": ci_l_sat.tolist(),
        "sat_ci_upper": ci_u_sat.tolist(),
        "sat_aic": aic_sat,
        "sat_durbin_watson": dw_sat,
        "sat_cooks_d": cooks_sat.tolist(),
        # Meta
        "preferred_model": preferred,
        "x_fine": x_fine.tolist(),
    }


# ============================================================================
#  time_trajectory
# ============================================================================

def time_trajectory(
    H: np.ndarray,
    gamma_phi: float,
    sink_rate: float,
    fluo_rate: float,
    t_span: Tuple[float, float] = (0.0, 200.0),
    n_pts: int = 2000,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Propagate Lindblad ODE in the time domain via solve_ivp.

    Correctly separates real and imaginary parts of the density matrix.

    Parameters
    ----------
    H : (N, N) Hamiltonian matrix.
    gamma_phi : dephasing rate.
    sink_rate : sink (charge recombination centre) rate kappa.
    fluo_rate : fluorescence (recombination) rate Gamma.
    t_span : (t_min, t_max) time interval.
    n_pts : number of time-evaluation points.

    Returns
    -------
    times : (n_pts,) time array.
    populations : (N, n_pts) diagonal elements rho_jj(t).
    cumulative_yield : (n_pts,) integrated sink yield eta(t).

    CRITICAL FIX: Original code had a dead-code bug where deta_correct
    was computed from the Liouvillian action but then immediately
    overwritten with a simpler (incorrect) expression. This function
    uses the real/imaginary separation consistently throughout.
    """
    N = H.shape[0]
    L_base, L_deph = liouvillian_parts(H, sink_rate, fluo_rate)
    L = L_base + gamma_phi * L_deph

    diag_indices = np.array([j + j * N for j in range(N)], dtype=int)
    sink_idx = (N - 1) + (N - 1) * N

    # Real-valued ODE: separate real and imaginary parts
    L_re = L.real
    L_im = L.imag

    def rhs(t, y):
        rho_re = y[: N * N]
        rho_im = y[N * N : 2 * N * N]

        drho_re = L_re @ rho_re - L_im @ rho_im
        drho_im = L_re @ rho_im + L_im @ rho_re
        deta = sink_rate * rho_re[sink_idx]
        return [*drho_re.tolist(), *drho_im.tolist(), deta]

    y0 = np.zeros(2 * N * N + 1, dtype=float)
    y0[0] = 1.0  # rho_11 real part = 1

    t_eval = np.linspace(t_span[0], t_span[1], n_pts)
    sol = solve_ivp(
        rhs, t_span, y0, method="RK45", t_eval=t_eval, rtol=1e-8, atol=1e-10
    )

    pops = np.array([sol.y[idx] for idx in diag_indices])  # (N, n_pts)
    cumulative = sol.y[2 * N * N]

    return sol.t, pops, cumulative


# ============================================================================
#  run_scaling
# ============================================================================

def run_scaling(
    n_max: int = 20,
    out_dir: Path | str | None = None,
) -> dict:
    """Run N-site chain scaling analysis with proper statistics.

    Uses the L_base + L_deph decomposition from enaqt.core for efficiency.
    Reports both linear fit AND saturation model with bootstrap confidence
    intervals, Durbin-Watson autocorrelation statistic, and Cook's distance.

    Parameters
    ----------
    n_max : maximum chain length (default 20).
    out_dir : directory for saving plots and JSON (default: auto-resolve).

    Returns
    -------
    dict with:
        - fit results (linear + saturation) with params, R^2, CIs, AIC
        - diagnostics (Durbin-Watson, Cook's D)
        - per-topology results
        - preferred model
    """
    if out_dir is None:
        out_dir = get_out_dir()
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    chain_sizes = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]
    chain_sizes = [n for n in chain_sizes if n <= n_max]
    if not chain_sizes:
        raise ValueError(f"n_max={n_max} is too small; need at least n=2")

    gp_arr = np.logspace(-3, 3, 200)
    results: dict = {"funnel": {}, "flat": {}, "disordered": {}}

    print(f"\n[nsite] Computing yields for {len(chain_sizes)} chain lengths x 3 topologies...")

    for N in chain_sizes:
        for topology, H in [
            ("funnel", hamiltonian_funnel(N, _DELTA)),
            ("flat", hamiltonian_flat(N, _DELTA)),
            ("disordered", hamiltonian_disordered(N, _DELTA, _DISORDER_SIGMA)),
        ]:
            gp, eta = optimal_dephasing(
                *liouvillian_parts(H, _KAPPA, _GAMMA),
                _KAPPA, N, n_points=200,
            )
            idx = int(np.argmax(eta))
            eta_peak = float(eta[idx])
            eta_zero = float(eta[0])
            eta_zeno = float(eta[-1])
            enhancement = eta_peak / eta_zero if eta_zero > 1e-16 else float("nan")

            results[topology][N] = {
                "N": N,
                "topology": topology,
                "gp": gp.tolist(),
                "eta": eta.tolist(),
                "gp_star": float(gp[idx]),
                "eta_peak": eta_peak,
                "eta_zero": eta_zero,
                "eta_zeno": eta_zeno,
                "enhancement": enhancement,
                "interior_peak": bool(0 < idx < len(eta) - 1),
            }

        r_f = results["funnel"][N]
        r_b = results["flat"][N]
        print(
            f"  N={N:2d}: funnel {r_f['enhancement']:.1f}x "
            f"(gp*={r_f['gp_star']:.3f})  |  "
            f"flat {r_b['enhancement']:.1f}x  |  "
            f"disordered {results['disordered'][N]['enhancement']:.1f}x"
        )

    # ---- Fit models ----
    Ns = np.array(chain_sizes)
    enh_funnel = np.array([results["funnel"][N]["enhancement"] for N in chain_sizes])
    fit_results = fit_models(Ns, enh_funnel, n_bootstrap=1000)

    # ---- Durbin-Watson interpretation ----
    dw_lin = fit_results["linear_durbin_watson"]
    dw_sat = fit_results["sat_durbin_watson"]
    fit_results["linear_dw_interpretation"] = (
        "positive autocorrelation" if dw_lin < 1.5 else
        "negative autocorrelation" if dw_lin > 2.5 else
        "no significant autocorrelation"
    )
    fit_results["sat_dw_interpretation"] = (
        "positive autocorrelation" if dw_sat < 1.5 else
        "negative autocorrelation" if dw_sat > 2.5 else
        "no significant autocorrelation"
    )

    # ---- Plot ----
    _plot_scaling(results, fit_results, chain_sizes, out_dir)

    # ---- JSON output ----
    output = {
        "experiment": "ENAQT N-Site Chain Scaling",
        "parameters": {
            "Delta": _DELTA,
            "kappa": _KAPPA,
            "Gamma": _GAMMA,
            "total_bias": _BIAS,
            "n_max": n_max,
        },
        "fit": fit_results,
        "funnel_results": {
            str(N): {k: v for k, v in results["funnel"][N].items() if k not in ("gp", "eta")}
            for N in chain_sizes
        },
        "flat_results": {
            str(N): {k: v for k, v in results["flat"][N].items() if k not in ("gp", "eta")}
            for N in chain_sizes
        },
        "disordered_results": {
            str(N): {k: v for k, v in results["disordered"][N].items() if k not in ("gp", "eta")}
            for N in chain_sizes
        },
    }
    save_json(output, out_dir, "enaqt_nsite_scaling.json")
    return output


def _plot_scaling(results: dict, fit: dict, chain_sizes: list, out_dir: Path) -> str:
    """Generate the 6-panel scaling figure."""
    fig = plt.figure(figsize=(18, 16))
    fig.suptitle(
        "ENAQT N-Site Chain Scaling  —  Saturation Model + Bootstrap 95% CIs",
        fontsize=15, fontweight="bold", y=0.99,
    )
    gs = gridspec.GridSpec(3, 2, hspace=0.48, wspace=0.38)
    gp_arr = np.logspace(-3, 3, 200)

    Ns = np.array(chain_sizes)
    n_sizes = len(chain_sizes)
    size_colors = plt.cm.turbo(np.linspace(0.05, 0.95, n_sizes))
    size_col_map = dict(zip(chain_sizes, size_colors))

    # ── (0,0): Bell curves ────────────────────────────────────────────
    ax0 = fig.add_subplot(gs[0, 0])
    for N in chain_sizes:
        r = results["funnel"][N]
        col = size_col_map[N]
        lw = 2.5 if N in (2, 5, 10, 20) else 1.2
        ax0.plot(
            r["gp"], r["eta"], lw=lw, color=col,
            label=f"N={N:2d}  ({r['enhancement']:.1f}x)"
            if N in (2, 5, 10, 15, 20) else None,
        )
        if r["interior_peak"]:
            ax0.plot(
                r["gp_star"], r["eta_peak"], "o", color=col, ms=6,
                markeredgecolor="white", markeredgewidth=0.8, zorder=5,
            )
    ax0.set_xscale("log")
    ax0.set_xlabel("Dephasing Rate  gamma_phi  [Delta]")
    ax0.set_ylabel("Transfer Yield  eta_inf")
    ax0.set_title("ENAQT Bell Curves  (Energy Funnel, bias = 5 Delta)")
    ax0.legend(title="N (enhancement)", fontsize=8, ncol=2)
    ax0.grid(True)
    ax0.text(
        0.03, 0.97,
        "Arrows show optimal gamma_phi* shifting left\n"
        "as N increases (longer chains need gentler noise)",
        transform=ax0.transAxes, fontsize=8.5, color="#444", va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
    )

    # ── (0,1): Enhancement vs N with BOTH models + CIs ────────────────
    ax1 = fig.add_subplot(gs[0, 1])
    enh_funnel = np.array([results["funnel"][N]["enhancement"] for N in chain_sizes])
    enh_flat = np.array([results["flat"][N]["enhancement"] for N in chain_sizes])
    enh_dis = np.array([results["disordered"][N]["enhancement"] for N in chain_sizes])

    ax1.plot(Ns, enh_funnel, "o-", color="#E74C3C", lw=2.5, ms=8,
             markeredgecolor="white", markeredgewidth=1.2,
             label="Energy funnel  (bias=5Delta)")
    ax1.plot(Ns, enh_dis, "s-", color="#9B59B6", lw=2.0, ms=7,
             markeredgecolor="white", markeredgewidth=1.0,
             label="Disordered  (sigma=2Delta)")
    ax1.plot(Ns, enh_flat, "^-", color="#2980B9", lw=2.0, ms=7,
             markeredgecolor="white", markeredgewidth=1.0,
             label="Flat chain  (no bias)")

    # Linear fit with CI
    x_fine = np.array(fit["x_fine"])
    lp = fit["linear_params"]
    ax1.plot(
        x_fine, lp["a"] * x_fine + lp["b"], "--", color="#E74C3C", lw=1.5,
        alpha=0.7, label=f"Linear: {lp['a']:.2f}N + {lp['b']:.2f}  (R\u00b2={fit['linear_r2']:.3f})",
    )
    ax1.fill_between(
        x_fine, fit["linear_ci_lower"], fit["linear_ci_upper"],
        color="#E74C3C", alpha=0.12,
    )

    # Saturation fit with CI
    sp = fit["sat_params"]
    sat_curve = sp["a"] * (x_fine ** sp["b"]) * np.exp(-sp["c"] * x_fine)
    ax1.plot(
        x_fine, sat_curve, "-.", color="#2ECC71", lw=1.8,
        alpha=0.8, label=f"Saturation: aN^b exp(-cN)  (R\u00b2={fit['sat_r2']:.3f})",
    )
    ax1.fill_between(
        x_fine, fit["sat_ci_lower"], fit["sat_ci_upper"],
        color="#2ECC71", alpha=0.12,
    )

    ax1.set_xlabel("Chain Length  N")
    ax1.set_ylabel("ENAQT Enhancement  (peak / zero-dephasing)")
    ax1.set_title(
        f"Enhancement vs N  |  Preferred: {fit['preferred_model']}  (AIC)\n"
        f"DW(linear)={fit['linear_durbin_watson']:.2f}, "
        f"DW(sat)={fit['sat_durbin_watson']:.2f}"
    )
    ax1.legend(fontsize=8)
    ax1.grid(True)

    # ── (1,0): Optimal gamma_phi* vs N ────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    gp_star_fun = [results["funnel"][N]["gp_star"] for N in chain_sizes]
    gp_star_flat = [results["flat"][N]["gp_star"] for N in chain_sizes]

    ax2.plot(Ns, gp_star_fun, "o-", color="#E74C3C", lw=2.2, ms=8,
             markeredgecolor="white", markeredgewidth=1.2, label="Energy funnel")
    ax2.plot(Ns, gp_star_flat, "^-", color="#2980B9", lw=2.0, ms=7,
             markeredgecolor="white", markeredgewidth=1.0, label="Flat chain")

    log_N = np.log(Ns)
    log_gp = np.log(np.array(gp_star_fun))
    fit2 = np.polyfit(log_N, log_gp, 1)
    Ns_fine = np.linspace(Ns.min(), Ns.max(), 100)
    ax2.plot(
        Ns_fine, np.exp(np.polyval(fit2, np.log(Ns_fine))), "--",
        color="#E74C3C", lw=1.3, alpha=0.7,
        label=f"Power law: ~ N^{fit2[0]:.2f}",
    )

    ax2.set_xlabel("Chain Length  N")
    ax2.set_ylabel("Optimal Dephasing  gamma_phi*  [Delta]")
    ax2.set_title("Longer Chains Need Gentler Noise")
    ax2.legend(fontsize=9)
    ax2.grid(True)
    ax2.set_yscale("log")

    # ── (1,1): Population flow at optimal gamma_phi (N=7) ─────────────
    ax3 = fig.add_subplot(gs[1, 1])
    N_demo = 7
    H_demo = hamiltonian_funnel(N_demo, _DELTA)
    gp_opt = results["funnel"][N_demo]["gp_star"]
    t, pops, eta_t = time_trajectory(
        H_demo, gp_opt, _KAPPA, _GAMMA, t_span=(0, 60.0)
    )

    site_colors = plt.cm.RdYlGn_r(np.linspace(0.05, 0.95, N_demo))
    for j in range(N_demo):
        lw = 2.5 if j in (0, N_demo - 1) else 1.3
        label = (
            f"Site {j + 1} (donor)" if j == 0 else
            f"Site {j + 1} (RC)" if j == N_demo - 1 else
            f"Site {j + 1}"
        )
        ax3.plot(t, pops[j], lw=lw, color=site_colors[j], label=label)
    ax3.plot(t, eta_t, "k-", lw=2.5, label=f"eta(t) yield (final={eta_t[-1]:.3f})")

    ax3.set_xlabel("Time  [Delta^-1]")
    ax3.set_ylabel("Site Population / Cumulative Yield")
    ax3.set_ylim(-0.02, 1.02)
    ax3.set_title(
        f"Population Flow  (N={N_demo}, gamma_phi*={gp_opt:.3f})\n"
        "Energy flows site 1 -> site 7 -> RC"
    )
    ax3.legend(fontsize=8, ncol=2)
    ax3.grid(True)

    # ── (2,0): Enhancement heat map ────────────────────────────────────
    ax4 = fig.add_subplot(gs[2, 0])

    Ns_heat = np.arange(2, min(n_max + 1, 21))
    gp_heat = np.logspace(-3, 3, 80)
    eta_zero_map = {}
    ETA_MAP = np.zeros((len(Ns_heat), len(gp_heat)))

    for i, Nh in enumerate(Ns_heat):
        H = hamiltonian_funnel(Nh, _DELTA)
        eta0 = analytical_yield(H, 1e-3, _KAPPA, _GAMMA)
        eta_zero_map[Nh] = eta0
        for j, gph in enumerate(gp_heat):
            ETA_MAP[i, j] = analytical_yield(H, gph, _KAPPA, _GAMMA)

    # Normalize each row by zero-dephasing value
    ENH_MAP = np.zeros_like(ETA_MAP)
    for i, Nh in enumerate(Ns_heat):
        e0 = max(eta_zero_map.get(Nh, 1e-8), 1e-8)
        ENH_MAP[i, :] = ETA_MAP[i, :] / e0

    im = ax4.pcolormesh(gp_heat, Ns_heat, ENH_MAP, cmap="hot", shading="auto")
    cb = plt.colorbar(im, ax=ax4, shrink=0.9)
    cb.set_label("ENAQT Enhancement  (peak / zero-deph)")

    ridge = [gp_heat[int(np.argmax(ENH_MAP[i, :]))] for i in range(len(Ns_heat))]
    ax4.plot(ridge, Ns_heat, "w--", lw=2.0, label="Optimal gamma_phi* (ridge)")

    ax4.set_xscale("log")
    ax4.set_xlabel("Dephasing Rate  gamma_phi  [Delta]")
    ax4.set_ylabel("Chain Length  N")
    ax4.set_title("Enhancement Heat Map  (Energy Funnel)\nWhite dashed = ENAQT ridge")
    ax4.legend(fontsize=9, loc="upper right")

    # ── (2,1): Diagnostics text panel ───────────────────────────────────
    ax5 = fig.add_subplot(gs[2, 1])
    ax5.set_xlim(0, 1)
    ax5.set_ylim(0, 1)
    ax5.axis("off")

    lp = fit["linear_params"]
    sp = fit["sat_params"]
    lin_cooks = np.array(fit["linear_cooks_d"])
    sat_cooks = np.array(fit["sat_cooks_d"])
    n_thresh = 4.0 / len(chain_sizes)

    diag_text = (
        "Model Comparison\n"
        "=" * 40 + "\n\n"
        f"Linear:        E(N) = {lp['a']:.3f} N + {lp['b']:.3f}\n"
        f"  R\u00b2 = {fit['linear_r2']:.4f}  |  "
        f"AIC = {fit['linear_aic']:.2f}\n"
        f"  Durbin-Watson = {fit['linear_durbin_watson']:.3f}\n"
        f"  ({fit['linear_dw_interpretation']})\n"
        f"  Cook's D > {n_thresh:.2f}: "
        f"{np.sum(lin_cooks > n_thresh)} of {len(chain_sizes)} points\n\n"
        f"Saturation:    E(N) = {sp['a']:.3f} * N^{sp['b']:.3f} * "
        f"exp(-{sp['c']:.3f}*N)\n"
        f"  R\u00b2 = {fit['sat_r2']:.4f}  |  "
        f"AIC = {fit['sat_aic']:.2f}\n"
        f"  Durbin-Watson = {fit['sat_durbin_watson']:.3f}\n"
        f"  ({fit['sat_dw_interpretation']})\n"
        f"  Cook's D > {n_thresh:.2f}: "
        f"{np.sum(sat_cooks > n_thresh)} of {len(chain_sizes)} points\n\n"
        f"Preferred model: {fit['preferred_model'].upper()} (lower AIC)\n"
        "=" * 40 + "\n\n"
        "Key Findings:\n"
        f"  N=2  (spin-boson):  {results['funnel'][2]['enhancement']:.1f}x\n"
    )
    if 7 in results["funnel"]:
        diag_text += (
            f"  N=7  (FMO-scale):   {results['funnel'][7]['enhancement']:.1f}x\n"
        )
    if 15 in results["funnel"]:
        diag_text += (
            f"  N=15 (microtubule?): {results['funnel'][15]['enhancement']:.1f}x\n"
        )
    if 20 in results["funnel"]:
        diag_text += (
            f"  N=20:               {results['funnel'][20]['enhancement']:.1f}x\n"
        )
    diag_text += (
        f"\nOptimal gp* ~ N^{fit2[0]:.2f} (power law)\n"
        "Shaded regions = 95% bootstrap CI"
    )

    ax5.text(0.05, 0.95, diag_text, transform=ax5.transAxes, fontsize=9,
             va="top", fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f8f8", alpha=0.95))

    out_path = out_dir / "enaqt_nsite_scaling.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {out_path}")
    plt.close()
    return str(out_path)


# ============================================================================
#  run_fmo_benchmark
# ============================================================================

def run_fmo_benchmark(
    out_dir: Path | str | None = None,
) -> dict:
    """Run FMO-7 benchmark with explicit caveats about Lindblad limitations.

    The FMO complex is a 7-site pigment-protein complex. This benchmark
    compares the actual FMO Hamiltonian against synthetic funnel/flat
    chains of the same size.

    IMPORTANT CAVEATS:
    - Lindblad dynamics provides only qualitative results for FMO.
      The full biological system involves non-Markovian environments,
      structured spectral densities, and finite-temperature effects
      that are not captured by a simple dephasing Lindblad model.
    - Any reported enhancement ratio is a model prediction, not
      experimentally validated.
    - This is benchmarking of the Lindblad computational framework,
      not biological validation of ENAQT in photosynthesis.

    Parameters
    ----------
    out_dir : directory for saving plot and JSON (default: auto-resolve).

    Returns
    -------
    dict with results + caveat strings.
    """
    if out_dir is None:
        out_dir = get_out_dir()
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    CAVEATS = [
        "Lindblad provides only qualitative results for FMO",
        "Reported enhancement is a model prediction, not experimentally validated",
        "This is benchmarking of the Lindblad framework, not biological validation",
    ]

    print("\n[nsite] Computing FMO-7 benchmark...")

    H_fmo = hamiltonian_fmo7()
    Lb_fmo, Ld_fmo = liouvillian_parts(H_fmo, _KAPPA, _GAMMA)
    gp_fmo, eta_fmo = optimal_dephasing(Lb_fmo, Ld_fmo, _KAPPA, 7, n_points=200)

    idx_fmo = int(np.argmax(eta_fmo))
    p_fmo = {
        "gp_star": float(gp_fmo[idx_fmo]),
        "eta_peak": float(eta_fmo[idx_fmo]),
        "eta_zero": float(eta_fmo[0]),
        "eta_zeno": float(eta_fmo[-1]),
        "enhancement": (
            float(eta_fmo[idx_fmo] / eta_fmo[0]) if eta_fmo[0] > 1e-16 else float("nan")
        ),
        "interior_peak": bool(0 < idx_fmo < len(eta_fmo) - 1),
    }

    # Also compute funnel N=7 for comparison
    H_fun7 = hamiltonian_funnel(7, _DELTA)
    Lb_fun, Ld_fun = liouvillian_parts(H_fun7, _KAPPA, _GAMMA)
    gp_fun, eta_fun = optimal_dephasing(Lb_fun, Ld_fun, _KAPPA, 7, n_points=200)

    idx_fun = int(np.argmax(eta_fun))
    p_fun7 = {
        "gp_star": float(gp_fun[idx_fun]),
        "eta_peak": float(eta_fun[idx_fun]),
        "eta_zero": float(eta_fun[0]),
        "enhancement": (
            float(eta_fun[idx_fun] / eta_fun[0]) if eta_fun[0] > 1e-16 else float("nan")
        ),
    }

    # Absolute yields
    abs_fmo_peak = p_fmo["eta_peak"]
    abs_fmo_zero = p_fmo["eta_zero"]
    abs_fun_peak = p_fun7["eta_peak"]
    abs_fun_zero = p_fun7["eta_zero"]

    print(f"  FMO-7:  enhancement={p_fmo['enhancement']:.2f}x, "
          f"peak_yield={abs_fmo_peak:.4f}, zero_yield={abs_fmo_zero:.4f}")
    print(f"  Funnel-7: enhancement={p_fun7['enhancement']:.2f}x, "
          f"peak_yield={abs_fun_peak:.4f}, zero_yield={abs_fun_zero:.4f}")

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.plot(
        gp_fmo, eta_fmo, lw=2.5, color="#F39C12",
        label=f"FMO-7 (actual, {p_fmo['enhancement']:.1f}x, "
              f"yield={abs_fmo_peak:.3f})",
    )
    ax.plot(
        gp_fun, eta_fun, lw=2.5, color="#E74C3C",
        label=f"Funnel N=7 (bias=5, {p_fun7['enhancement']:.1f}x, "
              f"yield={abs_fun_peak:.3f})",
    )

    # Mark peaks
    ax.plot(
        p_fmo["gp_star"], p_fmo["eta_peak"], "*", color="#F39C12", ms=14,
        zorder=5, markeredgecolor="white",
    )
    ax.plot(
        p_fun7["gp_star"], p_fun7["eta_peak"], "*", color="#E74C3C", ms=14,
        zorder=5, markeredgecolor="white",
    )

    # Caveat text box
    caveat_text = "\n".join(f"  {i + 1}. {c}" for i, c in enumerate(CAVEATS))
    ax.text(
        0.97, 0.97,
        f"CAVEATS:\n{caveat_text}",
        transform=ax.transAxes, fontsize=8.5, va="top", ha="right",
        color="#8B0000",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFF0F0", alpha=0.95),
    )

    ax.set_xscale("log")
    ax.set_xlabel("Dephasing Rate  gamma_phi  [Delta]")
    ax.set_ylabel("Transfer Yield  eta")
    ax.set_title(
        "FMO-7 Benchmark  —  Biological vs Synthetic Hamiltonian\n"
        "(FMO Hamiltonian in units of J_12 = 87.7 cm^-1)"
    )
    ax.legend(fontsize=9)
    ax.grid(True)

    out_path = out_dir / "enaqt_fmo_benchmark.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {out_path}")
    plt.close()

    output = {
        "experiment": "FMO-7 Lindblad Benchmark",
        "caveats": CAVEATS,
        "fmo7": {
            **p_fmo,
            "absolute_peak_yield": abs_fmo_peak,
            "absolute_zero_yield": abs_fmo_zero,
        },
        "funnel_n7": {
            **p_fun7,
            "absolute_peak_yield": abs_fun_peak,
            "absolute_zero_yield": abs_fun_zero,
        },
    }
    save_json(output, out_dir, "enaqt_fmo_benchmark.json")
    return output
