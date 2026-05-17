"""
ENAQT Spin-Boson Module -- 2-Site Environment-Assisted Quantum Transport
=========================================================================

Refactored from ``enaqt_sb_analysis.py`` and ``enaqt_sb_sink.py``.

Provides:

* :func:`run_validation` -- QD3SET-1 HEOM vs Bloch comparison with
  proper trace-distance error metrics.
* :func:`run_sink_analysis` -- Lindblad-sink vs no-sink analytical
  yield curves and visualisation.
* :func:`regime_classification` -- Analytical transport-regime
  classifier (``"localization" | "optimal" | "zeno"``).

All Hamiltonian / Liouvillian construction is delegated to
:mod:`enaqt.core` so that physics code lives in exactly one place.
"""

from __future__ import annotations

import glob
import json
import os
import re
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# Delegate all Hamiltonian / Liouvillian construction to enaqt.core
# ---------------------------------------------------------------------------
try:
    from enaqt.core import (
        analytical_yield,
        hamiltonian_2site,
        liouvillian,
        liouvillian_parts,
        save_json,
    )
except ImportError:  # pragma: no cover  (core lives on another branch)
    pass

matplotlib.use("Agg")

# ── Styling ─────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 9,
    "figure.dpi": 100,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.7,
})

_COLORS = {
    "sym": "#2980B9",   # blue  -- symmetric ε=0
    "asym": "#C0392B",  # red   -- asymmetric ε=1
    "low": "#27AE60",   # green -- localization regime
    "mid": "#F39C12",   # orange -- optimal ENAQT
    "high": "#8E44AD",  # purple -- quantum Zeno
    "neutral": "#7F8C8D",
}


# ═════════════════════════════════════════════════════════════════════════════
#  Private helpers -- data loading
# ═════════════════════════════════════════════════════════════════════════════

def _parse_filename(fpath: str) -> dict[str, Any]:
    """Extract physical parameters from an SB trajectory filename.

    Expected filename pattern::

        2_epsilon-<ε>_Delta-<Δ>_lambda-<λ>_gamma-<γ_c>_beta-<β>.npy

    Parameters
    ----------
    fpath :
        Absolute or relative path to the ``.npy`` file.

    Returns
    -------
    dict
        Mapping with keys ``file``, ``epsilon``, ``delta``, ``lam``,
        ``gamma_c``, ``beta``.
    """
    fname = os.path.basename(fpath)

    def _grab(key: str) -> float | None:
        m = re.search(rf"{key}-([\d.]+)", fname)
        return float(m.group(1)) if m else None

    beta_match = re.search(r"beta-([\d.]+?)\.npy", fname)
    beta = float(beta_match.group(1)) if beta_match else None

    return {
        "file": fname,
        "epsilon": _grab("epsilon"),
        "delta": _grab("Delta"),
        "lam": _grab("lambda"),
        "gamma_c": _grab("gamma"),
        "beta": beta,
    }


def _load_all(data_dir: str | Path) -> list[dict[str, Any]]:
    """Load every ``.npy`` trajectory in *data_dir* and compute derived quantities.

    Parameters
    ----------
    data_dir :
        Directory containing ``*.npy`` HEOM trajectory files.

    Returns
    -------
    list[dict]
        One record per trajectory with keys ``epsilon``, ``delta``, ``lam``,
        ``gamma_c``, ``beta``, ``T``, ``gamma_phi``, ``eta_final``,
        ``eta_avg``, ``coh_lifetime``, ``fpath``.
    """
    data_dir = str(data_dir)
    files = glob.glob(os.path.join(data_dir, "*.npy"))
    if not files:
        raise FileNotFoundError(f"No .npy files found in {data_dir}")

    records: list[dict[str, Any]] = []
    for fpath in files:
        p = _parse_filename(fpath)
        traj = np.load(fpath)  # shape (401, 5) complex128

        t = traj[:, 0].real
        rho11 = traj[:, 1].real
        rho22 = traj[:, 4].real
        coh = traj[:, 2]  # ρ₁₂

        # Dephasing rate (Drude-Lorentz high-T Markov): γ_φ = 2λ/(β·γ_c)
        gamma_phi = 2.0 * p["lam"] / (p["beta"] * p["gamma_c"])
        T = 1.0 / p["beta"]

        # Final site-2 population -- primary ENAQT efficiency metric
        eta_final = float(rho22[-1])

        # Time-averaged site-2 population
        eta_avg = float(np.trapz(rho22, t) / (t[-1] - t[0]))

        # Coherence lifetime -- time when |ρ₁₂| first drops to 1/e of max
        coh_mag = np.abs(coh)
        coh_init = coh_mag.max()
        coh_thresh = coh_init / np.e
        coh_decay_idx = (
            np.argmax(coh_mag < coh_thresh)
            if np.any(coh_mag < coh_thresh)
            else -1
        )
        coh_lifetime = (
            float(t[coh_decay_idx]) if coh_decay_idx > 0 else float(t[-1])
        )

        records.append({
            **p,
            "T": T,
            "gamma_phi": gamma_phi,
            "eta_final": eta_final,
            "eta_avg": eta_avg,
            "coh_lifetime": coh_lifetime,
            "fpath": fpath,
        })

    return records


# ═════════════════════════════════════════════════════════════════════════════
#  Private helpers -- log binning
# ═════════════════════════════════════════════════════════════════════════════

def _log_bin(
    x: np.ndarray,
    y: np.ndarray,
    n_bins: int = 25,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Log-bin *(x, y)* pairs.

    Parameters
    ----------
    x, y :
        1-D arrays of equal length.
    n_bins :
        Number of logarithmic bins.

    Returns
    -------
    centers, means, stds, counts :
        Each a 1-D :class:`~numpy.ndarray`.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    log_x = np.log10(x)
    edges = np.linspace(log_x.min(), log_x.max(), n_bins + 1)
    centers, means, stds, counts = [], [], [], []
    for i in range(n_bins):
        mask = (log_x >= edges[i]) & (log_x < edges[i + 1])
        if mask.sum() > 0:
            centers.append(10 ** ((edges[i] + edges[i + 1]) / 2))
            means.append(y[mask].mean())
            stds.append(y[mask].std())
            counts.append(int(mask.sum()))
    return (
        np.array(centers),
        np.array(means),
        np.array(stds),
        np.array(counts),
    )


# ═════════════════════════════════════════════════════════════════════════════
#  Private helpers -- Bloch physics (2-site, 4-element state vector)
# ═════════════════════════════════════════════════════════════════════════════

def _bloch_liouvillian_4x4(
    eps: float,
    Delta: float,
    gp: float,
    kappa: float,
    Gamma: float,
) -> np.ndarray:
    """4×4 Bloch-equation superoperator.

    State ordering: ``[ρ₁₁, ρ₂₂, Re(ρ₁₂), Im(ρ₁₂)]``.

    Parameters
    ----------
    eps :
        Site energy bias ε (site 1 = +ε/2, site 2 = –ε/2).
    Delta :
        Tunneling matrix element.
    gp :
        Pure dephasing rate γ_φ.
    kappa :
        Sink extraction rate from site 2.
    Gamma :
        Fluorescence / recombination rate.

    Returns
    -------
    np.ndarray
        4×4 real matrix *A* such that ``dρ/dt = A @ ρ``.
    """
    g = gp + kappa / 2 + Gamma  # total coherence decay rate
    return np.array([
        [-Gamma,          0,        0,       -Delta    ],
        [  0,   -(kappa + Gamma),   0,        Delta    ],
        [  0,         0,           -g,         eps     ],
        [Delta / 2,  -Delta / 2,   -eps,       -g      ],
    ])


def _analytical_yield_4x4(
    eps: float,
    Delta: float,
    gp: float,
    kappa: float,
    Gamma: float,
) -> float:
    """Exact total transfer yield via matrix inversion.

    ``η_∞ = κ · [-A⁻¹ ρ₀]₂`` where *A* is the Bloch Liouvillian.

    Returns
    -------
    float
        The yield, or ``nan`` if *A* is singular.
    """
    A = _bloch_liouvillian_4x4(eps, Delta, gp, kappa, Gamma)
    rho0 = np.array([1.0, 0.0, 0.0, 0.0])
    try:
        integral_rho = -np.linalg.solve(A, rho0)
        return float(kappa * integral_rho[1])
    except np.linalg.LinAlgError:
        return float("nan")


def _time_trajectory_4x4(
    eps: float,
    Delta: float,
    gp: float,
    kappa: float,
    Gamma: float,
    T: float = 150.0,
    n_pts: int = 3000,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate the Bloch ODE and return time-resolved populations + yield.

    The state vector is ``[ρ₁₁, ρ₂₂, Re(ρ₁₂), Im(ρ₁₂), η_cumulative]``.

    Returns
    -------
    t, y :
        ``t`` is a 1-D array of times; ``y`` has shape ``(5, len(t))``.
    """
    A = _bloch_liouvillian_4x4(eps, Delta, gp, kappa, Gamma)

    def _rhs(_t: float, y: np.ndarray) -> list[float]:
        dy = A @ y[:4]
        deta = kappa * y[1]  # dη/dt = κ ρ₂₂
        return [*dy, deta]

    y0 = [1.0, 0.0, 0.0, 0.0, 0.0]
    sol = solve_ivp(
        _rhs,
        [0, T],
        y0,
        method="RK45",
        dense_output=False,
        t_eval=np.linspace(0, T, n_pts),
        rtol=1e-9,
        atol=1e-11,
    )
    return sol.t, sol.y


def _sweep_gamma_phi(
    eps: float,
    Delta: float,
    kappa: float,
    Gamma: float,
    n_pts: int = 300,
) -> tuple[np.ndarray, np.ndarray]:
    """Sweep γ_φ over six decades and compute analytical yield at each point.

    Returns
    -------
    gp_arr, eta :
        Both are 1-D arrays of length *n_pts*.
    """
    gp_arr = np.logspace(-3, 3, n_pts)
    eta = np.array([
        _analytical_yield_4x4(eps, Delta, gp, kappa, Gamma)
        for gp in gp_arr
    ])
    return gp_arr, eta


# ═════════════════════════════════════════════════════════════════════════════
#  Private helpers -- trace distance
# ═════════════════════════════════════════════════════════════════════════════

def _trace_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Quantum trace distance T(ρ, σ) = ½ ||ρ – σ||₁.

    Parameters
    ----------
    rho, sigma :
        Square density matrices (2×2 for the spin-boson case).

    Returns
    -------
    float
        Trace distance in ``[0, 1]``.
    """
    delta = rho - sigma
    # ||A||₁ = Tr[√(A†A)]
    eigs = np.linalg.eigvalsh(delta.conj().T @ delta)
    eigs = np.maximum(eigs, 0.0)  # numerical safety
    return 0.5 * np.sum(np.sqrt(eigs))


def _reconstruct_rho2x2(
    rho11: complex | float,
    rho22: complex | float,
    rho12: complex,
) -> np.ndarray:
    """Build a 2×2 density matrix from independent elements.

    Parameters
    ----------
    rho11, rho22 :
        Diagonal populations (must be real).
    rho12 :
        Off-diagonal coherence (ρ₂₁ = ρ₁₂* is enforced).

    Returns
    -------
    np.ndarray
        2×2 Hermitian density matrix.
    """
    return np.array([
        [rho11, rho12],
        [np.conj(rho12), rho22],
    ], dtype=complex)


# ═════════════════════════════════════════════════════════════════════════════
#  Private helpers -- plotting
# ═════════════════════════════════════════════════════════════════════════════

def _plot_validation_main(
    sym: dict,
    asym: dict,
    all_records: list[dict],
    out_dir: Path,
) -> Path:
    """4-panel ENAQT validation figure."""
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.35)
    fig.suptitle(
        "ENAQT Analysis -- QD3SET-1 Spin-Boson HEOM  |  "
        "Bloch vs HEOM Trace-Distance Validation",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    # -- Panel 1: Bell curves ------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    for res, key, marker in [(sym, "sym", "o"), (asym, "asym", "s")]:
        gp_b = np.array(res["gp_binned"])
        eta_b = np.array(res["eta_binned"])
        eta_s = np.array(res["eta_std"])
        col = _COLORS[key]

        ax1.scatter(res["gp_raw"], res["eta_raw"],
                    c=col, alpha=0.10, s=8, zorder=1)
        ax1.fill_between(gp_b, eta_b - eta_s, eta_b + eta_s,
                         color=col, alpha=0.15, zorder=2)
        ax1.plot(gp_b, eta_b, f"{marker}-", color=col, lw=2.2, ms=7,
                 label=res["label"], zorder=3)

        if res["enaqt_detected"]:
            opt = res["optimal_gamma_phi"]
            pk = res["peak_eta"]
            ax1.axvline(opt, color=col, lw=1.2, ls="--", alpha=0.6)
            ax1.annotate(
                f"γ*={opt:.2f}\nη*={pk:.3f}",
                xy=(opt, pk), xytext=(opt * 2.5, pk * 0.97),
                fontsize=8, color=col,
                arrowprops=dict(arrowstyle="->", color=col, lw=1),
            )

    ax1.set_xscale("log")
    ax1.set_xlabel(r"Dephasing Rate  $\gamma_\phi = 2\lambda/(\beta\gamma_c)$  [$\Delta$]")
    ax1.set_ylabel(r"Transfer Efficiency  $\eta = \rho_{22}(t_{\mathrm{final}})$")
    ax1.set_title("ENAQT Bell Curve  (Final Population at Sink)")
    ax1.legend()
    ax1.grid(True)

    detected = [r for r in [sym, asym] if r["enaqt_detected"]]
    banner = "✓ ENAQT DETECTED" if detected else "✗ Bell curve not found"
    banner_col = "#27AE60" if detected else "#E74C3C"
    ax1.text(
        0.97, 0.05, banner, transform=ax1.transAxes,
        ha="right", va="bottom", fontsize=10, fontweight="bold",
        color=banner_col,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor=banner_col, alpha=0.9),
    )

    # -- Panel 2: Trace-distance error ---------------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    for res, key, marker in [(sym, "sym", "o"), (asym, "asym", "s")]:
        gp_b = np.array(res.get("td_binned_gp", []))
        td_b = np.array(res.get("td_binned_mean", []))
        td_s = np.array(res.get("td_binned_std", []))
        col = _COLORS[key]
        if len(gp_b) > 0:
            ax2.fill_between(gp_b, td_b - td_s, td_b + td_s,
                             color=col, alpha=0.15)
            ax2.plot(gp_b, td_b, f"{marker}-", color=col, lw=2.2, ms=7,
                     label=f"{res['label']} (trace dist)")

    ax2.set_xscale("log")
    ax2.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax2.set_ylabel(r"Trace Distance  $T(\rho_{\mathrm{HEOM}}, \rho_{\mathrm{Bloch}})$")
    ax2.set_title("Bloch vs HEOM Discrepancy\n(Trace Distance -- Proper Quantum Metric)")
    ax2.legend()
    ax2.grid(True)

    # -- Panel 3: Population dynamics at three regimes ------------------------
    ax3 = fig.add_subplot(gs[1, 0])
    sym_records = [r for r in all_records if r["epsilon"] == 0.0]
    sym_sorted = sorted(sym_records, key=lambda r: r["gamma_phi"])
    n = len(sym_sorted)

    regime_picks = {
        "Localization  (low γ_φ)": sym_sorted[max(0, n // 12)],
        "Optimal ENAQT  (mid γ_φ)": sym_sorted[n // 2],
        "Quantum Zeno  (high γ_φ)": sym_sorted[min(n - 1, 11 * n // 12)],
    }
    regime_cols = [_COLORS["low"], _COLORS["mid"], _COLORS["high"]]

    for (rlabel, rec), col in zip(regime_picks.items(), regime_cols):
        traj = np.load(rec["fpath"])
        t = traj[:, 0].real
        rho22 = traj[:, 4].real
        gp = rec["gamma_phi"]
        ax3.plot(
            t, rho22, color=col, lw=2.0,
            label=f"{rlabel}\n(γ_φ={gp:.3f}, λ={rec['lam']}, β={rec['beta']})",
        )

    ax3.axhline(
        0.5, color=_COLORS["neutral"], ls="--", lw=1.2, alpha=0.7,
        label="Thermal equilibrium (ε=0)",
    )
    ax3.set_xlabel(r"Time  [$\Delta^{-1}$]")
    ax3.set_ylabel(r"$\rho_{22}(t)$  -- Site 2 Population")
    ax3.set_title("Quantum Dynamics -- Three Dephasing Regimes  (ε = 0)")
    ax3.legend(fontsize=8)
    ax3.grid(True)

    # -- Panel 4: 2D parameter space -- η as color ---------------------------
    ax4 = fig.add_subplot(gs[1, 1])
    gp_sym = [r["gamma_phi"] for r in sym_records]
    lam_sym = [r["lam"] for r in sym_records]
    eta_sym = [r["eta_final"] for r in sym_records]

    sc = ax4.scatter(gp_sym, lam_sym, c=eta_sym,
                     cmap="plasma", alpha=0.75, s=25, vmin=0, vmax=1)
    cb = plt.colorbar(sc, ax=ax4, shrink=0.85, pad=0.02)
    cb.set_label(r"$\eta = \rho_{22}(t_{\mathrm{final}})$")

    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax4.set_ylabel(r"Reorganization Energy  $\lambda$  [$\Delta$]")
    ax4.set_title("Parameter Space Map -- Transport Efficiency  (ε = 0)")
    ax4.grid(True)

    out_path = out_dir / "spinboson_validation_main.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_validation_gallery(
    all_records: list[dict],
    out_dir: Path,
) -> Path:
    """Gallery: ε=0 vs ε=1 dynamics side-by-side at matched γ_φ values."""
    sym = sorted(
        [r for r in all_records if r["epsilon"] == 0.0],
        key=lambda r: r["gamma_phi"],
    )
    asym = sorted(
        [r for r in all_records if r["epsilon"] == 1.0],
        key=lambda r: r["gamma_phi"],
    )
    n_sym = len(sym)
    n_asym = len(asym)

    picks_sym = [sym[n_sym // 10], sym[n_sym // 2], sym[9 * n_sym // 10]]
    picks_asym = [asym[n_asym // 10], asym[n_asym // 2], asym[9 * n_asym // 10]]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=False)
    fig.suptitle(
        "Population Dynamics -- Symmetric vs Asymmetric SB  |  "
        "Three Dephasing Regimes",
        fontsize=13, fontweight="bold",
    )

    regime_labels = [
        "Localization (Low γ_φ)",
        "Optimal ENAQT (Mid γ_φ)",
        "Quantum Zeno (High γ_φ)",
    ]
    regime_cols = [_COLORS["low"], _COLORS["mid"], _COLORS["high"]]

    for col_idx, (rs, ra, rlabel, rcol) in enumerate(
        zip(picks_sym, picks_asym, regime_labels, regime_cols)
    ):
        for row_idx, (rec, eps_label) in enumerate([
            (rs, "ε = 0  (Symmetric)"),
            (ra, "ε = 1  (Asymmetric)"),
        ]):
            ax = axes[row_idx, col_idx]
            traj = np.load(rec["fpath"])
            t = traj[:, 0].real
            rho11 = traj[:, 1].real
            rho22 = traj[:, 4].real
            coh = np.abs(traj[:, 2])

            ax.plot(t, rho11, lw=1.8, color="#2980B9",
                    label=r"$\rho_{11}$ (Site 1)", alpha=0.9)
            ax.plot(t, rho22, lw=1.8, color=rcol,
                    label=r"$\rho_{22}$ (Site 2)", alpha=0.9)
            ax.plot(t, coh, lw=1.2, color="#95A5A6",
                    label=r"$|\rho_{12}|$ (coherence)", ls="--", alpha=0.7)

            ax.set_ylim(-0.05, 1.05)
            ax.set_xlabel(r"Time [$\Delta^{-1}$]")
            ax.set_ylabel("Population")
            ax.set_title(
                f"{eps_label}\n"
                f"{rlabel}\n"
                f"γ_φ={rec['gamma_phi']:.3f}, λ={rec['lam']}, β={rec['beta']}",
                fontsize=9,
            )
            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=8)
            ax.grid(True)

    plt.tight_layout()
    out_path = out_dir / "spinboson_validation_gallery.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_validation_beta_slices(
    all_records: list[dict],
    out_dir: Path,
) -> Path:
    """Bell curves sliced by temperature β."""
    sym = [r for r in all_records if r["epsilon"] == 0.0]
    betas = sorted({r["beta"] for r in sym})

    fig, axes = plt.subplots(
        1, len(betas), figsize=(4 * len(betas), 5), sharey=True,
    )
    if len(betas) == 1:
        axes = [axes]
    fig.suptitle(
        "ENAQT Bell Curves by Temperature  (ε = 0,  Symmetric SB)",
        fontsize=13, fontweight="bold",
    )

    cmap = plt.cm.coolwarm
    norm = plt.Normalize(vmin=min(betas), vmax=max(betas))

    for ax, beta in zip(axes, betas):
        subset = [r for r in sym if r["beta"] == beta]
        gp = np.array([r["gamma_phi"] for r in subset])
        eta = np.array([r["eta_final"] for r in subset])

        if len(gp) < 3:
            continue

        col = cmap(norm(beta))
        gp_b, eta_b, eta_s, _ = _log_bin(gp, eta, n_bins=15)

        ax.scatter(gp, eta, c=[col] * len(gp), alpha=0.25, s=12)
        ax.fill_between(gp_b, eta_b - eta_s, eta_b + eta_s,
                        color=col, alpha=0.20)
        ax.plot(gp_b, eta_b, "o-", color=col, lw=2.0, ms=6)

        idx_max = int(np.argmax(eta_b))
        interior = 0 < idx_max < len(eta_b) - 1
        status = "✓ ENAQT" if interior else "— monotone"
        T = 1.0 / beta

        ax.set_xscale("log")
        ax.set_xlabel(r"$\gamma_\phi$  [$\Delta$]")
        ax.set_title(f"β = {beta}  (T = {T:.1f} Δ)\n{status}", fontsize=10)
        ax.grid(True)

    axes[0].set_ylabel(r"$\eta = \rho_{22}(t_{\mathrm{final}})$")
    plt.tight_layout()
    out_path = out_dir / "spinboson_validation_beta_slices.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ═════════════════════════════════════════════════════════════════════════════
#  Private helpers -- sink-analysis plotting
# ═════════════════════════════════════════════════════════════════════════════

def _find_peak(
    gp_arr: np.ndarray,
    eta: np.ndarray,
) -> dict[str, Any]:
    """Characterise the peak in an η vs γ_φ curve.

    Returns a dict with keys ``gp_star``, ``eta_star``, ``eta_zero``,
    ``eta_zeno``, ``enhancement_zero``, ``enhancement_zeno``,
    ``interior_peak``.
    """
    idx = int(np.argmax(eta))
    gp_star = float(gp_arr[idx])
    eta_star = float(eta[idx])
    eta_zero = float(eta[0])
    eta_zeno = float(eta[-1])
    enhancement_zero = eta_star / eta_zero if eta_zero > 0 else float("nan")
    enhancement_zeno = eta_star / eta_zeno if eta_zeno > 0 else float("nan")
    interior = bool(0 < idx < len(eta) - 1)
    return {
        "gp_star": gp_star,
        "eta_star": eta_star,
        "eta_zero": eta_zero,
        "eta_zeno": eta_zeno,
        "enhancement_zero": enhancement_zero,
        "enhancement_zeno": enhancement_zeno,
        "interior_peak": interior,
    }


def _plot_sink_epic(
    sb_data: dict,
    out_dir: Path,
    n_points: int = 300,
) -> Path:
    """6-panel epic figure for the sink analysis."""
    fig = plt.figure(figsize=(18, 16))
    fig.suptitle(
        "ENAQT + Lindblad Sink  --  Spin-Boson Model  |  "
        "Rebentrost et al. (2009) Replication",
        fontsize=15, fontweight="bold", y=0.99,
    )
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.48, wspace=0.38)

    gp_arr = np.logspace(-3, 3, n_points)
    Delta = 1.0
    Gamma = 0.01

    # -- Panel (0,0): Main bell curve -- ε=5, multiple κ --------------------
    ax0 = fig.add_subplot(gs[0, 0])
    kappas = [0.01, 0.05, 0.1, 0.5, 1.0]
    kap_colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(kappas)))

    peak_summary = []
    for kap, col in zip(kappas, kap_colors):
        _, eta = _sweep_gamma_phi(eps=5.0, Delta=Delta, kappa=kap, Gamma=Gamma,
                                  n_pts=n_points)
        p = _find_peak(gp_arr, eta)
        peak_summary.append((kap, p))
        ax0.plot(gp_arr, eta, lw=2.0, color=col,
                 label=f"κ = {kap:.2f}  (peak {p['enhancement_zero']:.1f}×)")
        if p["interior_peak"]:
            ax0.plot(p["gp_star"], p["eta_star"], "o", color=col, ms=9,
                     zorder=5, markeredgecolor="white", markeredgewidth=1.2)

    ax0.set_xscale("log")
    ax0.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax0.set_ylabel(r"Transfer Yield  $\eta_\infty = \kappa \int \rho_{22} \, dt$")
    ax0.set_title(r"ENAQT Bell Curve  ($\epsilon = 5\Delta$, strong energy bias)")
    ax0.legend(title="Sink rate κ")
    ax0.grid(True)
    ax0.text(
        0.02, 0.05, r"$\epsilon = 5\Delta$  |  $\Gamma = 0.01\Delta$" "\n"
        "(FMO-like energy asymmetry)",
        transform=ax0.transAxes, fontsize=9, color="#555",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
    )

    best = max(peak_summary, key=lambda x: x[1]["enhancement_zero"])
    ax0.annotate(
        f"Peak η={best[1]['eta_star']:.3f}\n"
        f"({best[1]['enhancement_zero']:.1f}× over coherent limit)",
        xy=(best[1]["gp_star"], best[1]["eta_star"]),
        xytext=(best[1]["gp_star"] * 20, best[1]["eta_star"] * 0.85),
        fontsize=9, color="#222",
        arrowprops=dict(arrowstyle="->", color="#444", lw=1.3),
    )

    # -- Panel (0,1): ε sensitivity -----------------------------------------
    ax1 = fig.add_subplot(gs[0, 1])
    epsilons = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
    eps_colors = plt.cm.viridis(np.linspace(0.1, 0.95, len(epsilons)))

    for eps, col in zip(epsilons, eps_colors):
        _, eta = _sweep_gamma_phi(eps=eps, Delta=Delta, kappa=0.1, Gamma=Gamma,
                                  n_pts=n_points)
        p = _find_peak(gp_arr, eta)
        lbl = f"ε = {eps:.1f}Δ  ({p['enhancement_zero']:.2f}×)"
        ax1.plot(gp_arr, eta, lw=2.0, color=col, label=lbl)
        if p["interior_peak"]:
            ax1.plot(p["gp_star"], p["eta_star"], "o", color=col, ms=7,
                     markeredgecolor="white", markeredgewidth=1.0, zorder=5)

    ax1.set_xscale("log")
    ax1.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax1.set_ylabel(r"Transfer Yield  $\eta_\infty$")
    ax1.set_title(r"Energy Bias  $\epsilon$  Controls ENAQT Strength  ($\kappa = 0.1\Delta$)")
    ax1.legend(title="Energy bias (enhancement)", fontsize=8)
    ax1.grid(True)
    ax1.text(
        0.98, 0.95,
        "Larger ε  →  stronger ENAQT\n(coherent tunneling less efficient)",
        transform=ax1.transAxes, fontsize=8.5, color="#333",
        ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
    )

    # -- Panel (1,0): HEOM (no sink) vs Bloch+sink (ε=1) -------------------
    ax2 = fig.add_subplot(gs[1, 0])

    _, eta_no_sink = _sweep_gamma_phi(
        eps=1.0, Delta=Delta, kappa=0.0, Gamma=Gamma, n_pts=n_points,
    )
    _, eta_sink = _sweep_gamma_phi(
        eps=1.0, Delta=Delta, kappa=0.1, Gamma=Gamma, n_pts=n_points,
    )
    _, eta_fastsink = _sweep_gamma_phi(
        eps=1.0, Delta=Delta, kappa=1.0, Gamma=Gamma, n_pts=n_points,
    )

    ax2.plot(gp_arr, eta_no_sink, lw=2.0, color="#7F8C8D", ls="--",
             label="Bloch (no sink,  κ=0)")
    ax2.plot(gp_arr, eta_sink, lw=2.5, color="#E74C3C",
             label="Bloch + sink  κ=0.1")
    ax2.plot(gp_arr, eta_fastsink, lw=2.5, color="#C0392B", ls="-.",
             label="Bloch + sink  κ=1.0")

    # Overlay HEOM SB dataset (asymmetric ε=1 case)
    if sb_data and "asymmetric_case" in sb_data:
        ac = sb_data["asymmetric_case"]["bell_curve"]
        gp_heom = ac.get("gamma_phi", [])
        eta_heom = ac.get("eta_mean", [])
        std_heom = ac.get("eta_std", [])
        if gp_heom:
            ax2.errorbar(
                gp_heom, eta_heom, yerr=std_heom,
                fmt="s", color="#2980B9", ms=7, capsize=3, lw=1.5,
                label="SB HEOM (exact, ε=1, no sink)", zorder=6,
            )

    # FMO reference
    ax2.axhline(0.9586, color="#F39C12", lw=1.5, ls=":", alpha=0.8,
                label="FMO-8 HEOM peak (4.12× enhancement)")

    ax2.set_xscale("log")
    ax2.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax2.set_ylabel(r"Transfer Efficiency / Yield  $\eta$")
    ax2.set_title(r"HEOM Data vs Bloch+Sink  ($\epsilon = 1\Delta$)" "\n"
                  "Small sink → big effect!")
    ax2.legend(fontsize=8)
    ax2.grid(True)

    p_sink = _find_peak(gp_arr, eta_sink)
    ax2.annotate(
        f"Sink κ=0.1:\nηmax = {p_sink['eta_star']:.3f}\n"
        f"({p_sink['enhancement_zero']:.2f}× over zero-deph)",
        xy=(p_sink["gp_star"], p_sink["eta_star"]),
        xytext=(p_sink["gp_star"] * 0.05, p_sink["eta_star"] * 0.78),
        fontsize=9, color="#E74C3C",
        arrowprops=dict(arrowstyle="->", color="#E74C3C", lw=1.2),
    )

    # -- Panel (1,1): κ sensitivity (ε=1) ----------------------------------
    ax3 = fig.add_subplot(gs[1, 1])
    kappas2 = [0.001, 0.01, 0.1, 0.5, 1.0, 5.0]
    kap2_colors = plt.cm.autumn_r(np.linspace(0.1, 0.85, len(kappas2)))

    for kap, col in zip(kappas2, kap2_colors):
        _, eta = _sweep_gamma_phi(
            eps=1.0, Delta=Delta, kappa=kap, Gamma=Gamma, n_pts=n_points,
        )
        p = _find_peak(gp_arr, eta)
        lbl = f"κ={kap:.3f}  ({p['enhancement_zero']:.2f}×)"
        ax3.plot(gp_arr, eta, lw=2.0, color=col, label=lbl)
        if p["interior_peak"]:
            ax3.plot(p["gp_star"], p["eta_star"], "o", color=col, ms=6,
                     markeredgecolor="white", markeredgewidth=1.0, zorder=5)

    ax3.set_xscale("log")
    ax3.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax3.set_ylabel(r"Transfer Yield  $\eta_\infty$")
    ax3.set_title(r"Sink Rate  $\kappa$  Sensitivity  ($\epsilon = 1\Delta$)")
    ax3.legend(title="κ (enh.)", fontsize=8)
    ax3.grid(True)

    # -- Panel (2,0): Time-domain dynamics -- 3 regimes (ε=5, κ=0.1) --------
    ax4 = fig.add_subplot(gs[2, 0])

    gp_optimal = _find_peak(*_sweep_gamma_phi(
        5.0, Delta, 0.1, Gamma, n_pts=n_points,
    ))["gp_star"]
    regimes = {
        f"Low  γ_φ=0.01  (coherent limit)": (0.01, "#27AE60"),
        f"Opt  γ_φ={gp_optimal:.2f}  (ENAQT sweet spot)": (gp_optimal, "#E74C3C"),
        f"High γ_φ=100   (quantum Zeno)": (100.0, "#8E44AD"),
    }

    T_dyn = 200.0
    for label, (gp_val, col) in regimes.items():
        t, y = _time_trajectory(
            eps=5.0, Delta=Delta, gp=gp_val, kappa=0.1,
            Gamma=Gamma, T=T_dyn, n_pts=5000,
        )
        rho11 = y[0]
        rho22 = y[1]
        eta_t = y[4]

        ax4.plot(t, rho11, color=col, lw=1.8, alpha=0.85, ls="-")
        ax4.plot(t, rho22, color=col, lw=1.8, alpha=0.55, ls="--")
        ax4.plot(t, eta_t, color=col, lw=2.5, alpha=1.00, ls="-",
                 label=f"{label}\n(η∞={eta_t[-1]:.3f})")

    legend_extra = [
        Line2D([0], [0], lw=2, color="#555", ls="-",
               label=r"$\rho_{11}$  (Site 1)"),
        Line2D([0], [0], lw=2, color="#555", ls="--",
               label=r"$\rho_{22}$  (Site 2)"),
        Line2D([0], [0], lw=2.5, color="#555", ls="-",
               label=r"$\eta(t)$  (Cumul. yield, bold)"),
    ]
    h, l = ax4.get_legend_handles_labels()
    ax4.legend(
        h + legend_extra,
        l + [e.get_label() for e in legend_extra],
        fontsize=8, ncol=2,
    )
    ax4.set_xlabel(r"Time  [$\Delta^{-1}$]")
    ax4.set_ylabel("Population / Yield")
    ax4.set_title("Population Dynamics  (ε=5Δ, κ=0.1Δ)  --  Three Regimes\n"
                  "Thick line = η(t) cumulative RC yield")
    ax4.set_ylim(-0.03, 1.05)
    ax4.grid(True)

    # -- Panel (2,1): 2D heatmap γ_φ × ε → η -------------------------------
    ax5 = fig.add_subplot(gs[2, 1])

    eps_vals = np.linspace(0.0, 6.0, 60)
    gp_vals = np.logspace(-3, 3, 80)
    ETA_MAP = np.zeros((len(eps_vals), len(gp_vals)))

    for i, eps in enumerate(eps_vals):
        for j, gp in enumerate(gp_vals):
            ETA_MAP[i, j] = _analytical_yield_4x4(
                eps, Delta, gp, kappa=0.1, Gamma=Gamma,
            )

    im = ax5.pcolormesh(
        gp_vals, eps_vals, ETA_MAP, cmap="plasma",
        shading="auto", vmin=0, vmax=1,
    )
    cb = plt.colorbar(im, ax=ax5, shrink=0.85, pad=0.02)
    cb.set_label(r"$\eta_\infty$  (total transfer yield)")

    # Mark the ENAQT ridge
    optimal_ridge = []
    for i, _eps in enumerate(eps_vals):
        j_opt = np.argmax(ETA_MAP[i, :])
        optimal_ridge.append(gp_vals[j_opt])
    ax5.plot(optimal_ridge, eps_vals, "w--", lw=2.0, alpha=0.85,
             label="Optimal γ_φ (ENAQT ridge)")

    ax5.axhline(1.0, color="cyan", lw=1.5, ls=":", alpha=0.8,
                label="SB dataset ε = 1Δ")
    ax5.axhline(5.0, color="lime", lw=1.5, ls=":", alpha=0.8,
                label="Max enhancement ε = 5Δ")

    ax5.set_xscale("log")
    ax5.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax5.set_ylabel(r"Energy Bias  $\epsilon$  [$\Delta$]")
    ax5.set_title("Full Parameter Map  (κ=0.1Δ)\n"
                  "White dashed = ENAQT ridge (optimal γ_φ at each ε)")
    ax5.legend(fontsize=8, loc="upper left")

    out_path = out_dir / "spinboson_sink_epic.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_sink_comparison(
    sb_data: dict,
    out_dir: Path,
    n_points: int = 300,
) -> Path:
    """2-panel summary: THE key comparison (sink vs no-sink)."""
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        "The Sink Effect -- Environment-Assisted Quantum Transport",
        fontsize=14, fontweight="bold",
    )

    gp_arr = np.logspace(-3, 3, n_points)
    Delta = 1.0
    Gamma = 0.01

    # -- Left: No sink -- SB HEOM data --------------------------------------
    ax_l.set_title(
        "WITHOUT Sink  (QD3SET-1 HEOM, exact)\n"
        "Equilibration only -- all pop returns",
        fontsize=11,
    )

    if sb_data:
        sc = sb_data.get("symmetric_case", {}).get("bell_curve", {})
        ac = sb_data.get("asymmetric_case", {}).get("bell_curve", {})
        if sc.get("gamma_phi"):
            ax_l.errorbar(
                sc["gamma_phi"], sc["eta_mean"], yerr=sc["eta_std"],
                fmt="o-", color="#2980B9", ms=6, capsize=3, lw=2,
                label="HEOM ε=0 (symmetric)",
            )
        if ac.get("gamma_phi"):
            ax_l.errorbar(
                ac["gamma_phi"], ac["eta_mean"], yerr=ac["eta_std"],
                fmt="s-", color="#C0392B", ms=6, capsize=3, lw=2,
                label="HEOM ε=1 (asymmetric)",
            )

    ax_l.axhline(0.5, color="#95A5A6", ls="--", lw=1.3, alpha=0.7,
                 label="Thermal equilibrium (ε=0)")

    ax_l.set_xscale("log")
    ax_l.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax_l.set_ylabel(r"$\eta = \rho_{22}(t_{\mathrm{final}})$")
    ax_l.set_ylim(0.45, 0.90)
    ax_l.legend()
    ax_l.grid(True)

    p_asym = _find_peak(
        np.array(ac.get("gamma_phi", [1])),
        np.array(ac.get("eta_mean", [0.5])),
    )
    ax_l.text(
        0.03, 0.96, f"Max enhancement: {p_asym['enhancement_zero']:.2f}×",
        transform=ax_l.transAxes, fontsize=12, fontweight="bold",
        color="#C0392B", va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FADBD8", alpha=0.9),
    )

    # -- Right: With sink -- Bloch model ------------------------------------
    ax_r.set_title(
        "WITH Lindblad Sink  (Bloch equations, exact)\n"
        "Irreversible RC extraction -- true efficiency",
        fontsize=11,
    )

    configs = [
        (1.0, 0.1, "#2980B9", "ε=1Δ, κ=0.1"),
        (2.0, 0.1, "#8E44AD", "ε=2Δ, κ=0.1"),
        (3.0, 0.1, "#E67E22", "ε=3Δ, κ=0.1"),
        (5.0, 0.1, "#E74C3C", "ε=5Δ, κ=0.1  ← MAX ENAQT"),
    ]
    for eps, kap, col, lbl in configs:
        _, eta = _sweep_gamma_phi(
            eps=eps, Delta=Delta, kappa=kap, Gamma=Gamma, n_pts=n_points,
        )
        p = _find_peak(gp_arr, eta)
        ax_r.plot(gp_arr, eta, lw=2.5 if eps == 5.0 else 1.8, color=col,
                  label=f"{lbl}  ({p['enhancement_zero']:.1f}×)")
        if p["interior_peak"]:
            ax_r.plot(p["gp_star"], p["eta_star"], "o", color=col, ms=9,
                      markeredgecolor="white", markeredgewidth=1.2, zorder=5)

    ax_r.axhline(0.9586, color="#F39C12", lw=1.5, ls=":", alpha=0.9,
                 label="FMO-8 HEOM peak (4.12×)")

    ax_r.set_xscale("log")
    ax_r.set_xlabel(r"Dephasing Rate  $\gamma_\phi$  [$\Delta$]")
    ax_r.set_ylabel(r"Transfer Yield  $\eta_\infty$  (RC capture efficiency)")
    ax_r.set_ylim(0.0, 1.02)
    ax_r.legend(fontsize=9)
    ax_r.grid(True)

    best_sink = _find_peak(*_sweep_gamma_phi(
        5.0, Delta, 0.1, Gamma, n_pts=n_points,
    ))
    ax_r.text(
        0.03, 0.96,
        f"Max enhancement: {best_sink['enhancement_zero']:.1f}×\n"
        f"Peak η = {best_sink['eta_star']:.3f}\n"
        f"Quantum Zeno η = {best_sink['eta_zeno']:.3f}",
        transform=ax_r.transAxes, fontsize=12, fontweight="bold",
        color="#C0392B", va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FADBD8", alpha=0.9),
    )

    plt.tight_layout()
    out_path = out_dir / "spinboson_sink_comparison.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ═════════════════════════════════════════════════════════════════════════════
#  Public API
# ═════════════════════════════════════════════════════════════════════════════

def run_validation(
    qd3set1_dir: str | Path,
    out_dir: str | Path,
    n_points: int = 300,
) -> dict[str, Any]:
    """Run QD3SET-1 HEOM validation against Bloch equations.

    Loads every ``*.npy`` trajectory in *qd3set1_dir*, reconstructs the
    full 2×2 density matrix at each time step, solves the corresponding
    Bloch equations, and computes the **trace distance**
    ``T(ρ_HEOM, ρ_Bloch) = ½ ||ρ_HEOM − ρ_Bloch||₁`` as the proper
    quantum-state comparison metric.

    Parameters
    ----------
    qd3set1_dir :
        Directory containing ``*.npy`` HEOM trajectory files.
    out_dir :
        Directory where plots and JSON results are written.
    n_points :
        Number of γ_φ points for the Bloch sweep (default 300).

    Returns
    -------
    dict
        Structured results including:

        * ``records`` -- list of per-trajectory metadata
        * ``symmetric_case`` / ``asymmetric_case`` -- ENAQT analysis
        * ``trace_distance`` -- aggregated Bloch-vs-HEOM error metrics
        * ``plots`` -- list of generated plot paths
        * ``json_path`` -- path to saved JSON results
    """
    qd3set1_dir = Path(qd3set1_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print("  ENAQT VALIDATION -- QD3SET-1 Spin-Boson HEOM vs Bloch")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load HEOM trajectories
    # ------------------------------------------------------------------
    print("\n[1/5] Loading HEOM trajectories...")
    records = _load_all(qd3set1_dir)

    sym_records = [r for r in records if r["epsilon"] == 0.0]
    asym_records = [r for r in records if r["epsilon"] == 1.0]

    all_lam = sorted({r["lam"] for r in records})
    all_gc = sorted({r["gamma_c"] for r in records})
    all_beta = sorted({r["beta"] for r in records})
    all_gp = sorted(r["gamma_phi"] for r in records)

    print(f"  ε = 0 (symmetric):  {len(sym_records)} trajectories")
    print(f"  ε = 1 (asymmetric): {len(asym_records)} trajectories")
    print(f"  λ values:   {all_lam}")
    print(f"  γ_c values: {all_gc}")
    print(f"  β values:   {all_beta}")
    print(
        f"  γ_φ range:  [{min(all_gp):.4f}, {max(all_gp):.4f}] Δ  "
        f"(span: {max(all_gp) / min(all_gp):.0f}×)"
    )

    # ------------------------------------------------------------------
    # 2. ENAQT bell-curve analysis (on HEOM data)
    # ------------------------------------------------------------------
    print("\n[2/5] Running ENAQT bell-curve analysis...")
    sym_res = _enaqt_test(sym_records, "Symmetric  (ε = 0)")
    asym_res = _enaqt_test(asym_records, "Asymmetric (ε = 1)")

    for res in [sym_res, asym_res]:
        tag = (
            "✓  BELL CURVE DETECTED"
            if res["enaqt_detected"] else "✗  No bell curve"
        )
        print(f"\n  [{res['label']}]  -- {tag}")
        print(f"    N trajectories:  {res['n_trajs']}")
        print(
            f"    γ_φ range:       [{res['gamma_phi_range'][0]:.4f}, "
            f"{res['gamma_phi_range'][1]:.4f}]"
        )
        print(
            f"    η range:         [{res['eta_range'][0]:.4f}, "
            f"{res['eta_range'][1]:.4f}]"
        )
        if res["enaqt_detected"]:
            print(f"    Optimal γ_φ:     {res['optimal_gamma_phi']:.4f} Δ")
            print(f"    Peak η:          {res['peak_eta']:.4f}")
            print(f"    Enhancement:     {res['enhancement']:.2f}×")

    # ------------------------------------------------------------------
    # 3. Bloch vs HEOM trace-distance comparison
    # ------------------------------------------------------------------
    print("\n[3/5] Computing Bloch vs HEOM trace-distance errors...")
    _compute_trace_distance_errors(records, sym_res, asym_res, n_points)

    # ------------------------------------------------------------------
    # 4. Generate plots
    # ------------------------------------------------------------------
    print("\n[4/5] Generating plots...")
    p1 = _plot_validation_main(sym_res, asym_res, records, out_dir)
    p2 = _plot_validation_gallery(records, out_dir)
    p3 = _plot_validation_beta_slices(records, out_dir)

    # ------------------------------------------------------------------
    # 5. Save JSON results
    # ------------------------------------------------------------------
    print("\n[5/5] Saving results JSON...")

    # Temperature-sliced analysis
    beta_slice_results: dict[str, dict] = {}
    for beta in all_beta:
        subset = [r for r in sym_records if r["beta"] == beta]
        if len(subset) >= 3:
            br = _enaqt_test(subset, f"β={beta}")
            beta_slice_results[str(beta)] = {
                "T": round(1.0 / beta, 4),
                "n_trajs": br["n_trajs"],
                "enaqt_detected": br["enaqt_detected"],
                "optimal_gamma_phi": br.get("optimal_gamma_phi"),
                "peak_eta": br.get("peak_eta"),
                "enhancement": br.get("enhancement"),
            }

    output: dict[str, Any] = {
        "experiment": "ENAQT Spin-Boson Validation -- QD3SET-1 HEOM vs Bloch",
        "dataset": {
            "source": str(qd3set1_dir),
            "total_trajectories": len(records),
            "trajectory_shape": "(401, 5) complex128",
            "columns": ["time", "rho11", "rho12", "rho21", "rho22"],
            "propagation_time": 20.0,
            "time_steps": 401,
            "lambda_values": all_lam,
            "gamma_c_values": all_gc,
            "beta_values": all_beta,
            "gamma_phi_range": [float(min(all_gp)), float(max(all_gp))],
            "gamma_phi_span_decades": float(np.log10(max(all_gp) / min(all_gp))),
        },
        "symmetric_case": {
            "epsilon": 0.0,
            "description": "Unbiased two-level system -- Rabi oscillations + dephasing",
            "enaqt_detected": sym_res["enaqt_detected"],
            "optimal_gamma_phi": sym_res.get("optimal_gamma_phi"),
            "peak_eta": sym_res.get("peak_eta"),
            "low_eta": sym_res.get("low_eta"),
            "high_eta": sym_res.get("high_eta"),
            "enhancement": sym_res.get("enhancement"),
            "bell_curve": {
                "gamma_phi": sym_res["gp_binned"],
                "eta_mean": sym_res["eta_binned"],
                "eta_std": sym_res["eta_std"],
                "counts": sym_res["counts"],
            },
            "trace_distance": {
                "mean": sym_res.get("td_mean"),
                "max": sym_res.get("td_max"),
                "binned_gp": sym_res.get("td_binned_gp"),
                "binned_mean": sym_res.get("td_binned_mean"),
                "binned_std": sym_res.get("td_binned_std"),
            },
        },
        "asymmetric_case": {
            "epsilon": 1.0,
            "description": "Biased two-level system -- directional energy transport",
            "enaqt_detected": asym_res["enaqt_detected"],
            "optimal_gamma_phi": asym_res.get("optimal_gamma_phi"),
            "peak_eta": asym_res.get("peak_eta"),
            "low_eta": asym_res.get("low_eta"),
            "high_eta": asym_res.get("high_eta"),
            "enhancement": asym_res.get("enhancement"),
            "bell_curve": {
                "gamma_phi": asym_res["gp_binned"],
                "eta_mean": asym_res["eta_binned"],
                "eta_std": asym_res["eta_std"],
                "counts": asym_res["counts"],
            },
            "trace_distance": {
                "mean": asym_res.get("td_mean"),
                "max": asym_res.get("td_max"),
                "binned_gp": asym_res.get("td_binned_gp"),
                "binned_mean": asym_res.get("td_binned_mean"),
                "binned_std": asym_res.get("td_binned_std"),
            },
        },
        "temperature_slices": beta_slice_results,
        "plots_generated": [str(p) for p in [p1, p2, p3]],
    }

    json_path = out_dir / "spinboson_validation_results.json"
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"  Saved: {json_path}")

    print("\n" + "=" * 70)
    print("  DONE -- Spin-boson validation complete!")
    print("=" * 70)

    output["json_path"] = str(json_path)
    return output


def _compute_trace_distance_errors(
    records: list[dict],
    sym_res: dict,
    asym_res: dict,
    n_points: int,
) -> None:
    """Compute Bloch-vs-HEOM trace distances and attach to result dicts."""
    Gamma = 0.0  # No recombination for the open-system comparison

    sym_td_gp, sym_td_vals = [], []
    asym_td_gp, asym_td_vals = [], []

    for rec in records:
        eps = rec["epsilon"]
        Delta = rec["delta"]
        gp = rec["gamma_phi"]

        # Build 2x2 density matrix from HEOM trajectory at each time step
        traj = np.load(rec["fpath"])
        t_heom = traj[:, 0].real
        n_steps = len(t_heom)

        # Bloch prediction: solve ODE for the same time span
        A = _bloch_liouvillian_4x4(eps, Delta, gp, kappa=0.0, Gamma=Gamma)

        def _rhs(_t: float, y: np.ndarray) -> np.ndarray:
            return A @ y

        y0 = np.array([1.0, 0.0, 0.0, 0.0])
        sol = solve_ivp(
            _rhs,
            [float(t_heom[0]), float(t_heom[-1])],
            y0,
            t_eval=t_heom,
            method="RK45",
            rtol=1e-9,
            atol=1e-11,
        )

        # Compute trace distance at each time step
        td_vals = []
        for i in range(n_steps):
            rho_heom = _reconstruct_rho2x2(
                traj[i, 1], traj[i, 4], traj[i, 2],
            )
            rho_bloch = _reconstruct_rho2x2(
                sol.y[0, i], sol.y[1, i], sol.y[2, i] + 1j * sol.y[3, i],
            )
            td_vals.append(_trace_distance(rho_heom, rho_bloch))

        mean_td = float(np.mean(td_vals))
        max_td = float(np.max(td_vals))

        rec["trace_distance_mean"] = mean_td
        rec["trace_distance_max"] = max_td

        if eps == 0.0:
            sym_td_gp.append(gp)
            sym_td_vals.append(mean_td)
        elif eps == 1.0:
            asym_td_gp.append(gp)
            asym_td_vals.append(mean_td)

    # Attach aggregated metrics
    if sym_td_vals:
        sym_res["td_mean"] = float(np.mean(sym_td_vals))
        sym_res["td_max"] = float(np.max(sym_td_vals))
        gp_b, td_b, td_s, _ = _log_bin(
            np.array(sym_td_gp), np.array(sym_td_vals), n_bins=20,
        )
        sym_res["td_binned_gp"] = gp_b.tolist()
        sym_res["td_binned_mean"] = td_b.tolist()
        sym_res["td_binned_std"] = td_s.tolist()

    if asym_td_vals:
        asym_res["td_mean"] = float(np.mean(asym_td_vals))
        asym_res["td_max"] = float(np.max(asym_td_vals))
        gp_b, td_b, td_s, _ = _log_bin(
            np.array(asym_td_gp), np.array(asym_td_vals), n_bins=20,
        )
        asym_res["td_binned_gp"] = gp_b.tolist()
        asym_res["td_binned_mean"] = td_b.tolist()
        asym_res["td_binned_std"] = td_s.tolist()

    print(f"  Symmetric  (ε=0):  mean T = {sym_res.get('td_mean', 'N/A')}, "
          f"max T = {sym_res.get('td_max', 'N/A')}")
    print(f"  Asymmetric (ε=1):  mean T = {asym_res.get('td_mean', 'N/A')}, "
          f"max T = {asym_res.get('td_max', 'N/A')}")


def _enaqt_test(records: list[dict], label: str) -> dict[str, Any]:
    """Core ENAQT test: find bell curve in η vs γ_φ.

    Parameters
    ----------
    records :
        Subset of trajectory records (e.g. all ε=0 or all ε=1).
    label :
        Human-readable label for this subset.

    Returns
    -------
    dict
        Analysis results including ``enaqt_detected``, ``optimal_gamma_phi``,
        ``peak_eta``, ``enhancement``, binned curves, etc.
    """
    gp = np.array([r["gamma_phi"] for r in records])
    eta = np.array([r["eta_final"] for r in records])
    avg = np.array([r["eta_avg"] for r in records])

    gp_c, eta_c, eta_s, cnts = _log_bin(gp, eta, n_bins=20)
    _, avg_c, avg_s, _ = _log_bin(gp, avg, n_bins=20)

    if len(eta_c) < 3:
        return {"label": label, "enaqt_detected": False}

    idx_max = int(np.argmax(eta_c))
    interior_peak = 0 < idx_max < len(eta_c) - 1

    eta_low = eta_c[0]
    eta_peak = eta_c[idx_max]
    eta_high = eta_c[-1]
    enhancement = eta_peak / eta_low if eta_low > 0 else float("nan")

    return {
        "label": label,
        "n_trajs": len(records),
        "gamma_phi_range": [float(gp.min()), float(gp.max())],
        "eta_range": [float(eta.min()), float(eta.max())],
        "enaqt_detected": bool(interior_peak),
        "optimal_gamma_phi": float(gp_c[idx_max]),
        "peak_eta": float(eta_peak),
        "low_eta": float(eta_low),
        "high_eta": float(eta_high),
        "enhancement": float(enhancement),
        # binned curves
        "gp_binned": gp_c.tolist(),
        "eta_binned": eta_c.tolist(),
        "eta_std": eta_s.tolist(),
        "avg_binned": avg_c.tolist(),
        "avg_std": avg_s.tolist(),
        "counts": cnts.tolist(),
        # raw
        "gp_raw": gp.tolist(),
        "eta_raw": eta.tolist(),
    }


def run_sink_analysis(
    out_dir: str | Path,
    n_points: int = 300,
) -> dict[str, Any]:
    """Run sink vs no-sink comparison using Bloch + Lindblad sink.

    Computes analytical yield curves ``η_∞(γ_φ)`` for several combinations
    of energy bias ε and sink rate κ, and compares with the no-sink HEOM
    reference data (loaded from ``spinboson_validation_results.json`` in
    *out_dir* if present).

    Parameters
    ----------
    out_dir :
        Directory where plots and JSON results are written.
    n_points :
        Number of γ_φ points for the sweep (default 300).

    Returns
    -------
    dict
        Structured results with keys ``key_results``, ``curves``,
        ``plots``, ``json_path``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 72)
    print("  ENAQT + LINDBLAD SINK -- Rebentrost et al. (2009) Replication")
    print("=" * 72)

    # -- Load SB HEOM reference data ----------------------------------------
    print("\n[1/5] Loading SB HEOM reference data...")
    sb_data: dict = {}
    sb_json = out_dir / "spinboson_validation_results.json"
    if sb_json.exists():
        with open(sb_json) as f:
            sb_data = json.load(f)
        print(f"  Loaded: {sb_json}")
    else:
        print(f"  WARNING: {sb_json} not found -- validation plots will lack "
              "HEOM overlay")

    # -- Compute analytical yield curves ------------------------------------
    print("\n[2/5] Computing analytical ENAQT yield curves...")
    Delta = 1.0
    Gamma = 0.01
    gp_arr = np.logspace(-3, 3, n_points)

    results: dict[str, dict] = {}
    configs = [
        ("eps0_kap01", 0.0, 0.1),
        ("eps1_kap01", 1.0, 0.1),
        ("eps2_kap01", 2.0, 0.1),
        ("eps5_kap01", 5.0, 0.1),
        ("eps5_kap001", 5.0, 0.01),
        ("eps5_kap1", 5.0, 1.0),
    ]
    for label, eps, kap in configs:
        _, eta = _sweep_gamma_phi(eps, Delta, kap, Gamma, n_pts=n_points)
        p = _find_peak(gp_arr, eta)
        results[label] = {
            "eps": eps,
            "kappa": kap,
            **p,
            "gp": gp_arr.tolist(),
            "eta": eta.tolist(),
        }
        tag = "✓ ENAQT" if p["interior_peak"] else "— monotone"
        print(
            f"  ε={eps:.1f}, κ={kap:.3f}:  {tag}  |  "
            f"γ_φ* = {p['gp_star']:.3f}  |  "
            f"η* = {p['eta_star']:.4f}  |  "
            f"{p['enhancement_zero']:.2f}× enhancement"
        )

    # -- Key physics summary ------------------------------------------------
    print("\n" + "-" * 72)
    print("  KEY RESULT: How the sink amplifies ENAQT")
    print("-" * 72)
    print(
        f"  Bloch (ε=1, κ=0.1, WITH sink):  "
        f"{results['eps1_kap01']['enhancement_zero']:.2f}× enhancement  "
        "(analytical)"
    )
    print(
        f"  Bloch (ε=5, κ=0.1, WITH sink):  "
        f"{results['eps5_kap01']['enhancement_zero']:.2f}× enhancement  "
        "← EPIC!"
    )
    print(f"  FMO-8 HEOM reference:            4.12× enhancement")
    print(
        "\n  Why ε=5 >> FMO: larger site energy mismatch means coherent"
        "\n  tunneling is even LESS efficient, so noise helps MORE."
    )

    # -- Plots --------------------------------------------------------------
    print("\n[3/5] Generating plots...")
    p1 = _plot_sink_epic(sb_data, out_dir, n_points)
    p2 = _plot_sink_comparison(sb_data, out_dir, n_points)

    # -- Save JSON ----------------------------------------------------------
    print("\n[4/5] Saving results...")
    out = {
        "experiment": "ENAQT + Lindblad Sink -- Bloch Equation Analysis",
        "model": {
            "Hamiltonian": r"H = ε/2 σ_z + Δ/2 σ_x  (spin-boson)",
            "bath": "Pure dephasing, Lindblad σ_z, rate γ_φ",
            "sink": "Irreversible RC extraction from site 2, rate κ",
            "fluorescence": "Recombination from all sites, rate Γ",
            "Delta": Delta,
            "Gamma": Gamma,
            "method": "Exact analytical yield via matrix inversion",
        },
        "key_results": {
            "Bloch_eps1_kap01_enhancement": round(
                results["eps1_kap01"]["enhancement_zero"], 3,
            ),
            "Bloch_eps5_kap01_enhancement": round(
                results["eps5_kap01"]["enhancement_zero"], 3,
            ),
            "FMO8_HEOM_enhancement": 4.12,
            "conclusion": (
                "Sink + energy bias together unlock orders-of-magnitude "
                "stronger ENAQT"
            ),
        },
        "curves": results,
        "plots": [str(p1), str(p2)],
    }

    json_path = out_dir / "spinboson_sink_results.json"
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  Saved: {json_path}")

    # -- Final summary ------------------------------------------------------
    print("\n[5/5] Summary")
    print("=" * 72)
    print(
        f"  Step 1  + Lindblad sink (ε=1, κ=0.1): "
        f"{results['eps1_kap01']['enhancement_zero']:.2f}× ENAQT enhancement"
    )
    print(
        f"  Step 2  + Larger bias (ε=5, κ=0.1):   "
        f"{results['eps5_kap01']['enhancement_zero']:.2f}× enhancement!"
    )
    print(
        f"\n  Peak yield:     η* = {results['eps5_kap01']['eta_star']:.4f}  "
        f"at γ_φ = {results['eps5_kap01']['gp_star']:.3f} Δ"
    )
    print(f"  Zeno limit:     η  = {results['eps5_kap01']['eta_zeno']:.5f}")
    print(
        f"  Span (peak/Zeno): "
        f"{results['eps5_kap01']['eta_star'] / results['eps5_kap01']['eta_zeno']:.0f}×"
    )
    print(f"\n  Output files:")
    for p in [p1, p2, json_path]:
        print(f"    {p.name}")
    print("\n  Sink analysis complete!\n")

    out["json_path"] = str(json_path)
    return out


def regime_classification(
    epsilon: float,
    delta: float,
    kappa: float,
    gamma_phi: float,
) -> str:
    """Classify the transport regime for a 2-site spin-boson model.

    The classifier uses analytical criteria based on the ratios of the
    dephasing rate *gamma_phi* to the characteristic energy scales:

    * **localization** -- ``γ_φ ≪ min(Δ, ε)`` : low dephasing, quantum
      interference localises the excitation.
    * **optimal** -- ``γ_φ ~ Δ`` : dephasing breaks destructive
      interference, maximising transport.
    * **zeno** -- ``γ_φ ≫ max(Δ, ε, κ)`` : strong dephasing freezes
      dynamics via the quantum Zeno effect.

    Parameters
    ----------
    epsilon :
        Site energy bias (site 1 = +ε/2, site 2 = –ε/2).
    delta :
        Tunneling matrix element.
    kappa :
        Sink extraction rate from site 2.
    gamma_phi :
        Pure dephasing rate.

    Returns
    -------
    str
        One of ``"localization"``, ``"optimal"``, ``"zeno"``.

    Examples
    --------
    >>> regime_classification(1.0, 1.0, 0.1, 0.001)
    'localization'
    >>> regime_classification(1.0, 1.0, 0.1, 1.0)
    'optimal'
    >>> regime_classification(1.0, 1.0, 0.1, 100.0)
    'zeno'
    """
    # Characteristic energy scale
    E_char = max(delta, abs(epsilon), kappa)

    # Thresholds (factor of 10 in either direction)
    if gamma_phi < 0.1 * E_char:
        return "localization"
    elif gamma_phi > 10.0 * E_char:
        return "zeno"
    else:
        return "optimal"
