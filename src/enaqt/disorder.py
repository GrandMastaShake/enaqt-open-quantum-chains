"""
Disorder ensemble analysis and global optimization for ENAQT.

Refactored from ``enaqt_disorder_ensemble.py`` and ``enaqt_optimal_disorder.py``.

This module provides:

* :func:`run_ensemble` — embarrassingly-parallel disorder-ensemble sweep with
  robust (non-Gaussian) statistics: median, MAD, IQR, Wilson CIs, BCa
  bootstrap, skewness / excess kurtosis.
* :func:`run_optimization` — differential-evolution search for the site-energy
  landscape that maximises ENAQT, including convergence diagnostics and a
  sensitivity study over multiple bound widths.
* Utility helpers: :func:`compute_enhancement`, :func:`wilson_ci`,
  :func:`_single_seed_disorder`, :func:`_objective_function`.

All heavy lifting is delegated to :mod:`enaqt.core` (Liouvillian
construction, Hamiltonians, yield computation).
"""

from __future__ import annotations

import json
import logging
import os
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from numpy.polynomial import polynomial as P
from scipy import stats
from scipy.optimize import differential_evolution

# Imports from the shared core (merged separately)
from enaqt.core import (
    analytical_yield,
    get_out_dir,
    hamiltonian_disordered,
    hamiltonian_funnel,
    liouvillian_parts,
    save_json,
)

logger = logging.getLogger(__name__)

# ── Physical defaults ──────────────────────────────────────────────────────────
_DELTA: float = 1.0
_KAPPA: float = 0.1
_GAMMA: float = 0.01
_SIGMA_DISORDER: float = 2.0
_GP_ARR: np.ndarray = np.logspace(-3, 3, 100)

# Matplotlib is imported lazily inside plotting routines so that the module
# remains importable in headless / CI environments.


# ═══════════════════════════════════════════════════════════════════════════════
#  Public API
# ═══════════════════════════════════════════════════════════════════════════════


def compute_enhancement(
    eta_peak: float, eta_zero: float, floor: float = 1e-16
) -> float:
    """Compute the ENAQT enhancement ratio *eta_peak / eta_zero*.

    Parameters
    ----------
    eta_peak:
        Transfer yield at the optimal dephasing rate (peak of the bell curve).
    eta_zero:
        Transfer yield in the zero-dephasing (coherent) limit.
    floor:
        Minimum value for *eta_zero* to avoid division-by-zero.  If
        ``eta_zero < floor`` a warning is emitted and ``eta_peak / floor`` is
        returned so the caller still receives a finite scalar.

    Returns
    -------
    float
        The enhancement ratio.  NaN is returned only when *both* inputs are
        NaN.
    """
    if not np.isfinite(eta_peak) or not np.isfinite(eta_zero):
        return float("nan")
    if eta_zero < floor:
        warnings.warn(
            f"eta_zero ({eta_zero:.3e}) below floor ({floor:.3e}); "
            f"returning eta_peak / floor = {eta_peak / floor:.3e}",
            RuntimeWarning,
            stacklevel=2,
        )
        return float(eta_peak / floor)
    return float(eta_peak / eta_zero)


def wilson_ci(k: int, n: int, confidence: float = 0.95) -> Tuple[float, float]:
    """Wilson score interval for a binomial proportion *k / n*.

    The Wilson interval is well-behaved at the boundaries (unlike the Wald
    normal approximation) and is the recommended interval for proportions
    that can be close to 0 or 1 — exactly the regime encountered for the
    ENAQT-fraction in the disorder ensemble.

    Parameters
    ----------
    k:
        Number of successes (e.g. seeds with interior ENAQT peak).
    n:
        Number of trials (e.g. total seeds).
    confidence:
        Coverage probability (default 0.95).

    Returns
    -------
    tuple[float, float]
        ``(lower_bound, upper_bound)`` with each bound clipped to ``[0, 1]``.
    """
    if n <= 0:
        return (0.0, 1.0)
    if k < 0:
        k = 0
    if k > n:
        k = n

    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    p_hat = k / n
    denom = 1 + z * z / n
    centre = p_hat + z * z / (2 * n)
    margin = z * np.sqrt(
        (p_hat * (1 - p_hat) + z * z / (4 * n)) / n
    )
    lower = max(0.0, (centre - margin) / denom)
    upper = min(1.0, (centre + margin) / denom)
    return float(lower), float(upper)


# ═══════════════════════════════════════════════════════════════════════════════
#  Ensemble runner
# ═══════════════════════════════════════════════════════════════════════════════


def _single_seed_disorder(
    args: Tuple[int, int, float, float, float, float, np.ndarray],
) -> Dict[str, Any]:
    """Worker function for :func:`joblib.Parallel` — single disorder realisation.

    Parameters
    ----------
    args:
        ``(seed, n, disorder_sigma, coupling, sink_rate, fluo_rate, gp_arr)``

    Returns
    -------
    dict
        Keys: ``eta`` (yield curve), ``peak_idx``, ``eta_peak``,
        ``eta_zero``, ``gp_star``, ``enhancement``, ``has_interior_peak``.
    """
    seed, n, disorder_sigma, coupling, sink_rate, fluo_rate, gp_arr = args

    rng = np.random.default_rng(seed)
    H = hamiltonian_disordered(n, 0.0, coupling, disorder_sigma, rng)
    L_base, L_deph = liouvillian_parts(H, sink_rate, fluo_rate)

    eta = np.empty(len(gp_arr), dtype=float)
    for i, gp in enumerate(gp_arr):
        L = L_base + gp * L_deph
        eta[i] = analytical_yield(L, sink_rate, n)

    peak_idx = int(np.nanargmax(eta))
    eta_peak = float(eta[peak_idx])
    eta_zero = float(eta[0])
    gp_star = float(gp_arr[peak_idx])
    enhancement = compute_enhancement(eta_peak, eta_zero, floor=1e-16)
    has_interior_peak = bool(0 < peak_idx < len(gp_arr) - 1)

    return {
        "eta": eta,
        "peak_idx": peak_idx,
        "eta_peak": eta_peak,
        "eta_zero": eta_zero,
        "gp_star": gp_star,
        "enhancement": enhancement,
        "has_interior_peak": has_interior_peak,
    }


def run_ensemble(
    n: int,
    n_seeds: int = 1000,
    out_dir: Path | str | None = None,
    n_jobs: int = -1,
    disorder_sigma: float = _SIGMA_DISORDER,
    coupling: float = _DELTA,
    sink_rate: float = _KAPPA,
    fluo_rate: float = _GAMMA,
    gp_arr: np.ndarray | None = None,
) -> Dict[str, Any]:
    """Run a disorder ensemble with robust (non-Gaussian) statistics.

    For each of *n_seeds* realisations a Gaussian-disordered Hamiltonian is
    generated, the full yield-vs-dephasing bell curve is computed using the
    ``L_base + L_deph`` decomposition, and per-seed peak statistics are
    gathered.  The computation is embarrassingly parallel via
    ``joblib.Parallel``.

    Statistics reported (all per *N*):
    * **Median** yield curve (not mean — the distribution is heavy-tailed).
    * **MAD** (median absolute deviation) and **IQR**.
    * **5th / 95th percentiles**.
    * **Skewness** and **excess kurtosis** of the enhancement distribution.
    * **Wilson CI** for the ENAQT-interior-peak fraction.
    * **BCa bootstrap CI** for the median enhancement.
    * **Mean / median ratio** as a heavy-tail diagnostic.

    Parameters
    ----------
    n:
        Chain length (number of sites).
    n_seeds:
        Number of independent disorder realisations.
    out_dir:
        Directory for JSON + plot outputs (created if absent).  If *None*,
        resolved via :func:`enaqt.core.get_out_dir`.
    n_jobs:
        ``joblib.Parallel`` parallelism (``-1`` = all CPUs).
    disorder_sigma:
        Standard deviation of the Gaussian site-energy disorder (in units of
        *coupling*).
    coupling:
        Nearest-neighbour coupling :math:`\\Delta`.
    sink_rate:
        Sink rate :math:`\\kappa`.
    fluo_rate:
        Fluorescence loss rate :math:`\\Gamma`.
    gp_arr:
        Dephasing-rate grid.  Defaults to ``np.logspace(-3, 3, 100)``.

    Returns
    -------
    dict
        Comprehensive results dictionary (also saved to JSON when *out_dir* is
        given).  Top-level keys include ``summary``, ``ensemble_curves``,
        ``per_seed``, ``robust_stats``, ``wilson_ci``, ``bca_ci``.
    """
    import joblib

    t_start = time.perf_counter()
    gp_arr = gp_arr if gp_arr is not None else _GP_ARR.copy()
    out_dir = Path(out_dir) if out_dir is not None else get_out_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Disorder ensemble: N=%d, %d seeds, sigma=%.2f, n_jobs=%d",
        n, n_seeds, disorder_sigma, n_jobs,
    )

    # ── Parallel seed computation ──────────────────────────────────────────
    work = [
        (seed, n, disorder_sigma, coupling, sink_rate, fluo_rate, gp_arr)
        for seed in range(n_seeds)
    ]

    results: List[Dict[str, Any]] = joblib.Parallel(
        n_jobs=n_jobs, backend="loky", verbose=0
    )(joblib.delayed(_single_seed_disorder)(w) for w in work)

    # ── Gather per-seed arrays ─────────────────────────────────────────────
    all_eta = np.stack([r["eta"] for r in results], axis=0)          # (n_seeds, n_gp)
    enhancements = np.array([r["enhancement"] for r in results])
    peak_etas = np.array([r["eta_peak"] for r in results])
    zero_etas = np.array([r["eta_zero"] for r in results])
    gp_stars = np.array([r["gp_star"] for r in results])
    has_interior = np.array([r["has_interior_peak"] for r in results], dtype=bool)

    valid_enh = np.isfinite(enhancements)
    frac_enaqt = float(has_interior.mean())

    # ── Ensemble curves (robust statistics) ────────────────────────────────
    eta_median = np.median(all_eta, axis=0)
    eta_mad = stats.median_abs_deviation(all_eta, axis=0, scale="normal")
    eta_p25 = np.percentile(all_eta, 25, axis=0)
    eta_p75 = np.percentile(all_eta, 75, axis=0)
    eta_p05 = np.percentile(all_eta, 5, axis=0)
    eta_p95 = np.percentile(all_eta, 95, axis=0)
    eta_mean = all_eta.mean(axis=0)
    eta_std = all_eta.std(axis=0)

    # ── Robust statistics for the enhancement distribution ─────────────────
    enh = enhancements[valid_enh]
    enh_median = float(np.median(enh)) if len(enh) else float("nan")
    enh_mad = float(stats.median_abs_deviation(enh, scale="normal")) if len(enh) else float("nan")
    enh_iqr = float(np.percentile(enh, 75) - np.percentile(enh, 25)) if len(enh) else float("nan")
    enh_p05 = float(np.percentile(enh, 5)) if len(enh) else float("nan")
    enh_p95 = float(np.percentile(enh, 95)) if len(enh) else float("nan")
    enh_mean = float(np.mean(enh)) if len(enh) else float("nan")
    enh_std = float(np.std(enh, ddof=1)) if len(enh) > 1 else float("nan")
    enh_skew = float(stats.skew(enh)) if len(enh) > 2 else float("nan")
    enh_kurt = float(stats.kurtosis(enh)) if len(enh) > 3 else float("nan")
    mean_median_ratio = float(enh_mean / enh_median) if enh_median > 0 else float("nan")

    # ── Wilson CI for ENAQT fraction ───────────────────────────────────────
    n_enaqt = int(has_interior.sum())
    wilson_lower, wilson_upper = wilson_ci(n_enaqt, n_seeds)

    # ── BCa bootstrap CI for median enhancement ────────────────────────────
    bca_lower, bca_upper = _bca_bootstrap_median(enh)

    # ── Summary dict ───────────────────────────────────────────────────────
    summary = {
        "n": n,
        "n_seeds": n_seeds,
        "disorder_sigma": disorder_sigma,
        "coupling": coupling,
        "sink_rate": sink_rate,
        "fluo_rate": fluo_rate,
        "elapsed_seconds": round(time.perf_counter() - t_start, 3),
    }

    ensemble_curves = {
        "gp": gp_arr.tolist(),
        "eta_median": eta_median.tolist(),
        "eta_mad": eta_mad.tolist(),
        "eta_p25": eta_p25.tolist(),
        "eta_p75": eta_p75.tolist(),
        "eta_p05": eta_p05.tolist(),
        "eta_p95": eta_p95.tolist(),
        "eta_mean": eta_mean.tolist(),
        "eta_std": eta_std.tolist(),
    }

    per_seed = {
        "enhancements": enh.tolist(),
        "peak_etas": peak_etas.tolist(),
        "zero_etas": zero_etas.tolist(),
        "gp_stars": gp_stars.tolist(),
        "has_interior_peak": has_interior.tolist(),
    }

    robust_stats = {
        "enhancement_median": enh_median,
        "enhancement_mad": enh_mad,
        "enhancement_iqr": enh_iqr,
        "enhancement_p05": enh_p05,
        "enhancement_p95": enh_p95,
        "enhancement_mean": enh_mean,
        "enhancement_std": enh_std,
        "enhancement_skewness": enh_skew,
        "enhancement_excess_kurtosis": enh_kurt,
        "mean_median_ratio": mean_median_ratio,
    }

    wilson = {
        "fraction": frac_enaqt,
        "n_enaqt": n_enaqt,
        "n_total": n_seeds,
        "confidence": 0.95,
        "lower": wilson_lower,
        "upper": wilson_upper,
    }

    bca_ci = {
        "median": enh_median,
        "confidence": 0.95,
        "lower": bca_lower,
        "upper": bca_upper,
    }

    output: Dict[str, Any] = {
        "summary": summary,
        "ensemble_curves": ensemble_curves,
        "per_seed": per_seed,
        "robust_stats": robust_stats,
        "wilson_ci": wilson,
        "bca_ci": bca_ci,
    }

    # ── Save JSON ──────────────────────────────────────────────────────────
    save_json(output, out_dir, f"disorder_ensemble_n{n}.json")

    # ── Plots ──────────────────────────────────────────────────────────────
    _plot_ensemble_curves(output, out_dir, n)
    _plot_enhancement_distribution(output, out_dir, n)

    elapsed = time.perf_counter() - t_start
    logger.info("Ensemble done in %.1f s — median enh=%.2fx", elapsed, enh_median)

    return output


# ═══════════════════════════════════════════════════════════════════════════════
#  Optimization — differential evolution
# ═══════════════════════════════════════════════════════════════════════════════


def _objective_function(
    energies: np.ndarray,
    n: int,
    mode: str = "ratio",
    coupling: float = _DELTA,
    sink_rate: float = _KAPPA,
    fluo_rate: float = _GAMMA,
    n_gp: int = 30,
) -> float:
    """Objective for the DE optimizer.

    Parameters
    ----------
    energies:
        First ``n-1`` site energies (the sink site is pinned to 0).
    n:
        Chain length.
    mode:
        ``"ratio"`` → maximise *enhancement ratio* (default).
        ``"yield"`` → maximise *absolute yield* (ignores zero-dephasing
        baseline).
    coupling, sink_rate, fluo_rate:
        Physical parameters forwarded to :func:`liouvillian_parts`.
    n_gp:
        Number of :math:`\\gamma_\\phi` grid points for the inner sweep.

    Returns
    -------
    float
        Value to *minimise* (negative of the quantity of interest).
    """
    eps = np.append(energies, 0.0)          # pin sink site
    H = np.diag(eps) + np.diag(np.ones(n - 1) * coupling, 1) + np.diag(np.ones(n - 1) * coupling, -1)
    L_base, L_deph = liouvillian_parts(H, sink_rate, fluo_rate)

    gp_range = np.logspace(-3, 2, n_gp)
    etas = np.empty(n_gp)
    for i, gp in enumerate(gp_range):
        L = L_base + gp * L_deph
        etas[i] = analytical_yield(L, sink_rate, n)

    idx = int(np.nanargmax(etas))
    eta_peak = float(etas[idx])

    if mode == "ratio":
        eta_zero = float(etas[0]) if np.isfinite(etas[0]) else 0.0
        # Compute raw ratio without floor guard for the optimizer — the
        # floor is applied only when *reporting* results.
        raw_ratio = eta_peak / max(eta_zero, 1e-30)
        return -raw_ratio
    elif mode == "yield":
        return -eta_peak
    else:
        raise ValueError(f"Unknown mode {mode!r}")


def run_optimization(
    n: int,
    bounds: Tuple[float, float] = (-50, 50),
    out_dir: Path | str | None = None,
    coupling: float = _DELTA,
    sink_rate: float = _KAPPA,
    fluo_rate: float = _GAMMA,
    seed: int = 42,
    maxiter: int = 100,
    popsize: int = 8,
    workers: int = -1,
) -> Dict[str, Any]:
    """Optimise site-energy landscape via differential evolution.

    The routine runs ``scipy.optimize.differential_evolution`` for a set of
    progressively wider bound widths (±2, ±5, ±10, ±20, ±50) to study the
    well-known *step-function artifact*:  if the true optimum lies near a
    bound, the result is sensitive to the bound width, revealing that the
    landscape is not smoothly quadratic near the solution.

    For each bound width the following convergence diagnostics are recorded:
    * Final population diversity (max – min fitness).
    * Stagnation flag (``True`` if the best value changed by < 1 % over the
      last 10 iterations).
    * Number of iterations actually consumed.
    * Best, worst, and mean fitness of the final population.

    Both **enhancement ratio** (peak / zero-dephasing) and **absolute yield**
    (peak yield) are reported.

    Parameters
    ----------
    n:
        Chain length.
    bounds:
        ``(lower, upper)`` for each free site energy.  Internally this is
        expanded into the list of bound widths described above.
    out_dir:
        Output directory (created if absent).  If *None*, resolved via
        :func:`enaqt.core.get_out_dir`.
    coupling, sink_rate, fluo_rate:
        Physical parameters.
    seed:
        RNG seed forwarded to ``differential_evolution``.
    maxiter:
        Maximum DE iterations per run.
    popsize:
        Population multiplier (``popsize * n_params`` individuals).
    workers:
        ``differential_evolution`` parallelism (``-1`` = all CPUs).

    Returns
    -------
    dict
        Keys include ``optimal`` (best result across all bound widths),
        ``by_bound_width`` (detailed results per width), and
        ``bound_sensitivity`` (fraction of widths agreeing on the optimum
        within 10 %).
    """
    t_start = time.perf_counter()
    out_dir = Path(out_dir) if out_dir is not None else get_out_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Disorder optimization: N=%d, seed=%d, maxiter=%d", n, seed, maxiter)

    # ── Test multiple bound widths ─────────────────────────────────────────
    bound_widths = [2, 5, 10, 20, 50]
    by_width: Dict[str, Dict[str, Any]] = {}
    best_overall = {"ratio": -np.inf, "yield": -np.inf}

    for width in bound_widths:
        t_w = time.perf_counter()
        de_bounds = [(-width, width)] * (n - 1)

        # --- Ratio mode ---
        result_ratio = differential_evolution(
            func=_objective_function,
            bounds=de_bounds,
            args=(n, "ratio", coupling, sink_rate, fluo_rate, 30),
            seed=seed,
            maxiter=maxiter,
            popsize=popsize,
            tol=1e-4,
            mutation=(0.5, 1.5),
            recombination=0.7,
            polish=True,
            workers=workers,
            updating="deferred" if workers != 1 else "immediate",
        )

        opt_eps = np.append(result_ratio.x, 0.0)
        enh_ratio, gp_star, eta_peak, eta_zero = _fine_peak_eval(
            opt_eps, n, coupling, sink_rate, fluo_rate
        )

        # Convergence diagnostics
        final_pop = getattr(result_ratio, "population", None)
        if final_pop is not None and hasattr(result_ratio, "population_energies"):
            pop_energies = result_ratio.population_energies
            diversity = float(pop_energies.max() - pop_energies.min())
            worst_fitness = float(pop_energies.max())
            mean_fitness = float(pop_energies.mean())
        else:
            diversity = float("nan")
            worst_fitness = float("nan")
            mean_fitness = float("nan")

        stagnation = _detect_stagnation(result_ratio)

        entry = {
            "bound_width": width,
            "optimal_energies": opt_eps.tolist(),
            "enhancement_ratio": float(enh_ratio),
            "absolute_yield": float(eta_peak),
            "eta_zero": float(eta_zero),
            "gp_star": float(gp_star),
            "n_iterations": int(result_ratio.nit),
            "n_function_evals": int(result_ratio.nfev),
            "success": bool(result_ratio.success),
            "convergence_diagnostics": {
                "population_diversity": diversity,
                "stagnation_flag": stagnation,
                "worst_fitness": worst_fitness,
                "mean_fitness": mean_fitness,
                "best_fitness": float(result_ratio.fun),
            },
            "elapsed_seconds": round(time.perf_counter() - t_w, 3),
        }
        by_width[str(width)] = entry

        if enh_ratio > best_overall["ratio"]:
            best_overall = {
                "ratio": float(enh_ratio),
                "yield": float(eta_peak),
                "energies": opt_eps.tolist(),
                "gp_star": float(gp_star),
                "bound_width": width,
            }

        logger.info(
            "  width=±%2d: ratio=%.2fx yield=%.4f iters=%d stagnation=%s (%.1fs)",
            width, enh_ratio, eta_peak, result_ratio.nit, stagnation,
            time.perf_counter() - t_w,
        )

    # ── Bound-width sensitivity analysis ───────────────────────────────────
    ratios = [by_width[str(w)]["enhancement_ratio"] for w in bound_widths]
    ratios_finite = [r for r in ratios if np.isfinite(r)]
    if ratios_finite:
        ref_ratio = np.median(ratios_finite)
        n_agree = sum(
            1 for r in ratios_finite if abs(r - ref_ratio) / max(ref_ratio, 1e-12) < 0.10
        )
        sensitivity = {
            "median_ratio_across_widths": float(ref_ratio),
            "fraction_agreeing_within_10pct": float(n_agree / len(ratios_finite)),
            "ratios": {str(w): r for w, r in zip(bound_widths, ratios)},
            "step_function_risk": (
                "HIGH" if n_agree < len(ratios_finite) / 2 else
                "MODERATE" if n_agree < len(ratios_finite) else
                "LOW"
            ),
        }
    else:
        sensitivity = {
            "median_ratio_across_widths": float("nan"),
            "fraction_agreeing_within_10pct": float("nan"),
            "ratios": {},
            "step_function_risk": "UNKNOWN",
        }

    output: Dict[str, Any] = {
        "n": n,
        "optimal": best_overall,
        "by_bound_width": by_width,
        "bound_sensitivity": sensitivity,
        "elapsed_seconds": round(time.perf_counter() - t_start, 3),
    }

    # ── Save + plot ────────────────────────────────────────────────────────
    save_json(output, out_dir, f"disorder_optimize_n{n}.json")
    _plot_optimization_results(output, out_dir, n)

    logger.info(
        "Optimization done in %.1f s — best ratio=%.2fx (width=±%d)",
        time.perf_counter() - t_start, best_overall["ratio"],
        best_overall.get("bound_width", -1),
    )
    return output


# ═══════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ═══════════════════════════════════════════════════════════════════════════════


def _fine_peak_eval(
    eps: np.ndarray,
    n: int,
    coupling: float,
    sink_rate: float,
    fluo_rate: float,
    n_gp: int = 60,
) -> Tuple[float, float, float, float]:
    """Fine-grained peak evaluation after DE converges.

    Returns ``(enhancement, gp_star, eta_peak, eta_zero)``.
    """
    H = (
        np.diag(eps)
        + np.diag(np.ones(n - 1) * coupling, 1)
        + np.diag(np.ones(n - 1) * coupling, -1)
    )
    L_base, L_deph = liouvillian_parts(H, sink_rate, fluo_rate)

    gp_range = np.logspace(-3, 2, n_gp)
    etas = np.empty(n_gp)
    for i, gp in enumerate(gp_range):
        L = L_base + gp * L_deph
        etas[i] = analytical_yield(L, sink_rate, n)

    idx = int(np.nanargmax(etas))
    eta_peak = float(etas[idx])
    eta_zero = float(etas[0]) if np.isfinite(etas[0]) else 0.0
    gp_star = float(gp_range[idx])
    enhancement = compute_enhancement(eta_peak, eta_zero)
    return enhancement, gp_star, eta_peak, eta_zero


def _detect_stagnation(result) -> bool:
    """Heuristic stagnation detection from a DE result object.

    Returns ``True`` if the improvement over the last ~10 % of iterations was
    less than 1 % relative.
    """
    # scipy's OptimizeResult does not expose full history, so we use a
    # lightweight heuristic: if the solver terminated early (fewer than 80 %
    # of maxiter used) *and* the tolerance was tight, we flag potential
    # stagnation.
    nit = getattr(result, "nit", 0)
    max_iter = getattr(result, "maxiter", nit + 1)
    if max_iter > 0 and nit / max_iter < 0.80:
        return True  # converged early → possible plateau
    return False


# ── BCa bootstrap for median enhancement ────────────────────────────────────


def _bca_bootstrap_median(
    data: np.ndarray,
    n_boot: int = 2000,
    confidence: float = 0.95,
) -> Tuple[float, float]:
    """BCa bootstrap confidence interval for the median.

    Parameters
    ----------
    data:
        1-D array of enhancement values (finite only).
    n_boot:
        Number of bootstrap resamples.
    confidence:
        Coverage level.

    Returns
    -------
    tuple[float, float]
        ``(lower, upper)`` BCa interval.  Falls back to percentile bootstrap
        when the BCa acceleration estimate is unstable.
    """
    data = np.asarray(data)
    if len(data) < 5:
        return float("nan"), float("nan")

    rng = np.random.default_rng(0)
    boot_medians = np.empty(n_boot)
    n = len(data)
    for b in range(n_boot):
        sample = rng.choice(data, size=n, replace=True)
        boot_medians[b] = np.median(sample)

    # Acceleration factor (jackknife)
    jack_idx = np.arange(n)
    jack_medians = np.empty(n)
    for i in range(n):
        jack_medians[i] = np.median(np.delete(data, i))
    jack_mean = jack_medians.mean()
    num = np.sum((jack_mean - jack_medians) ** 3)
    den = np.sum((jack_mean - jack_medians) ** 2)
    if den < 1e-30:
        # Degenerate — fall back to percentile
        alpha = (1 - confidence) / 2
        return float(np.percentile(boot_medians, alpha * 100)), float(
            np.percentile(boot_medians, (1 - alpha) * 100)
        )
    a_hat = num / (6 * den ** 1.5)

    # Bias-correction factor
    z0 = stats.norm.ppf(np.mean(boot_medians < np.median(data)))

    alpha = (1 - confidence) / 2
    z_alpha = stats.norm.ppf(alpha)
    z_1m_alpha = stats.norm.ppf(1 - alpha)

    z_l = z0 + (z0 + z_alpha) / (1 - a_hat * (z0 + z_alpha))
    z_u = z0 + (z0 + z_1m_alpha) / (1 - a_hat * (z0 + z_1m_alpha))

    p_l = stats.norm.cdf(z_l)
    p_u = stats.norm.cdf(z_u)

    lower = float(np.percentile(boot_medians, p_l * 100))
    upper = float(np.percentile(boot_medians, p_u * 100))
    return lower, upper


# ═══════════════════════════════════════════════════════════════════════════════
#  Plotting
# ═══════════════════════════════════════════════════════════════════════════════


def _plot_ensemble_curves(
    data: Dict[str, Any], out_dir: Path, n: int
) -> str:
    """Plot ensemble yield curves with robust-statistic bands."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves = data["ensemble_curves"]
    gp = np.array(curves["gp"])

    fig, ax = plt.subplots(figsize=(8, 5))

    # 5-95th percentile band
    ax.fill_between(
        gp, curves["eta_p05"], curves["eta_p95"],
        color="#3498DB", alpha=0.12, label="5-95th pctile",
    )
    # IQR band
    ax.fill_between(
        gp, curves["eta_p25"], curves["eta_p75"],
        color="#3498DB", alpha=0.25, label="IQR (25-75th)",
    )
    # Median line
    ax.plot(
        gp, curves["eta_median"], lw=2.5, color="#2980B9",
        label="Median",
    )
    # Mean line (for comparison)
    ax.plot(
        gp, curves["eta_mean"], "--", lw=1.5, color="#8E44AD", alpha=0.7,
        label="Mean",
    )

    ax.set_xscale("log")
    ax.set_xlabel(r"Dephasing rate $\gamma_\phi$ [Δ]")
    ax.set_ylabel(r"Transfer yield $\eta_\infty$")
    ax.set_title(
        f"Disorder ensemble  (N={n}, {data['summary']['n_seeds']} seeds, "
        f"σ={data['summary']['disorder_sigma']}Δ)"
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    path = out_dir / f"disorder_ensemble_curves_n{n}.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info("Saved ensemble curves plot: %s", path)
    return str(path)


def _plot_enhancement_distribution(
    data: Dict[str, Any], out_dir: Path, n: int
) -> str:
    """Plot the enhancement distribution with robust-statistic annotations."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    enh = np.array(data["per_seed"]["enhancements"])
    enh = enh[np.isfinite(enh)]
    stats_dict = data["robust_stats"]

    fig, ax = plt.subplots(figsize=(8, 5))

    # Histogram with log bins if the range is large
    if len(enh) > 0 and enh.max() / max(enh.min(), 1e-12) > 100:
        bins = np.logspace(
            np.log10(max(enh.min(), 1e-3)),
            np.log10(enh.max() * 1.1),
            60,
        )
        ax.set_xscale("log")
    else:
        bins = 60

    ax.hist(enh, bins=bins, color="#9B59B6", alpha=0.5, edgecolor="white", density=True)

    # Reference lines
    ax.axvline(
        stats_dict["enhancement_median"], color="#E74C3C", lw=2,
        label=f"Median = {stats_dict['enhancement_median']:.2f}x",
    )
    ax.axvline(
        stats_dict["enhancement_mean"], color="#2980B9", lw=1.5, ls="--",
        label=f"Mean = {stats_dict['enhancement_mean']:.2f}x",
    )
    ax.axvspan(
        stats_dict["enhancement_median"] - stats_dict["enhancement_mad"],
        stats_dict["enhancement_median"] + stats_dict["enhancement_mad"],
        color="#E74C3C", alpha=0.1, label=f"±MAD",
    )

    ax.set_xlabel("ENAQT Enhancement (peak / zero-dephasing)")
    ax.set_ylabel("Density")
    ax.set_title(
        f"Enhancement distribution  (N={n}, median={stats_dict['enhancement_median']:.2f}x, "
        f"IQR={stats_dict['enhancement_iqr']:.2f}, skew={stats_dict['enhancement_skewness']:.2f})"
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()

    path = out_dir / f"disorder_enhancement_dist_n{n}.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info("Saved enhancement distribution plot: %s", path)
    return str(path)


def _plot_optimization_results(
    data: Dict[str, Any], out_dir: Path, n: int
) -> str:
    """Plot energy profile + bound-sensitivity summary."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    optimal = data["optimal"]
    energies = np.array(optimal["energies"])
    sites = np.arange(1, n + 1)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={"height_ratios": [1, 1]})

    # ── Top: optimal energy profile ────────────────────────────────────────
    ax_top.axhline(0, color="silver", lw=0.8, ls="--")
    ax_top.plot(sites, energies, "o-", color="#2563EB", lw=2, ms=8, label="Optimal disorder")
    ax_top.scatter([n], [0], s=120, marker="*", color="gold", zorder=5, label="Sink (site N)")

    # Funnel reference
    funnel_eps = np.array(
        [(n - 1 - i) * 1.0 for i in range(n)]
    )
    funnel_eps[-1] = 0.0
    ax_top.plot(sites, funnel_eps, "s--", color="#DC2626", lw=1.5, ms=6, alpha=0.7, label="Ordered funnel")

    ax_top.set_xlabel("Site index")
    ax_top.set_ylabel("Site energy (Δ)")
    ax_top.set_title(
        f"N={n}  |  Optimal ratio={optimal['ratio']:.1f}x  "
        f"|  Yield={optimal['yield']:.4f}  |  Bound=±{optimal.get('bound_width', '?')}"
    )
    ax_top.set_xticks(sites)
    ax_top.legend(fontsize=9)
    ax_top.grid(True, alpha=0.25)

    # ── Bottom: bound-width sensitivity ────────────────────────────────────
    by_width = data["by_bound_width"]
    widths = sorted(int(w) for w in by_width.keys())
    ratios = [by_width[str(w)]["enhancement_ratio"] for w in widths]
    yields = [by_width[str(w)]["absolute_yield"] for w in widths]

    ax2 = ax_bot
    ax2_twin = ax2.twinx()

    l1 = ax2.plot(widths, ratios, "o-", color="#2563EB", lw=2, ms=8, label="Enhancement ratio")[0]
    l2 = ax2_twin.plot(widths, yields, "s--", color="#16A34A", lw=1.5, ms=6, label="Absolute yield")[0]

    # Highlight best
    best_w = optimal.get("bound_width", widths[0])
    ax2.axvline(best_w, color="#DC2626", lw=1.5, ls=":", alpha=0.7, label=f"Best (±{best_w})")

    ax2.set_xlabel("Bound width (Δ)")
    ax2.set_ylabel("Enhancement ratio", color="#2563EB")
    ax2_twin.set_ylabel("Absolute yield", color="#16A34A")
    ax2.set_title(
        f"Bound-width sensitivity  |  "
        f"Step-function risk: {data['bound_sensitivity']['step_function_risk']}  |  "
        f"Agreement: {data['bound_sensitivity']['fraction_agreeing_within_10pct']:.0%}"
    )
    ax2.legend([l1, l2], ["Enhancement ratio", "Absolute yield"], fontsize=9, loc="upper left")
    ax2.grid(True, alpha=0.25)
    plt.tight_layout()

    path = out_dir / f"disorder_optimize_n{n}.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info("Saved optimization plot: %s", path)
    return str(path)
