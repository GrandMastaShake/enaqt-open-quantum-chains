"""
enaqt_optimal_disorder.py
─────────────────────────────────────────────────────────────────────────────
Optimal site-energy disorder for maximum ENAQT enhancement.

We already know (enaqt_disorder_ensemble.py) that random disorder amplifies
ENAQT -- sometimes dramatically. This script asks the natural next question:

    What is the BEST possible disorder configuration, and does it
    look like anything physically meaningful?

Method
------
Outer:  scipy.optimize.differential_evolution on N-1 free site energies
        (sink site ε_N ≡ 0 as reference; bounds ±5Δ per site)
Inner:  30-point log sweep of γ_φ to find peak enhancement for each candidate
Compare: optimal vs ordered funnel vs random ensemble (50 seeds, σ=2Δ)

Runtime: ~2-4 min for N = [3, 5, 7, 10, 12, 15]

Outputs
-------
enaqt_optimal_disorder.png   — energy profiles + scaling comparison
enaqt_optimal_results.json   — all numerical results
"""

import numpy as np
import matplotlib.pyplot as plt
import json
import time
from scipy.optimize import differential_evolution

# ── Physical constants ────────────────────────────────────────────────────────

J      = 1.0    # nearest-neighbour coupling (sets energy scale Δ ≡ J)
KAPPA  = 0.1    # irreversible sink rate (reaction centre)
GAMMA  = 0.01   # fluorescence loss (homogeneous)

# ── Core Lindblad machinery ───────────────────────────────────────────────────

def build_liouvillian_parts(eps, N):
    """
    Decompose Liouvillian into L(γ_φ) = L_base + γ_φ · L_deph.

    L_base  : coherent evolution + sink + fluorescence  [N²×N² complex]
    L_deph  : pure-dephasing superoperator              [N²×N² complex]

    Building them once and combining cheaply is what lets us sweep γ_φ
    fast inside the optimizer inner loop.
    """
    I = np.eye(N)
    H = np.diag(eps) + np.diag(np.ones(N - 1) * J, 1) + np.diag(np.ones(N - 1) * J, -1)

    # Coherent part: -i[H, ·] in vectorised form
    L_base = -1j * (np.kron(I, H) - np.kron(H.T, I))

    # Irreversible sink on the last site (reaction centre)
    PN = np.zeros((N, N)); PN[N - 1, N - 1] = 1.0
    L_base -= 0.5 * KAPPA * (np.kron(I, PN) + np.kron(PN, I))

    # Fluorescence loss (uniform)
    L_base -= GAMMA * np.eye(N * N)

    # Dephasing superoperator (sum of site projectors)
    L_deph = np.zeros((N * N, N * N), dtype=complex)
    for j in range(N):
        Pj = np.zeros((N, N)); Pj[j, j] = 1.0
        L_deph += np.kron(Pj, Pj) - 0.5 * np.kron(I, Pj) - 0.5 * np.kron(Pj, I)

    return L_base, L_deph


def compute_yield(L_base, L_deph, N, gamma_phi):
    """
    Analytical transport yield at dephasing rate γ_φ via Laplace transform:
        η = κ · [−L⁻¹ ρ₀]_{sink}
    Single matrix solve per γ_φ point.
    """
    L    = L_base + gamma_phi * L_deph
    rho0 = np.zeros(N * N); rho0[0] = 1.0          # excitation starts at site 1
    try:
        rho_ss = np.linalg.solve(L, -rho0)
        return max(0.0, KAPPA * rho_ss[(N - 1) * N + (N - 1)].real)
    except np.linalg.LinAlgError:
        return 0.0


def peak_enhancement(eps, N, n_gp=30):
    """
    Enhancement ratio η_peak / η_zero for site energies eps.

    Returns (enhancement, γ_φ*, η_peak, η_zero)
    """
    L_base, L_deph = build_liouvillian_parts(eps, N)

    eta_zero = compute_yield(L_base, L_deph, N, 1e-8)
    eta_zero = max(eta_zero, 1e-12)               # guard against zero-noise zeros

    gp_range = np.logspace(-3, 2, n_gp)
    etas     = np.array([compute_yield(L_base, L_deph, N, gp) for gp in gp_range])

    idx      = etas.argmax()
    return etas[idx] / eta_zero, gp_range[idx], etas[idx], eta_zero


# ── Optimizer wrapper ─────────────────────────────────────────────────────────

def neg_enhancement(free_eps, N):
    """Objective: −enhancement (we maximise by minimising the negative)."""
    eps = np.append(free_eps, 0.0)                # sink site pinned at 0
    enh, _, _, _ = peak_enhancement(eps, N)
    return -enh


def optimize_disorder(N, sigma_bound=5.0, seed=42, maxiter=100, popsize=8):
    """
    Global optimisation of site energies via differential evolution.

    Returns (optimal_eps, optimal_enhancement, γ_φ*)
    """
    t0     = time.time()
    bounds = [(-sigma_bound, sigma_bound)] * (N - 1)

    result = differential_evolution(
        neg_enhancement,
        bounds,
        args      = (N,),
        seed      = seed,
        maxiter   = maxiter,
        popsize   = popsize,
        tol       = 1e-4,
        mutation  = (0.5, 1.5),
        recombination = 0.7,
        polish    = True,          # L-BFGS-B refinement after DE
        workers   = 1,
    )

    opt_eps          = np.append(result.x, 0.0)
    enh, gp_star, _, _ = peak_enhancement(opt_eps, N, n_gp=60)   # fine sweep at end
    elapsed          = time.time() - t0
    print(f"    enhancement = {enh:>8.1f}x   gp* = {gp_star:.3f}D   ({elapsed:.1f}s)", flush=True)
    return opt_eps, enh, gp_star


# ── Baselines ─────────────────────────────────────────────────────────────────

def funnel_enhancement(N):
    """Linear energy funnel: ε_i = (N−1−i)·Δ, ε_N = 0."""
    eps = np.array([(N - 1 - i) * 1.0 for i in range(N)])
    eps[-1] = 0.0
    return peak_enhancement(eps, N, n_gp=60)


def random_ensemble_stats(N, n_seeds=50, sigma=2.0):
    """Median and mean enhancement over random disorder realisations."""
    enhs = []
    for s in range(n_seeds):
        rng = np.random.default_rng(s)
        eps = rng.normal(0, sigma, N); eps[-1] = 0.0
        enh, _, _, _ = peak_enhancement(eps, N)
        enhs.append(enh)
    return float(np.median(enhs)), float(np.mean(enhs))


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    chain_sizes = [3, 5, 7, 10, 12, 15]

    results       = {}
    optimal_eps_all = {}

    print("=" * 65)
    print("  Optimal Disorder Search — ENAQT Enhancement Maximisation")
    print("=" * 65)

    t_total = time.time()

    for N in chain_sizes:
        print(f"\n  N = {N}  ({N-1} free site-energy parameters)")

        # ── Optimal disorder ──────────────────────────────────────────────────
        print(f"    Optimising ...", end=" ", flush=True)
        opt_eps, opt_enh, opt_gp = optimize_disorder(N)

        # ── Baselines ─────────────────────────────────────────────────────────
        f_enh, f_gp, _, _ = funnel_enhancement(N)
        print(f"    Funnel:  {f_enh:.1f}x")

        med, mean = random_ensemble_stats(N, n_seeds=50, sigma=2.0)
        print(f"    Random (s=2D, 50 seeds):  median={med:.1f}x   mean={mean:.1f}x")

        results[N] = {
            "optimal_enhancement" : float(opt_enh),
            "optimal_gp_star"     : float(opt_gp),
            "optimal_eps"         : opt_eps.tolist(),
            "funnel_enhancement"  : float(f_enh),
            "funnel_gp_star"      : float(f_gp),
            "random_median"       : float(med),
            "random_mean"         : float(mean),
            "gain_vs_funnel"      : float(opt_enh / f_enh),
            "gain_vs_random_med"  : float(opt_enh / max(med, 1.0)),
        }
        optimal_eps_all[N] = opt_eps.tolist()

    print(f"\nTotal runtime: {time.time() - t_total:.1f}s")

    with open("enaqt_optimal_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Results saved -> enaqt_optimal_results.json")

    # ── Figure ────────────────────────────────────────────────────────────────

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(
        "Optimal Site-Energy Disorder for Maximum ENAQT Enhancement",
        fontsize=14, fontweight="bold", y=0.98
    )

    # Colour palette
    col_opt    = "#2563EB"   # blue   — optimal
    col_funnel = "#DC2626"   # red    — funnel
    col_med    = "#16A34A"   # green  — random median
    col_mean   = "#D97706"   # amber  — random mean

    # ── Row 1: optimal energy profiles for N = 5, 7, 10 ─────────────────────
    profile_Ns = [5, 7, 10]
    for col_idx, N in enumerate(profile_Ns):
        ax  = fig.add_subplot(2, 3, col_idx + 1)
        eps = np.array(results[N]["optimal_eps"])
        fun = np.array([(N - 1 - i) * 1.0 for i in range(N)]); fun[-1] = 0.0
        sites = np.arange(1, N + 1)

        ax.axhline(0, color="silver", lw=0.8, ls="--")
        ax.plot(sites, fun, "s--", color=col_funnel, lw=1.5, ms=7,
                alpha=0.75, label="Ordered funnel")
        ax.plot(sites, eps, "o-",  color=col_opt,    lw=2,   ms=9,
                label="Optimal disorder")
        ax.scatter([N], [0], s=120, marker="*", color="gold",
                   zorder=5, label="Sink (site N)")

        ax.set_xlabel("Site index", fontsize=11)
        ax.set_ylabel("Site energy (Δ)", fontsize=11)
        ax.set_title(
            f"N = {N}   |   {results[N]['optimal_enhancement']:.0f}× optimal  "
            f"vs  {results[N]['funnel_enhancement']:.0f}× funnel",
            fontsize=10
        )
        ax.legend(fontsize=8, loc="upper right")
        ax.set_xticks(sites)
        ax.grid(True, alpha=0.25)

    # ── Row 2 left: enhancement scaling (log y) ───────────────────────────────
    ax = fig.add_subplot(2, 3, 4)
    Ns_plot  = sorted(results.keys())
    opt_enhs = [results[N]["optimal_enhancement"]  for N in Ns_plot]
    fun_enhs = [results[N]["funnel_enhancement"]    for N in Ns_plot]
    med_enhs = [results[N]["random_median"]         for N in Ns_plot]
    mea_enhs = [results[N]["random_mean"]           for N in Ns_plot]

    ax.plot(Ns_plot, opt_enhs,  "o-",  color=col_opt,    lw=2,   ms=8,  label="Optimal disorder")
    ax.plot(Ns_plot, fun_enhs,  "s--", color=col_funnel, lw=1.5, ms=6,  label="Ordered funnel")
    ax.plot(Ns_plot, med_enhs,  "^:",  color=col_med,    lw=1.5, ms=6,  label="Random median (σ=2Δ)")
    ax.plot(Ns_plot, mea_enhs,  "v:",  color=col_mean,   lw=1.5, ms=6,  alpha=0.7, label="Random mean (σ=2Δ)")
    ax.set_yscale("log")
    ax.set_xlabel("Chain length N", fontsize=11)
    ax.set_ylabel("Enhancement ratio (log)", fontsize=11)
    ax.set_title("Enhancement scaling comparison", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25, which="both")

    # ── Row 2 middle: optimal γ_φ* scaling ───────────────────────────────────
    ax = fig.add_subplot(2, 3, 5)
    gp_opt = [results[N]["optimal_gp_star"]  for N in Ns_plot]
    gp_fun = [results[N]["funnel_gp_star"]   for N in Ns_plot]

    ax.plot(Ns_plot, gp_opt, "D-",  color=col_opt,    lw=2,   ms=8,  label="Optimal disorder")
    ax.plot(Ns_plot, gp_fun, "s--", color=col_funnel, lw=1.5, ms=6,  label="Ordered funnel")

    # Power-law fit to optimal
    log_N  = np.log(Ns_plot)
    log_gp = np.log(gp_opt)
    slope, intercept = np.polyfit(log_N, log_gp, 1)
    N_fit  = np.linspace(min(Ns_plot), max(Ns_plot), 200)
    ax.plot(N_fit, np.exp(intercept) * N_fit ** slope, "k--", lw=1.5,
            label=f"gp* ~ N^{slope:.2f}")

    ax.set_xlabel("Chain length N", fontsize=11)
    ax.set_ylabel("γ_φ* (Δ)", fontsize=11)
    ax.set_title("Optimal dephasing rate vs N", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    # ── Row 2 right: how much better is optimal vs baselines? ─────────────────
    ax = fig.add_subplot(2, 3, 6)
    gain_f = [results[N]["gain_vs_funnel"]     for N in Ns_plot]
    gain_r = [results[N]["gain_vs_random_med"] for N in Ns_plot]

    ax.plot(Ns_plot, gain_f, "o-",  color=col_opt,  lw=2,   ms=8,  label="Optimal / funnel")
    ax.plot(Ns_plot, gain_r, "s--", color=col_med,  lw=1.5, ms=6,  label="Optimal / random median")
    ax.axhline(1, color="silver", lw=1, ls="-")
    ax.set_xlabel("Chain length N", fontsize=11)
    ax.set_ylabel("Gain factor", fontsize=11)
    ax.set_title("Advantage of optimal over baselines", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("enaqt_optimal_disorder.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Figure saved -> enaqt_optimal_disorder.png")

    # ── Terminal summary ──────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print(f"  {'N':>3}  {'Optimal':>10}  {'Funnel':>10}  {'Rand.Med':>10}  {'vsF':>10}")
    print("  " + "-" * 55)
    for N in Ns_plot:
        r = results[N]
        print(f"  {N:>3}  {r['optimal_enhancement']:>9.1f}×  "
              f"{r['funnel_enhancement']:>9.1f}×  "
              f"{r['random_median']:>9.1f}×  "
              f"{r['gain_vs_funnel']:>9.1f}×")
    print("=" * 65)
