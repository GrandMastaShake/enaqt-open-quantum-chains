"""
enaqt_2d_lattice.py
─────────────────────────────────────────────────────────────────────────────
ENAQT on 2D square lattices: does the physics survive beyond 1D chains?

We built a complete picture of 1D ENAQT (chains, scaling laws, disorder
universality, optimal step-function). The natural next question:

    Does ENAQT work in 2D? How does enhancement scale with grid size L?
    Does disorder universality survive? What does the bell curve look like?

Setup
-----
L×L square lattice. Sites indexed k = i*L + j  (row i, column j).
Nearest-neighbour coupling J along both axes (isotropic).
Source: site (0,0) — top-left corner. Excitation starts here.
Sink:   site (L-1,L-1) — bottom-right corner. RC extraction.

Three topologies compared:
  Funnel   : diagonal energy gradient, ε_{ij} = (2-i-j)/(2(L-1)) * ε_tot
             (high energy at source, zero at sink)
  Flat     : all sites degenerate
  Disorder : Gaussian random site energies, σ = 2Δ

For disorder universality: 50 seeds per L, L ∈ {2,3,4,5,6}

Outputs
-------
enaqt_2d_results.json        — all numerical results
enaqt_2d_lattice.png         — main figure
"""

import numpy as np
import matplotlib.pyplot as plt
import json
import time
from itertools import product

# ── Constants ─────────────────────────────────────────────────────────────────

J      = 1.0    # nearest-neighbour coupling (Δ ≡ J)
KAPPA  = 0.1    # sink rate
GAMMA  = 0.01   # fluorescence loss
EPS_TOT = 5.0   # total diagonal energy bias

# ── Lattice construction ───────────────────────────────────────────────────────

def site_index(i, j, L):
    return i * L + j

def build_hamiltonian_2d(eps_grid, L):
    """
    Build Hamiltonian for L×L lattice.
    eps_grid[i,j] = on-site energy at row i, col j.
    Coupling J between all nearest neighbours (up/down/left/right).
    """
    N = L * L
    H = np.diag(eps_grid.flatten().astype(complex))
    for i, j in product(range(L), range(L)):
        k = site_index(i, j, L)
        if j + 1 < L:
            m = site_index(i, j + 1, L)
            H[k, m] = J;  H[m, k] = J
        if i + 1 < L:
            m = site_index(i + 1, j, L)
            H[k, m] = J;  H[m, k] = J
    return H

def funnel_energies(L):
    """Diagonal gradient: high at (0,0), zero at (L-1,L-1)."""
    eps = np.zeros((L, L))
    for i, j in product(range(L), range(L)):
        frac = (i + j) / (2 * (L - 1)) if L > 1 else 0.0
        eps[i, j] = EPS_TOT * (1.0 - frac)
    eps[L-1, L-1] = 0.0   # sink pinned to reference
    return eps

def flat_energies(L):
    return np.zeros((L, L))

def random_energies(L, sigma, rng):
    eps = rng.normal(0, sigma, (L, L))
    eps[L-1, L-1] = 0.0
    return eps

# ── Liouvillian machinery ──────────────────────────────────────────────────────

def build_liouvillian_parts_2d(eps_grid, L):
    """
    L(γ_φ) = L_base + γ_φ · L_deph
    Source: site (0,0). Sink: site (L-1,L-1).
    """
    N    = L * L
    I    = np.eye(N)
    H    = build_hamiltonian_2d(eps_grid, L)

    # Coherent part
    L_base = -1j * (np.kron(I, H) - np.kron(H.T, I))

    # Sink on bottom-right corner
    sink_k = site_index(L-1, L-1, L)
    Psink  = np.zeros((N, N)); Psink[sink_k, sink_k] = 1.0
    L_base -= 0.5 * KAPPA * (np.kron(I, Psink) + np.kron(Psink, I))

    # Fluorescence
    L_base -= GAMMA * np.eye(N * N)

    # Dephasing superoperator
    L_deph = np.zeros((N * N, N * N), dtype=complex)
    for k in range(N):
        Pk = np.zeros((N, N)); Pk[k, k] = 1.0
        L_deph += np.kron(Pk, Pk) - 0.5*np.kron(I, Pk) - 0.5*np.kron(Pk, I)

    return L_base, L_deph, sink_k

def compute_yield_2d(L_base, L_deph, N_sites, sink_k, gamma_phi):
    """Analytical yield via Laplace transform."""
    Lop  = L_base + gamma_phi * L_deph
    rho0 = np.zeros(N_sites * N_sites); rho0[0] = 1.0  # excitation at (0,0)
    try:
        rho_ss = np.linalg.solve(Lop, -rho0)
        return max(0.0, KAPPA * rho_ss[sink_k * N_sites + sink_k].real)
    except np.linalg.LinAlgError:
        return 0.0

def bell_curve(eps_grid, L, n_gp=80):
    """Compute η(γ_φ) and return (gp_array, eta_array, enhancement, gp_star)."""
    N_sites = L * L
    L_base, L_deph, sink_k = build_liouvillian_parts_2d(eps_grid, L)

    gp_range = np.logspace(-3, 2, n_gp)
    etas     = np.array([compute_yield_2d(L_base, L_deph, N_sites, sink_k, gp)
                         for gp in gp_range])

    eta_zero = compute_yield_2d(L_base, L_deph, N_sites, sink_k, 1e-8)
    eta_zero = max(eta_zero, 1e-12)

    idx      = etas.argmax()
    enh      = etas[idx] / eta_zero
    return gp_range, etas, float(enh), float(gp_range[idx]), float(etas[idx]), float(eta_zero)

# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    grid_sizes = [2, 3, 4, 5, 6]
    N_SEEDS    = 50
    SIGMA      = 2.0

    results    = {}
    bell_curves = {}

    print("=" * 65)
    print("  ENAQT on 2D Square Lattices")
    print("=" * 65)

    t0 = time.time()

    for L in grid_sizes:
        N_sites = L * L
        print(f"\n  L={L}  ({N_sites} sites, {N_sites**2}x{N_sites**2} Liouvillian)")

        # ── Funnel ────────────────────────────────────────────────────────────
        eps_f   = funnel_energies(L)
        gp_arr, eta_arr, enh_f, gps_f, peak_f, zero_f = bell_curve(eps_f, L)
        print(f"    Funnel:  {enh_f:.2f}x   (gp*={gps_f:.3f}D,  peak_eta={peak_f:.3f})")
        bell_curves[f"funnel_{L}"] = (gp_arr.tolist(), eta_arr.tolist())

        # ── Flat ──────────────────────────────────────────────────────────────
        eps_flat = flat_energies(L)
        _, _, enh_flat, _, _, _ = bell_curve(eps_flat, L)
        print(f"    Flat:    {enh_flat:.2f}x")

        # ── Disorder ensemble ─────────────────────────────────────────────────
        enhs_dis = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(seed)
            eps_d = random_energies(L, SIGMA, rng)
            _, _, enh_d, _, _, _ = bell_curve(eps_d, L, n_gp=50)
            enhs_dis.append(enh_d)

        med_d  = float(np.median(enhs_dis))
        mean_d = float(np.mean(enhs_dis))
        std_d  = float(np.std(enhs_dis))
        frac   = float(np.mean([e > 1.01 for e in enhs_dis]))

        print(f"    Disorder ({N_SEEDS} seeds): median={med_d:.1f}x  "
              f"mean={mean_d:.1f}x  ENAQT={frac*100:.0f}%")

        results[L] = {
            "N_sites"            : N_sites,
            "funnel_enhancement" : enh_f,
            "funnel_gp_star"     : gps_f,
            "funnel_peak_eta"    : peak_f,
            "funnel_zero_eta"    : zero_f,
            "flat_enhancement"   : enh_flat,
            "disorder_median"    : med_d,
            "disorder_mean"      : mean_d,
            "disorder_std"       : std_d,
            "disorder_enaqt_frac": frac,
            "disorder_all"       : enhs_dis,
        }

    print(f"\n  Total runtime: {time.time()-t0:.1f}s")

    with open("enaqt_2d_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("  Saved -> enaqt_2d_results.json")

    # ── Figure ────────────────────────────────────────────────────────────────

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle("ENAQT on 2D Square Lattices  (source: corner (0,0)  ->  sink: corner (L-1,L-1))",
                 fontsize=13, fontweight="bold", y=0.98)

    col_funnel = "#DC2626"
    col_dis    = "#2563EB"
    col_flat   = "#9CA3AF"

    Ls      = sorted(results.keys())
    fun_e   = [results[L]["funnel_enhancement"]  for L in Ls]
    med_e   = [results[L]["disorder_median"]      for L in Ls]
    mean_e  = [results[L]["disorder_mean"]        for L in Ls]
    flat_e  = [results[L]["flat_enhancement"]     for L in Ls]
    gp_star = [results[L]["funnel_gp_star"]       for L in Ls]
    frac_e  = [results[L]["disorder_enaqt_frac"]*100 for L in Ls]
    Nsites  = [results[L]["N_sites"]              for L in Ls]

    # ── Row 1: bell curves for funnel, L=2,3,4,5,6 ───────────────────────────
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(Ls)))
    ax1 = fig.add_subplot(2, 3, 1)
    for i, L in enumerate(Ls):
        gp_arr = bell_curves[f"funnel_{L}"][0]
        et_arr = bell_curves[f"funnel_{L}"][1]
        ax1.semilogx(gp_arr, et_arr, color=colors[i], lw=2,
                     label=f"L={L} ({L*L} sites)")
    ax1.set_xlabel("Dephasing rate gp (D)", fontsize=11)
    ax1.set_ylabel("Transfer yield eta", fontsize=11)
    ax1.set_title("Bell curves: 2D funnel topologies", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.25)

    # ── Row 1: enhancement vs N_sites (log-log) ───────────────────────────────
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.loglog(Nsites, fun_e,  "o-",  color=col_funnel, lw=2, ms=8, label="2D funnel")
    ax2.loglog(Nsites, med_e,  "s--", color=col_dis,    lw=1.5, ms=7, label="2D disorder median")
    ax2.loglog(Nsites, mean_e, "^:",  color=col_dis,    lw=1.5, ms=7, alpha=0.6, label="2D disorder mean")

    # Power-law fit for funnel
    log_n = np.log(Nsites)
    log_e = np.log(fun_e)
    slope_2d, intercept_2d = np.polyfit(log_n, log_e, 1)
    n_fit = np.logspace(np.log10(min(Nsites)), np.log10(max(Nsites)), 100)
    ax2.loglog(n_fit, np.exp(intercept_2d) * n_fit**slope_2d, "k--", lw=1.5,
               label=f"~ N_sites^{slope_2d:.2f}")

    ax2.set_xlabel("Number of sites N = L^2 (log)", fontsize=11)
    ax2.set_ylabel("Enhancement (log)", fontsize=11)
    ax2.set_title("Enhancement scaling with system size", fontsize=11)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.25, which="both")

    # ── Row 1: compare 2D vs 1D scaling ──────────────────────────────────────
    ax3 = fig.add_subplot(2, 3, 3)
    # 1D funnel data from our existing results
    Ns_1d    = [2, 3, 5, 7, 10, 15]
    enhs_1d  = [2.8, 6.3, 14.3, 22.8, 32.4, 37.9]
    ax3.plot(Ns_1d,  enhs_1d, "o-",  color="#16A34A", lw=2, ms=8, label="1D chain (fixed 5D bias)")
    ax3.plot(Nsites, fun_e,   "s--", color=col_funnel, lw=2, ms=8, label="2D lattice (diag 5D bias)")
    ax3.set_xlabel("Chain length N  /  Sites L^2", fontsize=11)
    ax3.set_ylabel("Enhancement", fontsize=11)
    ax3.set_title("1D chain vs 2D lattice comparison", fontsize=11)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.25)

    # ── Row 2: gp* vs L ──────────────────────────────────────────────────────
    ax4 = fig.add_subplot(2, 3, 4)
    ax4.plot(Ls, gp_star, "D-", color=col_funnel, lw=2, ms=9)
    slope_gp, ic_gp = np.polyfit(np.log(Ls), np.log(gp_star), 1)
    L_fit = np.linspace(min(Ls), max(Ls), 100)
    ax4.plot(L_fit, np.exp(ic_gp)*L_fit**slope_gp, "k--", lw=1.5,
             label=f"gp* ~ L^{slope_gp:.2f}")
    ax4.set_xlabel("Grid size L", fontsize=11)
    ax4.set_ylabel("Optimal dephasing gp* (D)", fontsize=11)
    ax4.set_title("Optimal dephasing rate vs grid size", fontsize=11)
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.25)

    # ── Row 2: disorder universality fraction ─────────────────────────────────
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.bar(Ls, frac_e, color=col_dis, alpha=0.7, edgecolor="navy")
    ax5.axhline(95, color="tomato", lw=1.5, ls="--", label="1D threshold (95-100%)")
    ax5.set_xlabel("Grid size L  (L^2 sites)", fontsize=11)
    ax5.set_ylabel("ENAQT fraction (%)", fontsize=11)
    ax5.set_title("Disorder universality in 2D", fontsize=11)
    ax5.set_ylim(0, 105)
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.25, axis="y")

    # ── Row 2: disorder distributions ────────────────────────────────────────
    ax6 = fig.add_subplot(2, 3, 6)
    for i, L in enumerate(Ls):
        data = [x for x in results[L]["disorder_all"] if x > 0]
        log_data = np.log10(np.array(data) + 1)
        ax6.violinplot([log_data], positions=[L], widths=0.6,
                       showmedians=True, showextrema=True)
    ax6.set_xlabel("Grid size L", fontsize=11)
    ax6.set_ylabel("log10(enhancement + 1)", fontsize=11)
    ax6.set_title("Enhancement distributions (2D disorder, 50 seeds)", fontsize=11)
    ax6.grid(True, alpha=0.25)
    ax6.set_xticks(Ls)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("enaqt_2d_lattice.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Saved -> enaqt_2d_lattice.png")

    # ── Terminal summary ──────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print(f"  {'L':>3}  {'Sites':>6}  {'Funnel':>10}  {'Dis.Med':>10}  {'ENAQT%':>8}  {'gp*':>6}")
    print("  " + "-" * 55)
    for L in Ls:
        r = results[L]
        print(f"  {L:>3}  {r['N_sites']:>6}  "
              f"{r['funnel_enhancement']:>9.1f}x  "
              f"{r['disorder_median']:>9.1f}x  "
              f"{r['disorder_enaqt_frac']*100:>7.0f}%  "
              f"{r['funnel_gp_star']:>6.3f}")
    print("=" * 65)

    # ── Scaling law summary ───────────────────────────────────────────────────
    print(f"\n  2D funnel enhancement ~ N_sites^{slope_2d:.2f}  (log-log fit)")
    print(f"  2D optimal dephasing  ~ L^{slope_gp:.2f}           (power law fit)")
    print(f"\n  Recall 1D: enhancement ~ N^1.0  (linear),  gp* ~ N^-1.24")
