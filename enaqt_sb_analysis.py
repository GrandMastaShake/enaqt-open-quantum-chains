#!/usr/bin/env python3
"""
ENAQT Analysis — QD3SET-1 Spin-Boson Model
==========================================
Environment-Assisted Quantum Transport: testing the Rebentrost et al. (2009)
prediction on exact HEOM spin-boson trajectories from the SB (1) dataset.

Physics:
    Spin-boson Hamiltonian: H = ε/2 σ_z + Δ/2 σ_x
    Bath: Drude-Lorentz spectral density J(ω) = 2λγω/(ω²+γ²)
    Effective dephasing rate: γ_φ = 2λ / (β · γ_c)

    ENAQT prediction (Rebentrost 2009):
        η vs γ_φ should show a non-monotonic bell curve:
          - Low γ_φ  → quantum interference traps excitation (low efficiency)
          - Optimal γ_φ → noise breaks interference → peak transport
          - High γ_φ → quantum Zeno freezes dynamics (low efficiency)

Data:
    SB (1)/SB/data/*.npy  — 1000 trajectories, shape (401, 5) complex128
    Columns: [time, ρ₁₁, ρ₁₂, ρ₂₁, ρ₂₂]
    Filenames: 2_epsilon-ε_Delta-Δ_lambda-λ_gamma-γ_beta-β.npy

Author: Kimi Agent / Ember Professional
"""

import sys, os, glob, re, json
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR = r"C:\Users\alexa\Desktop\Death_Star\Ember\Professional\tauNOW\SB (1)\SB\data"
OUT_DIR  = r"C:\Users\alexa\Desktop\Death_Star\Ember\Professional\tauNOW\Kimi_Agent_ENAQT"

# ── Styling ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'legend.fontsize': 9,
    'figure.dpi': 100,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.7,
})

COLORS = {
    'sym':  '#2980B9',   # blue  — symmetric ε=0
    'asym': '#C0392B',   # red   — asymmetric ε=1
    'low':  '#27AE60',   # green — localization regime
    'mid':  '#F39C12',   # orange — optimal ENAQT
    'high': '#8E44AD',   # purple — quantum Zeno
    'neutral': '#7F8C8D',
}


# ═══════════════════════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def parse_filename(fpath: str) -> dict:
    """Extract physical parameters from SB trajectory filename."""
    fname = os.path.basename(fpath)
    def grab(key):
        m = re.search(rf'{key}-([\d.]+)', fname)
        return float(m.group(1)) if m else None
    return {
        'file':    fname,
        'epsilon': grab('epsilon'),
        'delta':   grab('Delta'),
        'lam':     grab('lambda'),
        'gamma_c': grab('gamma'),
        'beta':    float(re.search(r'beta-([\d.]+?)\.npy', fname).group(1)),
    }


def load_all(data_dir: str) -> list:
    """Load all .npy trajectories and compute derived quantities."""
    files = glob.glob(os.path.join(data_dir, '*.npy'))
    if not files:
        raise FileNotFoundError(f"No .npy files found in {data_dir}")

    records = []
    for fpath in files:
        p = parse_filename(fpath)
        traj = np.load(fpath)        # (401, 5) complex128

        t     = traj[:, 0].real
        rho11 = traj[:, 1].real
        rho22 = traj[:, 4].real
        coh   = traj[:, 2]           # ρ₁₂ — off-diagonal coherence

        # Dephasing rate (Drude-Lorentz high-T Markov): γ_φ = 2λ/(β·γ_c)
        gamma_phi = 2.0 * p['lam'] / (p['beta'] * p['gamma_c'])

        # Temperature in units of Δ (ℏ=k_B=Δ=1)
        T = 1.0 / p['beta']

        # Final site-2 population — primary ENAQT efficiency metric
        eta_final = float(rho22[-1])

        # Time-averaged site-2 population — captures relaxation speed
        eta_avg = float(np.trapz(rho22, t) / (t[-1] - t[0]))

        # Coherence lifetime — time when |ρ₁₂| first drops to 1/e of initial
        coh_mag = np.abs(coh)
        coh_init = coh_mag.max()
        coh_thresh = coh_init / np.e
        coh_decay_idx = np.argmax(coh_mag < coh_thresh) if np.any(coh_mag < coh_thresh) else -1
        coh_lifetime = float(t[coh_decay_idx]) if coh_decay_idx > 0 else float(t[-1])

        records.append({
            **p,
            'T':            T,
            'gamma_phi':    gamma_phi,
            'eta_final':    eta_final,
            'eta_avg':      eta_avg,
            'coh_lifetime': coh_lifetime,
            'fpath':        fpath,
        })

    print(f"  Loaded {len(records)} trajectories from {data_dir}")
    return records


# ═══════════════════════════════════════════════════════════════════════════════
#  ENAQT ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def log_bin(x, y, n_bins=25):
    """Log-bin (x, y) pairs → (bin_centers, means, stds)."""
    x = np.array(x)
    y = np.array(y)
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
    return (np.array(centers), np.array(means),
            np.array(stds), np.array(counts))


def enaqt_test(records: list, label: str) -> dict:
    """Core ENAQT test: find bell curve in η vs γ_φ."""
    gp  = np.array([r['gamma_phi']  for r in records])
    eta = np.array([r['eta_final']  for r in records])
    avg = np.array([r['eta_avg']    for r in records])

    gp_c, eta_c, eta_s, cnts = log_bin(gp, eta, n_bins=20)
    _,    avg_c, avg_s, _    = log_bin(gp, avg, n_bins=20)

    # Bell-curve test: peak is in the interior (not at endpoints)
    if len(eta_c) < 3:
        return {'label': label, 'enaqt_detected': False}

    idx_max = int(np.argmax(eta_c))
    interior_peak = 0 < idx_max < len(eta_c) - 1

    # Enhancement vs low-dephasing baseline
    eta_low  = eta_c[0]
    eta_peak = eta_c[idx_max]
    eta_high = eta_c[-1]
    enhancement = eta_peak / eta_low if eta_low > 0 else float('nan')

    # Regime boundaries (half-peak width)
    half_peak = (eta_peak + eta_low) / 2
    left_idx  = np.searchsorted(eta_c[:idx_max+1], half_peak)
    right_arr = eta_c[idx_max:]
    right_idx = idx_max + np.searchsorted(-right_arr, -half_peak)

    return {
        'label':              label,
        'n_trajs':            len(records),
        'gamma_phi_range':    [float(gp.min()), float(gp.max())],
        'eta_range':          [float(eta.min()), float(eta.max())],
        'enaqt_detected':     bool(interior_peak),
        'optimal_gamma_phi':  float(gp_c[idx_max]),
        'peak_eta':           float(eta_peak),
        'low_eta':            float(eta_low),
        'high_eta':           float(eta_high),
        'enhancement':        float(enhancement),
        # binned curves
        'gp_binned':          gp_c.tolist(),
        'eta_binned':         eta_c.tolist(),
        'eta_std':            eta_s.tolist(),
        'avg_binned':         avg_c.tolist(),
        'avg_std':            avg_s.tolist(),
        'counts':             cnts.tolist(),
        # raw
        'gp_raw':             gp.tolist(),
        'eta_raw':            eta.tolist(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════════

def plot_main(sym, asym, all_records, out_dir):
    """4-panel ENAQT figure."""
    fig = plt.figure(figsize=(16, 12))
    gs  = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.35)
    fig.suptitle(
        'ENAQT Analysis — QD3SET-1 Spin-Boson HEOM  |  1000 Exact Trajectories',
        fontsize=14, fontweight='bold', y=0.98,
    )

    # ── Panel 1: Bell curves ───────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])

    for res, key, marker in [(sym, 'sym', 'o'), (asym, 'asym', 's')]:
        gp_b  = np.array(res['gp_binned'])
        eta_b = np.array(res['eta_binned'])
        eta_s = np.array(res['eta_std'])
        col   = COLORS[key]

        ax1.scatter(res['gp_raw'], res['eta_raw'],
                    c=col, alpha=0.10, s=8, zorder=1)
        ax1.fill_between(gp_b, eta_b - eta_s, eta_b + eta_s,
                         color=col, alpha=0.15, zorder=2)
        ax1.plot(gp_b, eta_b, f'{marker}-', color=col, lw=2.2, ms=7,
                 label=res['label'], zorder=3)

        if res['enaqt_detected']:
            opt = res['optimal_gamma_phi']
            pk  = res['peak_eta']
            ax1.axvline(opt, color=col, lw=1.2, ls='--', alpha=0.6)
            ax1.annotate(f'γ*={opt:.2f}\nη*={pk:.3f}',
                         xy=(opt, pk), xytext=(opt * 2.5, pk * 0.97),
                         fontsize=8, color=col,
                         arrowprops=dict(arrowstyle='->', color=col, lw=1))

    ax1.set_xscale('log')
    ax1.set_xlabel('Dephasing Rate  γ_φ = 2λ/(βγ_c)  [Δ]')
    ax1.set_ylabel('Transfer Efficiency  η = ρ₂₂(t_final)')
    ax1.set_title('ENAQT Bell Curve  (Final Population at Sink)')
    ax1.legend()
    ax1.grid(True)

    detected = [r for r in [sym, asym] if r['enaqt_detected']]
    banner = '✓ ENAQT DETECTED' if detected else '✗ Bell curve not found'
    banner_col = '#27AE60' if detected else '#E74C3C'
    ax1.text(0.97, 0.05, banner, transform=ax1.transAxes,
             ha='right', va='bottom', fontsize=10, fontweight='bold',
             color=banner_col,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor=banner_col, alpha=0.9))

    # ── Panel 2: Time-averaged efficiency ─────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])

    for res, key, marker in [(sym, 'sym', 'o'), (asym, 'asym', 's')]:
        gp_b  = np.array(res['gp_binned'])
        avg_b = np.array(res['avg_binned'])
        avg_s = np.array(res['avg_std'])
        col   = COLORS[key]

        ax2.fill_between(gp_b, avg_b - avg_s, avg_b + avg_s,
                         color=col, alpha=0.15)
        ax2.plot(gp_b, avg_b, f'{marker}-', color=col, lw=2.2, ms=7,
                 label=res['label'])

    ax2.set_xscale('log')
    ax2.set_xlabel('Dephasing Rate  γ_φ  [Δ]')
    ax2.set_ylabel('Time-Averaged ρ₂₂')
    ax2.set_title('Integrated Transport Efficiency\n(∫ρ₂₂ dt / T_total)')
    ax2.legend()
    ax2.grid(True)

    # ── Panel 3: Population dynamics at three regimes ─────────────────────────
    ax3 = fig.add_subplot(gs[1, 0])

    sym_records = [r for r in all_records if r['epsilon'] == 0.0]
    sym_sorted  = sorted(sym_records, key=lambda r: r['gamma_phi'])
    n           = len(sym_sorted)

    regime_picks = {
        'Localization  (low γ_φ)':  sym_sorted[max(0, n // 12)],
        'Optimal ENAQT  (mid γ_φ)': sym_sorted[n // 2],
        'Quantum Zeno  (high γ_φ)': sym_sorted[min(n-1, 11*n//12)],
    }
    regime_cols = [COLORS['low'], COLORS['mid'], COLORS['high']]

    for (rlabel, rec), col in zip(regime_picks.items(), regime_cols):
        traj  = np.load(rec['fpath'])
        t     = traj[:, 0].real
        rho22 = traj[:, 4].real
        gp    = rec['gamma_phi']
        ax3.plot(t, rho22, color=col, lw=2.0,
                 label=f'{rlabel}\n(γ_φ={gp:.3f}, λ={rec["lam"]}, β={rec["beta"]})')

    ax3.axhline(0.5, color=COLORS['neutral'], ls='--', lw=1.2, alpha=0.7,
                label='Thermal equilibrium (ε=0)')
    ax3.set_xlabel('Time  [Δ⁻¹]')
    ax3.set_ylabel('ρ₂₂(t)  — Site 2 Population')
    ax3.set_title('Quantum Dynamics — Three Dephasing Regimes  (ε = 0)')
    ax3.legend(fontsize=8)
    ax3.grid(True)

    # ── Panel 4: 2D parameter space — η as color ──────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])

    gp_sym  = [r['gamma_phi'] for r in sym_records]
    lam_sym = [r['lam']       for r in sym_records]
    eta_sym = [r['eta_final'] for r in sym_records]

    sc = ax4.scatter(gp_sym, lam_sym, c=eta_sym,
                     cmap='plasma', alpha=0.75, s=25,
                     vmin=0, vmax=1)
    cb = plt.colorbar(sc, ax=ax4, shrink=0.85, pad=0.02)
    cb.set_label('η = ρ₂₂(t_final)')

    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('Dephasing Rate  γ_φ  [Δ]')
    ax4.set_ylabel('Reorganization Energy  λ  [Δ]')
    ax4.set_title('Parameter Space Map — Transport Efficiency  (ε = 0)')
    ax4.grid(True)

    out_path = os.path.join(out_dir, 'enaqt_sb_analysis.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {out_path}")
    plt.close()
    return out_path


def plot_dynamics_gallery(all_records, out_dir):
    """Gallery: ε=0 vs ε=1 dynamics side-by-side at matched γ_φ values."""
    sym  = sorted([r for r in all_records if r['epsilon'] == 0.0],
                  key=lambda r: r['gamma_phi'])
    asym = sorted([r for r in all_records if r['epsilon'] == 1.0],
                  key=lambda r: r['gamma_phi'])
    n_sym  = len(sym)
    n_asym = len(asym)

    picks_sym  = [sym[n_sym // 10], sym[n_sym // 2], sym[9 * n_sym // 10]]
    picks_asym = [asym[n_asym // 10], asym[n_asym // 2], asym[9 * n_asym // 10]]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=False)
    fig.suptitle('Population Dynamics — Symmetric vs Asymmetric SB  |  Three Dephasing Regimes',
                 fontsize=13, fontweight='bold')

    regime_labels = ['Localization (Low γ_φ)', 'Optimal ENAQT (Mid γ_φ)', 'Quantum Zeno (High γ_φ)']
    regime_cols   = [COLORS['low'], COLORS['mid'], COLORS['high']]

    for col_idx, (rs, ra, rlabel, rcol) in enumerate(
            zip(picks_sym, picks_asym, regime_labels, regime_cols)):

        for row_idx, (rec, eps_label) in enumerate([(rs, 'ε = 0  (Symmetric)'),
                                                     (ra, 'ε = 1  (Asymmetric)')]):
            ax = axes[row_idx, col_idx]
            traj  = np.load(rec['fpath'])
            t     = traj[:, 0].real
            rho11 = traj[:, 1].real
            rho22 = traj[:, 4].real
            coh   = np.abs(traj[:, 2])

            ax.plot(t, rho11, lw=1.8, color='#2980B9', label='ρ₁₁ (Site 1)', alpha=0.9)
            ax.plot(t, rho22, lw=1.8, color=rcol,      label='ρ₂₂ (Site 2)', alpha=0.9)
            ax.plot(t, coh,   lw=1.2, color='#95A5A6', label='|ρ₁₂| (coherence)', ls='--', alpha=0.7)

            ax.set_ylim(-0.05, 1.05)
            ax.set_xlabel('Time [Δ⁻¹]')
            ax.set_ylabel('Population')
            ax.set_title(
                f'{eps_label}\n'
                f'{rlabel}\n'
                f'γ_φ={rec["gamma_phi"]:.3f}, λ={rec["lam"]}, β={rec["beta"]}',
                fontsize=9,
            )
            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=8)
            ax.grid(True)

    plt.tight_layout()
    out_path = os.path.join(out_dir, 'enaqt_sb_dynamics_gallery.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {out_path}")
    plt.close()
    return out_path


def plot_beta_slices(all_records, out_dir):
    """Bell curves sliced by temperature β — shows ENAQT at each temperature."""
    sym = [r for r in all_records if r['epsilon'] == 0.0]
    betas = sorted(set(r['beta'] for r in sym))

    fig, axes = plt.subplots(1, len(betas), figsize=(4 * len(betas), 5), sharey=True)
    if len(betas) == 1:
        axes = [axes]
    fig.suptitle('ENAQT Bell Curves by Temperature  (ε = 0,  Symmetric SB)',
                 fontsize=13, fontweight='bold')

    cmap = plt.cm.coolwarm
    norm = plt.Normalize(vmin=min(betas), vmax=max(betas))

    for ax, beta in zip(axes, betas):
        subset = [r for r in sym if r['beta'] == beta]
        gp  = np.array([r['gamma_phi'] for r in subset])
        eta = np.array([r['eta_final'] for r in subset])

        if len(gp) < 3:
            continue

        col = cmap(norm(beta))
        gp_b, eta_b, eta_s, _ = log_bin(gp, eta, n_bins=15)

        ax.scatter(gp, eta, c=[col]*len(gp), alpha=0.25, s=12)
        ax.fill_between(gp_b, eta_b - eta_s, eta_b + eta_s,
                        color=col, alpha=0.20)
        ax.plot(gp_b, eta_b, 'o-', color=col, lw=2.0, ms=6)

        idx_max = int(np.argmax(eta_b))
        interior = 0 < idx_max < len(eta_b) - 1
        status = '✓ ENAQT' if interior else '— monotone'
        T = 1.0 / beta

        ax.set_xscale('log')
        ax.set_xlabel('γ_φ  [Δ]')
        ax.set_title(f'β = {beta}  (T = {T:.1f} Δ)\n{status}', fontsize=10)
        ax.grid(True)

    axes[0].set_ylabel('η = ρ₂₂(t_final)')
    plt.tight_layout()
    out_path = os.path.join(out_dir, 'enaqt_sb_beta_slices.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {out_path}")
    plt.close()
    return out_path


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("\n" + "═" * 70)
    print("  ENAQT ANALYSIS — QD3SET-1 Spin-Boson HEOM Trajectories")
    print("═" * 70)

    # ── Load ──────────────────────────────────────────────────────────────────
    print("\n[1/5] Loading trajectories...")
    records = load_all(DATA_DIR)

    sym_records  = [r for r in records if r['epsilon'] == 0.0]
    asym_records = [r for r in records if r['epsilon'] == 1.0]

    # Parameter survey
    all_lam  = sorted(set(r['lam']     for r in records))
    all_gc   = sorted(set(r['gamma_c'] for r in records))
    all_beta = sorted(set(r['beta']    for r in records))
    all_gp   = sorted(r['gamma_phi'] for r in records)

    print(f"  ε = 0 (symmetric):  {len(sym_records)} trajectories")
    print(f"  ε = 1 (asymmetric): {len(asym_records)} trajectories")
    print(f"  λ values:   {all_lam}")
    print(f"  γ_c values: {all_gc}")
    print(f"  β values:   {all_beta}")
    print(f"  γ_φ range:  [{min(all_gp):.4f}, {max(all_gp):.4f}]  (span: {max(all_gp)/min(all_gp):.0f}×)")

    # ── ENAQT analysis ────────────────────────────────────────────────────────
    print("\n[2/5] Running ENAQT bell-curve analysis...")
    sym_res  = enaqt_test(sym_records,  'Symmetric  (ε = 0)')
    asym_res = enaqt_test(asym_records, 'Asymmetric (ε = 1)')

    for res in [sym_res, asym_res]:
        tag = '✓  BELL CURVE DETECTED' if res['enaqt_detected'] else '✗  No bell curve'
        print(f"\n  [{res['label']}]  — {tag}")
        print(f"    N trajectories:  {res['n_trajs']}")
        print(f"    γ_φ range:       [{res['gamma_phi_range'][0]:.4f}, {res['gamma_phi_range'][1]:.4f}]")
        print(f"    η range:         [{res['eta_range'][0]:.4f}, {res['eta_range'][1]:.4f}]")
        if res['enaqt_detected']:
            enh = res['enhancement']
            print(f"    Optimal γ_φ:     {res['optimal_gamma_phi']:.4f} Δ")
            print(f"    Peak η:          {res['peak_eta']:.4f}")
            print(f"    Low-deph η:      {res['low_eta']:.4f}")
            print(f"    High-deph η:     {res['high_eta']:.4f}")
            print(f"    Enhancement:     {enh:.2f}× over localization baseline")

    # ── Plots ─────────────────────────────────────────────────────────────────
    print("\n[3/5] Generating plots...")
    p1 = plot_main(sym_res, asym_res, records, OUT_DIR)
    p2 = plot_dynamics_gallery(records, OUT_DIR)
    p3 = plot_beta_slices(records, OUT_DIR)

    # ── Save JSON ─────────────────────────────────────────────────────────────
    print("\n[4/5] Saving results JSON...")

    # Build temperature-sliced analysis
    beta_slice_results = {}
    for beta in all_beta:
        subset = [r for r in sym_records if r['beta'] == beta]
        if len(subset) >= 3:
            br = enaqt_test(subset, f'β={beta}')
            beta_slice_results[str(beta)] = {
                'T': round(1.0 / beta, 4),
                'n_trajs': br['n_trajs'],
                'enaqt_detected': br['enaqt_detected'],
                'optimal_gamma_phi': br.get('optimal_gamma_phi'),
                'peak_eta': br.get('peak_eta'),
                'enhancement': br.get('enhancement'),
            }

    output = {
        'experiment': 'ENAQT Spin-Boson Analysis — QD3SET-1 HEOM',
        'date': '2026-04-29',
        'dataset': {
            'source': 'SB (1)/SB/data/',
            'total_trajectories': len(records),
            'trajectory_shape': '(401, 5) complex128',
            'columns': ['time', 'rho11', 'rho12', 'rho21', 'rho22'],
            'propagation_time': 20.0,
            'time_steps': 401,
            'lambda_values': all_lam,
            'gamma_c_values': all_gc,
            'beta_values': all_beta,
            'gamma_phi_range': [float(min(all_gp)), float(max(all_gp))],
            'gamma_phi_span_decades': float(np.log10(max(all_gp) / min(all_gp))),
        },
        'symmetric_case': {
            'epsilon': 0.0,
            'description': 'Unbiased two-level system — Rabi oscillations + dephasing',
            'enaqt_detected': sym_res['enaqt_detected'],
            'optimal_gamma_phi': sym_res.get('optimal_gamma_phi'),
            'peak_eta': sym_res.get('peak_eta'),
            'low_eta': sym_res.get('low_eta'),
            'high_eta': sym_res.get('high_eta'),
            'enhancement': sym_res.get('enhancement'),
            'bell_curve': {
                'gamma_phi': sym_res['gp_binned'],
                'eta_mean': sym_res['eta_binned'],
                'eta_std': sym_res['eta_std'],
                'counts': sym_res['counts'],
            },
        },
        'asymmetric_case': {
            'epsilon': 1.0,
            'description': 'Biased two-level system — directional energy transport',
            'enaqt_detected': asym_res['enaqt_detected'],
            'optimal_gamma_phi': asym_res.get('optimal_gamma_phi'),
            'peak_eta': asym_res.get('peak_eta'),
            'low_eta': asym_res.get('low_eta'),
            'high_eta': asym_res.get('high_eta'),
            'enhancement': asym_res.get('enhancement'),
            'bell_curve': {
                'gamma_phi': asym_res['gp_binned'],
                'eta_mean': asym_res['eta_binned'],
                'eta_std': asym_res['eta_std'],
                'counts': asym_res['counts'],
            },
        },
        'temperature_slices': beta_slice_results,
        'plots_generated': [p1, p2, p3],
    }

    json_path = os.path.join(OUT_DIR, 'enaqt_sb_results.json')
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"  Saved: {json_path}")

    # ── Final summary ─────────────────────────────────────────────────────────
    print("\n[5/5] Summary")
    print("═" * 70)
    print(f"  Trajectories analyzed:   {len(records)}")
    print(f"  Dephasing range:         {min(all_gp):.4f} → {max(all_gp):.4f} Δ  "
          f"({np.log10(max(all_gp)/min(all_gp)):.1f} decades)")

    sym_tag  = '✓ ENAQT DETECTED' if sym_res['enaqt_detected']  else '✗ not detected'
    asym_tag = '✓ ENAQT DETECTED' if asym_res['enaqt_detected'] else '✗ not detected'
    print(f"\n  Symmetric  (ε=0):  {sym_tag}")
    if sym_res['enaqt_detected']:
        print(f"    Optimal γ_φ = {sym_res['optimal_gamma_phi']:.3f} Δ  |  "
              f"η* = {sym_res['peak_eta']:.4f}  |  "
              f"{sym_res['enhancement']:.2f}× enhancement")

    print(f"\n  Asymmetric (ε=1):  {asym_tag}")
    if asym_res['enaqt_detected']:
        print(f"    Optimal γ_φ = {asym_res['optimal_gamma_phi']:.3f} Δ  |  "
              f"η* = {asym_res['peak_eta']:.4f}  |  "
              f"{asym_res['enhancement']:.2f}× enhancement")

    print(f"\n  Output files:")
    for p in [p1, p2, p3, json_path]:
        print(f"    {os.path.basename(p)}")

    print("\n  DONE — Epic ENAQT test complete!\n")
    return output


if __name__ == '__main__':
    main()
