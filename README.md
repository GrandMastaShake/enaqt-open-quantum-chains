# ENAQT: Environment-Assisted Quantum Transport in Open Quantum Chains

**Validation, Scaling Laws, and Disorder Universality**

A computational study demonstrating that environmental noise constructively enhances quantum
transport in open chains — and that structural disorder amplifies rather than destroys this effect.

---

## Overview

This repository contains four analysis scripts and a paper draft exploring ENAQT across
increasing system complexity:

| Script | What it does | Key result |
|--------|-------------|------------|
| `enaqt_sb_analysis.py` | Validates Lindblad framework against 1,000 exact HEOM trajectories (QD3SET-1) | 1.27× enhancement, machine-precision agreement |
| `enaqt_sb_sink.py` | Adds irreversible Lindblad sink (reaction center model) | 7.20× enhancement at ε=5Δ, κ=0.1 |
| `enaqt_nsite_chain.py` | N-site chain scaling, N=2–20, three topologies | Enhancement ~2.12N; γ_φ* ~N^−1.24; FMO-7 achieves 32.1× |
| `enaqt_disorder_ensemble.py` | 100-seed disorder ensemble + σ sweep | ENAQT universal (95–100%); disorder median 244× at N=15 |

Paper draft (submission-ready): `PAPER_ENAQT_DRAFT.md`

---

## Requirements

```bash
pip install numpy scipy matplotlib
```

Python 3.10+. No GPU needed — all computations run on CPU via `numpy.linalg.solve`.

The QD3SET-1 dataset (spin-boson HEOM subset) is required for `enaqt_sb_analysis.py`:

> Ullah et al. (2023), *Frontiers in Physics* **11**, 1223973.  
> Download: https://doi.org/10.25452/figshare.plus.c.6389553

Place the `.npy` trajectory files in `../SB/data/` relative to this directory.

---

## Running the Analysis

Run scripts in order — each builds on the previous:

```bash
# 1. HEOM validation (requires QD3SET-1 data)
python enaqt_sb_analysis.py

# 2. Lindblad sink analysis (no external data needed)
python enaqt_sb_sink.py

# 3. N-site chain scaling (no external data needed)
python enaqt_nsite_chain.py

# 4. Disorder ensemble (no external data needed, ~60s runtime)
python enaqt_disorder_ensemble.py
```

Each script writes PNG figures and a JSON results file to the working directory.

---

## Key Results

### Ordered energy funnel (bias = 5Δ, κ = 0.1Δ, Γ = 0.01Δ)

| N | Enhancement | Optimal γ_φ* [Δ] |
|---|-------------|-----------------|
| 2 | 2.8× | 5.11 |
| 5 | 14.3× | 1.28 |
| 7 | 22.8× | 0.84 |
| 10 | 32.4× | 0.55 |
| 15 | 37.9× | 0.37 |
| **FMO-7** | **32.1×** | **1.57** |

Scaling laws: enhancement ≈ 2.12N (linear, R²≈0.97); γ_φ* ≈ C·N^(−1.24)

### Disorder ensemble (σ = 2Δ, 100 seeds)

| N | ENAQT fraction | Median enhancement | Mean enhancement |
|---|---------------|-------------------|-----------------|
| 7 | 100% | 8.9× | 97.0× |
| 10 | 100% | 35.7× | 594.8× |
| 15 | 100% | **244.5×** | **6916×** |

Disorder amplifies ENAQT via Anderson localization. Mean grows as σ^5–6 with disorder strength.

---

## Outputs

```
enaqt_sb_analysis.png          — HEOM bell curves (ε=0 and ε=1)
enaqt_sb_dynamics_gallery.png  — 12 sample HEOM trajectories
enaqt_sb_beta_slices.png       — Temperature-sliced transport map
enaqt_sb_sink.png              — Lindblad yield vs ε sweep (with sink)
enaqt_sink_vs_nosink.png       — Sink vs no-sink comparison
enaqt_nsite_scaling.png        — Enhancement and γ_φ* scaling with N
enaqt_nsite_dynamics.png       — Population dynamics at three regimes
enaqt_disorder_ensemble.png    — 6-panel disorder ensemble summary
enaqt_disorder_paper_figure.png — 2-panel publication figure
```

---

## Citation

If you use this code or results, please cite:

```
Harper, A. (2026). Environment-Assisted Quantum Transport in Open Quantum Chains:
Validation, Scaling Laws, and Disorder Universality.
[Preprint / New Journal of Physics submission]
```

And the underlying dataset:

```
Ullah, A., Herrera Rodriguez, L. E., Dral, P. O., & Kananenka, A. A. (2023).
QD3SET-1: A quantum dissipative dynamics dataset.
Frontiers in Physics, 11, 1223973. https://doi.org/10.3389/fphy.2023.1223973
```

And the original ENAQT theory paper:

```
Rebentrost, P., Mohseni, M., Kassal, I., Lloyd, S., & Aspuru-Guzik, A. (2009).
Environment-assisted quantum transport.
New Journal of Physics, 11, 033003. https://doi.org/10.1088/1367-2630/11/3/033003
```
