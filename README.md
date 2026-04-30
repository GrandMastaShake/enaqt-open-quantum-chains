# Environment-Assisted Quantum Transport in Open Quantum Chains

*Noise helps. Disorder helps more. And it works almost universally.*

---

## The idea

In quantum systems, environmental noise is usually the enemy — it destroys
coherence and scrambles information. But in the right conditions, noise does
something counterintuitive: it **helps** excitation get from A to B.

This is Environment-Assisted Quantum Transport (ENAQT). It was proposed in 2009
to explain how photosynthesis achieves near-perfect efficiency despite operating
in a warm, wet, noisy cell. Too little noise → quantum interference traps the
excitation. Too much → the system freezes (quantum Zeno effect). At just the
right level, noise opens transport channels that pure quantum mechanics closes.

This repository documents a systematic computational study of that effect.

---

## What we found

**1. The framework is exact.**
We validated an analytical Lindblad model against 1,000 numerically exact
simulations from the QD3SET-1 database. Agreement is at machine precision
(error < 10⁻¹⁵). The math is right.

**2. The sink is what makes it real.**
Without an irreversible "reaction center" that permanently captures arriving
excitation, ENAQT is a subtle 1.27× effect. Add the sink and it jumps to
**7.20×** — the same physics, but now you're measuring actual transfer
efficiency instead of just thermalization rate.

**3. Longer chains are better, in a specific way.**
In a linear energy funnel, enhancement grows almost linearly with chain length
(~2.12× per site added) up to about N = 15. But longer chains also need
*gentler* noise — the optimal dephasing rate falls as N⁻¹·²⁴. The actual
FMO-7 photosynthetic complex achieves **32.1×** at a dephasing rate that falls
squarely within the biological room-temperature window. Evolution found the
sweet spot.

**4. Disorder makes it stronger — not weaker.**
This one surprised us. We ran 100 random disorder configurations for each chain
length. ENAQT appeared in **95–100% of all of them**. And at large N, random
disorder *outperforms* the carefully designed energy funnel:

| Chain length | Ordered funnel | Disordered median | Disordered mean |
|:---:|:---:|:---:|:---:|
| N = 7  | 22.8× | 8.9×   | 97×    |
| N = 10 | 32.4× | 35.7×  | 595×   |
| N = 15 | 37.9× | **244×** | **6,916×** |

The reason: disorder creates Anderson localization — it traps excitation in the
coherent limit — making the noise-assisted route look even more miraculous by
comparison. The mean grows as σ⁵⁻⁶ with disorder strength. Structural
heterogeneity isn't the enemy of quantum transport. It's a resource.

---

## What's in this repo

| File | What it does |
|------|-------------|
| `enaqt_sb_analysis.py` | Loads 1,000 HEOM trajectories, validates the Lindblad model |
| `enaqt_sb_sink.py` | Adds the reaction center sink, sweeps energy bias and sink rate |
| `enaqt_nsite_chain.py` | Scales to N = 2–20 sites, fits scaling laws, benchmarks FMO-7 |
| `enaqt_disorder_ensemble.py` | 100-seed disorder ensemble + disorder strength sweep (~60s) |
| `PAPER_ENAQT_DRAFT.md` | Full paper draft (submission-ready) |
| `main.tex` + `references.bib` | LaTeX source for journal/preprint submission |

Each script is self-contained and writes figures + a JSON results file when run.

---

## Quick start

```bash
pip install numpy scipy matplotlib
python enaqt_sb_sink.py        # no external data needed
python enaqt_nsite_chain.py    # no external data needed
python enaqt_disorder_ensemble.py  # no external data needed, ~60s
```

The first script (`enaqt_sb_analysis.py`) requires the QD3SET-1 spin-boson
dataset — download from
[figshare](https://doi.org/10.25452/figshare.plus.c.6389553) and place the
`.npy` files in `../SB/data/`.

---

## Reproducing the paper figures

| Figure | Script | Output file |
|--------|--------|-------------|
| HEOM bell curves | `enaqt_sb_analysis.py` | `enaqt_sb_analysis.png` |
| Sink vs. no-sink | `enaqt_sb_sink.py` | `enaqt_sink_vs_nosink.png` |
| N-site scaling | `enaqt_nsite_chain.py` | `enaqt_nsite_scaling.png` |
| Disorder ensemble | `enaqt_disorder_ensemble.py` | `enaqt_disorder_paper_figure.png` |

---

## Citation

Preprint available on bioRxiv (link forthcoming).

```bibtex
@article{harper2026enaqt,
  author  = {Harper, Alexander},
  title   = {Environment-Assisted Quantum Transport in Open Quantum Chains:
             Validation, Scaling Laws, and Disorder Universality},
  year    = {2026},
  note    = {bioRxiv preprint}
}
```

This work uses the QD3SET-1 database:

> Ullah et al. (2023). *Frontiers in Physics* **11**, 1223973.
> https://doi.org/10.3389/fphy.2023.1223973

And builds on the original ENAQT theory:

> Rebentrost et al. (2009). *New Journal of Physics* **11**, 033003.
> https://doi.org/10.1088/1367-2630/11/3/033003
