"""Quick figure regeneration from saved JSON — no recomputation needed."""
import numpy as np
import matplotlib.pyplot as plt
import json

with open("enaqt_optimal_results.json") as f:
    results = {int(k): v for k, v in json.load(f).items()}

col_opt    = "#2563EB"
col_funnel = "#DC2626"
col_med    = "#16A34A"
col_mean   = "#D97706"

fig = plt.figure(figsize=(16, 10))
fig.suptitle("Optimal Site-Energy Disorder for Maximum ENAQT Enhancement",
             fontsize=14, fontweight="bold", y=0.98)

# Row 1: energy profiles for N = 5, 7, 10
for col_idx, N in enumerate([5, 7, 10]):
    ax    = fig.add_subplot(2, 3, col_idx + 1)
    eps   = np.array(results[N]["optimal_eps"])
    fun   = np.array([(N - 1 - i) * 1.0 for i in range(N)]); fun[-1] = 0.0
    sites = np.arange(1, N + 1)

    ax.axhline(0, color="silver", lw=0.8, ls="--")
    ax.plot(sites, fun, "s--", color=col_funnel, lw=1.5, ms=7, alpha=0.75, label="Ordered funnel")
    ax.plot(sites, eps, "o-",  color=col_opt,    lw=2,   ms=9, label="Optimal disorder")
    ax.scatter([N], [0], s=120, marker="*", color="gold", zorder=5, label="Sink (site N)")

    ax.set_xlabel("Site index", fontsize=11)
    ax.set_ylabel("Site energy (D)", fontsize=11)
    ax.set_title(
        f"N={N}  |  {results[N]['optimal_enhancement']:.2e}x optimal  "
        f"vs  {results[N]['funnel_enhancement']:.2e}x funnel",
        fontsize=9
    )
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xticks(sites)
    ax.grid(True, alpha=0.25)

# Row 2 left: enhancement scaling
ax      = fig.add_subplot(2, 3, 4)
Ns      = sorted(results.keys())
opt_e   = [results[N]["optimal_enhancement"]  for N in Ns]
fun_e   = [results[N]["funnel_enhancement"]    for N in Ns]
med_e   = [results[N]["random_median"]         for N in Ns]
mean_e  = [results[N]["random_mean"]           for N in Ns]

ax.plot(Ns, opt_e,  "o-",  color=col_opt,    lw=2,   ms=8, label="Optimal disorder")
ax.plot(Ns, fun_e,  "s--", color=col_funnel, lw=1.5, ms=6, label="Ordered funnel")
ax.plot(Ns, med_e,  "^:",  color=col_med,    lw=1.5, ms=6, label="Random median (s=2D)")
ax.plot(Ns, mean_e, "v:",  color=col_mean,   lw=1.5, ms=6, alpha=0.7, label="Random mean (s=2D)")
ax.set_yscale("log")
ax.set_xlabel("Chain length N", fontsize=11)
ax.set_ylabel("Enhancement ratio (log)", fontsize=11)
ax.set_title("Enhancement scaling comparison", fontsize=11)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.25, which="both")

# Row 2 middle: optimal gp* scaling
ax     = fig.add_subplot(2, 3, 5)
gp_opt = [results[N]["optimal_gp_star"]  for N in Ns]
gp_fun = [results[N]["funnel_gp_star"]   for N in Ns]

ax.plot(Ns, gp_opt, "D-",  color=col_opt,    lw=2,   ms=8, label="Optimal disorder")
ax.plot(Ns, gp_fun, "s--", color=col_funnel, lw=1.5, ms=6, label="Ordered funnel")

slope, intercept = np.polyfit(np.log(Ns), np.log(gp_opt), 1)
N_fit = np.linspace(min(Ns), max(Ns), 200)
ax.plot(N_fit, np.exp(intercept) * N_fit**slope, "k--", lw=1.5,
        label=f"gp* ~ N^{slope:.2f}")

ax.set_xlabel("Chain length N", fontsize=11)
ax.set_ylabel("gp* (D)", fontsize=11)
ax.set_title("Optimal dephasing rate vs N", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.25)

# Row 2 right: gain over baselines
ax     = fig.add_subplot(2, 3, 6)
gain_f = [results[N]["gain_vs_funnel"]     for N in Ns]
gain_r = [results[N]["gain_vs_random_med"] for N in Ns]

ax.plot(Ns, gain_f, "o-",  color=col_opt, lw=2,   ms=8, label="Optimal / funnel")
ax.plot(Ns, gain_r, "s--", color=col_med, lw=1.5, ms=6, label="Optimal / random median")
ax.axhline(1, color="silver", lw=1)
ax.set_yscale("log")
ax.set_xlabel("Chain length N", fontsize=11)
ax.set_ylabel("Gain factor (log)", fontsize=11)
ax.set_title("Advantage of optimal over baselines", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.25, which="both")

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("enaqt_optimal_disorder.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved -> enaqt_optimal_disorder.png")

# Terminal summary
print("\n" + "="*65)
print(f"  {'N':>3}  {'Optimal':>16}  {'Funnel':>16}  {'Opt>Funnel?':>12}")
print("  " + "-"*55)
for N in Ns:
    r = results[N]
    winner = "OPT" if r["optimal_enhancement"] > r["funnel_enhancement"] else "FUNNEL"
    print(f"  {N:>3}  {r['optimal_enhancement']:>16.3e}x  "
          f"{r['funnel_enhancement']:>16.3e}x  {winner:>12}")
print("="*65)
