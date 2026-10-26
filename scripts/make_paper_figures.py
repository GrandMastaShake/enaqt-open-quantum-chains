"""Draw the four main-text figures of main.tex (v5) from paper_numbers.json.

Usage:  python -I scripts/make_paper_figures.py [paper_numbers.json] [out_dir]
Writes fig1_two_site.png, fig2_funnel.png, fig3_disorder.png, fig4_optimal.png.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_numbers import GP300, funnel, two_site, yields  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUES = ["#b9d3f2", "#8db7ea", "#5f99df", "#2a78d6", "#1f5fae", "#164784", "#0e2f59"]

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "savefig.dpi": 300,
    }
)


def tag(ax, letter):
    ax.text(
        -0.14, 1.04, f"({letter})", transform=ax.transAxes, fontweight="bold", color=INK
    )


def fig_two_site(d, out):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8), constrained_layout=True)
    for eps, col in zip([0.0, 1.0, 2.0, 5.0], [BLUES[1], BLUES[3], BLUES[5], ORANGE]):
        e = yields(two_site(eps), GP300)
        a.semilogx(GP300, e, color=col, label=rf"$\varepsilon={eps:g}\Delta$")
        i = int(np.argmax(e))
        if 0 < i < len(e) - 1:
            a.plot(GP300[i], e[i], "o", ms=5, color=col, mec="white", mew=1)
    a.set_xlabel(r"dephasing rate $\gamma_\phi$ [$\Delta$]")
    a.set_ylabel(r"yield $\eta$")
    a.set_ylim(0, 1)
    a.legend(loc="upper right", fontsize=7, ncol=2)
    tag(a, "a")
    for k, col in zip([0.01, 0.1, 1.0], [BLUE, ORANGE, AQUA]):
        e = yields(two_site(5.0), GP300, kappa=k)
        b.semilogx(GP300, e, color=col, label=rf"$\kappa={k:g}\Delta$")
    b.set_xlabel(r"dephasing rate $\gamma_\phi$ [$\Delta$]")
    b.set_ylabel(r"yield $\eta$ ($\varepsilon=5\Delta$)")
    b.set_ylim(0, 1)
    b.legend(loc="upper left", fontsize=8)
    tag(b, "b")
    fig.savefig(out / "fig1_two_site.png")
    plt.close(fig)


def fig_funnel(d, out):
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6), constrained_layout=True)
    a, b, c = axs
    for n, col in zip([2, 3, 5, 7, 10, 15, 20], BLUES):
        a.semilogx(GP300, yields(funnel(n), GP300), color=col, lw=1.6, label=f"N={n}")
    a.set_xlabel(r"$\gamma_\phi$ [$\Delta$]")
    a.set_ylabel(r"yield $\eta$")
    a.legend(fontsize=7, ncol=2, loc="upper left")
    tag(a, "a")
    ns = np.array(sorted(int(k) for k in d["funnel"]))
    enh = np.array([d["funnel"][str(n)]["enh"] for n in ns])
    gps = np.array([d["funnel"][str(n)]["gp_star"] for n in ns])
    f = d["funnel_fits"]
    enh0 = np.array([d["funnel"][str(n)]["enh_true0"] for n in ns])
    b.semilogy(
        ns,
        enh,
        "o",
        ms=4.5,
        color=BLUE,
        mec="white",
        mew=0.8,
        label=r"baseline $\gamma_\phi=10^{-3}\Delta$",
    )
    b.semilogy(
        ns,
        enh0,
        "^",
        ms=4.5,
        color=ORANGE,
        mec="white",
        mew=0.8,
        label=r"baseline $\gamma_\phi=0$",
    )
    x = np.linspace(2, 10, 50)
    b.semilogy(
        x,
        f["lin_le10_slope"] * x + f["lin_le10_intercept"],
        "--",
        color=MUTED,
        lw=1.2,
        label=r"linear fit $N\leq10$",
    )
    b.set_xlabel("chain length N")
    b.set_ylabel(r"enhancement $\eta_{peak}/\eta_{zero}$")
    b.legend(fontsize=7, loc="upper left")
    tag(b, "b")
    c.loglog(ns, gps, "o", ms=4.5, color=BLUE, mec="white", mew=0.8)
    c.loglog(
        ns,
        f["gp_star_prefactor"] * ns ** f["gp_star_power"],
        "--",
        color=MUTED,
        lw=1.2,
        label=rf"$\propto N^{{{f['gp_star_power']:.2f}}}$",
    )
    c.set_xticks([2, 5, 10, 20])
    c.set_yticks([0.3, 1, 3])
    for axis in (c.xaxis, c.yaxis):
        axis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    c.set_xlabel("chain length N")
    c.set_ylabel(r"optimal $\gamma_\phi^*$ [$\Delta$]")
    c.legend(fontsize=8)
    tag(c, "c")
    fig.savefig(out / "fig2_funnel.png")
    plt.close(fig)


def fig_disorder(d, out):
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6), constrained_layout=True)
    a, b, c = axs
    ns = np.array(sorted(int(k) for k in d["disorder"]))
    q = {
        k: np.array([d["disorder"][str(n)][k] for n in ns])
        for k in [
            "median",
            "p25",
            "p75",
            "median_eta_peak",
            "median_eta0",
            "p25_eta_peak",
            "p75_eta_peak",
            "p25_eta0",
            "p75_eta0",
        ]
    }
    ref = np.array([d["funnel"][str(n)]["enh"] for n in ns])
    a.fill_between(ns, q["p25"], q["p75"], color=ORANGE, alpha=0.18, lw=0)
    a.semilogy(ns, q["median"], "o-", ms=4, color=ORANGE, label="disorder median (IQR)")
    a.semilogy(ns, ref, "s-", ms=4, color=BLUE, label="ordered funnel")
    a.set_xlabel("chain length N")
    a.set_ylabel(r"$\eta_{peak}/\eta_{zero}$")
    a.legend(fontsize=7, loc="upper left")
    tag(a, "a")
    b.fill_between(
        ns, q["p25_eta_peak"], q["p75_eta_peak"], color=ORANGE, alpha=0.18, lw=0
    )
    b.semilogy(
        ns, q["median_eta_peak"], "o-", ms=4, color=ORANGE, label=r"$\eta_{peak}$"
    )
    b.fill_between(ns, q["p25_eta0"], q["p75_eta0"], color=AQUA, alpha=0.18, lw=0)
    b.semilogy(ns, q["median_eta0"], "^-", ms=4, color=AQUA, label=r"$\eta_{zero}$")
    b.set_xlabel("chain length N")
    b.set_ylabel(r"absolute yield (median, IQR)")
    b.legend(fontsize=7, loc="lower left")
    tag(b, "b")
    sig = np.array(sorted(float(k) for k in d["sigma_sweep_n7"]))
    s = {
        k: np.array([d["sigma_sweep_n7"][str(x)][k] for x in sig])
        for k in ["median", "p25", "p75"]
    }
    c.fill_between(sig, s["p25"], s["p75"], color=ORANGE, alpha=0.18, lw=0)
    c.semilogy(sig, s["median"], "o-", ms=4, color=ORANGE, label="median (IQR)")
    fit = d["sigma_fit"]
    xs = np.linspace(1, 5, 50)
    c.semilogy(
        xs,
        np.exp(fit["log_median_intercept"] + fit["log_median_slope"] * xs),
        "--",
        color=MUTED,
        lw=1.2,
        label=rf"$\propto e^{{{fit['log_median_slope']:.2f}\sigma/\Delta}}$",
    )
    c.set_xlabel(r"disorder strength $\sigma$ [$\Delta$]  (N=7)")
    c.set_ylabel(r"$\eta_{peak}/\eta_{zero}$")
    c.legend(fontsize=7, loc="upper left")
    tag(c, "c")
    fig.savefig(out / "fig3_disorder.png")
    plt.close(fig)


def fig_optimal(d, out):
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7), constrained_layout=True)
    a, b = axs
    for n, col in zip(["5", "7", "10"], [BLUES[2], BLUES[4], BLUES[6]]):
        eps = d["optimal"][n]["eps"]
        a.plot(
            range(1, len(eps) + 1),
            eps,
            "o-",
            ms=4,
            color=col,
            lw=1.4,
            label=f"optimised N={n}",
        )
    n = 10
    a.plot(
        range(1, n + 1),
        [(n - 1 - i) * 1.0 for i in range(n)],
        "s--",
        ms=3.5,
        color=ORANGE,
        lw=1.2,
        label="funnel N=10",
    )
    a.set_xlabel("site index (sink = last)")
    a.set_ylabel(r"site energy [$\Delta$]")
    a.legend(fontsize=7, loc="upper right")
    tag(a, "a")
    ns = np.array(sorted(int(k) for k in d["optimal"]))
    get = lambda k: np.array([d["optimal"][str(x)][k] for x in ns])  # noqa: E731
    b.semilogy(
        ns,
        get("opt_eta_peak"),
        "o-",
        ms=4,
        color=BLUE,
        label=r"optimised $\eta_{peak}$",
    )
    b.semilogy(
        ns,
        get("opt_eta0"),
        "o--",
        ms=4,
        color=BLUE,
        mfc="white",
        label=r"optimised $\eta_{zero}$",
    )
    b.semilogy(
        ns,
        get("funnel_eta_peak"),
        "s-",
        ms=4,
        color=ORANGE,
        label=r"funnel $\eta_{peak}$",
    )
    b.semilogy(
        ns,
        get("funnel_eta0"),
        "s--",
        ms=4,
        color=ORANGE,
        mfc="white",
        label=r"funnel $\eta_{zero}$",
    )
    b.set_xlabel("chain length N")
    b.set_ylabel("absolute yield")
    b.legend(fontsize=7, loc="lower left", ncol=2)
    tag(b, "b")
    fig.savefig(out / "fig4_optimal.png")
    plt.close(fig)


if __name__ == "__main__":
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "paper_numbers.json"
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    data = json.loads(src.read_text())
    fig_two_site(data, out_dir)
    fig_funnel(data, out_dir)
    fig_disorder(data, out_dir)
    fig_optimal(data, out_dir)
    print("figures written to", out_dir)
