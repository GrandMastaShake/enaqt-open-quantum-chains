"""Recompute every number quoted in main.tex (v5) from scratch.

Self-contained: numpy only (plus the optimiser output enaqt_optimal_results.json
and the 2D results enaqt_2d_results.json for the two sections that reuse them).

Model (Sec. II of the paper)
    d rho/dt = -i[H, rho] + gamma_phi * D[rho] - (kappa/2){P_sink, rho} - Gamma rho
    yield eta = kappa * [ -L^{-1} vec(rho0) ]_{sink, sink}

Conventions stated in the paper:
    * N-site chains use nearest-neighbour coupling Delta = 1.
    * The two-site "Bloch" model uses H = (eps/2) sigma_z + (Delta/2) sigma_x,
      i.e. inter-site coupling Delta/2.
    * eta_zero is the yield at the lowest grid point gamma_phi = 1e-3 (1D),
      as in the original scripts; the true gamma_phi = 0 value is also reported.

Usage:  python -I scripts/paper_numbers.py  [out.json]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
KAPPA, GAMMA = 0.1, 0.01
GP300 = np.logspace(-3, 3, 300)
GP100 = np.logspace(-3, 3, 100)


# ----------------------------------------------------------------------------
# Liouvillian and yield
# ----------------------------------------------------------------------------
def liouvillian_parts(
    H: np.ndarray, sink: int, kappa: float = KAPPA, gamma: float = GAMMA
):
    """Return (L_base, L_deph) in the column-stack convention."""
    n = H.shape[0]
    eye = np.eye(n)
    l_base = -1j * (np.kron(eye, H) - np.kron(H.T, eye))
    p_s = np.zeros((n, n))
    p_s[sink, sink] = 1.0
    l_base = l_base - 0.5 * kappa * (np.kron(eye, p_s) + np.kron(p_s, eye))
    l_base = l_base - gamma * np.eye(n * n)
    l_deph = np.zeros((n * n, n * n), dtype=complex)
    for j in range(n):
        p_j = np.zeros((n, n))
        p_j[j, j] = 1.0
        l_deph += np.kron(p_j, p_j) - 0.5 * np.kron(eye, p_j) - 0.5 * np.kron(p_j, eye)
    return l_base, l_deph


def yields(H, gps, sink=None, src=0, kappa=KAPPA, gamma=GAMMA) -> np.ndarray:
    n = H.shape[0]
    sink = n - 1 if sink is None else sink
    l_base, l_deph = liouvillian_parts(H, sink, kappa, gamma)
    rho0 = np.zeros(n * n, dtype=complex)
    rho0[src + src * n] = 1.0
    idx = sink + sink * n
    return np.array(
        [kappa * (-np.linalg.solve(l_base + g * l_deph, rho0))[idx].real for g in gps]
    )


def chain(eps, coupling: float = 1.0) -> np.ndarray:
    eps = np.asarray(eps, dtype=float)
    H = np.diag(eps)
    for j in range(len(eps) - 1):
        H[j, j + 1] = H[j + 1, j] = coupling
    return H


def funnel(n: int, total: float = 5.0) -> np.ndarray:
    return chain([total / 2 - j * total / (n - 1) for j in range(n)])


def two_site(eps: float) -> np.ndarray:
    return np.array([[eps / 2, 0.5], [0.5, -eps / 2]])


def summarize(gps, eta, eta_true0) -> dict:
    i = int(np.argmax(eta))
    return dict(
        gp_star=float(gps[i]),
        eta_peak=float(eta[i]),
        eta0_grid=float(eta[0]),
        eta0_true=float(eta_true0),
        enh=float(eta[i] / eta[0]),
        enh_true0=float(eta[i] / eta_true0),
        eta_zeno=float(eta[-1]),
        interior=bool(0 < i < len(eta) - 1),
    )


def linfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    a, b = np.polyfit(x, y, 1)
    r2 = 1 - np.sum((y - (a * x + b)) ** 2) / np.sum((y - y.mean()) ** 2)
    return float(a), float(b), float(r2)


# ----------------------------------------------------------------------------
# FMO (Adolphs & Renger 2006 site energies / couplings, cm^-1)
# ----------------------------------------------------------------------------
FMO_CM = np.array(
    [
        [12445, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9],
        [-87.7, 12520, 30.8, 8.2, 0.7, 11.8, 4.3],
        [5.5, 30.8, 12205, -53.5, -2.2, -9.6, 6.0],
        [-5.9, 8.2, -53.5, 12335, -70.7, -17.0, -63.3],
        [6.7, 0.7, -2.2, -70.7, 12490, 81.1, -1.3],
        [-13.7, 11.8, -9.6, -17.0, 81.1, 12640, 39.7],
        [-9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 12450],
    ],
    dtype=float,
)


def fmo_hamiltonian() -> np.ndarray:
    H = FMO_CM / 87.7
    return H - np.eye(7) * H[0, 0]


def main(out_path: Path) -> dict:
    t0 = time.time()
    out: dict = {
        "conventions": {"kappa": KAPPA, "Gamma": GAMMA, "gp_grid_1d": [1e-3, 1e3]}
    }

    # 1. two-site model, energy-bias sweep and sink-rate sweep -----------------
    out["two_site_eps"] = {}
    for eps in [0.0, 0.5, 1.0, 2.0, 5.0]:
        H = two_site(eps)
        out["two_site_eps"][str(eps)] = summarize(
            GP300, yields(H, GP300), yields(H, [0.0])[0]
        )
    out["two_site_kappa_eps5"] = {}
    for k in [0.01, 0.1, 1.0]:
        H = two_site(5.0)
        e = yields(H, GP300, kappa=k)
        s = summarize(GP300, e, yields(H, [0.0], kappa=k)[0])
        s["peak_over_zeno"] = float(e.max() / e[-1])
        out["two_site_kappa_eps5"][str(k)] = s
    # N = 2 chain with coupling Delta (the funnel convention) at eps_tot = 5
    out["two_site_coupling_delta_eps5"] = summarize(
        GP300, yields(funnel(2), GP300), yields(funnel(2), [0.0])[0]
    )

    # 2. funnel and flat chains ---------------------------------------------
    ns = list(range(2, 21))
    out["funnel"] = {}
    for n in ns:
        H = funnel(n)
        out["funnel"][str(n)] = summarize(GP300, yields(H, GP300), yields(H, [0.0])[0])
    out["flat"] = {}
    for n in [3, 7, 15]:
        H = chain([0.0] * n)
        out["flat"][str(n)] = summarize(GP300, yields(H, GP300), yields(H, [0.0])[0])
    enh = np.array([out["funnel"][str(n)]["enh"] for n in ns])
    gps = np.array([out["funnel"][str(n)]["gp_star"] for n in ns])
    n_arr = np.array(ns, float)
    m = n_arr <= 10
    a, b, r2 = linfit(n_arr[m], enh[m])
    a_all, b_all, r2_all = linfit(n_arr, enh)
    p, c, r2p = linfit(np.log(n_arr), np.log(gps))
    out["funnel_fits"] = dict(
        lin_le10_slope=a,
        lin_le10_intercept=b,
        lin_le10_r2=r2,
        lin_all_slope=a_all,
        lin_all_r2=r2_all,
        gp_star_power=p,
        gp_star_prefactor=float(np.exp(c)),
        gp_star_r2=r2p,
        plateau_n12_20=[float(enh[n_arr >= 12].min()), float(enh[n_arr >= 12].max())],
    )

    # 3. FMO ------------------------------------------------------------------
    Hf = fmo_hamiltonian()
    out["fmo"] = {}
    for label, sink in [("sink_bchl3", 2), ("sink_bchl7", 6)]:
        out["fmo"][label] = summarize(
            GP300, yields(Hf, GP300, sink=sink), yields(Hf, [0.0], sink=sink)[0]
        )
    # Ishizaki-Fleming high-T estimate of the pure-dephasing rate, in Delta_FMO units
    lam, wc_inv_fs, kT = 35.0, 50.0, 208.5  # cm^-1, fs, cm^-1 (300 K)
    gamma_c_cm = 1.0 / (wc_inv_fs * 1e-15 * 2.99792458e10 * 2 * np.pi)
    out["fmo"]["dephasing_estimate_delta_units"] = float(
        2 * lam * kT / gamma_c_cm / 87.7
    )
    out["fmo"]["delta_fmo_period_fs"] = float(1e15 / (87.7 * 2.99792458e10 * 2 * np.pi))

    # 4. disorder ensemble ----------------------------------------------------
    def dis_h(n, sigma, seed):
        rng = np.random.default_rng(seed)
        return chain(rng.normal(0.0, sigma, n))

    out["disorder"] = {}
    total_interior = 0
    total_real = 0
    for n in [2, 3, 4, 5, 6, 7, 8, 10, 12, 15]:
        rows = []
        for s in range(100):
            H = dis_h(n, 2.0, s)
            rows.append(summarize(GP100, yields(H, GP100), yields(H, [0.0])[0]))
        eg = np.array([r["enh"] for r in rows])
        ep = np.array([r["eta_peak"] for r in rows])
        e0 = np.array([r["eta0_grid"] for r in rows])
        enaqt = np.array([r["interior"] and r["enh"] > 1.01 for r in rows])
        total_interior += int(enaqt.sum())
        total_real += len(rows)
        out["disorder"][str(n)] = dict(
            frac_enaqt=float(enaqt.mean()),
            median=float(np.median(eg)),
            p25=float(np.percentile(eg, 25)),
            p75=float(np.percentile(eg, 75)),
            mean=float(eg.mean()),
            max=float(eg.max()),
            median_true0=float(np.median([r["enh_true0"] for r in rows])),
            median_eta_peak=float(np.median(ep)),
            median_eta0=float(np.median(e0)),
            p25_eta_peak=float(np.percentile(ep, 25)),
            p75_eta_peak=float(np.percentile(ep, 75)),
            p25_eta0=float(np.percentile(e0, 25)),
            p75_eta0=float(np.percentile(e0, 75)),
        )
    out["disorder_total"] = dict(enaqt=total_interior, realizations=total_real)

    out["sigma_sweep_n7"] = {}
    sig = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
    for sg in sig:
        eg, ep, e0 = [], [], []
        for s in range(50):
            e = yields(dis_h(7, sg, s), GP100)
            eg.append(e.max() / e[0])
            ep.append(e.max())
            e0.append(e[0])
        eg = np.array(eg)
        out["sigma_sweep_n7"][str(sg)] = dict(
            median=float(np.median(eg)),
            p25=float(np.percentile(eg, 25)),
            p75=float(np.percentile(eg, 75)),
            mean=float(eg.mean()),
            median_eta_peak=float(np.median(ep)),
            median_eta0=float(np.median(e0)),
        )
    sg_arr = np.array(sig)
    med = np.array([out["sigma_sweep_n7"][str(s)]["median"] for s in sig])
    mask = sg_arr >= 1.0
    a, b, r2 = linfit(sg_arr[mask], np.log(med[mask]))
    out["sigma_fit"] = dict(log_median_slope=a, log_median_intercept=b, r2=r2)

    # 5. optimal-disorder configurations: absolute yields -----------------------
    opt = json.loads((ROOT / "enaqt_optimal_results.json").read_text())
    gp30 = np.logspace(-3, 2, 30)
    out["optimal"] = {}
    for key, v in opt.items():
        n = int(key)
        Ho = chain(v["optimal_eps"])
        Hl = chain([(n - 1 - i) * 1.0 for i in range(n)])
        eo, el = yields(Ho, gp30), yields(Hl, gp30)
        out["optimal"][key] = dict(
            eps=[float(x) for x in v["optimal_eps"]],
            opt_eta_peak=float(eo.max()),
            opt_eta0=float(yields(Ho, [1e-8])[0]),
            funnel_eta_peak=float(el.max()),
            funnel_eta0=float(yields(Hl, [1e-8])[0]),
            opt_gp_star=float(gp30[int(np.argmax(eo))]),
        )

    # 6. 2D lattice (reused from enaqt_2d_lattice.py output) ----------------------
    out["lattice_2d"] = json.loads((ROOT / "enaqt_2d_results.json").read_text())
    for v in out["lattice_2d"].values():
        v.pop("disorder_all", None)

    out["runtime_s"] = time.time() - t0
    out_path.write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "paper_numbers.json"
    main(target)
    print("wrote", target)
