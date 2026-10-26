"""
Pin the manuscript's key analytical numbers (fast, deterministic).

All numbers use the exact steady-state yield
``eta_inf = kappa * [-L^-1 rho_0]_NN`` from :mod:`enaqt.core`, with

    Gamma (fluorescence) = 0.01,  kappa (sink) = 0.1,  Delta = 1,

and the dephasing rate swept over ``np.logspace(-3, 3, n_points)``.
The enhancement is ``max(eta) / eta(gamma_phi = 1e-3)``.

Grid sizes follow the scripts that produced the manuscript numbers:
300 points for the 2-site Bloch model (``enaqt_sb_sink.py``) and 200 points
for the N-site chains and FMO-7 (``enaqt_nsite_chain.py``).  The location of
the optimum ``gamma_phi*`` is a grid point, so it is compared with a looser
tolerance than the enhancement itself.
"""

import numpy as np
import pytest

from enaqt.core import (
    analytical_yield,
    hamiltonian_2site,
    hamiltonian_flat,
    hamiltonian_fmo7,
    hamiltonian_funnel,
    liouvillian,
    liouvillian_parts,
    optimal_dephasing,
)
from enaqt.spinboson import _analytical_yield_4x4

GAMMA = 0.01
KAPPA = 0.1
DELTA = 1.0
FUNNEL_BIAS = 5.0


def _sweep(H, n_points):
    """Return (gamma_phi, eta) over logspace(-3, 3, n_points)."""
    n = H.shape[0]
    L_base, L_deph = liouvillian_parts(H, KAPPA, GAMMA)
    return optimal_dephasing(L_base, L_deph, KAPPA, n, n_points=n_points)


def _peak(H, n_points):
    """Return (enhancement, gamma_phi_star, eta_peak, peak_index)."""
    gp, eta = _sweep(H, n_points)
    idx = int(np.argmax(eta))
    return float(eta[idx] / eta[0]), float(gp[idx]), float(eta[idx]), idx


# ---------------------------------------------------------------------------
#  2-site Bloch-convention model  H = (eps/2) sigma_z + (Delta/2) sigma_x
# ---------------------------------------------------------------------------


class TestTwoSiteBloch:
    def test_eps5_headline(self):
        enh, gp_star, eta_peak, idx = _peak(hamiltonian_2site(5.0, DELTA), 300)
        assert enh == pytest.approx(7.20, rel=1e-2)
        assert gp_star == pytest.approx(4.92, rel=2e-2)
        assert eta_peak == pytest.approx(0.704, rel=1e-2)
        assert 0 < idx < 299, "optimum must be an interior peak"

    @pytest.mark.parametrize("eps, expected", [(1.0, 1.26), (2.0, 2.05)])
    def test_bias_dependence(self, eps, expected):
        enh, _, _, _ = _peak(hamiltonian_2site(eps, DELTA), 300)
        assert enh == pytest.approx(expected, rel=1e-2)

    def test_enhancement_grows_with_bias(self):
        enh = [_peak(hamiltonian_2site(e, DELTA), 300)[0] for e in (1, 2, 5)]
        assert enh[0] < enh[1] < enh[2]

    def test_core_matches_4x4_bloch_formula(self):
        """N=2 Liouvillian yield == closed 4x4 Bloch yield to 1e-12."""
        max_diff = 0.0
        for eps in (0.0, 1.0, 2.0, 5.0):
            for gp in np.logspace(-3, 3, 25):
                for kappa, gamma in ((KAPPA, GAMMA), (1.0, 0.0), (0.01, 0.1)):
                    H = hamiltonian_2site(eps, DELTA)
                    L = liouvillian(H, gp, kappa, gamma)
                    a = analytical_yield(L, kappa, 2)
                    b = _analytical_yield_4x4(eps, DELTA, gp, kappa, gamma)
                    max_diff = max(max_diff, abs(a - b))
        assert max_diff < 1e-12

    def test_funnel_n2_is_bloch_with_double_delta(self):
        """Convention check: the N-site funnel couples sites with J = Delta,
        whereas the Bloch model uses Delta/2, so funnel(N=2) == 2site(eps, 2J)."""
        np.testing.assert_allclose(
            hamiltonian_funnel(2, FUNNEL_BIAS, 1.0),
            hamiltonian_2site(FUNNEL_BIAS, 2.0),
        )


# ---------------------------------------------------------------------------
#  N-site chains (enaqt_nsite_chain.py conventions: bias 5 Delta, J = Delta)
# ---------------------------------------------------------------------------


class TestNSiteFunnel:
    @pytest.mark.parametrize("n, expected", [(7, 22.8), (10, 32.4), (15, 37.9)])
    def test_funnel_enhancement(self, n, expected):
        H = hamiltonian_funnel(n, FUNNEL_BIAS, DELTA)
        enh, _, _, idx = _peak(H, 200)
        assert enh == pytest.approx(expected, rel=1e-2)
        assert 0 < idx < 199

    def test_optimal_dephasing_decreases_with_length(self):
        gps = [
            _peak(hamiltonian_funnel(n, FUNNEL_BIAS, DELTA), 200)[1]
            for n in (7, 10, 15)
        ]
        assert gps[0] > gps[1] > gps[2]

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_flat_chain_has_no_enaqt(self, n):
        enh, _, _, idx = _peak(hamiltonian_flat(n, DELTA), 200)
        assert idx == 0, "flat chain: yield is maximal in the coherent limit"
        assert enh == pytest.approx(1.0, abs=1e-12)


class TestFMO7:
    def test_fmo7_enhancement(self):
        enh, gp_star, _, idx = _peak(hamiltonian_fmo7(), 200)
        assert enh == pytest.approx(32.1, rel=1e-2)
        assert gp_star == pytest.approx(1.57, rel=2e-2)
        assert 0 < idx < 199

    def test_fmo7_units(self):
        H = hamiltonian_fmo7()
        assert H[0, 0] == 0.0
        assert H[0, 1] == pytest.approx(-1.0)  # J_12 = 87.7 cm^-1 is the unit
        np.testing.assert_allclose(H, H.T)


def test_sweep_matches_full_liouvillian():
    """L_base + gp * L_deph decomposition reproduces the full Liouvillian."""
    H = hamiltonian_funnel(4, FUNNEL_BIAS, DELTA)
    gp, eta = _sweep(H, 20)
    for g, e in zip(gp, eta):
        L = liouvillian(H, g, KAPPA, GAMMA)
        assert analytical_yield(L, KAPPA, 4) == pytest.approx(e, rel=1e-10, abs=1e-14)
