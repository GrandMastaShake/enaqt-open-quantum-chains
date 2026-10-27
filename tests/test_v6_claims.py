"""
Pin the qualitative claims added in manuscript v6 (fast, deterministic).

* Mismatch rule: in a linear funnel the optimal dephasing rate equals the
  nearest-neighbour mismatch ``bias / (N - 1)``.
* Bias symmetry: on a chain the yield is unchanged by ``eps_j -> -eps_j``
  (uphill and downhill funnels give the same yield).
* Flat chains show no ENAQT and beat the funnel's peak yield.
* Two-site resonance: ``gamma_phi*`` tracks the bias.
"""

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from enaqt.core import (
    analytical_yield,
    hamiltonian_2site,
    hamiltonian_flat,
    hamiltonian_funnel,
    liouvillian_parts,
)

GAMMA = 0.01
KAPPA = 0.1
BIAS = 5.0


def _eta(H, gp, kappa=KAPPA, gamma=GAMMA):
    L_base, L_deph = liouvillian_parts(H, kappa, gamma)
    return analytical_yield(L_base + gp * L_deph, kappa, H.shape[0])


def _gp_star(H, kappa=KAPPA, gamma=GAMMA):
    L_base, L_deph = liouvillian_parts(H, kappa, gamma)
    n = H.shape[0]
    res = minimize_scalar(
        lambda x: -analytical_yield(L_base + 10.0**x * L_deph, kappa, n),
        bounds=(-3.0, 3.0),
        method="bounded",
        options={"xatol": 1e-6},
    )
    return 10.0**res.x, -res.fun


@pytest.mark.parametrize("n", [3, 5, 10, 15])
def test_mismatch_rule(n):
    gp_star, _ = _gp_star(hamiltonian_funnel(n, BIAS, 1.0))
    assert gp_star / (BIAS / (n - 1)) == pytest.approx(1.0, abs=0.03)


@pytest.mark.parametrize("n", [3, 6, 9])
@pytest.mark.parametrize("gp", [0.0, 0.3, 3.0])
def test_bias_symmetry(n, gp):
    H = hamiltonian_funnel(n, BIAS, 1.0)
    H_up = H.copy()
    np.fill_diagonal(H_up, -np.diag(H))
    assert _eta(H_up, gp) == pytest.approx(_eta(H, gp), rel=1e-10)


@pytest.mark.parametrize("n", [3, 5, 7])
@pytest.mark.parametrize("kappa", [0.01, 1.0, 100.0])
@pytest.mark.parametrize("gamma", [1e-3, 1e-1])
def test_flat_chain_no_enaqt(n, kappa, gamma):
    H = hamiltonian_flat(n, 1.0)
    eta = [_eta(H, gp, kappa, gamma) for gp in np.logspace(-3, 3, 61)]
    assert np.all(np.diff(eta) <= 1e-12)


@pytest.mark.parametrize("n", [3, 7, 15])
def test_flat_beats_funnel(n):
    flat = _eta(hamiltonian_flat(n, 1.0), 0.0)
    _, funnel_peak = _gp_star(hamiltonian_funnel(n, BIAS, 1.0))
    assert flat > funnel_peak


@pytest.mark.parametrize("eps", [0.5, 1.0, 2.0, 5.0])
def test_two_site_resonance(eps):
    gp_star, _ = _gp_star(hamiltonian_2site(eps, 1.0))
    assert gp_star / eps == pytest.approx(1.0, abs=0.15)
