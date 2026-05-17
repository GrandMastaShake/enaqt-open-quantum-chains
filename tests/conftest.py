"""
Shared pytest fixtures for the ENAQT test suite.

These fixtures provide standard Hamiltonians and Liouvillians used across
the 5 test groups, ensuring consistent test inputs and avoiding duplication.
"""

import numpy as np
import pytest

# Physical constants used throughout the test suite
DELTA = 1.0     # tunneling matrix element (energy unit)
KAPPA = 0.1     # sink rate
GAMMA = 0.01    # fluorescence recombination rate


@pytest.fixture
def fixture_2site_H():
    """Standard 2-site spin-boson Hamiltonian with epsilon=1, Delta=1.

    H = epsilon/2 * sigma_z + Delta/2 * sigma_x
      = [[+epsilon/2, Delta/2],
         [Delta/2, -epsilon/2]]

    With epsilon=1, Delta=1:
      H = [[+0.5, 0.5],
           [0.5, -0.5]]
    """
    from enaqt.core import hamiltonian_2site
    return hamiltonian_2site(epsilon=1.0, delta=DELTA)


@pytest.fixture
def fixture_3site_H():
    """Standard 3-site chain Hamiltonian with bias=3, coupling=1.

    Energy funnel: site 1 highest (+1.5), site 3 lowest (-1.5).
    Nearest-neighbor coupling = Delta = 1.0.
    """
    from enaqt.core import hamiltonian_funnel
    return hamiltonian_funnel(n=3, bias=3.0, coupling=DELTA)


@pytest.fixture
def fixture_liouvillian_2site():
    """Pre-built Liouvillian for a 2-site system (N^2 = 4 dimensional).

    Parameters: epsilon=1, Delta=1, gamma_phi=0.1, kappa=0.1, Gamma=0.01.
    Used for testing Liouvillian properties without repeated construction.
    """
    from enaqt.core import hamiltonian_2site, liouvillian
    H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
    return liouvillian(H, dephasing_rate=0.1, sink_rate=KAPPA, fluo_rate=GAMMA)


@pytest.fixture
def fixture_liouvillian_2site_closed():
    """Closed-system Liouvillian (no sink, no fluorescence) for conservation tests.

    With kappa=Gamma=0, the Liouvillian has a zero eigenvalue and trace is preserved.
    """
    from enaqt.core import hamiltonian_2site, liouvillian
    H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
    return liouvillian(H, dephasing_rate=0.1, sink_rate=0.0, fluo_rate=0.0)


@pytest.fixture
def rng():
    """Reproducible random number generator for disorder tests."""
    return np.random.default_rng(42)
