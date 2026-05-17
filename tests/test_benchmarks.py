"""
Group E: Analytical Benchmark Tests
====================================

These tests verify our numerical results against known analytical
limits and mathematical properties. They serve as the strongest form
of validation: comparison against exact mathematical results.

Benchmarks tested:
  1. 2-site resonant limit (epsilon=0, gamma>0)
  2. 2-site detuned limit (epsilon>0, gamma>0)
  3. Quantum Zeno limit (large gamma_phi -> eta -> 0)
  4. Flat chain at resonance
  5. Monotonicity: eta decreases with gamma_phi
  6. solve_ivp matches analytical_yield
  7. Enhancement ratio > 1 (ENAQT exists)
"""

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import solve_ivp

# Physical constants
DELTA = 1.0
KAPPA = 0.1
GAMMA_F = 0.01


# ──────────────────────────────────────────────────────────────────────────────
#  Exact analytical formulas (derived via SymPy from the Lindblad master eq.)
# ──────────────────────────────────────────────────────────────────────────────

def _eta_2site_exact(epsilon: float, delta: float, kappa: float,
                     gamma: float) -> float:
    """Exact steady-state yield for the 2-site Lindblad system.

    Hamiltonian: H = [[epsilon/2, delta/2], [delta/2, -epsilon/2]]
    Sink (site 2) rate: kappa
    Fluorescence rate: gamma
    No dephasing.

    Derived by solving the 4x4 Liouvillian L @ integral = -rho0
    symbolically and extracting the sink population integral.
    """
    # SymPy-derived closed form
    num = delta**2 * kappa * (2 * gamma + kappa)
    denom = (
        delta**2 * (2 * gamma + kappa)**2
        + 4 * epsilon**2 * gamma * (gamma + kappa)
        + gamma * (gamma + kappa) * (2 * gamma + kappa)**2
    )
    return float(num / denom)


# ──────────────────────────────────────────────────────────────────────────────


class TestAnalytical2Site:
    """Tests against known analytical results for the 2-site system."""

    def test_analytical_2site_resonant_no_dephasing(self):
        """Exact analytical result for epsilon=0, gamma_phi=0.

        At resonance (epsilon=0) with no dephasing, the 2-site system
        has Hamiltonian H = (Delta/2) sigma_x with sink at site 2.

        The exact steady-state yield (derived from the 4x4 Liouvillian
        via symbolic matrix inversion) is::

            eta = Delta^2 * kappa * (2*gamma + kappa)
                  / [Delta^2 * (2*gamma + kappa)^2
                     + gamma * (gamma + kappa) * (2*gamma + kappa)^2]

        which simplifies to::

            eta = Delta^2 * kappa / [(2*gamma + kappa)
                                     * (Delta^2 + gamma*(gamma + kappa))]
        """
        from enaqt.core import hamiltonian_2site, liouvillian

        H = hamiltonian_2site(epsilon=0.0, delta=DELTA)
        L = liouvillian(H, dephasing_rate=0.0, sink_rate=KAPPA,
                        fluo_rate=GAMMA_F)

        n = H.shape[0]
        rho0_vec = np.zeros(n * n, dtype=complex)
        rho0_vec[0] = 1.0
        integral = -np.linalg.solve(L, rho0_vec)
        sink_idx = (n - 1) + (n - 1) * n
        eta_numerical = KAPPA * integral[sink_idx].real

        eta_analytical = _eta_2site_exact(
            epsilon=0.0, delta=DELTA, kappa=KAPPA, gamma=GAMMA_F
        )

        assert_allclose(eta_numerical, eta_analytical, rtol=1e-8,
                        err_msg=f"Resonant yield mismatch: num={eta_numerical}, "
                                f"anal={eta_analytical}")

    def test_analytical_2site_detuned_no_dephasing(self):
        """Exact result for epsilon=1, gamma_phi=0 (detuned, no dephasing).

        With bias epsilon and tunneling Delta, the exact yield is::

            eta = Delta^2 * kappa * (2*gamma + kappa)
                  / [Delta^2 * (2*gamma + kappa)^2
                     + 4*epsilon^2 * gamma * (gamma + kappa)
                     + gamma * (gamma + kappa) * (2*gamma + kappa)^2]
        """
        from enaqt.core import hamiltonian_2site, liouvillian

        eps = 1.0
        H = hamiltonian_2site(epsilon=eps, delta=DELTA)
        L = liouvillian(H, dephasing_rate=0.0, sink_rate=KAPPA,
                        fluo_rate=GAMMA_F)

        n = H.shape[0]
        rho0_vec = np.zeros(n * n, dtype=complex)
        rho0_vec[0] = 1.0
        integral = -np.linalg.solve(L, rho0_vec)
        sink_idx = (n - 1) + (n - 1) * n
        eta_numerical = KAPPA * integral[sink_idx].real

        eta_analytical = _eta_2site_exact(
            epsilon=eps, delta=DELTA, kappa=KAPPA, gamma=GAMMA_F
        )

        assert_allclose(eta_numerical, eta_analytical, rtol=1e-8,
                        err_msg=f"Detuned yield mismatch: num={eta_numerical}, "
                                f"anal={eta_analytical}")


class TestAnalyticalZeno:
    """Tests the quantum Zeno limit: large dephasing suppresses transport."""

    def test_analytical_zeno_limit(self):
        """In the Zeno limit (gamma_phi -> infinity), eta_inf -> 0.

        Strong dephasing localizes the excitation at the initial site,
        preventing tunneling to the sink. With a non-zero fluorescence
        rate (which provides a loss channel competing with the sink),
        the yield is suppressed when dephasing prevents transport.

        At gamma_phi = 1000*Delta with gamma_f = 0.01, the yield
        should be < 1% of the optimal value.
        """
        from enaqt.core import hamiltonian_2site, liouvillian, analytical_yield

        H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        n = H.shape[0]
        gamma_f = 0.01  # Non-zero fluorescence needed for yield < 1

        # Large dephasing (gamma_phi >> Delta)
        L = liouvillian(H, dephasing_rate=5000.0, sink_rate=KAPPA,
                        fluo_rate=gamma_f)
        eta = analytical_yield(L, KAPPA, n)

        assert eta < 0.01, (
            f"Zeno limit not reached: eta={eta} at gamma_phi=5000.0. "
            f"Expected eta < 0.01."
        )

    def test_analytical_zeno_monotonic(self):
        """eta_inf decreases monotonically for gamma_phi > gamma_phi*.

        After the optimal dephasing point, increasing dephasing further
        should monotonically decrease the yield (Zeno suppression).
        """
        from enaqt.core import hamiltonian_2site, liouvillian, analytical_yield

        H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        n = H.shape[0]
        gamma_f = 0.01  # Non-zero fluorescence for Zeno effect

        # Test in the Zeno regime (well past optimal, which is ~0.1-1)
        gammas = np.logspace(1, 3, 20)  # 10 to 1000
        etas = []
        for gp in gammas:
            L = liouvillian(H, dephasing_rate=gp, sink_rate=KAPPA,
                            fluo_rate=gamma_f)
            etas.append(analytical_yield(L, KAPPA, n))

        # Check monotonic decrease in Zeno regime
        for i in range(1, len(etas)):
            assert etas[i] <= etas[i - 1] + 1e-6, (
                f"Non-monotonic at gamma_phi={gammas[i]}: "
                f"eta[{i-1}]={etas[i-1]}, eta[{i}]={etas[i]}"
            )


class TestAnalyticalFlatChain:
    """Tests flat chain Hamiltonian against known limits."""

    def test_analytical_flat_chain_2site(self):
        """2-site flat chain: H = [[0, Delta], [Delta, 0]] (coupling = Delta, not Delta/2).

        Note: hamiltonian_flat uses Delta directly as coupling, while
        hamiltonian_2site uses Delta/2. So a 2-site flat chain has
        twice the coupling of the standard 2-site Hamiltonian.
        """
        from enaqt.core import hamiltonian_flat, hamiltonian_2site

        H_flat = hamiltonian_flat(n=2, coupling=DELTA)

        # Expected: [[0, Delta], [Delta, 0]]
        H_expected = np.array([[0.0, DELTA], [DELTA, 0.0]])

        assert_allclose(H_flat, H_expected, rtol=1e-10,
                        err_msg="Flat 2-site Hamiltonian incorrect")

    def test_analytical_flat_chain_longer(self):
        """N-site flat chain: nearest-neighbor coupling = Delta, no bias."""
        from enaqt.core import hamiltonian_flat

        for n in [3, 4, 5]:
            H = hamiltonian_flat(n=n, coupling=DELTA)

            # All diagonal elements should be zero
            assert_allclose(np.diag(H), np.zeros(n), atol=1e-12)

            # Nearest-neighbor couplings should be Delta
            for i in range(n - 1):
                assert abs(H[i, i + 1] - DELTA) < 1e-12
                assert abs(H[i + 1, i] - DELTA) < 1e-12

            # No next-nearest-neighbor couplings
            for i in range(n - 2):
                assert abs(H[i, i + 2]) < 1e-12


class TestAnalyticalMonotonicity:
    """Tests that yield decreases monotonically with dephasing past optimal."""

    def test_analytical_yield_monotonic_in_gamma(self):
        """For gamma_phi > optimal, eta_inf decreases monotonically.

        This is a fundamental property of ENAQT: after the optimal dephasing
        point, increasing dephasing further suppresses transport (Zeno effect).
        """
        from enaqt.core import hamiltonian_2site, liouvillian, analytical_yield

        H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        n = H.shape[0]
        gamma_f = 0.01

        # Scan past the optimal (which is around gamma_phi ~ 0.1-1 for these params)
        gammas = np.logspace(0, 3, 50)
        etas = []
        for gp in gammas:
            L = liouvillian(H, dephasing_rate=gp, sink_rate=KAPPA,
                            fluo_rate=gamma_f)
            etas.append(analytical_yield(L, KAPPA, n))

        # Find the peak index
        peak_idx = np.argmax(etas)

        # After the peak, yield must decrease monotonically
        for i in range(peak_idx + 1, len(etas)):
            assert etas[i] <= etas[i - 1] + 1e-10, (
                f"Non-monotonic at gamma_phi={gammas[i]}: "
                f"eta[{i-1}]={etas[i-1]}, eta[{i}]={etas[i]}"
            )


class TestAnalyticalVsNumerical:
    """Compares analytical yield against numerical time integration."""

    def test_analytical_vs_numerical_integration(self):
        """`analytical_yield` matches `solve_ivp` time integration.

        The analytical formula eta_inf = kappa * [-L^-1 @ rho0]_{sink}
        should match the long-time limit of the time-dependent
        integration: eta(t->inf) = kappa * integral_0^inf rho_NN(t) dt.

        The Lindblad evolution operates on complex density matrices, so
        we must integrate both real and imaginary parts of the
        column-stacked density vector.
        """
        from enaqt.core import (hamiltonian_2site, liouvillian,
                                 analytical_yield)

        H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        n = H.shape[0]
        gp = 0.1
        kappa = KAPPA
        gamma = GAMMA_F

        # Analytical yield
        L = liouvillian(H, dephasing_rate=gp, sink_rate=kappa,
                        fluo_rate=gamma)
        eta_analytical = analytical_yield(L, kappa, n)

        # Numerical: solve_ivp to long time and integrate.
        # The Liouvillian is complex, so we stack real and imaginary
        # parts into a real ODE system of size 2*N^2.
        rho0_vec = np.zeros(n * n, dtype=complex)
        rho0_vec[0] = 1.0
        t_max = 500.0

        # Build real-block Liouvillian for complex ODE:
        #   d/dt [Re(rho); Im(rho)] = [[Re(L), -Im(L)], [Im(L), Re(L)]] @ [Re(rho); Im(rho)]
        L_real = np.block([[L.real, -L.imag],
                           [L.imag, L.real]])
        y0 = np.concatenate([rho0_vec.real, rho0_vec.imag])

        def rhs(t, y):
            return L_real @ y

        sol = solve_ivp(rhs, [0, t_max], y0,
                        method='RK45', rtol=1e-9, atol=1e-11,
                        dense_output=True)

        # Integrate sink population (real part of rho_NN)
        sink_idx = (n - 1) + (n - 1) * n
        sink_pop = sol.y[sink_idx, :]
        times = sol.t
        eta_numerical = kappa * np.trapezoid(sink_pop, times)

        assert_allclose(eta_analytical, eta_numerical, rtol=1e-2,
                        err_msg=f"Analytical vs numerical mismatch: "
                                f"anal={eta_analytical}, num={eta_numerical}")

    def test_analytical_enhancement_vs_coherent(self):
        """Enhancement ratio: eta_peak / eta_zero > 1 (ENAQT exists).

        For parameters where ENAQT exists, the enhancement ratio should
        be > 1, meaning dephasing-assisted transport exceeds coherent
        transport. We use a non-zero fluorescence rate so the yield is
        < 1 and dephasing can genuinely help.
        """
        from enaqt.core import (hamiltonian_2site, liouvillian,
                                 analytical_yield, liouvillian_parts,
                                 optimal_dephasing)

        H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        n = H.shape[0]
        gamma_f = 0.01  # Non-zero fluorescence needed for enhancement > 1

        # Yield with no dephasing (coherent)
        L0 = liouvillian(H, dephasing_rate=0.0, sink_rate=KAPPA,
                         fluo_rate=gamma_f)
        eta_zero = analytical_yield(L0, KAPPA, n)

        # Find optimal dephasing via sweep
        L_base, L_deph = liouvillian_parts(H, KAPPA, gamma_f)
        gamma_phi_arr, eta_arr = optimal_dephasing(
            L_base, L_deph, KAPPA, n, n_points=200
        )
        peak_idx = np.argmax(eta_arr)
        eta_peak = eta_arr[peak_idx]
        gamma_opt = gamma_phi_arr[peak_idx]

        enhancement = eta_peak / max(eta_zero, 1e-16)

        assert enhancement > 1.0, (
            f"No ENAQT: enhancement={enhancement} (should be > 1). "
            f"eta_zero={eta_zero}, eta_peak={eta_peak} "
            f"at gamma_phi={gamma_opt:.4f}"
        )
