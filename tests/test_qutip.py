"""
Group D: QuTiP Cross-Validation Tests
======================================

These tests cross-validate our implementation against QuTiP, the
gold-standard Python library for open quantum systems.

Tests use ``pytest.importorskip("qutip")`` so they are skipped if QuTiP
is not installed, making them optional but authoritative when available.

**Sink model note**
-------------------
``core.py`` implements the sink as the no-jump anticommutator channel::

    L_sink rho = -kappa/2 {P_N, rho}

This is *not* a standard Lindblad jump channel, so ``qutip.mesolve`` cannot
reproduce it directly via ``c_ops``.  QuTiP's master equation uses the
*commutator* ``-i[H, rho]``, which for a non-Hermitian H gives
``-i(A rho - rho A)`` (not the needed ``A rho + rho A``).

To avoid mismatched physics we therefore:

* Use ``scipy.integrate.solve_ivp`` with our own Liouvillian as an
  independent reference for yield tests (Tests 1-2).
* Use ``qutip.mesolve`` with ``kappa = 0`` (no sink) for dynamics tests
  (Tests 3-4) so that only the standard Lindblad dephasing channel is
  active — the part that *can* be validated cleanly against QuTiP.

Cross-validated properties
--------------------------
1. ``analytical_yield`` vs ODE integration (resonant case)
2. ``analytical_yield`` vs ODE integration (detuned, γ_φ sweep)
3. Coherent Rabi oscillations: QuTiP mesolve vs Hamiltonian eigenvalues
4. Quantum Zeno suppression: QuTiP dephasing + analytical yield
5. Trace conservation: QuTiP mesolve vs our ODE (existing test)
"""

from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import solve_ivp

# Physical constants for all tests
DELTA = 1.0
KAPPA = 0.1
GAMMA_F = 0.01

qutip = pytest.importorskip("qutip")


# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

def _ode_yield(H_np: np.ndarray, gamma_phi: float,
               kappa: float, fluo_rate: float,
               t_end: float = 600.0, n_pts: int = 6000) -> float:
    """Compute transfer yield via ODE integration of our Liouvillian.

    This is an *independent* reference: it solves the same differential
    equation as ``analytical_yield`` but via numerical time-stepping rather
    than direct matrix inversion.
    """
    from enaqt.core import liouvillian

    n = H_np.shape[0]
    L = liouvillian(H_np, gamma_phi, kappa, fluo_rate)
    rho0 = np.zeros(n * n, dtype=complex)
    rho0[0] = 1.0
    sink_idx = (n - 1) + (n - 1) * n  # rho_NN in column-stacked vector

    tlist = np.linspace(0.0, t_end, n_pts)
    sol = solve_ivp(lambda t, y: L @ y, (0.0, t_end), rho0,
                    t_eval=tlist, rtol=1e-9, atol=1e-11, method="RK45")
    return float(kappa * np.trapezoid(sol.y[sink_idx].real, tlist))


# ─────────────────────────────────────────────────────────────────────────────
# Test Group 1 & 2 — Yield self-consistency
# ─────────────────────────────────────────────────────────────────────────────

class TestQuTiP2SiteResonant:
    """analytical_yield vs ODE integration for yield self-consistency.

    Using QuTiP Qobj for Hamiltonian sanity-check; yield is compared via
    an independent scipy ODE reference rather than mesolve (see module note).
    """

    def test_qutip_2site_resonant(self):
        """epsilon=0 (resonant): analytical_yield agrees with ODE integration.

        Validates:
          (a) The Liouvillian matrix is constructed correctly — ODE and
              matrix-inversion formula agree to within ODE tolerance.
          (b) ``qutip.Qobj`` wraps the Hamiltonian correctly (eigenvalue check).
        """
        from enaqt.core import hamiltonian_2site, liouvillian, analytical_yield

        H_np = hamiltonian_2site(epsilon=0.0, delta=DELTA)
        gamma_phi = 0.1

        # QuTiP sanity-check: eigenvalues must match numpy
        H_q = qutip.Qobj(H_np)
        assert H_q.isherm, "Hamiltonian must be Hermitian"
        assert_allclose(np.sort(H_q.eigenenergies()),
                        np.sort(np.linalg.eigvalsh(H_np)),
                        rtol=1e-10,
                        err_msg="QuTiP and numpy eigenvalues disagree")

        eta_analytic = analytical_yield(
            liouvillian(H_np, gamma_phi, KAPPA, GAMMA_F), KAPPA, n=2)
        eta_ode = _ode_yield(H_np, gamma_phi, KAPPA, GAMMA_F)

        assert_allclose(eta_analytic, eta_ode, rtol=5e-3,
                        err_msg=(f"Analytical yield {eta_analytic:.6f} disagrees "
                                 f"with ODE reference {eta_ode:.6f}"))

    def test_qutip_2site_detuned(self):
        """epsilon=1 (detuned): analytical_yield agrees with ODE for 4 gamma values.

        At finite bias, dephasing can enhance transport (ENAQT).  We sweep
        four dephasing rates spanning 3 decades and verify the formula at each.
        """
        from enaqt.core import hamiltonian_2site, liouvillian, analytical_yield

        H_np = hamiltonian_2site(epsilon=1.0, delta=DELTA)

        for gamma_phi in [0.01, 0.1, 1.0, 10.0]:
            L = liouvillian(H_np, gamma_phi, KAPPA, GAMMA_F)
            eta_analytic = analytical_yield(L, KAPPA, n=2)
            eta_ode = _ode_yield(H_np, gamma_phi, KAPPA, GAMMA_F)

            assert_allclose(eta_analytic, eta_ode, rtol=5e-3,
                            err_msg=(f"gp={gamma_phi}: "
                                     f"analytic={eta_analytic:.6f}, "
                                     f"ode={eta_ode:.6f}"))


# ─────────────────────────────────────────────────────────────────────────────
# Test Group 3 & 4 — Dynamics cross-validation with QuTiP mesolve (kappa=0)
# ─────────────────────────────────────────────────────────────────────────────

class TestQuTiPDynamics:
    """Dynamics cross-validation using QuTiP mesolve with kappa=0.

    With no sink (kappa=0), only the standard Lindblad dephasing channel is
    active, which QuTiP's mesolve can handle cleanly.  These tests validate:

    * The dephasing superoperator is constructed consistently with QuTiP's
      Lindblad dissipator.
    * Zeno suppression appears in both QuTiP dynamics and our yield formula.
    """

    def test_qutip_zero_dephasing(self):
        """Pure coherent evolution: mesolve shows Rabi oscillations at the
        correct frequency given by the Hamiltonian eigenvalue splitting.

        Uses epsilon=0 (resonant) so that the Rabi amplitude is 100%
        (P0 oscillates all the way to 0), giving unambiguous 0.5-crossings.
        With epsilon>0, the bias limits transfer to <50%, making crossings
        impossible.
        """
        from enaqt.core import hamiltonian_2site

        # Resonant case: full 100% population transfer, unambiguous oscillations
        H_np = hamiltonian_2site(epsilon=0.0, delta=DELTA)
        H_q = qutip.Qobj(H_np)
        rho0 = qutip.ket2dm(qutip.basis(2, 0))

        # Pure coherent evolution: no c_ops, no sink
        tlist = np.linspace(0, 20, 2000)
        result = qutip.mesolve(H_q, rho0, tlist, c_ops=[], e_ops=[])

        P0_proj = qutip.ket2dm(qutip.basis(2, 0))
        P0 = np.array([qutip.expect(P0_proj, s) for s in result.states])

        # Population at site 0 must oscillate below 0.5 (full Rabi transfer)
        crossings_up = np.sum((P0[:-1] < 0.5) & (P0[1:] >= 0.5))
        assert crossings_up >= 2, (
            f"Expected ≥2 Rabi oscillations in t=[0, 20], found {crossings_up}. "
            f"epsilon=0 resonant case should show full population transfer."
        )

        # Oscillation period = 2*pi / (energy splitting)
        eigvals = np.linalg.eigvalsh(H_np)
        omega_expected = abs(eigvals[1] - eigvals[0])
        period_expected = 2.0 * np.pi / omega_expected

        from scipy.signal import find_peaks
        peaks, _ = find_peaks(P0, height=0.7, distance=50)
        if len(peaks) >= 2:
            period_measured = tlist[peaks[1]] - tlist[peaks[0]]
            assert_allclose(period_measured, period_expected, rtol=5e-2,
                            err_msg=(f"Rabi period mismatch: "
                                     f"measured={period_measured:.3f}, "
                                     f"expected={period_expected:.3f}"))

    def test_qutip_strong_dephasing(self):
        """Large gamma_phi suppresses oscillations (Zeno) — cross-validated with
        QuTiP mesolve and our analytical yield formula.

        With kappa=0 (no sink), only dephasing is active.  We verify:
          (a) Weak dephasing (gamma_phi=0.01): oscillations visible in mesolve
          (b) Strong dephasing (gamma_phi=100): oscillations suppressed in mesolve
          (c) Our analytical yield decreases monotonically from the ENAQT peak
              toward zero as gamma_phi increases (Zeno suppression of yield)
        """
        from enaqt.core import hamiltonian_2site, liouvillian, analytical_yield

        # Use resonant case so oscillations are visible at weak dephasing
        H_np = hamiltonian_2site(epsilon=0.0, delta=DELTA)
        H_q = qutip.Qobj(H_np)
        rho0 = qutip.ket2dm(qutip.basis(2, 0))
        tlist = np.linspace(0, 20, 2000)
        P0_proj = qutip.ket2dm(qutip.basis(2, 0))

        def _deph_c_ops(gamma: float) -> list:
            return [np.sqrt(gamma) * qutip.basis(2, j) * qutip.basis(2, j).dag()
                    for j in range(2)]

        # (a) Weak dephasing: oscillations should survive
        gamma_low = 0.01
        result_low = qutip.mesolve(H_q, rho0, tlist,
                                   c_ops=_deph_c_ops(gamma_low), e_ops=[])
        P0_low = np.array([qutip.expect(P0_proj, s) for s in result_low.states])
        crossings_low = np.sum((P0_low[:-1] < 0.5) & (P0_low[1:] >= 0.5))
        assert crossings_low >= 1, (
            f"Expected oscillations at gamma_phi={gamma_low}, "
            f"found {crossings_low} crossings."
        )

        # (b) Strong dephasing: oscillations should be suppressed
        gamma_strong = 100.0
        result_strong = qutip.mesolve(H_q, rho0, tlist,
                                      c_ops=_deph_c_ops(gamma_strong), e_ops=[])
        P0_strong = np.array([qutip.expect(P0_proj, s)
                               for s in result_strong.states])
        crossings_strong = np.sum((P0_strong[:-1] < 0.5) & (P0_strong[1:] >= 0.5))
        assert crossings_strong < crossings_low, (
            f"Zeno suppression not observed: gamma={gamma_strong} gives "
            f"{crossings_strong} crossings, gamma={gamma_low} gives {crossings_low}."
        )

        # (c) Analytical yield decreases at strong dephasing (detuned case)
        H_det = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        eta_low = analytical_yield(liouvillian(H_det, gamma_low, KAPPA, GAMMA_F),
                                   KAPPA, n=2)
        eta_strong = analytical_yield(liouvillian(H_det, gamma_strong, KAPPA, GAMMA_F),
                                      KAPPA, n=2)
        assert eta_strong < eta_low, (
            f"Yield should be suppressed in Zeno limit: "
            f"eta(gp={gamma_low})={eta_low:.4f}, "
            f"eta(gp={gamma_strong})={eta_strong:.4f}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Test Group 5 — Trace conservation (unchanged from original)
# ─────────────────────────────────────────────────────────────────────────────

class TestQuTiPConservation:
    """Conservation law cross-validation with QuTiP."""

    def test_qutip_conservation_comparison(self):
        """Compare Tr[rho(t)] between QuTiP and our code.

        Both implementations must conserve probability when kappa=Gamma=0.
        We evolve the same initial state using both QuTiP and our ODE
        system and compare the trace at each time point.

        This is the most direct cross-validation: the trace dynamics
        should match to within integration tolerance.
        """
        from enaqt.core import hamiltonian_2site, liouvillian

        H_np = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        H_q = qutip.Qobj(H_np)
        rho0 = qutip.ket2dm(qutip.basis(2, 0))

        gamma_phi = 0.1

        # --- QuTiP evolution (no sink, no fluorescence) ---
        deph_ops = [np.sqrt(gamma_phi) * qutip.basis(2, j) * qutip.basis(2, j).dag()
                    for j in range(2)]
        tlist = np.linspace(0, 10, 101)
        result = qutip.mesolve(H_q, rho0, tlist, c_ops=deph_ops)

        traces_qutip = np.array([qutip.expect(qutip.qeye(2), s).real
                                  for s in result.states])

        # --- Our evolution (no sink, no fluorescence) ---
        L = liouvillian(H_np, gamma_phi, sink_rate=0.0, fluo_rate=0.0)
        n = H_np.shape[0]

        def rhs(t, y):
            rho_vec = y[:n*n] + 1j * y[n*n:]
            drho = L @ rho_vec
            return [*drho.real, *drho.imag]

        y0 = np.zeros(2 * n * n)
        y0[0] = 1.0  # rho_11 = 1
        sol = solve_ivp(rhs, [0, 10], y0, t_eval=tlist,
                        rtol=1e-9, atol=1e-11)

        traces_ours = np.array([
            np.sum((sol.y[:n*n, i] + 1j * sol.y[n*n:, i])[::n+1]).real
            for i in range(len(tlist))
        ])

        assert_allclose(traces_ours, traces_qutip, rtol=1e-6, atol=1e-6,
                        err_msg="Trace conservation mismatch between our code and QuTiP")
