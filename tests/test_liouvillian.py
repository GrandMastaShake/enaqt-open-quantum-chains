"""
Group B: Liouvillian Superoperator Tests
=========================================

These tests verify the Lindblad Liouvillian L that generates the open
quantum system dynamics: d(rho_vec)/dt = L @ rho_vec.

Physical properties tested:
  - Correct N^2 x N^2 dimensions
  - Linearity in the dephasing rate gamma_phi
  - Additive decomposition: L = L_base + gamma_phi * L_deph
  - All eigenvalues have non-positive real part (dissipative dynamics)
  - Known analytical limits for simple Hamiltonians
"""

import numpy as np
from numpy.testing import assert_allclose
import pytest


# Physical constants
DELTA = 1.0
KAPPA = 0.1
GAMMA_F = 0.01


class TestLiouvillianSize:
    """Tests for Liouvillian matrix dimensions."""

    def test_liouvillian_size(self):
        """L must be N^2 x N^2 for an N-site Hamiltonian.

        The Liouvillian acts on the vectorized density matrix via the
        column-stack isomorphism: vec(rho) is an N^2-dimensional vector,
        so L must be N^2 x N^2.

        This follows from the vec() isomorphism: vec(A X B) = (B^T kron A) vec(X).
        """
        from enaqt.core import hamiltonian_nsite, liouvillian

        for n in [2, 3, 4, 5, 7]:
            H = hamiltonian_nsite(n=n, bias=5.0, coupling=DELTA)
            L = liouvillian(H, dephasing_rate=0.1, sink_rate=KAPPA,
                            fluo_rate=GAMMA_F)
            expected_size = n * n
            assert L.shape == (expected_size, expected_size), (
                f"Expected shape ({expected_size},{expected_size}), got {L.shape} for n={n}"
            )


class TestLiouvillianLinearity:
    """Tests for linearity of L in the dephasing rate."""

    def test_liouvillian_linear_in_dephasing(self):
        """L(gamma_phi) = L(0) + gamma_phi * dL/dgamma (linear in gamma_phi).

        The dephasing Lindblad term is proportional to gamma_phi, so the
        full Liouvillian must be an affine function: L(gamma) = L(0) + gamma * L_deph.

        We verify this by checking that L evaluated at three different gamma
        values satisfies the linear relation exactly.
        """
        from enaqt.core import hamiltonian_2site, liouvillian

        H = hamiltonian_2site(epsilon=1.0, delta=DELTA)

        gamma_values = [0.0, 0.1, 0.5, 1.0, 2.0, 10.0]
        L_zero = liouvillian(H, dephasing_rate=0.0,
                             sink_rate=KAPPA, fluo_rate=GAMMA_F)
        L_unit = liouvillian(H, dephasing_rate=1.0,
                             sink_rate=KAPPA, fluo_rate=GAMMA_F)

        # L_deph = L(1) - L(0) (the dephasing-only contribution at unit rate)
        L_deph = L_unit - L_zero

        for gp in gamma_values:
            L_gp = liouvillian(H, dephasing_rate=gp,
                               sink_rate=KAPPA, fluo_rate=GAMMA_F)
            L_expected = L_zero + gp * L_deph
            assert_allclose(L_gp, L_expected, rtol=1e-10,
                            err_msg=f"L not linear in gamma_phi at gp={gp}")

    def test_liouvillian_parts_reconstruct(self):
        """L_base + gamma_phi * L_deph must equal the full L.

        liouvillian_parts returns (L_base, L_deph) such that:
          L_full = L_base + gamma_phi * L_deph

        This decomposition is the key optimization that enables 2-14x
        speedup for gamma_phi sweeps, since L_base and L_deph are
        precomputed once and only a vectorized addition is needed for
        each gamma_phi value.
        """
        from enaqt.core import hamiltonian_nsite, liouvillian, liouvillian_parts

        for n in [2, 3, 5]:
            H = hamiltonian_nsite(n=n, bias=3.0, coupling=DELTA)
            L_base, L_deph = liouvillian_parts(H, sink_rate=KAPPA,
                                               fluo_rate=GAMMA_F)

            for gp in [0.0, 0.01, 0.1, 1.0, 10.0]:
                L_full = liouvillian(H, dephasing_rate=gp, sink_rate=KAPPA,
                                     fluo_rate=GAMMA_F)
                L_reconstructed = L_base + gp * L_deph
                assert_allclose(L_full, L_reconstructed, rtol=1e-10,
                                err_msg=f"Reconstruction failed for n={n}, gp={gp}")


class TestLiouvillianEigenvalues:
    """Tests for Liouvillian spectral properties."""

    def test_liouvillian_eigenvalues_nonpositive(self):
        """All eigenvalues must have non-positive real part: Re(lambda) <= 0.

        This is the defining property of a dissipative (completely positive
        trace-preserving) quantum dynamical semigroup. The Liouvillian is
        the generator of a CPTP map, so its spectrum must lie in the left
        half-plane.

        A positive real part would imply exponential growth of some state,
        violating complete positivity. A zero eigenvalue corresponds to the
        steady state (attractor).

        Reference: Gorini, Kossakowski, Sudarshan & Lindblad (1976).
        """
        from enaqt.core import hamiltonian_nsite, liouvillian

        for n in [2, 3, 4, 5]:
            H = hamiltonian_nsite(n=n, bias=3.0, coupling=DELTA)
            for gp in [0.0, 0.1, 1.0, 10.0]:
                L = liouvillian(H, dephasing_rate=gp, sink_rate=KAPPA,
                                fluo_rate=GAMMA_F)
                eigs = np.linalg.eigvals(L)
                max_real = np.max(eigs.real)

                assert max_real <= 1e-10, (
                    f"Found eigenvalue with positive real part: {max_real} "
                    f"for n={n}, gp={gp}. Spectrum violates dissipativity."
                )

    def test_liouvillian_zero_eigenvalue_exists(self):
        """At least one eigenvalue has Re(lambda) approx 0 (steady state exists).

        For a CPTP map, the Liouvillian always has at least one zero
        eigenvalue corresponding to the invariant steady state. With a
        sink, this is slightly shifted negative, but in the no-sink limit
        it must be exactly zero.
        """
        from enaqt.core import hamiltonian_2site, liouvillian

        H = hamiltonian_2site(epsilon=1.0, delta=DELTA)
        L = liouvillian(H, dephasing_rate=0.1, sink_rate=0.0,
                        fluo_rate=0.0)  # Must be 0 for zero eigenvalue
        eigs = np.linalg.eigvals(L)
        max_real = np.max(eigs.real)
        assert abs(max_real) < 1e-8, (
            f"No zero eigenvalue found (max Re(lambda) = {max_real}). "
            f"Liouvillian without sink should have at least one zero eigenvalue."
        )


class TestLiouvillianFlatChain:
    """Tests against known analytical limits for flat chains."""

    def test_liouvillian_flat_chain_limit(self):
        """Flat chain (all energies equal) at gamma=0, kappa=0: unitary evolution.

        For a flat chain with no dephasing and no sink, the dynamics is
        purely Hamiltonian: L = -i(I kron H - H^T kron I). The eigenvalues
        of L are then purely imaginary: lambda = -i(E_i - E_j) for all i,j.

        This means all eigenvalues have exactly zero real part when there
        is no dissipation, corresponding to purely oscillatory (unitary)
        dynamics with no decay.

        We verify that Re(lambda) = 0 for all eigenvalues of a flat chain
        Liouvillian with gamma_phi = kappa = Gamma = 0.
        """
        from enaqt.core import hamiltonian_flat, liouvillian

        for n in [2, 3, 4, 5]:
            H = hamiltonian_flat(n=n, coupling=DELTA)
            L = liouvillian(H, dephasing_rate=0.0, sink_rate=0.0,
                            fluo_rate=0.0)
            eigs = np.linalg.eigvals(L)

            # All eigenvalues must be purely imaginary (zero real part)
            assert_allclose(eigs.real, np.zeros_like(eigs.real), atol=1e-10,
                            err_msg=f"Flat chain n={n} has non-zero real eigenvalue parts "
                                    f"when gamma=kappa=Gamma=0")

    def test_liouvillian_flat_chain_with_sink(self):
        """Flat chain with sink and no dephasing: coherent transport.

        Adding a sink introduces negative real parts to eigenvalues,
        corresponding to irreversible population loss. The trace of the
        Liouvillian equals the sum of eigenvalues and must be negative
        (population is lost to the sink).
        """
        from enaqt.core import hamiltonian_flat, liouvillian

        H = hamiltonian_flat(n=3, coupling=DELTA)
        L = liouvillian(H, dephasing_rate=0.0, sink_rate=KAPPA,
                        fluo_rate=0.0)
        eigs = np.linalg.eigvals(L)

        # Trace = sum of eigenvalues
        trace = np.trace(L)
        eig_sum = np.sum(eigs)
        assert_allclose(trace, eig_sum, rtol=1e-10)

        # With sink, trace must be negative (population lost)
        assert trace.real < 0, (
            f"Trace of L with sink must be negative, got {trace}"
        )

    def test_liouvillian_consistent_across_topologies(self):
        """Same physical parameters on different topologies give different L.

        The Liouvillian should reflect the Hamiltonian structure, so
        different topologies (funnel, flat, disordered) should produce
        different Liouvillians even with identical rates.
        """
        from enaqt.core import (hamiltonian_funnel, hamiltonian_flat,
                                 hamiltonian_disordered, liouvillian)

        rng = np.random.default_rng(42)
        H_funnel = hamiltonian_funnel(n=5, bias=5.0, coupling=DELTA)
        H_flat = hamiltonian_flat(n=5, coupling=DELTA)
        H_dis = hamiltonian_disordered(n=5, bias=5.0, coupling=DELTA,
                                        disorder_sigma=2.0, rng=rng)

        L_funnel = liouvillian(H_funnel, dephasing_rate=0.1,
                               sink_rate=KAPPA, fluo_rate=GAMMA_F)
        L_flat = liouvillian(H_flat, dephasing_rate=0.1,
                             sink_rate=KAPPA, fluo_rate=GAMMA_F)
        L_dis = liouvillian(H_dis, dephasing_rate=0.1,
                            sink_rate=KAPPA, fluo_rate=GAMMA_F)

        # Funnel != flat
        assert not np.allclose(L_funnel, L_flat, atol=1e-10), (
            "Funnel and flat topologies produced identical Liouvillians"
        )

        # Funnel != disordered
        assert not np.allclose(L_funnel, L_dis, atol=1e-10), (
            "Funnel and disordered topologies produced identical Liouvillians"
        )
