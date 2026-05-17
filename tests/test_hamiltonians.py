"""
Group A: Hamiltonian Construction Tests
========================================

These tests verify that all Hamiltonian generators produce physically
correct matrices: Hermiticity, correct structure, expected energy
gradients, and reproducibility under fixed random seeds.
"""

import numpy as np
from numpy.testing import assert_allclose
import pytest


class TestHamiltonian2Site:
    """Tests for the 2-site spin-boson Hamiltonian."""

    def test_hamiltonian_2site_hermitian(self):
        """H must be Hermitian: H = H^dagger (H.conj().T == H).

        Physics: All Hamiltonians in quantum mechanics must be Hermitian
        because they represent observable energies. Non-Hermitian H would
        give complex eigenvalues, violating the spectral theorem.
        """
        from enaqt.core import hamiltonian_2site

        # Test a range of physically relevant parameters
        params = [(0.0, 1.0), (1.0, 1.0), (5.0, 1.0), (1.0, 0.5), (0.0, 0.0)]
        for eps, delta in params:
            H = hamiltonian_2site(epsilon=eps, delta=delta)
            assert_allclose(H, H.conj().T, rtol=1e-10,
                            err_msg=f"H not Hermitian for eps={eps}, delta={delta}")

    def test_hamiltonian_2site_values(self):
        """Verify explicit matrix elements for known parameter sets.

        For H = eps/2 * sigma_z + delta/2 * sigma_x:
          H[0,0] = +eps/2, H[1,1] = -eps/2, H[0,1] = H[1,0] = delta/2
        """
        from enaqt.core import hamiltonian_2site

        H = hamiltonian_2site(epsilon=2.0, delta=1.0)
        assert_allclose(H[0, 0], 1.0, rtol=1e-10)   # +eps/2 = +1
        assert_allclose(H[1, 1], -1.0, rtol=1e-10)  # -eps/2 = -1
        assert_allclose(H[0, 1], 0.5, rtol=1e-10)   # +delta/2 = +0.5
        assert_allclose(H[1, 0], 0.5, rtol=1e-10)   # symmetric

        H_res = hamiltonian_2site(epsilon=0.0, delta=2.0)
        assert_allclose(H_res[0, 0], 0.0, rtol=1e-10)
        assert_allclose(H_res[0, 1], 1.0, rtol=1e-10)  # delta/2 = 1


class TestHamiltonianNSite:
    """Tests for N-site chain Hamiltonians."""

    def test_hamiltonian_nsite_structure(self):
        """N-site chain must have correct tridiagonal structure.

        Only nearest-neighbor couplings H[j, j+1] and H[j+1, j] should be
        non-zero off the diagonal. All other off-diagonal elements must be
        exactly zero. This reflects the physical chain topology where each
        site only connects to its immediate neighbors.
        """
        from enaqt.core import hamiltonian_nsite

        for n in [2, 3, 5, 7, 10]:
            H = hamiltonian_nsite(n=n, bias=5.0, coupling=1.0)
            assert H.shape == (n, n), f"Shape mismatch for n={n}"

            # Check Hermiticity
            assert_allclose(H, H.conj().T, rtol=1e-10,
                            err_msg=f"H not Hermitian for n={n}")

            # Check tridiagonal: all elements with |i-j| > 1 must be zero
            for i in range(n):
                for j in range(n):
                    if abs(i - j) > 1:
                        assert H[i, j] == 0, (
                            f"Non-zero off-tridiagonal H[{i},{j}]={H[i,j]} for n={n}"
                        )

            # Check nearest-neighbor couplings are non-zero
            for j in range(n - 1):
                assert H[j, j + 1] != 0, f"Missing coupling at {j},{j+1} for n={n}"

    def test_hamiltonian_funnel_gradient(self):
        """Linear energy gradient: E[i] - E[i-1] = -bias/(n-1) for all i > 0.

        The energy funnel creates a linear gradient from site 1 (highest)
        to site n (lowest). The energy drop between adjacent sites must be
        uniform: Delta_E = bias / (n-1).

        Physics: This mimics the energy landscape in photosynthetic FMO
        complexes where energy decreases monotonically toward the reaction
        center, guiding exciton transport downhill.
        """
        from enaqt.core import hamiltonian_funnel

        for n in [3, 5, 7, 10]:
            for bias in [1.0, 3.0, 5.0, 10.0]:
                H = hamiltonian_funnel(n=n, bias=bias, coupling=1.0)
                diag = np.diag(H).real

                # Expected energy differences (site i+1 minus site i)
                expected_drop = bias / (n - 1)

                for i in range(1, n):
                    actual_drop = diag[i - 1] - diag[i]
                    assert_allclose(actual_drop, expected_drop, rtol=1e-10,
                                    err_msg=f"Gradient mismatch at i={i}, n={n}, bias={bias}")

    def test_hamiltonian_funnel_first_last(self):
        """Energy funnel: first site = +bias/2, last site = -bias/2."""
        from enaqt.core import hamiltonian_funnel

        H = hamiltonian_funnel(n=7, bias=5.0, coupling=1.0)
        diag = np.diag(H).real

        assert_allclose(diag[0], 2.5, rtol=1e-10)   # +bias/2
        assert_allclose(diag[-1], -2.5, rtol=1e-10)  # -bias/2


class TestHamiltonianFMO7:
    """Tests for the FMO-7 Hamiltonian against published reference data."""

    def test_hamiltonian_fmo7_matches_reference(self):
        """FMO-7 Hamiltonian must match Adolphs & Renger 2006 values.

        Reference data (in cm^-1, from Adolphs & Renger, Biophysical J. 91, 2778):
          H[0,0]=12445  H[0,1]=-87.7  H[0,2]=5.5    H[0,3]=-5.9
          H[1,1]=12520  H[1,2]=30.8   H[1,3]=8.2    H[1,4]=0.7
          H[2,2]=12205  H[2,3]=-53.5  H[2,4]=-2.2   H[2,5]=-9.6
          H[3,3]=12335  H[3,4]=-70.7  H[3,5]=-17.0  H[3,6]=-63.0
          H[4,4]=12490  H[4,5]=81.1   H[4,6]=-1.3
          H[5,5]=12640  H[5,6]=39.7
          H[6,6]=12450

        After conversion to Delta units and shifting site 1 to zero energy,
        all 28 unique upper-triangular elements (including diagonal) must
        match the reference implementation.
        """
        from enaqt.core import hamiltonian_fmo7

        H = hamiltonian_fmo7()
        assert H.shape == (7, 7), f"Expected (7,7), got {H.shape}"

        # Reference Hamiltonian in cm^-1 from Adolphs & Renger 2006
        H_ref_cm = np.array([
            [12445, -87.7,   5.5,  -5.9,   6.7, -13.7,  -9.9],
            [-87.7, 12520,  30.8,   8.2,   0.7,  11.8,   4.3],
            [  5.5,  30.8, 12205, -53.5,  -2.2,  -9.6,   6.0],
            [ -5.9,   8.2, -53.5, 12335, -70.7, -17.0, -63.0],
            [  6.7,   0.7,  -2.2, -70.7, 12490,  81.1,  -1.3],
            [-13.7,  11.8,  -9.6, -17.0,  81.1, 12640,  39.7],
            [ -9.9,   4.3,   6.0, -63.0,  -1.3,  39.7, 12450],
        ], dtype=float)

        # Convert to same units: divide by 87.7, shift site 1 to zero
        H_ref = H_ref_cm / 87.7
        H_ref -= np.eye(7) * H_ref[0, 0]

        # Verify all 49 elements (28 unique due to symmetry)
        assert_allclose(H, H_ref, rtol=1e-10,
                        err_msg="FMO-7 Hamiltonian does not match reference data")

        # Also verify Hermiticity
        assert_allclose(H, H.T, rtol=1e-10,
                        err_msg="FMO-7 Hamiltonian not symmetric")


class TestHamiltonianDisordered:
    """Tests for disordered Hamiltonian with Gaussian random site energies."""

    def test_hamiltonian_disordered_reproducible(self):
        """Same rng seed must produce exactly the same Hamiltonian.

        Disorder realizations must be deterministic for reproducible
        ensemble averages. Two calls with the same Generator must yield
        identical H matrices element-for-element.

        Physics: In Anderson localization studies, disorder ensemble
        statistics require many realizations. Reproducibility is essential
        for verifying ensemble-averaged results.
        """
        from enaqt.core import hamiltonian_disordered

        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(42)

        H1 = hamiltonian_disordered(n=5, bias=5.0, coupling=1.0,
                                     disorder_sigma=2.0, rng=rng1)
        H2 = hamiltonian_disordered(n=5, bias=5.0, coupling=1.0,
                                     disorder_sigma=2.0, rng=rng2)

        assert_allclose(H1, H2, rtol=1e-10,
                        err_msg="Same RNG seed produced different Hamiltonians")

    def test_hamiltonian_disordered_different_seeds(self):
        """Different rng seeds must produce different diagonal elements.

        The off-diagonal (coupling) should remain the same, but site
        energies should differ between different random seeds.
        """
        from enaqt.core import hamiltonian_disordered

        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(12345)

        H1 = hamiltonian_disordered(n=5, bias=5.0, coupling=1.0,
                                     disorder_sigma=2.0, rng=rng1)
        H2 = hamiltonian_disordered(n=5, bias=5.0, coupling=1.0,
                                     disorder_sigma=2.0, rng=rng2)

        # Off-diagonal (couplings) should be identical
        off_diag_1 = H1 - np.diag(np.diag(H1))
        off_diag_2 = H2 - np.diag(np.diag(H2))
        assert_allclose(off_diag_1, off_diag_2, rtol=1e-10)

        # Diagonal (site energies) should differ
        diag_1 = np.diag(H1)
        diag_2 = np.diag(H2)
        assert not np.allclose(diag_1, diag_2, rtol=1e-10), (
            "Different seeds produced identical disorder realizations"
        )

    def test_hamiltonian_disordered_hermitian(self):
        """Disordered Hamiltonian must remain Hermitian for any realization."""
        from enaqt.core import hamiltonian_disordered

        for seed in [1, 42, 99, 256]:
            rng = np.random.default_rng(seed)
            H = hamiltonian_disordered(n=7, bias=5.0, coupling=1.0,
                                        disorder_sigma=1.5, rng=rng)
            assert_allclose(H, H.conj().T, rtol=1e-10,
                            err_msg=f"Disordered H not Hermitian for seed={seed}")

    def test_hamiltonian_disordered_tridiagonal(self):
        """Disordered Hamiltonian must remain tridiagonal (no long-range couplings)."""
        from enaqt.core import hamiltonian_disordered

        rng = np.random.default_rng(42)
        H = hamiltonian_disordered(n=7, bias=5.0, coupling=1.0,
                                    disorder_sigma=2.0, rng=rng)
        for i in range(7):
            for j in range(7):
                if abs(i - j) > 1:
                    assert H[i, j] == 0, (
                        f"Non-tridiagonal element H[{i},{j}]={H[i,j]}"
                    )
