"""
Smoke tests for the analysis entry points (tiny sizes, temporary output dirs).

These exercise the full code paths of :mod:`enaqt.nsite`,
:mod:`enaqt.disorder`, :mod:`enaqt.spinboson` and :mod:`enaqt.cli` --
including plotting and JSON output -- and check a few physical invariants
of the results.  Each test runs in a few seconds at most.
"""

import json

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from enaqt import cli, disorder, nsite, spinboson
from enaqt.core import (
    hamiltonian_disordered,
    hamiltonian_funnel,
    liouvillian_parts,
    optimal_dephasing,
)

KAPPA = 0.1
GAMMA = 0.01

# ---------------------------------------------------------------------------
#  nsite
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def scaling_result(tmp_path_factory):
    out = tmp_path_factory.mktemp("scaling")
    return out, nsite.run_scaling(n_max=7, out_dir=out)


class TestNSite:
    def test_run_scaling_outputs(self, scaling_result):
        out, res = scaling_result
        assert (out / "enaqt_nsite_scaling.json").exists()
        assert (out / "enaqt_nsite_scaling.png").exists()
        assert sorted(res["funnel_results"], key=int) == ["2", "3", "4", "5", "6", "7"]
        assert res["parameters"]["total_bias"] == 5.0
        assert res["fit"]["preferred_model"] in ("linear", "saturation")
        assert len(res["fit"]["linear_cooks_d"]) == 6

    def test_run_scaling_uses_bias_5_funnel(self, scaling_result):
        """Regression: the funnel must use bias=5, coupling=1 (paper value)."""
        _, res = scaling_result
        assert res["funnel_results"]["7"]["enhancement"] == pytest.approx(
            22.8, rel=1e-2
        )

    def test_run_scaling_flat_and_disordered(self, scaling_result):
        _, res = scaling_result
        for r in res["flat_results"].values():
            assert r["enhancement"] == pytest.approx(1.0)
            assert not r["interior_peak"]
        # Disordered chain = zero-mean Anderson model, seed 42 per N
        H = hamiltonian_disordered(7, 0.0, 1.0, 2.0, np.random.default_rng(42))
        gp, eta = optimal_dephasing(*liouvillian_parts(H, KAPPA, GAMMA), KAPPA, 7)
        assert res["disordered_results"]["7"]["enhancement"] == pytest.approx(
            float(eta.max() / eta[0])
        )

    def test_run_scaling_rejects_tiny_n_max(self, tmp_path):
        with pytest.raises(ValueError):
            nsite.run_scaling(n_max=1, out_dir=tmp_path)

    def test_run_fmo_benchmark(self, tmp_path):
        res = nsite.run_fmo_benchmark(out_dir=tmp_path)
        assert res["fmo7"]["enhancement"] == pytest.approx(32.1, rel=1e-2)
        assert res["funnel_n7"]["enhancement"] == pytest.approx(22.8, rel=1e-2)
        assert len(res["caveats"]) == 3
        assert json.loads((tmp_path / "enaqt_fmo_benchmark.json").read_text())

    def test_time_trajectory_conserves_probability(self):
        H = hamiltonian_funnel(3, 5.0, 1.0)
        t, pops, eta_t = nsite.time_trajectory(
            H, 1.0, KAPPA, 0.0, t_span=(0.0, 20.0), n_pts=50
        )
        assert pops.shape == (3, 50)
        # Without fluorescence, population lost from the chain goes to the sink
        np.testing.assert_allclose(pops.sum(axis=0) + eta_t, 1.0, atol=1e-6)
        assert np.all(np.diff(eta_t) >= -1e-10)

    def test_fit_models_linear_data(self):
        x = np.arange(2, 12, dtype=float)
        fit = nsite.fit_models(x, 3.0 * x + 1.0, n_bootstrap=20)
        assert fit["linear_params"]["a"] == pytest.approx(3.0, rel=1e-6)
        assert fit["linear_params"]["b"] == pytest.approx(1.0, rel=1e-6)
        assert fit["linear_r2"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
#  disorder
# ---------------------------------------------------------------------------


class TestDisorder:
    def test_run_ensemble_small(self, tmp_path):
        res = disorder.run_ensemble(
            n=3, n_seeds=8, out_dir=tmp_path, n_jobs=1, gp_arr=np.logspace(-3, 3, 30)
        )
        assert res["summary"]["n_seeds"] == 8
        assert len(res["per_seed"]["enhancements"]) == 8
        assert all(e >= 1.0 - 1e-12 for e in res["per_seed"]["enhancements"])
        assert 0.0 <= res["wilson_ci"]["lower"] <= res["wilson_ci"]["fraction"]
        assert res["wilson_ci"]["fraction"] <= res["wilson_ci"]["upper"] <= 1.0
        assert (tmp_path / "disorder_ensemble_n3.json").exists()

    def test_run_optimization_small(self, tmp_path):
        res = disorder.run_optimization(
            3, out_dir=tmp_path, workers=1, maxiter=3, popsize=4
        )
        assert res["n"] == 3
        assert set(res["by_bound_width"]) == {"2", "5", "10", "20", "50"}
        assert res["optimal"]["ratio"] >= 1.0
        assert (tmp_path / "disorder_optimize_n3.json").exists()

    def test_compute_enhancement_floor(self):
        assert disorder.compute_enhancement(0.5, 0.25) == pytest.approx(2.0)
        with pytest.warns(RuntimeWarning):
            assert np.isfinite(disorder.compute_enhancement(0.5, 0.0))

    def test_wilson_ci(self):
        lo, hi = disorder.wilson_ci(50, 100)
        assert lo < 0.5 < hi
        assert disorder.wilson_ci(0, 10)[0] == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
#  spinboson
# ---------------------------------------------------------------------------


def _write_bloch_dataset(data_dir):
    """Write QD3SET-1-style .npy files generated by the Bloch model itself.

    Columns: time, rho11, rho12, rho21, rho22 (complex), t in [0, 20].
    gamma_phi = 2 lambda / (beta gamma_c) with gamma_c = 1.
    """
    t = np.linspace(0.0, 20.0, 401)
    for eps in (0.0, 1.0):
        for beta in (1.0, 2.0):
            for lam in np.logspace(-2.5, 1.5, 10):
                gp = 2.0 * lam / beta
                A = spinboson._bloch_liouvillian_4x4(eps, 1.0, gp, 0.0, 0.0)
                sol = solve_ivp(
                    lambda _t, y: A @ y,
                    (0.0, 20.0),
                    [1.0, 0.0, 0.0, 0.0],
                    t_eval=t,
                    rtol=1e-10,
                    atol=1e-12,
                )
                r12 = sol.y[2] + 1j * sol.y[3]
                arr = np.stack(
                    [t + 0j, sol.y[0] + 0j, r12, np.conj(r12), sol.y[1] + 0j], axis=1
                )
                name = (
                    f"2_epsilon-{eps}_Delta-1.0_lambda-{lam:.6f}"
                    f"_gamma-1.0_beta-{beta}.npy"
                )
                np.save(data_dir / name, arr)


@pytest.fixture(scope="module")
def validation_result(tmp_path_factory):
    data = tmp_path_factory.mktemp("qd3set1")
    out = tmp_path_factory.mktemp("validation")
    _write_bloch_dataset(data)
    return out, spinboson.run_validation(data, out, n_points=40)


class TestSpinBoson:
    def test_run_validation_self_consistent(self, validation_result):
        """Bloch-generated 'HEOM' data must match the Bloch model (T ~ 0)."""
        out, res = validation_result
        assert res["dataset"]["total_trajectories"] == 40
        assert res["symmetric_case"]["trace_distance"]["max"] < 1e-4
        assert res["asymmetric_case"]["trace_distance"]["max"] < 1e-4
        assert (out / "spinboson_validation_results.json").exists()

    def test_run_sink_analysis(self, tmp_path, validation_result):
        # Reuse the validation JSON so the HEOM-overlay code path also runs.
        out, _ = validation_result
        res = spinboson.run_sink_analysis(out, n_points=40)
        key = res["key_results"]
        assert key["Bloch_eps5_kap01_enhancement"] == pytest.approx(7.2, rel=5e-2)
        assert key["Bloch_eps1_kap01_enhancement"] == pytest.approx(1.26, rel=5e-2)
        assert (out / "spinboson_sink_results.json").exists()

    def test_run_sink_analysis_without_heom_data(self, tmp_path):
        """Regression: must not crash when no validation JSON is present."""
        res = spinboson.run_sink_analysis(tmp_path, n_points=20)
        assert len(res["plots"]) == 2

    def test_load_all_requires_files(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            spinboson._load_all(tmp_path)

    def test_parse_filename(self):
        p = spinboson._parse_filename(
            "2_epsilon-1.0_Delta-1.0_lambda-0.5_gamma-2.0_beta-0.25.npy"
        )
        assert (p["epsilon"], p["delta"], p["lam"]) == (1.0, 1.0, 0.5)
        assert (p["gamma_c"], p["beta"]) == (2.0, 0.25)

    def test_trace_distance(self):
        rho = np.array([[1.0, 0.0], [0.0, 0.0]])
        sigma = np.array([[0.0, 0.0], [0.0, 1.0]])
        assert spinboson._trace_distance(rho, rho) == pytest.approx(0.0)
        assert spinboson._trace_distance(rho, sigma) == pytest.approx(1.0)

    @pytest.mark.parametrize(
        "gp, regime", [(0.001, "localization"), (1.0, "optimal"), (100.0, "zeno")]
    )
    def test_regime_classification(self, gp, regime):
        assert spinboson.regime_classification(1.0, 1.0, 0.1, gp) == regime


# ---------------------------------------------------------------------------
#  cli
# ---------------------------------------------------------------------------


class TestCLI:
    def test_no_command_prints_help(self, capsys):
        assert cli.main([]) == 0
        assert "usage" in capsys.readouterr().out.lower()

    def test_fmo(self, tmp_path, capsys):
        assert cli.main(["fmo", "--out-dir", str(tmp_path)]) == 0
        assert "Enhancement : 32.1x" in capsys.readouterr().out

    def test_disorder(self, tmp_path):
        argv = ["disorder", "--n", "3", "--n-seeds", "4", "--n-jobs", "1"]
        assert cli.main(argv + ["--out-dir", str(tmp_path)]) == 0
        assert (tmp_path / "disorder_ensemble_n3.json").exists()

    def test_validate_missing_dir(self, tmp_path):
        missing = tmp_path / "nope"
        argv = ["validate", "--qd3set1-dir", str(missing), "--out-dir", str(tmp_path)]
        assert cli.main(argv) == 1

    def test_optimize_bad_bounds(self, tmp_path):
        argv = ["optimize", "--bounds-min", "5", "--bounds-max", "1"]
        assert cli.main(argv + ["--out-dir", str(tmp_path)]) == 1

    def test_invalid_n_jobs(self):
        with pytest.raises(SystemExit):
            cli.main(["disorder", "--n-jobs", "0"])

    def test_out_dir_env(self, tmp_path, monkeypatch):
        monkeypatch.setenv("ENAQT_OUT_DIR", str(tmp_path / "env"))
        assert cli.get_out_dir() == (tmp_path / "env").resolve()
        assert cli.get_out_dir(str(tmp_path / "arg")) == (tmp_path / "arg").resolve()

    def test_summaries_match_result_schemas(self, scaling_result, capsys):
        """Regression: printers used keys the analyses never produce and
        crashed formatting the 'N/A' placeholder with a numeric spec."""
        _, res = scaling_result
        cli._print_scaling_summary(res)
        cli._print_validate_summary({})
        cli._print_disorder_summary({})
        cli._print_optimize_summary({"optimal": {"ratio": None}})
        cli._print_fmo_summary({})
        out = capsys.readouterr().out
        assert "N= 7: enhancement=22.8x" in out
        assert "N/A" in out

    def test_fmt(self):
        assert cli._fmt(1.23456, ".2f") == "1.23"
        assert cli._fmt(None, ".2f") == "N/A"
        assert cli._fmt("x", ".2f") == "N/A"
