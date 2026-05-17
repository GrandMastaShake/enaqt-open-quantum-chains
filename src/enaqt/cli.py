"""
ENAQT Command-Line Interface
==============================

argparse-based CLI with subcommands for running ENAQT analyses::

    enaqt validate --qd3set1-dir PATH --out-dir PATH [--n-points 300]
    enaqt scaling --n-max 20 --out-dir PATH
    enaqt disorder --n 7 --n-seeds 1000 --out-dir PATH [--n-jobs -1]
    enaqt optimize --n 7 --bounds-min -50 --bounds-max 50 --out-dir PATH
    enaqt fmo --out-dir PATH
    enaqt all --out-dir PATH [--n-max 20] [--n-seeds 1000]

Each subcommand calls the appropriate module function, prints a results summary,
and returns exit code 0 on success / 1 on error.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import traceback
from pathlib import Path
from typing import Callable, Sequence, Tuple

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logger = logging.getLogger("enaqt")


def _configure_logging(verbose: bool) -> None:
    """Set up logging level and format based on --verbose flag."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# Output directory helper
# ---------------------------------------------------------------------------

DEFAULT_OUT_DIR = Path("./results")
ENV_OUT_DIR_VAR = "ENAQT_OUT_DIR"


def get_out_dir(cli_value: str | None = None) -> Path:
    """Resolve output directory from CLI arg > ENV > default.

    Resolution order:
      1. ``cli_value`` if explicitly provided (non-None).
      2. ``$ENAQT_OUT_DIR`` environment variable if set and non-empty.
      3. Fall back to ``./results``.

    The directory is created automatically if it does not exist.
    """
    if cli_value is not None:
        path = Path(cli_value)
    else:
        env_val = os.environ.get(ENV_OUT_DIR_VAR, "").strip()
        path = Path(env_val) if env_val else DEFAULT_OUT_DIR

    path = path.expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# Input validation helpers
# ---------------------------------------------------------------------------

def _validate_positive_int(value: str, name: str) -> int:
    """Validate that *value* is a positive integer."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{name} must be an integer, got {value!r}")
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{name} must be > 0, got {ivalue}")
    return ivalue


def _validate_bounds_min(value: str) -> float:
    """Validate bounds-min as a float."""
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"--bounds-min must be a number, got {value!r}")
    return fvalue


def _validate_bounds_max(value: str, min_val: float) -> float:
    """Validate bounds-max as a float greater than bounds-min."""
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"--bounds-max must be a number, got {value!r}")
    if fvalue <= min_val:
        raise argparse.ArgumentTypeError(
            f"--bounds-max ({fvalue}) must be greater than --bounds-min ({min_val})"
        )
    return fvalue


def _validate_n_jobs(value: str) -> int:
    """Validate --n-jobs (positive integer or -1 for all cores)."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"--n-jobs must be an integer, got {value!r}")
    if ivalue == 0 or ivalue < -1:
        raise argparse.ArgumentTypeError(
            f"--n-jobs must be > 0 or -1 (all cores), got {ivalue}"
        )
    return ivalue


# ---------------------------------------------------------------------------
# Lazily import analysis modules so the CLI is fast for --help
# ---------------------------------------------------------------------------

_LAZY_IMPORTS_DONE = False
_run_validation: Callable | None = None
_run_sink_analysis: Callable | None = None
_run_scaling: Callable | None = None
_run_fmo_benchmark: Callable | None = None
_run_ensemble: Callable | None = None
_run_optimization: Callable | None = None


def _ensure_imports() -> None:
    """Attempt to import the analysis sub-modules once."""
    global _LAZY_IMPORTS_DONE
    global _run_validation, _run_sink_analysis
    global _run_scaling, _run_fmo_benchmark
    global _run_ensemble, _run_optimization

    if _LAZY_IMPORTS_DONE:
        return

    # spinboson
    try:
        from enaqt.spinboson import run_validation, run_sink_analysis

        _run_validation = run_validation
        _run_sink_analysis = run_sink_analysis
    except ImportError:
        logger.debug("enaqt.spinboson not available (feature branch not merged)")

    # nsite
    try:
        from enaqt.nsite import run_scaling, run_fmo_benchmark

        _run_scaling = run_scaling
        _run_fmo_benchmark = run_fmo_benchmark
    except ImportError:
        logger.debug("enaqt.nsite not available (feature branch not merged)")

    # disorder
    try:
        from enaqt.disorder import run_ensemble, run_optimization

        _run_ensemble = run_ensemble
        _run_optimization = run_optimization
    except ImportError:
        logger.debug("enaqt.disorder not available (feature branch not merged)")

    _LAZY_IMPORTS_DONE = True


def _require_module(func_name: str, func: Callable | None) -> Callable:
    """Raise a clean error if a module function is not available."""
    if func is None:
        raise RuntimeError(
            f"The '{func_name}' analysis module is not available. "
            "Make sure the corresponding feature branch is installed."
        )
    return func


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog="enaqt",
        description=(
            "ENAQT — Environment-Assisted Quantum Transport analysis toolkit. "
            "Run validation, scaling, disorder ensemble, optimization, and FMO "
            "benchmarks from the command line."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  enaqt --help\n"
            "  enaqt validate --qd3set1-dir ./SB/data --out-dir ./results\n"
            "  enaqt scaling --n-max 20 --verbose\n"
            "  enaqt disorder --n 7 --n-seeds 1000 --n-jobs 4\n"
            "  enaqt optimize --n 7 --bounds-min -50 --bounds-max 50\n"
            "  enaqt fmo --out-dir ./results\n"
            "  enaqt all --out-dir ./results --n-max 20 --n-seeds 1000\n\n"
            "Environment variable:\n"
            "  ENAQT_OUT_DIR   default output directory (used when --out-dir is omitted)\n"
        ),
    )
    subparsers = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="Available analysis commands",
        help="Analysis type to run",
    )

    # -- shared helpers -------------------------------------------------------
    def add_out_dir(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--out-dir",
            type=str,
            default=None,
            help=(
                "Output directory for results (default: resolve from "
                "$ENAQT_OUT_DIR or ./results/)"
            ),
        )

    def add_verbose(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "-v",
            "--verbose",
            action="store_true",
            default=False,
            help="Enable verbose/debug logging",
        )

    # -- validate --------------------------------------------------------------
    p_validate = subparsers.add_parser(
        "validate",
        help="QD3SET-1 spin-boson validation against HEOM trajectories",
        description="Run QD3SET-1 spin-boson validation with HEOM comparison.",
    )
    add_out_dir(p_validate)
    add_verbose(p_validate)
    p_validate.add_argument(
        "--qd3set1-dir",
        type=str,
        required=True,
        help="Path to the QD3SET-1 SB(1)/SB/data/ directory containing .npy files",
    )
    p_validate.add_argument(
        "--n-points",
        type=lambda v: _validate_positive_int(v, "--n-points"),
        default=300,
        help="Number of gamma_phi sweep points (default: 300)",
    )

    # -- scaling ---------------------------------------------------------------
    p_scaling = subparsers.add_parser(
        "scaling",
        help="N-site chain scaling analysis",
        description="Run ENAQT scaling analysis for chain lengths N=2..N_MAX.",
    )
    add_out_dir(p_scaling)
    add_verbose(p_scaling)
    p_scaling.add_argument(
        "--n-max",
        type=lambda v: _validate_positive_int(v, "--n-max"),
        default=20,
        help="Maximum chain length N (default: 20)",
    )

    # -- disorder --------------------------------------------------------------
    p_disorder = subparsers.add_parser(
        "disorder",
        help="Disorder ensemble statistical analysis",
        description="Run disorder ensemble with Gaussian random site energies.",
    )
    add_out_dir(p_disorder)
    add_verbose(p_disorder)
    p_disorder.add_argument(
        "--n",
        type=lambda v: _validate_positive_int(v, "--n"),
        default=7,
        help="Chain length N (default: 7)",
    )
    p_disorder.add_argument(
        "--n-seeds",
        type=lambda v: _validate_positive_int(v, "--n-seeds"),
        default=1000,
        help="Number of disorder realizations (default: 1000)",
    )
    p_disorder.add_argument(
        "--n-jobs",
        type=_validate_n_jobs,
        default=-1,
        help=(
            "Number of parallel jobs (-1 = all cores, default: -1)"
        ),
    )

    # -- optimize --------------------------------------------------------------
    p_optimize = subparsers.add_parser(
        "optimize",
        help="Optimal site-energy landscape optimization",
        description="Run global optimization of site energies for maximum ENAQT.",
    )
    add_out_dir(p_optimize)
    add_verbose(p_optimize)
    p_optimize.add_argument(
        "--n",
        type=lambda v: _validate_positive_int(v, "--n"),
        default=7,
        help="Chain length N (default: 7)",
    )
    p_optimize.add_argument(
        "--bounds-min",
        type=_validate_bounds_min,
        default=-50.0,
        help="Lower bound for site energies (default: -50)",
    )
    p_optimize.add_argument(
        "--bounds-max",
        type=float,
        default=50.0,
        help="Upper bound for site energies (default: 50)",
    )

    # -- fmo -------------------------------------------------------------------
    p_fmo = subparsers.add_parser(
        "fmo",
        help="FMO-7 photosynthetic complex benchmark",
        description="Run FMO-7 benchmark analysis.",
    )
    add_out_dir(p_fmo)
    add_verbose(p_fmo)

    # -- all -------------------------------------------------------------------
    p_all = subparsers.add_parser(
        "all",
        help="Run the full analysis pipeline",
        description=(
            "Run the complete ENAQT pipeline: validate, scaling, disorder, "
            "optimize, and FMO benchmark. Each step is executed sequentially."
        ),
    )
    add_out_dir(p_all)
    add_verbose(p_all)
    p_all.add_argument(
        "--n-max",
        type=lambda v: _validate_positive_int(v, "--n-max"),
        default=20,
        help="Maximum chain length for scaling (default: 20)",
    )
    p_all.add_argument(
        "--n-seeds",
        type=lambda v: _validate_positive_int(v, "--n-seeds"),
        default=1000,
        help="Number of disorder seeds (default: 1000)",
    )
    p_all.add_argument(
        "--n-jobs",
        type=_validate_n_jobs,
        default=-1,
        help="Number of parallel jobs for disorder step (default: -1 = all cores)",
    )
    p_all.add_argument(
        "--qd3set1-dir",
        type=str,
        default=None,
        help=(
            "Path to QD3SET-1 data directory (optional; if omitted, "
            "the validation step is skipped)"
        ),
    )
    p_all.add_argument(
        "--bounds-min",
        type=_validate_bounds_min,
        default=-50.0,
        help="Lower bound for site energies in optimize step (default: -50)",
    )
    p_all.add_argument(
        "--bounds-max",
        type=float,
        default=50.0,
        help="Upper bound for site energies in optimize step (default: 50)",
    )

    return parser


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _handle_validate(args: argparse.Namespace) -> int:
    """Handle the ``validate`` subcommand."""
    out_dir = get_out_dir(args.out_dir)
    qd3set1_dir = Path(args.qd3set1_dir).expanduser().resolve()
    if not qd3set1_dir.exists():
        logger.error("QD3SET-1 directory does not exist: %s", qd3set1_dir)
        return 1
    if not qd3set1_dir.is_dir():
        logger.error("QD3SET-1 path is not a directory: %s", qd3set1_dir)
        return 1

    logger.info("Running QD3SET-1 validation")
    logger.info("  Data directory : %s", qd3set1_dir)
    logger.info("  Output directory: %s", out_dir)
    logger.info("  n_points      : %d", args.n_points)

    func = _require_module("run_validation", _run_validation)
    result = func(qd3set1_dir, out_dir, n_points=args.n_points)

    _print_validate_summary(result)
    return 0


def _handle_scaling(args: argparse.Namespace) -> int:
    """Handle the ``scaling`` subcommand."""
    out_dir = get_out_dir(args.out_dir)
    logger.info("Running N-site scaling analysis")
    logger.info("  Output directory: %s", out_dir)
    logger.info("  n_max           : %d", args.n_max)

    func = _require_module("run_scaling", _run_scaling)
    result = func(n_max=args.n_max, out_dir=out_dir)

    _print_scaling_summary(result)
    return 0


def _handle_disorder(args: argparse.Namespace) -> int:
    """Handle the ``disorder`` subcommand."""
    out_dir = get_out_dir(args.out_dir)
    logger.info("Running disorder ensemble")
    logger.info("  Output directory: %s", out_dir)
    logger.info("  N               : %d", args.n)
    logger.info("  n_seeds         : %d", args.n_seeds)
    logger.info("  n_jobs          : %d", args.n_jobs)

    func = _require_module("run_ensemble", _run_ensemble)
    result = func(n=args.n, n_seeds=args.n_seeds, out_dir=out_dir, n_jobs=args.n_jobs)

    _print_disorder_summary(result)
    return 0


def _handle_optimize(args: argparse.Namespace) -> int:
    """Handle the ``optimize`` subcommand."""
    out_dir = get_out_dir(args.out_dir)
    if args.bounds_max <= args.bounds_min:
        logger.error(
            "--bounds-max (%.1f) must be greater than --bounds-min (%.1f)",
            args.bounds_max,
            args.bounds_min,
        )
        return 1

    logger.info("Running optimal disorder optimization")
    logger.info("  Output directory: %s", out_dir)
    logger.info("  N               : %d", args.n)
    logger.info("  bounds          : [%.1f, %.1f]", args.bounds_min, args.bounds_max)

    func = _require_module("run_optimization", _run_optimization)
    bounds: Tuple[float, float] = (args.bounds_min, args.bounds_max)
    result = func(n=args.n, bounds=bounds, out_dir=out_dir)

    _print_optimize_summary(result)
    return 0


def _handle_fmo(args: argparse.Namespace) -> int:
    """Handle the ``fmo`` subcommand."""
    out_dir = get_out_dir(args.out_dir)
    logger.info("Running FMO-7 benchmark")
    logger.info("  Output directory: %s", out_dir)

    func = _require_module("run_fmo_benchmark", _run_fmo_benchmark)
    result = func(out_dir=out_dir)

    _print_fmo_summary(result)
    return 0


def _handle_all(args: argparse.Namespace) -> int:
    """Handle the ``all`` subcommand — run the full pipeline sequentially."""
    out_dir = get_out_dir(args.out_dir)
    overall_success = True
    results: dict[str, dict] = {}

    logger.info("=" * 60)
    logger.info("Running full ENAQT pipeline")
    logger.info("  Output directory: %s", out_dir)
    logger.info("  n_max           : %d", args.n_max)
    logger.info("  n_seeds         : %d", args.n_seeds)
    logger.info("=" * 60)

    # -- Step 1: validate (only if data directory provided) -------------------
    if args.qd3set1_dir is not None:
        logger.info("\n[Pipeline 1/5] QD3SET-1 validation")
        try:
            qd3set1_dir = Path(args.qd3set1_dir).expanduser().resolve()
            func = _require_module("run_validation", _run_validation)
            results["validate"] = func(qd3set1_dir, out_dir, n_points=300)
            _print_validate_summary(results["validate"])
        except Exception as exc:
            logger.error("Validation step failed: %s", exc)
            overall_success = False
    else:
        logger.info("[Pipeline 1/5] Validation skipped (no --qd3set1-dir)")

    # -- Step 2: scaling ------------------------------------------------------
    logger.info("\n[Pipeline 2/5] N-site scaling")
    try:
        func = _require_module("run_scaling", _run_scaling)
        results["scaling"] = func(n_max=args.n_max, out_dir=out_dir)
        _print_scaling_summary(results["scaling"])
    except Exception as exc:
        logger.error("Scaling step failed: %s", exc)
        overall_success = False

    # -- Step 3: disorder -----------------------------------------------------
    logger.info("\n[Pipeline 3/5] Disorder ensemble")
    try:
        func = _require_module("run_ensemble", _run_ensemble)
        results["disorder"] = func(
            n=7, n_seeds=args.n_seeds, out_dir=out_dir, n_jobs=args.n_jobs
        )
        _print_disorder_summary(results["disorder"])
    except Exception as exc:
        logger.error("Disorder step failed: %s", exc)
        overall_success = False

    # -- Step 4: optimize -----------------------------------------------------
    logger.info("\n[Pipeline 4/5] Optimal disorder optimization")
    try:
        func = _require_module("run_optimization", _run_optimization)
        bounds: Tuple[float, float] = (args.bounds_min, args.bounds_max)
        results["optimize"] = func(n=7, bounds=bounds, out_dir=out_dir)
        _print_optimize_summary(results["optimize"])
    except Exception as exc:
        logger.error("Optimization step failed: %s", exc)
        overall_success = False

    # -- Step 5: fmo ----------------------------------------------------------
    logger.info("\n[Pipeline 5/5] FMO-7 benchmark")
    try:
        func = _require_module("run_fmo_benchmark", _run_fmo_benchmark)
        results["fmo"] = func(out_dir=out_dir)
        _print_fmo_summary(results["fmo"])
    except Exception as exc:
        logger.error("FMO step failed: %s", exc)
        overall_success = False

    # -- Overall summary ------------------------------------------------------
    logger.info("\n" + "=" * 60)
    logger.info("Pipeline complete")
    logger.info("  Steps succeeded: %d / 5", sum(1 for k in results if results[k]))
    logger.info("  Output directory: %s", out_dir)
    logger.info("=" * 60)

    return 0 if overall_success else 1


# ---------------------------------------------------------------------------
# Summary printers
# ---------------------------------------------------------------------------

def _print_validate_summary(result: dict | None) -> None:
    """Print human-readable summary for validation results."""
    if result is None:
        print("  No validation results.")
        return
    sym = result.get("symmetric_case", {})
    asym = result.get("asymmetric_case", {})
    print("\n  --- Validation Results ---")
    print(f"  Trajectories analyzed: {result.get('total_trajectories', 'N/A')}")
    print(f"  Symmetric  (ε=0):  detected={sym.get('enaqt_detected', False)}, "
          f"peak η={sym.get('peak_eta', 'N/A'):.4f}, "
          f"enhancement={sym.get('enhancement', 'N/A'):.2f}x")
    print(f"  Asymmetric (ε=1):  detected={asym.get('enaqt_detected', False)}, "
          f"peak η={asym.get('peak_eta', 'N/A'):.4f}, "
          f"enhancement={asym.get('enhancement', 'N/A'):.2f}x")


def _print_scaling_summary(result: dict | None) -> None:
    """Print human-readable summary for scaling results."""
    if result is None:
        print("  No scaling results.")
        return
    scaling = result.get("scaling_laws", {})
    print("\n  --- Scaling Results ---")
    print(f"  Enhancement slope  : {scaling.get('enhancement_linear_slope', 'N/A'):.3f}")
    print(f"  gp* power-law exp  : {scaling.get('gp_star_power_law_exponent', 'N/A'):.3f}")
    funnel = result.get("funnel_results", {})
    for N in sorted(funnel.keys(), key=lambda x: int(x)):
        r = funnel[N]
        print(f"  N={N:>2}: enhancement={r.get('enhancement', 'N/A'):.1f}x, "
              f"gp*={r.get('gp_star', 'N/A'):.3f}")


def _print_disorder_summary(result: dict | None) -> None:
    """Print human-readable summary for disorder ensemble results."""
    if result is None:
        print("  No disorder results.")
        return
    print("\n  --- Disorder Ensemble Results ---")
    print(f"  N                 : {result.get('N', 'N/A')}")
    print(f"  Seeds             : {result.get('n_seeds', 'N/A')}")
    print(f"  Mean enhancement  : {result.get('ensemble_enhancement_mean', 'N/A'):.1f}x")
    print(f"  Median enhancement: {result.get('ensemble_enhancement_med', 'N/A'):.1f}x")
    print(f"  Std enhancement   : {result.get('ensemble_enhancement_std', 'N/A'):.1f}x")
    print(f"  Interior peak frac: {result.get('frac_interior_peak', 'N/A')}")


def _print_optimize_summary(result: dict | None) -> None:
    """Print human-readable summary for optimization results."""
    if result is None:
        print("  No optimization results.")
        return
    print("\n  --- Optimization Results ---")
    for N, r in sorted(result.items(), key=lambda kv: int(kv[0]) if str(kv[0]).isdigit() else 0):
        if isinstance(r, dict):
            print(f"  N={N:>2}: optimal={r.get('optimal_enhancement', 'N/A'):.1f}x, "
                  f"funnel={r.get('funnel_enhancement', 'N/A'):.1f}x, "
                  f"gain_vs_funnel={r.get('gain_vs_funnel', 'N/A'):.1f}x")


def _print_fmo_summary(result: dict | None) -> None:
    """Print human-readable summary for FMO benchmark results."""
    if result is None:
        print("  No FMO results.")
        return
    fmo = result.get("fmo7_result", result)
    print("\n  --- FMO-7 Benchmark Results ---")
    print(f"  Enhancement : {fmo.get('enhancement', fmo.get('FMO7_enhancement', 'N/A')):.1f}x")
    print(f"  Peak η      : {fmo.get('eta_peak', 'N/A'):.4f}")
    print(f"  Optimal gp* : {fmo.get('gp_star', 'N/A'):.3f}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    """Main CLI entry point.

    Parameters
    ----------
    argv : list[str] | None
        Command-line arguments (defaults to ``sys.argv[1:]``).

    Returns
    -------
    int
        Exit code: 0 on success, 1 on error.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    _configure_logging(args.verbose)
    _ensure_imports()

    handlers: dict[str, Callable[[argparse.Namespace], int]] = {
        "validate": _handle_validate,
        "scaling": _handle_scaling,
        "disorder": _handle_disorder,
        "optimize": _handle_optimize,
        "fmo": _handle_fmo,
        "all": _handle_all,
    }

    handler = handlers.get(args.command, lambda _: (parser.print_help(), 0)[1])

    try:
        return handler(args)
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 130
    except RuntimeError as exc:
        logger.error("%s", exc)
        return 1
    except Exception as exc:
        logger.error("Unexpected error: %s", exc)
        if args.verbose:
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
