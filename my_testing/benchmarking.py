"""
benchmark.py
============
Runs the covering pipeline probe-and-prune experiment across three TBN system
types (linear cascade, binary tree, DNA cascade), sweeping over module sizes and
covering strengths t in {3, 4, 5, 6, 7, 8}.

For each (system, size, t) cell:
  - Probe-only flag ON  → records best projected total runtime for each t
  - Probe-only flag OFF → records best projected total AND actual full runtime

Results are written to:
  /Users/archit/Projects/Hilbert Basis Algorithm/my_testing/logs/benchmarking/

One master JSON results file is produced alongside per-run log files, making
it straightforward to load into a plotting script for grouped bar charts.

Usage
-----
  # Probe-only (fast, estimates only)
  python benchmark.py --probe

  # Full runs (probe then enumerate on best k for each t)
  python benchmark.py

  # Specific systems / sizes only
  python benchmark.py --systems cascade --cascade-sizes 5 6 7

  # Include base case (k=n) in the sweep
  python benchmark.py --include-base

  # Use greedy fallback when online covering fetch fails
  python benchmark.py --fallback-greedy

  # Save Hilbert basis vectors to disk (full mode only)
  python benchmark.py --save
"""

import argparse
import json
import os
import sys
import time
from datetime import datetime

# ---- Import pipeline functions directly ----
# covering_pipeline.py must be on the Python path (same directory is fine).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from claude_pipeline_new import (
    load_monomers,
    get_all_unique_domains,
    run_covering_sweep,
    cleanup_normaliz_files,
    start_input_listener,
    K_MAX,
)

# -------------------------
# System registry
# -------------------------

BASE_TBN_DIR = "/Users/archit/Projects/Hilbert Basis Algorithm/example-tbns"
LOG_BASE_DIR = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/logs/benchmarking"

SYSTEMS = {
    "cascade": {
        "label":    "Linear Cascade",
        "path_fn":  lambda m: os.path.join(
            BASE_TBN_DIR, f"monomers_cascade_n{m}.txt"
        ),
        "sizes":    [5, 6, 7, 8, 9],
        "size_key": "m",
    },
    "binary": {
        "label":    "Binary Tree",
        "path_fn":  lambda d: os.path.join(
            BASE_TBN_DIR, f"monomers_binary_tree_d{d}.txt"
        ),
        "sizes":    [3, 4],
        "size_key": "d",
    },
    "dna": {
        "label":    "DNA Cascade",
        "path_fn":  lambda m: os.path.join(
            BASE_TBN_DIR, f"monomers_dna_tbn_depth{m}.txt"
        ),
        "sizes":    [4, 5, 6, 7],
        "size_key": "m",
    },
}

T_VALUES = [3, 4, 5, 6, 7, 8]

# -------------------------
# Helpers
# -------------------------

def safe_float(x):
    return round(x, 6) if x is not None else None


def make_run_id(system_key, size, t):
    return f"{system_key}_size{size}_t{t}"


# -------------------------
# Single (system, size, t) cell
# -------------------------

def run_cell(
    system_key, size, t,
    mode, probe_only, include_base,
    fallback_greedy, save, save_dir,
    log,
):
    """
    Run one (system, size, t) cell.
    Returns a dict with timing results, or None on failure.
    """
    spec     = SYSTEMS[system_key]
    path     = spec["path_fn"](size)
    run_id   = make_run_id(system_key, size, t)

    log.write(
        f"\n{'#'*70}\n"
        f"# Cell: system={system_key} ({spec['label']}), "
        f"{spec['size_key']}={size}, t={t}\n"
        f"# Monomer file: {path}\n"
        f"{'#'*70}\n"
    )
    log.flush()

    if not os.path.exists(path):
        msg = f"  SKIPPED: monomer file not found: {path}"
        print(msg)
        log.write(msg + "\n")
        log.flush()
        return None

    all_monomers = load_monomers(path)
    all_domains  = get_all_unique_domains(all_monomers)
    n_monomers   = len(all_monomers)
    n_domains    = len(all_domains)
    n            = n_monomers if mode == "monomer" else n_domains

    log.write(
        f"  n_monomers={n_monomers}, n_domains={n_domains}, "
        f"n_for_covering={n}\n"
    )

    # Sanity: need at least one valid k = t+1 .. min(K_MAX, n)
    if t + 1 > min(K_MAX, n):
        msg = (
            f"  SKIPPED: t={t} leaves no valid k "
            f"(need t+1 <= min(K_MAX={K_MAX}, n={n}))"
        )
        print(msg)
        log.write(msg + "\n")
        log.flush()
        return None

    # La Jolla limits check (warn only — fallback_greedy handles it)
    if not fallback_greedy and (n >= 100 or t > 8):
        msg = (
            f"  WARNING: n={n} or t={t} may exceed La Jolla limits "
            f"(n<100, t<=8). Consider --fallback-greedy."
        )
        print(msg)
        log.write(msg + "\n")

    cell_start = time.time()

    try:
        best_k, best_projected, full_result = run_covering_sweep(
            t               = t,
            all_monomers    = all_monomers,
            mode            = mode,
            log             = log,
            fallback_greedy = fallback_greedy,
            probe_only      = probe_only,
            include_base    = include_base,
            save            = save,
            save_dir        = os.path.join(save_dir, run_id) if save_dir else None,
        )
    except Exception as e:
        msg = f"  ERROR in run_covering_sweep: {e}"
        print(msg)
        log.write(msg + "\n")
        log.flush()
        return None

    cell_elapsed = time.time() - cell_start

    result = {
        "run_id":          run_id,
        "system":          system_key,
        "system_label":    spec["label"],
        "size_key":        spec["size_key"],
        "size":            size,
        "t":               t,
        "mode":            mode,
        "n_monomers":      n_monomers,
        "n_domains":       n_domains,
        "best_k":          best_k,
        "projected_total": safe_float(best_projected),
        "cell_wall_time":  round(cell_elapsed, 3),
        # Full-run fields (None when probe_only=True)
        "full_wall_time":        None,
        "full_normaliz_time":    None,
        "full_overhead_time":    None,
        "full_unique_vectors":   None,
    }

    if full_result is not None:
        result.update({
            "full_wall_time":      safe_float(full_result["total_wall_time"]),
            "full_normaliz_time":  safe_float(full_result["total_normaliz_time"]),
            "full_overhead_time":  safe_float(full_result["overhead_time"]),
            "full_unique_vectors": full_result["unique_vectors"],
        })

    log.write(
        f"\n  Cell result: best_k={best_k}, "
        f"projected={safe_float(best_projected)}, "
        f"cell_wall={cell_elapsed:.2f}s"
    )
    if full_result:
        log.write(
            f", full_normaliz={result['full_normaliz_time']}s, "
            f"unique_vectors={result['full_unique_vectors']}"
        )
    log.write("\n")
    log.flush()

    return result


# -------------------------
# Main benchmark loop
# -------------------------

def run_benchmark(args):
    os.makedirs(LOG_BASE_DIR, exist_ok=True)
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path   = os.path.join(
        LOG_BASE_DIR,
        f"benchmark_{timestamp}_{'probe' if args.probe else 'full'}.log"
    )
    json_path  = os.path.join(
        LOG_BASE_DIR,
        f"benchmark_{timestamp}_{'probe' if args.probe else 'full'}_results.json"
    )

    # Determine which systems and sizes to run
    systems_to_run = {}
    for key, spec in SYSTEMS.items():
        if args.systems and key not in args.systems:
            continue
        # Allow per-system size overrides from CLI
        size_override = getattr(args, f"{key}_sizes", None)
        sizes = size_override if size_override else spec["sizes"]
        systems_to_run[key] = sizes

    t_values = args.t_values

    total_cells = sum(
        len(sizes) * len(t_values)
        for sizes in systems_to_run.values()
    )

    print(f"\n{'='*70}")
    print(f"Benchmark run: {timestamp}")
    print(f"Mode: {'probe-only' if args.probe else 'probe + full enumeration'}")
    print(f"Systems: {list(systems_to_run.keys())}")
    print(f"t values: {t_values}")
    print(f"Total cells: {total_cells}")
    print(f"Log: {log_path}")
    print(f"Results JSON: {json_path}")
    print(f"{'='*70}\n")

    start_input_listener()
    all_results = []
    cell_num    = 0

    with open(log_path, "w") as log:
        log.write(
            f"Benchmark — Pareto-Optimal Polymer Enumeration\n"
            f"Started: {datetime.now()}\n"
            f"Mode: {'probe-only' if args.probe else 'probe + full'}\n"
            f"Systems: {list(systems_to_run.keys())}\n"
            f"t_values: {t_values}\n"
            f"include_base: {args.include_base}\n"
            f"fallback_greedy: {args.fallback_greedy}\n"
            f"covering_mode: {args.mode}\n"
            f"total_cells: {total_cells}\n"
            + "=" * 70 + "\n"
        )

        for system_key, sizes in systems_to_run.items():
            for size in sizes:
                for t in t_values:
                    cell_num += 1
                    print(
                        f"\n[{cell_num}/{total_cells}] "
                        f"system={system_key}, "
                        f"{SYSTEMS[system_key]['size_key']}={size}, t={t}"
                    )

                    cleanup_normaliz_files()

                    result = run_cell(
                        system_key      = system_key,
                        size            = size,
                        t               = t,
                        mode            = args.mode,
                        probe_only      = args.probe,
                        include_base    = args.include_base,
                        fallback_greedy = args.fallback_greedy,
                        save            = args.save,
                        save_dir        = args.save_dir if args.save else None,
                        log             = log,
                    )

                    if result is not None:
                        all_results.append(result)
                        # Checkpoint after every cell
                        _write_json(all_results, json_path)

                        projected_str = (
                            f"{result['projected_total']:.3f}s"
                            if result['projected_total'] is not None
                            else "N/A"
                        )
                        full_str = (
                            f"  full_normaliz={result['full_normaliz_time']:.3f}s"
                            if result['full_normaliz_time'] is not None
                            else ""
                        )
                        print(
                            f"  → best_k={result['best_k']}, "
                            f"projected={projected_str}"
                            f"{full_str}"
                        )

        log.write(
            f"\n{'='*70}\n"
            f"Benchmark complete: {datetime.now()}\n"
            f"Total cells attempted: {cell_num}\n"
            f"Cells with results:    {len(all_results)}\n"
        )

    _write_json(all_results, json_path)

    print(f"\n{'='*70}")
    print(f"Benchmark complete. {len(all_results)}/{total_cells} cells succeeded.")
    print(f"Results: {json_path}")
    print(f"Log:     {log_path}")

    _print_summary_table(all_results, t_values)

    return all_results


def _write_json(results, path):
    with open(path, "w") as f:
        json.dump(
            {
                "generated": datetime.now().isoformat(),
                "results":   results,
            },
            f,
            indent=2,
        )


def _print_summary_table(results, t_values):
    """Print a compact summary table: rows = (system, size), cols = t values."""
    if not results:
        return

    # Group by (system, size)
    from collections import defaultdict
    table = defaultdict(dict)
    for r in results:
        row_key = (r["system_label"], r["size_key"], r["size"])
        table[row_key][r["t"]] = r

    col_w = 22
    header = f"{'System / Size':<30}" + "".join(
        f"{'t=' + str(t):<{col_w}}" for t in t_values
    )
    print(f"\n{'='*len(header)}")
    print("SUMMARY (projected total / full normaliz time)")
    print(header)
    print("-" * len(header))

    for (label, size_key, size), t_map in sorted(table.items()):
        row = f"{label} {size_key}={size:<4}"
        row = f"{row:<30}"
        for t in t_values:
            if t not in t_map:
                cell = "—"
            else:
                r    = t_map[t]
                proj = f"{r['projected_total']:.2f}s" if r["projected_total"] else "?"
                full = f"/{r['full_normaliz_time']:.2f}s" if r["full_normaliz_time"] else ""
                cell = f"k={r['best_k']} {proj}{full}"
            row += f"{cell:<{col_w}}"
        print(row)

    print("=" * len(header))


# -------------------------
# CLI
# -------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark covering pipeline across TBN system types, sizes, and "
            "covering strengths t in {3..8}."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--probe", action="store_true",
        help="Probe-only mode: report best projected runtime per (system, size, t) "
             "without running full enumeration.",
    )
    parser.add_argument(
        "--mode", choices=["monomer", "domain"], default="monomer",
        help="Covering universe: monomers or domain types. Default: monomer.",
    )
    parser.add_argument(
        "--systems", nargs="+", choices=list(SYSTEMS.keys()),
        default=None,
        help="Which systems to benchmark. Default: all.",
    )
    parser.add_argument(
        "--t-values", nargs="+", type=int, default=T_VALUES,
        dest="t_values", metavar="T",
        help=f"Covering strength values to sweep. Default: {T_VALUES}.",
    )
    parser.add_argument(
        "--cascade-sizes", nargs="+", type=int, default=None,
        dest="cascade_sizes", metavar="M",
        help="Override module sizes for linear cascade. Default: 5 6 7 8 9.",
    )
    parser.add_argument(
        "--binary-sizes", nargs="+", type=int, default=None,
        dest="binary_sizes", metavar="D",
        help="Override depth values for binary tree. Default: 3 4.",
    )
    parser.add_argument(
        "--dna-sizes", nargs="+", type=int, default=None,
        dest="dna_sizes", metavar="M",
        help="Override module sizes for DNA cascade. Default: 4 5 6 7.",
    )
    parser.add_argument(
        "--include-base", action="store_true", dest="include_base",
        help="Also probe/run k=n (full system) for each (system, size, t).",
    )
    parser.add_argument(
        "--fallback-greedy", action="store_true", dest="fallback_greedy",
        help="Compute covering locally when online fetch fails.",
    )
    parser.add_argument(
        "--save", action="store_true",
        help="Save Hilbert basis vectors to disk (full mode only).",
    )
    parser.add_argument(
        "--save-dir", type=str,
        default="/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/hilbert_benchmark_output",
        dest="save_dir", metavar="PATH",
        help="Root directory for saved Hilbert basis files.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_benchmark(args)