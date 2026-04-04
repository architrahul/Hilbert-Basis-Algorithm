"""
fill_missing.py
===============
Fills in missing benchmark data points by running the covering pipeline
directly at k=25 for each (system, size, t) cell that was skipped or
returned k=None in the main benchmark run.

For each cell, runs:
  1. Probe at k=25 → projected total runtime
  2. Full enumeration at k=25 → actual normaliz runtime

Results are appended to a JSON file and a summary table is printed at the end.

Usage
-----
  python fill_missing.py                    # run all missing cells
  python fill_missing.py --probe            # probe-only (no full enumeration)
  python fill_missing.py --fallback-greedy  # use greedy if online fetch fails
  python fill_missing.py --save             # save Hilbert basis vectors
"""

import json
import os
import sys
import time
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from claude_pipeline_new import (
    load_monomers,
    get_all_unique_domains,
    load_covering_blocks,
    probe_k,
    full_run_k,
    cleanup_normaliz_files,
    start_input_listener,
    save_polymer_vectors,
)

import argparse

LOG_BASE_DIR = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/logs/benchmarking"
BASE_TBN_DIR = "/Users/archit/Projects/Hilbert Basis Algorithm/example-tbns"

FIXED_K = 25

# -------------------------
# All missing cells
# (system_key, size, t, monomer_path)
# -------------------------

def monomer_path(system_key, size):
    if system_key == "cascade":
        return os.path.join(BASE_TBN_DIR, f"monomers_cascade_n{size}.txt")
    elif system_key == "binary":
        return os.path.join(BASE_TBN_DIR, f"monomers_binary_tree_d{size}.txt")
    elif system_key == "dna":
        return os.path.join(BASE_TBN_DIR, f"monomers_dna_tbn_depth{size}.txt")
    raise ValueError(f"Unknown system: {system_key}")


SYSTEM_LABELS = {
    "cascade": "Linear Cascade",
    "binary":  "Binary Tree",
    "dna":     "DNA Cascade",
}

SIZE_KEYS = {
    "cascade": "m",
    "binary":  "d",
    "dna":     "m",
}

# Every cell that showed — or k=None in the summary table
MISSING_CELLS = [
    # Binary Tree d=3, t=3
    ("binary",  3, 3),
    # Binary Tree d=4 (all t values)
    ("binary",  4, 3),
    ("binary",  4, 4),
    ("binary",  4, 5),
    ("binary",  4, 6),
    ("binary",  4, 7),
    ("binary",  4, 8),
    # DNA Cascade m=4, t=3..7
    ("dna",     4, 3),
    ("dna",     4, 4),
    ("dna",     4, 5),
    ("dna",     4, 6),
    ("dna",     4, 7),
    # DNA Cascade m=5, t=3..5
    ("dna",     5, 3),
    ("dna",     5, 4),
    ("dna",     5, 5),
    # DNA Cascade m=6, t=3..4
    ("dna",     6, 3),
    ("dna",     6, 4),
    # DNA Cascade m=7, t=3
    ("dna",     7, 3),
    # Linear Cascade m=5, t=3..5
    ("cascade", 5, 3),
    ("cascade", 5, 4),
    ("cascade", 5, 5),
    # Linear Cascade m=6, t=3..4
    ("cascade", 6, 3),
    ("cascade", 6, 4),
    # Linear Cascade m=7, t=3
    ("cascade", 7, 3),
    # Linear Cascade m=8, t=3
    ("cascade", 8, 3),
    # Linear Cascade m=9, t=3
    ("cascade", 9, 3),
]


# -------------------------
# Run one cell at fixed k=25
# -------------------------

def run_cell_fixed_k(system_key, size, t, mode, probe_only,
                     fallback_greedy, save, save_dir, log):
    path = monomer_path(system_key, size)
    label = SYSTEM_LABELS[system_key]
    size_key = SIZE_KEYS[system_key]

    log.write(
        f"\n{'#'*70}\n"
        f"# {label} {size_key}={size}, t={t}, k={FIXED_K}\n"
        f"# {path}\n"
        f"{'#'*70}\n"
    )
    log.flush()

    if not os.path.exists(path):
        msg = f"  SKIPPED: file not found: {path}"
        print(msg); log.write(msg + "\n"); log.flush()
        return None

    all_monomers = load_monomers(path)
    all_domains  = get_all_unique_domains(all_monomers)
    n_monomers   = len(all_monomers)
    n_domains    = len(all_domains)
    n            = n_monomers if mode == "monomer" else n_domains

    log.write(f"  n_monomers={n_monomers}, n_domains={n_domains}, n={n}\n")

    if FIXED_K > n:
        msg = f"  SKIPPED: k={FIXED_K} > n={n}"
        print(msg); log.write(msg + "\n"); log.flush()
        return None

    # Load covering blocks for C(n, k=25, t)
    try:
        blocks = load_covering_blocks(n, FIXED_K, t, fallback_greedy=fallback_greedy)
    except RuntimeError as e:
        msg = f"  SKIPPED: covering unavailable: {e}"
        print(msg); log.write(msg + "\n"); log.flush()
        return None

    num_blocks = len(blocks)
    log.write(f"  C({n},{FIXED_K},{t}): {num_blocks} blocks\n")
    log.flush()

    cell_start = time.time()

    # ---- Probe ----
    projected, probe_times, _ = probe_k(
        FIXED_K, blocks, n, all_monomers, mode, all_domains, n_monomers,
        best_projected=None,  # no pruning — we always want this estimate
        log=log,
    )

    result = {
        "run_id":          f"{system_key}_size{size}_t{t}",
        "system":          system_key,
        "system_label":    label,
        "size_key":        size_key,
        "size":            size,
        "t":               t,
        "mode":            mode,
        "n_monomers":      n_monomers,
        "n_domains":       n_domains,
        "best_k":          FIXED_K,
        "projected_total": round(projected, 6) if projected is not None else None,
        "cell_wall_time":  None,
        "full_wall_time":        None,
        "full_normaliz_time":    None,
        "full_overhead_time":    None,
        "full_unique_vectors":   None,
    }

    if probe_only or projected is None:
        result["cell_wall_time"] = round(time.time() - cell_start, 3)
        proj_str = f"{projected:.3f}s" if projected else "N/A"
        print(f"  → k={FIXED_K}, projected={proj_str}")
        log.write(f"  Probe-only result: projected={proj_str}\n")
        log.flush()
        return result

    # ---- Full enumeration ----
    print(f"  Running full enumeration at k={FIXED_K} ...")
    full_result, _ = full_run_k(
        FIXED_K, blocks, n, all_monomers, mode, all_domains, n_monomers, log
    )

    result["cell_wall_time"] = round(time.time() - cell_start, 3)

    if full_result is not None:
        result.update({
            "full_wall_time":      round(full_result["total_wall_time"],      6),
            "full_normaliz_time":  round(full_result["total_normaliz_time"],  6),
            "full_overhead_time":  round(full_result["overhead_time"],        6),
            "full_unique_vectors": full_result["unique_vectors"],
        })

        if save and save_dir:
            os.makedirs(save_dir, exist_ok=True)
            out = os.path.join(save_dir, f"hilbert_{system_key}_size{size}_k{FIXED_K}_t{t}_{mode}.txt")
            save_polymer_vectors(
                full_result["vectors"], out,
                n_monomers=n_monomers,
                comment=f"k={FIXED_K}, t={t}, mode={mode}"
            )

        print(
            f"  → k={FIXED_K}, projected={result['projected_total']:.3f}s, "
            f"full_normaliz={result['full_normaliz_time']:.3f}s, "
            f"unique_vectors={result['full_unique_vectors']}"
        )

    return result


# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser(description="Fill missing benchmark cells at k=25.")
    parser.add_argument("--probe", action="store_true",
                        help="Probe-only: no full enumeration.")
    parser.add_argument("--mode", choices=["monomer", "domain"], default="monomer")
    parser.add_argument("--fallback-greedy", action="store_true", dest="fallback_greedy")
    parser.add_argument("--save", action="store_true")
    parser.add_argument(
        "--save-dir", type=str, dest="save_dir",
        default="/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/hilbert_benchmark_output",
    )
    args = parser.parse_args()

    os.makedirs(LOG_BASE_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path  = os.path.join(LOG_BASE_DIR, f"fill_missing_{timestamp}.log")
    json_path = os.path.join(LOG_BASE_DIR, f"fill_missing_{timestamp}_results.json")

    print(f"\n{'='*70}")
    print(f"Fill missing cells — k={FIXED_K} forced")
    print(f"Cells to run: {len(MISSING_CELLS)}")
    print(f"Mode: {'probe-only' if args.probe else 'probe + full'}")
    print(f"Log:  {log_path}")
    print(f"{'='*70}\n")

    start_input_listener()
    all_results = []

    with open(log_path, "w") as log:
        log.write(
            f"Fill Missing Cells — k={FIXED_K} forced\n"
            f"Started: {datetime.now()}\n"
            f"mode={args.mode}, probe_only={args.probe}, "
            f"fallback_greedy={args.fallback_greedy}\n"
            f"Total cells: {len(MISSING_CELLS)}\n"
            + "=" * 70 + "\n"
        )

        for i, (system_key, size, t) in enumerate(MISSING_CELLS, 1):
            label    = SYSTEM_LABELS[system_key]
            size_key = SIZE_KEYS[system_key]
            print(f"\n[{i}/{len(MISSING_CELLS)}] {label} {size_key}={size}, t={t}")

            cleanup_normaliz_files()

            result = run_cell_fixed_k(
                system_key      = system_key,
                size            = size,
                t               = t,
                mode            = args.mode,
                probe_only      = args.probe,
                fallback_greedy = args.fallback_greedy,
                save            = args.save,
                save_dir        = args.save_dir,
                log             = log,
            )

            if result is not None:
                all_results.append(result)
                with open(json_path, "w") as f:
                    json.dump({"generated": datetime.now().isoformat(),
                               "results": all_results}, f, indent=2)

        log.write(
            f"\n{'='*70}\n"
            f"Done: {datetime.now()}\n"
            f"Cells succeeded: {len(all_results)}/{len(MISSING_CELLS)}\n"
        )

    # ---- Summary table ----
    print(f"\n{'='*70}")
    print(f"Completed {len(all_results)}/{len(MISSING_CELLS)} cells.")
    print(f"Results: {json_path}\n")

    col_w = 30
    header = f"{'Cell':<35}" + f"{'projected':>{col_w}}" + f"{'full_normaliz':>{col_w}}"
    print(header)
    print("-" * len(header))
    for r in all_results:
        cell  = f"{r['system_label']} {r['size_key']}={r['size']}, t={r['t']}"
        proj  = f"{r['projected_total']:.3f}s"  if r['projected_total']  else "N/A"
        full  = f"{r['full_normaliz_time']:.3f}s" if r['full_normaliz_time'] else "—"
        print(f"{cell:<35}{proj:>{col_w}}{full:>{col_w}}")
    print("=" * len(header))

    cleanup_normaliz_files()


if __name__ == "__main__":
    main()