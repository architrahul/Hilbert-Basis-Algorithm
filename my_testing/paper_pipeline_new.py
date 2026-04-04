import subprocess
import time
from datetime import datetime
import os
import sys
import argparse
import itertools
import random
from collections import OrderedDict
import requests
from bs4 import BeautifulSoup
import threading
from export_polymers import save_polymer_vectors

# -------------------------
# Config
# -------------------------

DEFAULT_MONOMER_FILE = "/Users/archit/Projects/Hilbert Basis Algorithm/example-tbns/monomers_cascade_n8.txt"
PYTHON_SCRIPT        = "monomers_to_normaliz.py"
NORMALIZ_EXE         = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz"

TMP_MONOMERS = "tmp_monomers.txt"
TMP_EQS      = "eqs.in"
PROBE_LIMIT  = 100
K_MAX        = 25

# -------------------------
# Utility functions
# -------------------------

def cleanup_normaliz_files():
    for pattern in ["eqs.out", "eqs.gen", "eqs.inv", "eqs.cst", "eqs.typ", "eqs.egn"]:
        if os.path.exists(pattern):
            os.remove(pattern)


def read_hilbert_basis(base_filename="eqs"):
    output_file = f"{base_filename}.out"
    if not os.path.exists(output_file):
        return []

    vectors = []
    with open(output_file, "r") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if "Hilbert basis elements:" in line:
            start = i + 1
            break
    if start is None:
        return []

    for line in lines[start:]:
        line = line.strip()
        if not line or line.startswith("***") or "extreme rays" in line.lower():
            break
        if line.startswith("#"):
            continue
        try:
            vectors.append([int(x) for x in line.split()])
        except ValueError:
            pass

    return vectors


def expand_vector_to_full_space(reduced_vector, selected_indices, n):
    """Monomer mode: maps reduced monomer-subset vector back to full monomer space."""
    full_vector = [0] * n
    for i, idx in enumerate(selected_indices):
        if i < len(reduced_vector):
            full_vector[idx] = reduced_vector[i]
    return tuple(full_vector)


def expand_vector_to_full_monomer_space(reduced_vector, filtered_monomers_indices, n_monomers):
    """Domain mode: maps filtered-monomer vector back to full monomer space."""
    monomer_part = reduced_vector[:len(filtered_monomers_indices)]
    full_vector  = [0] * n_monomers
    for i, idx in enumerate(filtered_monomers_indices):
        if i < len(monomer_part):
            full_vector[idx] = monomer_part[i]
    return tuple(full_vector)


_skip_current = False

def _listen_for_skip():
    global _skip_current
    while True:
        if input().strip().lower() == 's':
            _skip_current = True

def start_input_listener():
    threading.Thread(target=_listen_for_skip, daemon=True).start()

def check_and_clear_skip() -> bool:
    global _skip_current
    if _skip_current:
        _skip_current = False
        return True
    return False


# -------------------------
# Domain helpers
# -------------------------

def get_domains_from_monomer(monomer: str) -> set:
    return {domain.rstrip("*") for domain in monomer.split()}


def get_all_unique_domains(monomers: list) -> list:
    seen = OrderedDict()
    for monomer in monomers:
        for domain in get_domains_from_monomer(monomer):
            seen[domain] = None
    return list(seen.keys())


def filter_monomers_by_domains(monomers: list, selected_domains: list) -> list:
    """Return only monomers whose domains are all within the selected set."""
    selected_set = set(selected_domains)
    return [m for m in monomers if get_domains_from_monomer(m).issubset(selected_set)]


# -------------------------
# Covering design fetch / compute
# -------------------------

def fetch_covering_online(v: int, k: int, t: int) -> list:
    url = f"https://ljcr.dmgordon.org/show_cover.php?v={v}&k={k}&t={t}"
    r   = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=30)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "lxml")
    pre  = soup.find("pre")
    if pre is None:
        raise RuntimeError(f"No covering stored online for C({v},{k},{t})")

    blocks = []
    for line in pre.text.strip().splitlines():
        row = [int(x) for x in line.split()]
        if row:
            blocks.append(row)

    if not blocks:
        raise RuntimeError(f"No covering stored online for C({v},{k},{t})")
    return blocks


def compute_covering_greedy(v: int, k: int, t: int) -> list:
    print(f"  Computing greedy covering C({v},{k},{t}) ...")
    universe  = list(range(1, v + 1))
    uncovered = set(itertools.combinations(universe, t))
    blocks    = []
    max_cov   = len(list(itertools.combinations(range(k), t)))

    while uncovered:
        best_block = next(itertools.combinations(universe, k))
        best_count = -1
        for block in itertools.combinations(universe, k):
            count = len(set(itertools.combinations(block, t)) & uncovered)
            if count > best_count:
                best_count = count
                best_block = block
            if best_count == max_cov:
                break
        blocks.append(list(best_block))
        uncovered -= set(itertools.combinations(best_block, t))

    print(f"  Greedy covering complete: {len(blocks)} blocks.")
    return blocks


def load_covering_blocks(v: int, k: int, t: int, fallback_greedy: bool = False) -> list:
    """
    Returns list of blocks for C(v,k,t).
    Raises RuntimeError if unavailable and fallback_greedy is False.
    When k == v the single full block is returned directly.
    """
    if k == v:
        return [list(range(1, v + 1))]

    try:
        return fetch_covering_online(v, k, t)
    except (RuntimeError, requests.RequestException) as e:
        if not fallback_greedy:
            raise RuntimeError(
                f"Online fetch failed for C({v},{k},{t}): {e}. "
                f"Use --fallback-greedy to compute locally."
            ) from e
        print(f"  Online fetch failed ({e}). Falling back to greedy.")
        return compute_covering_greedy(v, k, t)


# -------------------------
# Single-block Normaliz runner
# -------------------------

def run_normaliz_on_subset(subset_monomers: list) -> tuple:
    """
    Write monomers, invoke monomers_to_normaliz.py, run Normaliz.
    Returns (elapsed_seconds, hilbert_vectors_raw).
    """
    with open(TMP_MONOMERS, "w") as f:
        for m in subset_monomers:
            f.write(m + "\n")

    subprocess.run(
        ["python", PYTHON_SCRIPT],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )

    t0 = time.time()
    subprocess.run(
        [NORMALIZ_EXE, "-d", "-N", TMP_EQS],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    return time.time() - t0, read_hilbert_basis()


# -------------------------
# Probe phase for one k value
# -------------------------

def probe_k(k, t, blocks, n, all_monomers, mode, all_domains, n_monomers,
            best_projected, log):
    """
    Run up to PROBE_LIMIT randomly sampled blocks for the given k.

    Early-abort rule: if cumulative probe wall time exceeds best_projected at
    any point during the probe, this k is pruned immediately.

    Returns:
        projected_total  – avg_block_time * num_blocks, or None if pruned/empty
        probe_times      – list of individual block times recorded
        num_blocks       – total number of blocks in the covering design
    """
    num_blocks  = len(blocks)
    probe_size  = min(PROBE_LIMIT, num_blocks)
    sample_idxs = random.sample(range(num_blocks), probe_size)

    probe_times = []
    log.write(f"    Probing k={k}: {num_blocks} blocks total, sampling {probe_size}\n")
    log.flush()

    for i, idx in enumerate(sample_idxs):
        if check_and_clear_skip():
            log.write(f"    k={k}: SKIPPED by user during probe\n")
            log.flush()
            return None, probe_times, num_blocks

        block = blocks[idx]
        cleanup_normaliz_files()

        if mode == "monomer":
            selected_indices = [x - 1 for x in block]
            subset           = [all_monomers[j] for j in selected_indices]
        else:
            selected_domains = [all_domains[x - 1] for x in block]
            subset           = filter_monomers_by_domains(all_monomers, selected_domains)
            if not subset:
                continue

        elapsed, _ = run_normaliz_on_subset(subset)
        probe_times.append(elapsed)
        cumulative = sum(probe_times)

        # Early abort: cumulative probe time already exceeds best projected total
        if best_projected is not None and cumulative > best_projected:
            log.write(
                f"    k={k}: PROBE PRUNED at iteration {i+1}/{probe_size} "
                f"(cumulative={cumulative:.3f}s > best_projected={best_projected:.3f}s)\n"
            )
            log.flush()
            print(
                f"  [t={t}] k={k}: pruned at block {i+1}/{probe_size} "
                f"(cumulative={cumulative:.3f}s > best={best_projected:.3f}s)"
            )
            return None, probe_times, num_blocks

    if not probe_times:
        return None, probe_times, num_blocks

    avg_time        = sum(probe_times) / len(probe_times)
    projected_total = avg_time * num_blocks

    log.write(
        f"    k={k}: probe complete — "
        f"probe_iters={len(probe_times)}, avg={avg_time:.4f}s, "
        f"num_blocks={num_blocks}, projected_total={projected_total:.3f}s\n"
    )
    log.flush()
    return projected_total, probe_times, num_blocks


# -------------------------
# Full enumeration for one k value
# -------------------------

def full_run_k(k, blocks, n, all_monomers, mode, all_domains, n_monomers, log):
    """
    Run Normaliz on every block for the given k, collecting all Hilbert vectors.
    Returns (result_dict, total_normaliz_time).
    """
    num_blocks          = len(blocks)
    all_hilbert_vectors = set()
    times               = []
    wall_start          = time.time()

    log.write(f"    Full run k={k}: {num_blocks} blocks\n")
    log.flush()

    for idx, block in enumerate(blocks):
        if check_and_clear_skip():
            log.write(f"    k={k}: ABORTED by user during full run\n")
            log.flush()
            return None, float("inf")

        cleanup_normaliz_files()

        if mode == "monomer":
            selected_indices = [x - 1 for x in block]
            subset           = [all_monomers[j] for j in selected_indices]
        else:
            selected_domains = [all_domains[x - 1] for x in block]
            subset           = filter_monomers_by_domains(all_monomers, selected_domains)
            if not subset:
                continue
            selected_indices = [
                i for i, m in enumerate(all_monomers) if m in subset
            ]

        elapsed, raw_vectors = run_normaliz_on_subset(subset)
        times.append(elapsed)

        for rv in raw_vectors:
            if mode == "monomer":
                all_hilbert_vectors.add(
                    expand_vector_to_full_space(rv, selected_indices, n)
                )
            else:
                all_hilbert_vectors.add(
                    expand_vector_to_full_monomer_space(rv, selected_indices, n_monomers)
                )

        if (idx + 1) % 50 == 0 or (idx + 1) == num_blocks:
            print(f"  k={k}: {idx+1}/{num_blocks} blocks done "
                  f"({sum(times):.1f}s, {len(all_hilbert_vectors)} vectors)")

    wall_time     = time.time() - wall_start
    normaliz_time = sum(times)
    overhead      = wall_time - normaliz_time
    avg_normaliz  = normaliz_time / num_blocks if num_blocks else 0.0

    log.write(
        f"    k={k}: FULL RUN COMPLETE — "
        f"wall={wall_time:.3f}s  normaliz={normaliz_time:.3f}s  "
        f"overhead={overhead:.3f}s  avg/block={avg_normaliz:.4f}s  "
        f"unique_vectors={len(all_hilbert_vectors)}\n"
    )
    log.flush()

    return {
        "k":                           k,
        "num_blocks":                  num_blocks,
        "total_wall_time":             wall_time,
        "total_normaliz_time":         normaliz_time,
        "overhead_time":               overhead,
        "avg_normaliz_time_per_block": avg_normaliz,
        "unique_vectors":              len(all_hilbert_vectors),
        "vectors":                     all_hilbert_vectors,
    }, normaliz_time


## -------------------------
# Main covering sweep for one t value
# -------------------------

def run_covering_sweep(
    t,
    all_monomers,
    mode,
    log,
    fallback_greedy = False,
    probe_only      = False,
    include_base    = False,
    save            = False,
    save_dir        = None,
):
    """
    Sweep k from t+1 up to min(K_MAX, n), probing each k with probe-and-prune.
    Skips any k where no covering design is available (no crash).

    probe_only=True  → returns (best_k, best_projected, None)
    probe_only=False → runs full enumeration on best k;
                       returns (best_k, best_projected, full_result_dict)

    include_base=True additionally probes/runs k=n.
    If include_base=True, k=n is always run FIRST to establish the initial upper bound.
    """
    n_monomers  = len(all_monomers)
    all_domains = get_all_unique_domains(all_monomers)
    n_domains   = len(all_domains)
    n           = n_monomers if mode == "monomer" else n_domains

    # Candidate k values excluding base case
    k_values = list(range(t + 1, min(K_MAX, n) + 1))
    if n in k_values:
        k_values.remove(n)

    log.write(
        f"\n{'='*60}\n"
        f"  Covering sweep: t={t}, mode={mode}, n={n}\n"
        f"  include_base={include_base}\n"
        f"  k_values={([n] if include_base else []) + k_values}\n"
        f"{'='*60}\n"
    )
    log.flush()

    best_projected = None
    best_k         = None
    best_blocks    = None
    probe_summary  = {}  # k -> projected_total or None

    # ------------------------------------------------------------
    # STEP 1: Run base case FIRST if requested
    # ------------------------------------------------------------
    if include_base:
        print(f"\n  [t={t}] Probing base case k={n} first ...")

        try:
            base_blocks = load_covering_blocks(n, n, t, fallback_greedy=fallback_greedy)
        except RuntimeError as e:
            print(f"  Base case k={n} unavailable: {e}")
            log.write(f"    k={n}: SKIPPED (base case unavailable: {e})\n")
            log.flush()
            probe_summary[n] = None
        else:
            projected, probe_times, num_blocks = probe_k(
                n, t, base_blocks, n, all_monomers, mode, all_domains, n_monomers,
                None, log
            )
            probe_summary[n] = projected

            if projected is not None:
                best_projected = projected
                best_k         = n
                best_blocks    = base_blocks
                print(f"  [t={t}] Base case k={n}: projected={projected:.3f}s")

                # --- NEW: full run + save for base case ---
                if not probe_only:
                    print(f"\n  [t={t}] Running full enumeration on base case k={n} ...")
                    base_result, _ = full_run_k(
                        n, base_blocks, n, all_monomers, mode, all_domains, n_monomers, log
                    )
                    if base_result is not None and save and save_dir:
                        os.makedirs(save_dir, exist_ok=True)
                        output_path = os.path.join(
                            save_dir, f"hilbert_k{n}_t{t}_{mode}.txt"
                        )
                        save_polymer_vectors(
                            base_result["vectors"],
                            output_path,
                            n_monomers=n_monomers,
                            comment=f"covering mode, k={n}, t={t}, mode={mode} (base case)"
                        )
                        print(f"  Saved base case vectors → {output_path}")

    # ------------------------------------------------------------
    # STEP 2: Probe remaining k values using current best for pruning
    # ------------------------------------------------------------
    for k in k_values:
        print(f"\n  [t={t}] Probing k={k} ...")

        try:
            blocks = load_covering_blocks(n, k, t, fallback_greedy=fallback_greedy)
        except RuntimeError as e:
            print(f"  Skipping k={k}: {e}")
            log.write(f"    k={k}: SKIPPED (covering unavailable: {e})\n")
            log.flush()
            probe_summary[k] = None
            continue

        projected, probe_times, num_blocks = probe_k(
            k, t, blocks, n, all_monomers, mode, all_domains, n_monomers,
            best_projected, log
        )

        probe_summary[k] = projected

        if projected is None:
            continue

        if best_projected is None:
            print(f"  [t={t}] k={k}: projected={projected:.3f}s  (first estimate)")
        else:
            print(f"  [t={t}] k={k}: projected={projected:.3f}s  "
                  f"(best so far: {best_projected:.3f}s)")

        if best_projected is None or projected < best_projected:
            best_projected = projected
            best_k         = k
            best_blocks    = blocks

    # ---- Probe sweep summary ----
    summary_str = {
        k: (f"{v:.3f}s" if v is not None else "pruned/skipped")
        for k, v in probe_summary.items()
    }
    log.write(
        f"\n  t={t}: Probe sweep complete.\n"
        f"    Probe summary (projected totals): {summary_str}\n"
    )

    if best_k is None:
        log.write(f"  t={t}: No valid k found.\n")
        log.flush()
        print(f"  [t={t}] No valid k found.")
        return None, None, None

    log.write(f"    Best k={best_k}, projected={best_projected:.3f}s\n")
    log.flush()
    print(f"\n  [t={t}] Best k={best_k}, projected total={best_projected:.3f}s")

    if probe_only:
        return best_k, best_projected, None

    # ---- Full enumeration on best k ----
    print(f"\n  [t={t}] Running full enumeration on k={best_k} ...")
    full_result, _ = full_run_k(
        best_k, best_blocks, n, all_monomers, mode, all_domains, n_monomers, log
    )

    if full_result is not None and save and save_dir:
        os.makedirs(save_dir, exist_ok=True)
        output_path = os.path.join(
            save_dir, f"hilbert_k{best_k}_t{t}_{mode}.txt"
        )
        save_polymer_vectors(
            full_result["vectors"],
            output_path,
            n_monomers=n_monomers,
            comment=f"covering mode, k={best_k}, t={t}, mode={mode}"
        )

    return best_k, best_projected, full_result


# -------------------------
# Naive enumeration
# -------------------------

def run_naive(k, all_monomers, mode, log):
    """
    Exhaustive: runs Normaliz on every k-subset of the universe.
    Returns (result_dict, total_normaliz_time).
    """
    n_monomers  = len(all_monomers)
    all_domains = get_all_unique_domains(all_monomers)
    n_domains   = len(all_domains)
    n           = n_monomers if mode == "monomer" else n_domains

    from math import comb
    num_subsets = comb(n, k)

    print(f"\n[naive / {mode} mode] k={k}, running all {num_subsets} subsets")
    log.write(f"\nnaive {mode} k={k}: {num_subsets} subsets\n")
    log.flush()

    all_hilbert_vectors = set()
    times      = []
    wall_start = time.time()

    for idx, block in enumerate(itertools.combinations(range(n), k)):
        if check_and_clear_skip():
            print("Naive run aborted by user.")
            log.write("  ABORTED by user.\n")
            log.flush()
            return None, float("inf")

        cleanup_normaliz_files()
        selected_indices = list(block)

        if mode == "monomer":
            subset = [all_monomers[i] for i in selected_indices]
        else:
            selected_domains = [all_domains[i] for i in selected_indices]
            subset           = filter_monomers_by_domains(all_monomers, selected_domains)
            if not subset:
                continue
            selected_indices = [
                i for i, m in enumerate(all_monomers) if m in subset
            ]

        elapsed, raw_vectors = run_normaliz_on_subset(subset)
        times.append(elapsed)

        for rv in raw_vectors:
            if mode == "monomer":
                all_hilbert_vectors.add(
                    expand_vector_to_full_space(rv, selected_indices, n)
                )
            else:
                all_hilbert_vectors.add(
                    expand_vector_to_full_monomer_space(rv, selected_indices, n_monomers)
                )

        if (idx + 1) % 50 == 0 or (idx + 1) == num_subsets:
            print(f"  Progress: {idx+1}/{num_subsets} "
                  f"({sum(times):.1f}s, {len(all_hilbert_vectors)} vectors)")

    wall_time     = time.time() - wall_start
    normaliz_time = sum(times)
    overhead      = wall_time - normaliz_time
    avg_normaliz  = normaliz_time / num_subsets if num_subsets else 0.0

    log.write(
        f"\nk={k}: FULL RUN COMPLETE — "
        f"wall={wall_time:.3f}s  normaliz={normaliz_time:.3f}s  "
        f"overhead={overhead:.3f}s  avg/block={avg_normaliz:.4f}s  "
        f"unique_vectors={len(all_hilbert_vectors)}\n"
    )
    log.flush()

    return {
        "k":                           k,
        "num_blocks":                  num_subsets,
        "total_wall_time":             wall_time,
        "total_normaliz_time":         normaliz_time,
        "overhead_time":               overhead,
        "avg_normaliz_time_per_block": avg_normaliz,
        "unique_vectors":              len(all_hilbert_vectors),
        "vectors":                     all_hilbert_vectors,
    }, normaliz_time


# -------------------------
# Helpers for CLI
# -------------------------

def load_monomers(monomer_file: str) -> list:
    monomers = []
    with open(monomer_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                line = line.split(":", 1)[1].strip()
            if line:
                monomers.append(line)
    return monomers


def validate_args(args, n_monomers, n_domains):
    errors   = []
    warnings = []
    n        = n_monomers if args.mode == "monomer" else n_domains

    if args.strategy == "covering":
        if args.t < 1:
            errors.append(f"--t must be >= 1 (got {args.t}).")
        if args.t >= min(K_MAX, n):
            errors.append(
                f"--t must be < min(K_MAX={K_MAX}, n={n})={min(K_MAX,n)} "
                f"so that at least one valid k = t+1 .. {min(K_MAX,n)} exists."
            )
        if not args.fallback_greedy and (n >= 100 or args.t > 8):
            issues = []
            if n >= 100:   issues.append(f"n={n} >= 100")
            if args.t > 8: issues.append(f"t={args.t} > 8")
            errors.append(
                f"La Jolla repository limits: n < 100, t <= 8. "
                f"Parameters exceed limits: {', '.join(issues)}. Add --fallback-greedy."
            )
    else:  # naive
        for flag, val, default, name in [
            (args.t,            args.t,            5,     "--t"),
            (args.include_base, args.include_base, False, "--include-base"),
            (args.probe,        args.probe,        False, "--probe"),
        ]:
            if val != default:
                warnings.append(f"{name} has no effect in naive mode.")
        if not (1 <= args.k_naive <= (n_monomers if args.mode == "monomer" else n_domains)):
            errors.append(
                f"--k-naive must be in [1, n] "
                f"(got {args.k_naive}, n={n_monomers if args.mode == 'monomer' else n_domains})."
            )

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  [!] {w}")
    if errors:
        print("\nErrors:")
        for i, e in enumerate(errors, 1):
            print(f"  [{i}] {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Pareto-Optimal Polymer Enumeration via Hilbert Basis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Covering strategy: probe sweep for t=5, pick best k\n"
            "  python covering_pipeline.py --t 5\n\n"
            "  # Probe-only (no full enumeration)\n"
            "  python covering_pipeline.py --t 5 --probe\n\n"
            "  # Domain mode, include full-system base case, save output\n"
            "  python covering_pipeline.py --t 4 --mode domain --include-base --save\n\n"
            "  # Naive exhaustive enumeration at k=8\n"
            "  python covering_pipeline.py --strategy naive --k-naive 8\n\n"
            "Type 's' + Enter during a run to skip the current k or abort naive."
        )
    )

    parser.add_argument("--strategy", choices=["covering", "naive"],
                        default="covering")
    parser.add_argument("--t", type=int, default=5, metavar="INT",
                        help="Covering strength. Sweeps k from t+1 to min(25, n). Default: 5.")
    parser.add_argument("--mode", choices=["monomer", "domain"], default="monomer")
    parser.add_argument("--monomer-file", type=str, default=DEFAULT_MONOMER_FILE,
                        dest="monomer_file")
    parser.add_argument("--include-base", action="store_true", dest="include_base",
                        help="Also probe/run k=n (full system).")
    parser.add_argument("--fallback-greedy", action="store_true", dest="fallback_greedy",
                        help="Compute covering locally if online fetch fails.")
    parser.add_argument("--probe", action="store_true",
                        help="Probe-only: estimate runtimes without full enumeration.")
    parser.add_argument("--save", action="store_true",
                        help="Save Hilbert basis vectors to disk.")
    parser.add_argument("--save-dir", type=str, default="./hilbert_output",
                        dest="save_dir", metavar="PATH")
    parser.add_argument("--k-naive", type=int, default=8, dest="k_naive",
                        metavar="INT",
                        help="k for naive exhaustive enumeration. Default: 8.")

    args = parser.parse_args()

    if args.save:
        os.makedirs(args.save_dir, exist_ok=True)

    cleanup_normaliz_files()

    all_monomers = load_monomers(args.monomer_file)
    n_monomers   = len(all_monomers)
    all_domains  = get_all_unique_domains(all_monomers)
    n_domains    = len(all_domains)

    print(f"\nLoaded {n_monomers} monomers, {n_domains} unique domain types")
    print(f"Strategy: {args.strategy}  |  Mode: {args.mode}")

    validate_args(args, n_monomers, n_domains)

    os.makedirs("logs", exist_ok=True)
    start_input_listener()

    # ---- Covering ----
    if args.strategy == "covering":
        log_file = (
            f"logs/log_covering_{args.mode}_n{n_monomers}_t{args.t}"
            f"{'_probe' if args.probe else ''}.txt"
        )
        with open(log_file, "a") as log:
            log.write(
                f"Covering Strategy — Pareto-Optimal Polymer Enumeration\n"
                f"Started: {datetime.now()}\n"
                f"Mode: {args.mode}  |  t={args.t}  |  probe_only={args.probe}\n"
                f"n_monomers={n_monomers}  n_domains={n_domains}\n"
                f"include_base={args.include_base}  fallback_greedy={args.fallback_greedy}\n"
                + "=" * 70 + "\n"
            )
            try:
                best_k, best_projected, full_result = run_covering_sweep(
                    t               = args.t,
                    all_monomers    = all_monomers,
                    mode            = args.mode,
                    log             = log,
                    fallback_greedy = args.fallback_greedy,
                    probe_only      = args.probe,
                    include_base    = args.include_base,
                    save            = args.save,
                    save_dir        = args.save_dir,
                )
            except KeyboardInterrupt:
                print("\nInterrupted.")
                log.write("\nInterrupted by user.\n")
                log.flush()
                cleanup_normaliz_files()
                return

        if best_k is not None:
            print(f"\nBest k={best_k}, projected total={best_projected:.3f}s")
            if full_result:
                print(
                    f"Full run: wall={full_result['total_wall_time']:.3f}s  "
                    f"normaliz={full_result['total_normaliz_time']:.3f}s  "
                    f"unique_vectors={full_result['unique_vectors']}"
                )

    # ---- Naive ----
    else:
        k        = args.k_naive
        log_file = f"logs/log_naive_{args.mode}_n{n_monomers}_k{k}.txt"
        with open(log_file, "a") as log:
            log.write(
                f"Naive Enumeration — Pareto-Optimal Polymer Enumeration\n"
                f"Started: {datetime.now()}\n"
                f"Mode: {args.mode}  |  k={k}\n"
                f"n_monomers={n_monomers}  n_domains={n_domains}\n"
                + "=" * 70 + "\n"
            )
            try:
                result, _ = run_naive(k, all_monomers, args.mode, log)
            except KeyboardInterrupt:
                print("\nInterrupted.")
                log.write("\nInterrupted by user.\n")
                log.flush()
                cleanup_normaliz_files()
                return

            if result is not None:
                print(
                    f"\nNaive k={k}: normaliz={result['total_normaliz_time']:.3f}s  "
                    f"unique_vectors={result['unique_vectors']}"
                )
                if args.save:
                    save_polymer_vectors(
                        result["vectors"],
                        os.path.join(args.save_dir, f"hilbert_naive_k{k}_{args.mode}.txt"),
                        n_monomers=n_monomers,
                        comment=f"naive mode, k={k}"
                    )

    cleanup_normaliz_files()







def main2():
    """
    Leakage analysis: runs full Pareto-optimal enumeration (monomer mode, k=25, t=5)
    for both the full (M) and incomplete (M') 8-module linear cascade systems.

    Calls load_covering_blocks + full_run_k directly (bypassing the probe sweep)
    since k=25 and t=5 are fixed by experimental design.

    Output saved to:
        /Users/archit/Projects/Hilbert Basis Algorithm/my_testing/logs/leakage_analysis/
    """
    SAVE_DIR = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/logs/leakage_analysis"
    MONOMER_FILES = {
        "full":       "/Users/archit/Projects/Hilbert Basis Algorithm/example-tbns/monomers_cascade_n8.txt",
        "incomplete": "/Users/archit/Projects/Hilbert Basis Algorithm/example-tbns/monomers_cascade_n8_incomplete.txt",
    }
    T    = 5
    K    = 25
    MODE = "monomer"  # <-- changed

    os.makedirs(SAVE_DIR, exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    cleanup_normaliz_files()
    start_input_listener()

    for label, monomer_file in MONOMER_FILES.items():
        print(f"\n{'='*70}")
        print(f"  System  : {label.upper()}  ({monomer_file})")
        print(f"  Mode    : {MODE}  |  t={T}  |  k={K}")
        print(f"{'='*70}")

        all_monomers = load_monomers(monomer_file)
        n_monomers   = len(all_monomers)
        all_domains  = get_all_unique_domains(all_monomers)
        n_domains    = len(all_domains)
        n            = n_monomers  # <-- changed: monomer mode uses n_monomers

        print(f"  Loaded {n_monomers} monomers, {n_domains} unique domain types")

        if K > n:
            print(f"  WARNING: k={K} > n={n} for '{label}' system — skipping.")
            continue

        log_file = os.path.join(
            SAVE_DIR,
            f"log_leakage_{MODE}_{label}_n{n_monomers}_t{T}_k{K}.txt"
        )

        with open(log_file, "a") as log:
            log.write(
                f"Leakage Analysis — Pareto-Optimal Polymer Enumeration\n"
                f"Started    : {datetime.now()}\n"
                f"System     : {label}  ({monomer_file})\n"
                f"Mode       : {MODE}  |  t={T}  |  k={K}\n"
                f"n_monomers={n_monomers}  n_domains={n_domains}\n"
                + "=" * 70 + "\n"
            )
            log.flush()

            # --- Load covering design for fixed k, t ---
            try:
                blocks = load_covering_blocks(n, K, T, fallback_greedy=False)
            except RuntimeError as e:
                msg = f"  [{label}] Could not load covering blocks: {e}\n"
                print(msg)
                log.write(msg)
                log.flush()
                continue

            print(f"  [{label}] Loaded {len(blocks)} covering blocks. Running full enumeration ...")

            # --- Full enumeration (no probe, direct run) ---
            try:
                result, _ = full_run_k(
                    K, blocks, n, all_monomers, MODE, all_domains, n_monomers, log
                )
            except KeyboardInterrupt:
                print("\n  Interrupted.")
                log.write("\nInterrupted by user.\n")
                log.flush()
                cleanup_normaliz_files()
                return

            if result is None:
                msg = f"  [{label}] full_run_k returned None (possibly skipped/timed out).\n"
                print(msg)
                log.write(msg)
                log.flush()
                continue

            # --- Summary ---
            summary = (
                f"\n  [{label}] k={K}: "
                f"wall={result['total_wall_time']:.3f}s  "
                f"normaliz={result['total_normaliz_time']:.3f}s  "
                f"unique_vectors={result['unique_vectors']}\n"
            )
            print(summary)
            log.write(summary)
            log.flush()

            # --- Save vectors ---
            out_path = os.path.join(
                SAVE_DIR,
                f"hilbert_k{K}_t{T}_{MODE}_{label}.txt"
            )
            save_polymer_vectors(
                result["vectors"],
                out_path,
                n_monomers=n_monomers,
                comment=f"leakage analysis, covering mode, k={K}, t={T}, system={label}"
            )
            print(f"  Saved vectors → {out_path}")

    cleanup_normaliz_files()
    print("\n  Leakage analysis complete.")

if __name__ == "__main__":
    main()