import subprocess
import time
from datetime import datetime
import os
import sys
import argparse
import itertools
from collections import OrderedDict
import requests
from bs4 import BeautifulSoup
import threading

# Config

monomer_file = "monomers.txt"
python_script = "monomers_to_normaliz.py"
normaliz_exe = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz"
# log_file = f"log_{n_for_covering}_{k_start}_{t}.txt" ; defined in main


tmp_monomers = "tmp_monomers.txt"
tmp_eqs = "eqs.in"


PROBE_LIMIT = 100

"""
FLAGS
-----
--mode monomer (default)
    Covering design is over monomer subsets. Each block is a k-subset of
    monomers; Normaliz is run directly on those monomers.

--mode domain
    Covering design is over binding-site types (domains). Each block selects
    k domain types; monomers are filtered to those whose domains are all
    within the selected set before running Normaliz.

--include-base (NOT IN USE RN)
    Also run the base case k=n (all monomers or all domains), which gives
    the exact full Hilbert basis at the cost of the longest runtime.

--fallback-greedy
    If a covering design C(n,k,t) is not found in the online database,
    compute one locally using a greedy set-cover algorithm instead of
    skipping that k value. Warning: can be slow for large n and k.

--tolerance [float]  (default: 1.0)
    Pruning leeway multiplier. A k value is only pruned if:
        estimated_total_time > best_time_so_far * tolerance
    E.g. --tolerance 1.2 allows up to 20% slack before pruning,
    which helps avoid premature pruning due to random timing spikes.

Interative  (during run)
----------------------------------
s   Skip the current k value immediately (does not update best time).
"""

# Utility functions

def cleanup_normaliz_files():
    patterns = ["eqs.out", "eqs.gen", "eqs.inv", "eqs.cst", "eqs.typ", "eqs.egn"]
    for pattern in patterns:
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
            vec = [int(x) for x in line.split()]
            vectors.append(vec)
        except ValueError:
            pass

    return vectors


def expand_vector_to_full_space(reduced_vector, selected_indices, n):
    """Used in monomer mode: maps reduced monomer-subset vector back to full monomer space."""
    full_vector = [0] * n
    for i, idx in enumerate(selected_indices):
        if i < len(reduced_vector):
            full_vector[idx] = reduced_vector[i]
    return tuple(full_vector)


def expand_vector_to_full_monomer_space(reduced_vector, filtered_monomers_indices, n_monomers):
    """Used in domain mode: maps filtered-monomer vector back to full monomer space."""
    monomer_part = reduced_vector[:len(filtered_monomers_indices)]
    full_vector = [0] * n_monomers
    for i, idx in enumerate(filtered_monomers_indices):
        if i < len(monomer_part):
            full_vector[idx] = monomer_part[i]
    return tuple(full_vector)

_skip_current = False

def _listen_for_skip():
    global _skip_current
    while True:
        cmd = input()
        if cmd.strip().lower() == 's':
            _skip_current = True

def start_input_listener():
    t = threading.Thread(target=_listen_for_skip, daemon=True)
    t.start()

def check_and_clear_skip() -> bool:
    global _skip_current
    if _skip_current:
        _skip_current = False
        return True
    return False
# -------------------------
# Domain helpers (for domain mode)
# -------------------------

def get_domains_from_monomer(monomer: str) -> set[str]:
    return {domain.rstrip("*") for domain in monomer.split()}


def get_all_unique_domains(monomers: list[str]) -> list[str]:
    all_domains: dict[str, None] = OrderedDict()
    for monomer in monomers:
        for domain in get_domains_from_monomer(monomer):
            all_domains[domain] = None
    return list(all_domains.keys())


def filter_monomers_by_domains(monomers: list[str], selected_domains: list[str]) -> list[str]:
    """Return only monomers whose domains are all within the selected set."""
    selected_set = set(selected_domains)
    return [m for m in monomers if get_domains_from_monomer(m).issubset(selected_set)]


# -------------------------
# Online covering fetch
# -------------------------

def fetch_covering_online(v: int, k: int, t: int) -> list[list[int]]:
    url = f"https://ljcr.dmgordon.org/show_cover.php?v={v}&k={k}&t={t}"
    print(f"Fetching covering C({v},{k},{t})...")

    r = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=30)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "lxml")
    pre = soup.find("pre")

    if pre is None:
        raise RuntimeError(f"No covering stored online for C({v},{k},{t})")

    blocks = []
    for line in pre.text.strip().splitlines():
        row = [int(x) for x in line.split()]
        if row:
            blocks.append(row)

    if not blocks:
        raise RuntimeError(f"No covering stored online for C({v},{k},{t})")

    print(f"Retrieved {len(blocks)} blocks.")
    return blocks


# -------------------------
# Greedy covering design
# -------------------------

def compute_covering_greedy(v: int, k: int, t: int) -> list[list[int]]:
    print(f"Computing greedy covering C({v},{k},{t}) ...")

    universe = list(range(1, v + 1))
    uncovered = set(itertools.combinations(universe, t))
    blocks = []
    max_possible = len(list(itertools.combinations(range(k), t)))

    while uncovered:
        best_block: tuple[int, ...] = next(itertools.combinations(universe, k))
        best_count = -1

        for block in itertools.combinations(universe, k):
            count = len(set(itertools.combinations(block, t)) & uncovered)
            if count > best_count:
                best_count = count
                best_block = block
            if best_count == max_possible:
                break

        blocks.append(list(best_block))
        uncovered -= set(itertools.combinations(best_block, t))

    print(f"Greedy covering complete: {len(blocks)} blocks.")
    return blocks


# -------------------------
# Unified covering loader
# -------------------------

def load_covering_blocks(v: int, k: int, t: int, fallback_greedy: bool = False):
    if k == v:
        yield list(range(1, v + 1))
        return

    try:
        for block in fetch_covering_online(v, k, t):
            yield block
        return
    except (RuntimeError, requests.RequestException) as e:
        if not fallback_greedy:
            raise RuntimeError(
                f"Online fetch failed for C({v},{k},{t}): {e}\n"
                f"  Tip: run with --fallback-greedy to compute locally."
            ) from e
        print(f"Online fetch failed ({e}). Falling back to greedy computation.")

    for block in compute_covering_greedy(v, k, t):
        yield block


# -------------------------
# Core pipeline: MONOMER MODE
# -------------------------

def run_for_k_value_monomer(k, all_monomers, n, t, min_normaliz_time, log, fallback_greedy, tolerance):
    """
    Covering design over monomers: each block is a k-subset of monomer indices.
    Runs Normaliz on each subset of monomers directly.
    """
    print(f"\n{'='*70}")
    print(f"[monomer mode] Testing k={k}")
    print(f"{'='*70}\n")

    try:
        blocks = list(load_covering_blocks(n, k, t, fallback_greedy=fallback_greedy))
    except (RuntimeError, requests.RequestException) as e:
        print(f"Skipping k={k}: {e}")
        log.write(f"\nk={k}: SKIPPED (covering unavailable: {e})\n")
        log.flush()
        return None, min_normaliz_time

    num_subsets = len(blocks)
    print(f"Using covering design: {num_subsets} subsets\n")

    all_hilbert_vectors: set[tuple[int, ...]] = set()
    times = []
    probe_size = min(PROBE_LIMIT, num_subsets)

    for phase, index_range in [("PROBE", range(probe_size)),
                                ("FULL",  range(probe_size, num_subsets))]:

        if phase == "FULL":
            probe_time = sum(times)
            avg_time = probe_time / probe_size if probe_size else 0
            estimated_total = avg_time * num_subsets
            print(f"Probe complete: avg={avg_time:.3f}s, estimated total={estimated_total:.2f}s")

            if estimated_total > min_normaliz_time * tolerance:
                print(f"PRUNED k={k}: estimated {estimated_total:.2f}s > best {min_normaliz_time:.2f}s * tolerance {tolerance}")
                log.write(f"\nk={k}: PRUNED\n"
                          f"  Probe time: {probe_time:.3f}s for {probe_size} subsets\n"
                          f"  Estimated: {estimated_total:.2f}s vs best {min_normaliz_time:.2f}s * tolerance {tolerance}\n")
                log.flush()
                return None, min_normaliz_time

        for idx in index_range:
            if check_and_clear_skip():
                print(f"Skipping k={k} by user request.")
                log.write(f"\nk={k}: SKIPPED by user\n")
                log.flush()
                return None, min_normaliz_time
            block = blocks[idx]
            cleanup_normaliz_files()

            selected_indices = [x - 1 for x in block]
            subset = [all_monomers[i] for i in selected_indices]

            with open(tmp_monomers, "w") as f:
                for m in subset:
                    f.write(m + "\n")

            subprocess.run(["python", python_script],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

            start = time.time()
            subprocess.run([normaliz_exe, "-d", "-N", tmp_eqs],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            times.append(time.time() - start)

            for rv in read_hilbert_basis():
                all_hilbert_vectors.add(expand_vector_to_full_space(rv, selected_indices, n))

    return _finish_run(k, num_subsets, times, all_hilbert_vectors, log)


# -------------------------
# Core pipeline: DOMAIN MODE
# -------------------------

def run_for_k_value_domain(k, all_monomers, all_domains, n_domains, n_monomers,
                            t, min_normaliz_time, log, fallback_greedy, tolerance):
    """
    Covering design over binding site types (domains): each block is a k-subset
    of domain indices. Monomers are filtered to only those whose domains are all
    within the selected subset, then Normaliz is run on those monomers.
    """
    print(f"\n{'='*70}")
    print(f"[domain mode] Testing k={k}")
    print(f"{'='*70}\n")

    try:
        blocks = list(load_covering_blocks(n_domains, k, t, fallback_greedy=fallback_greedy))
    except (RuntimeError, requests.RequestException) as e:
        print(f"Skipping k={k}: {e}")
        log.write(f"\nk={k}: SKIPPED (covering unavailable: {e})\n")
        log.flush()
        return None, min_normaliz_time

    num_subsets = len(blocks)
    print(f"Using covering design: {num_subsets} subsets\n")

    all_hilbert_vectors: set[tuple[int, ...]] = set()
    times = []
    probe_size = min(PROBE_LIMIT, num_subsets)

    for phase, index_range in [("PROBE", range(probe_size)),
                                ("FULL",  range(probe_size, num_subsets))]:

        if phase == "FULL":
            probe_time = sum(times)
            avg_time = probe_time / probe_size if probe_size else 0
            estimated_total = avg_time * num_subsets
            print(f"Probe complete: avg={avg_time:.3f}s, estimated total={estimated_total:.2f}s")

            if estimated_total > min_normaliz_time * tolerance:
                print(f"PRUNED k={k}: estimated {estimated_total:.2f}s > best {min_normaliz_time:.2f}s * tolerance {tolerance}")
                log.write(f"\nk={k}: PRUNED\n"
                          f"  Probe time: {probe_time:.3f}s for {probe_size} subsets\n"
                          f"  Estimated: {estimated_total:.2f}s vs best {min_normaliz_time:.2f}s * tolerance {tolerance}\n")
                log.flush()
                return None, min_normaliz_time

        for idx in index_range:
            if check_and_clear_skip():
                print(f"Skipping k={k} by user request.")
                log.write(f"\nk={k}: SKIPPED by user\n")
                log.flush()
                return None, min_normaliz_time
            block = blocks[idx]
            cleanup_normaliz_files()

            # Convert 1-indexed block to domain names
            selected_domains = [all_domains[i - 1] for i in block]

            filtered_monomers = filter_monomers_by_domains(all_monomers, selected_domains)
            if not filtered_monomers:
                continue

            filtered_indices = [
                i for i, m in enumerate(all_monomers) if m in filtered_monomers
            ]

            with open(tmp_monomers, "w") as f:
                for m in filtered_monomers:
                    f.write(m + "\n")

            subprocess.run(["python", python_script],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

            start = time.time()
            subprocess.run([normaliz_exe, "-d", "-N", tmp_eqs],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            times.append(time.time() - start)

            for rv in read_hilbert_basis():
                all_hilbert_vectors.add(
                    expand_vector_to_full_monomer_space(rv, filtered_indices, n_monomers)
                )

    return _finish_run(k, num_subsets, times, all_hilbert_vectors, log)


# -------------------------
# Shared finish helper
# -------------------------

def _finish_run(k, num_subsets, times, all_hilbert_vectors, log):
    total_time = sum(times)
    avg = (total_time / num_subsets) if num_subsets else 0.0

    log.write(f"\nk={k}: FULL RUN COMPLETE\n"
              f"  Number of subsets: {num_subsets}\n"
              f"  Total Normaliz time: {total_time:.3f}s\n"
              f"  Avg Normaliz time per subset: {avg:.3f}s\n"
              f"  Unique Hilbert basis vectors: {len(all_hilbert_vectors)}\n")
    log.flush()

    print(f"\nk = {k}: total time={total_time:.3f}s, "
          f"avg time/subset={avg:.3f}s, unique vectors={len(all_hilbert_vectors)}")

    return {
        "k": k,
        "num_subsets": num_subsets,
        "total_normaliz_time": total_time,
        "avg_normaliz_time_per_subset": avg,
        "unique_vectors": len(all_hilbert_vectors)
    }, total_time


# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser(description="Vary-k Hilbert basis pipeline")
    parser.add_argument(
        "--fallback-greedy",
        action="store_true",
        help="Compute covering locally if not found online."
    )
    parser.add_argument(
        "--mode",
        choices=["monomer", "domain"],
        default="monomer",
        help=(
            "'monomer' (default): covering design over monomer subsets. "
            "'domain': covering design over binding-site-type subsets; "
            "monomers are filtered to those using only the selected domains."
        )
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0,
        help="Pruning tolerance multiplier. Only prune if estimated_time > best_time * tolerance. "
            "Default 1.0 (no leeway). E.g. 1.2 allows 20%% slack before pruning."
    )
    args = parser.parse_args()

    cleanup_normaliz_files()

    all_monomers = []
    with open(monomer_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                line = line.split(":", 1)[1].strip()
            if line:
                all_monomers.append(line)

    n_monomers = len(all_monomers)
    all_domains = get_all_unique_domains(all_monomers)
    n_domains = len(all_domains)

    print(f"\nDetected {n_monomers} monomers, {n_domains} unique binding site types")
    print(f"Mode: {args.mode}")
    print(f"Greedy fallback: {'ENABLED' if args.fallback_greedy else 'DISABLED'}")
    print(f"Note: Online covering lookup is limited to n < 100, k ≤ 25, t ≤ 8\n")

    # n for the covering design depends on mode
    n_for_covering = n_monomers if args.mode == "monomer" else n_domains

    t = int(input("Enter value of t (default=5): ") or 5)
    k_start = int(input(f"Enter starting k (≤ {n_for_covering})\nk decreases from starting value to t\n(default 25)): ").strip() or 25)
    include_base = input(f"Include base case k={n_for_covering}? (y/n, default=n): ").strip().lower() or "n"
    start_input_listener()

    k_values = list(range(k_start, t, -1))
    if include_base == "y" and n_for_covering not in k_values:
        k_values = [n_for_covering] + k_values

    log_file = f"log_{n_for_covering}_{k_start}_{t}.txt"

    min_total_time = float("inf")
    results = []

    with open(log_file, "a") as log:
        log.write(f"Vary k Benchmark\n"
                  f"Started: {datetime.now()}\n"
                  f"Mode: {args.mode}\n"
                  f"n_monomers={n_monomers}, n_domains={n_domains}, t={t}\n"
                  f"k values: {k_values}\n"
                  f"Greedy fallback: {args.fallback_greedy}\n"
                  + "="*70 + "\n")

        try:
            for k in k_values:
                if args.mode == "monomer":
                    result, new_time = run_for_k_value_monomer(
                        k, all_monomers, n_monomers, t, min_total_time, log, args.fallback_greedy, args.tolerance
                    )
                else:
                    result, new_time = run_for_k_value_domain(
                        k, all_monomers, all_domains, n_domains, n_monomers,
                        t, min_total_time, log, args.fallback_greedy, args.tolerance
                    )

                if result is None:
                    continue

                results.append(result)
                if new_time < min_total_time:
                    min_total_time = new_time

        except KeyboardInterrupt:
            log.flush()

    if results:
        best = min(results, key=lambda r: r["total_normaliz_time"])
        print(f"\nOptimal k = {best['k']} ({best['total_normaliz_time']:.2f}s)")
    else:
        print("\nNo k completed fully.")


if __name__ == "__main__":
    main()