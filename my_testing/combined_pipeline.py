import subprocess
import time
from datetime import datetime
import os
import sys
import argparse
import itertools
from collections import OrderedDict
from unittest import result
import requests
from bs4 import BeautifulSoup
import threading
from export_polymers import save_polymer_vectors
import numpy as np

# Config

monomer_file = "/Users/archit/Projects/Hilbert Basis Algorithm/example-tbns/cascade_n10.txt"
python_script = "monomers_to_normaliz.py"
normaliz_exe = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz"

save = True
save_dir = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/pareto_optimal_set_cascade10"
os.makedirs(save_dir, exist_ok=True)


tmp_monomers = "tmp_monomers.txt"
tmp_eqs = "eqs.in"

PROBE_LIMIT = 100

"""
FLAGS
-----
--strategy covering (default)
    Use the covering design strategy. A covering design C(n, k, t) is fetched
    or computed so that every t-element subset of the universe appears in at
    least one block; Normaliz is run only on those blocks. The algorithm sweeps
    k from k-start down to t+1 with probe-and-prune to find the fastest k.

--strategy naive
    Exhaustive subset enumeration. Normaliz is run on every k-subset of the
    universe (C(n,k,k) blocks in covering-design terms), for a single fixed k
    given by --k-start. This recovers P*_k exactly with no approximation from
    the covering design, at the cost of C(n,k) Normaliz invocations.
    --t and --tolerance are ignored in naive mode.

--mode monomer (default)
    The universe for subset selection is the set of monomers. Each block is a
    k-subset of monomers; Normaliz is run directly on those monomers.

--mode domain
    The universe for subset selection is the set of binding-site types (domains).
    Each block selects k domain types; monomers are filtered to those whose
    domains all lie within the selected set before running Normaliz.

--t [int]  (default: 5)  [covering strategy only]
    The covering strength: every t-element subset of the universe is guaranteed
    to appear in at least one block. Must satisfy 1 <= t < k-start.

--k-start [int]  (default: 25)
    Covering strategy: the largest value of k to sweep from (down to t+1).
    Naive strategy: the single fixed k value to use.
    Must satisfy k-start <= n in both cases, and k-start > t in covering mode.

--include-base  [covering strategy only]
    Also run the base case k=n (all monomers or all domain types), recovering
    the exact full Hilbert basis P* at the cost of the longest possible runtime.

--fallback-greedy  [covering strategy only]
    If a covering design C(n,k,t) is not found in the La Jolla Covering
    Repository, compute one locally using a greedy set-cover algorithm.
    Warning: can be slow for large n and k.
    Note: online lookup is limited to n < 100, k <= 25, t <= 8.

--tolerance [float]  (default: 1.0)  [covering strategy only]
    Pruning leeway multiplier. A k value is only pruned if:
        estimated_total_time > best_time_so_far * tolerance
    Must be >= 1.0. E.g. --tolerance 1.2 allows up to 20% slack before
    pruning, which helps avoid premature pruning due to random timing spikes.

Interactive (during run)
------------------------
s   Skip the current k value (covering) or abort the naive run immediately.
    Does not update best time.
"""

# -------------------------
# Utility functions
# -------------------------

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
    print(f"Fetching covering C({v},{k},{t}) from La Jolla Covering Repository...")

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

# Fetch covering local

import os
import csv


def fetch_covering_local(v: int, k: int, t: int, base_dir: str) -> list[list[int]]:
    """
    Load a covering design C(v,k,t) from a local CSV file.

    Expected file format:
        {base_dir}/C_{v}_k_{t}/C{v}_{k}_{t}_actual_cover.csv

    Returns:
        blocks: list of blocks, where each block is a list[int]
    """

    path = os.path.join(
        base_dir,
        f"C_{v}_k_{t}",
        f"C{v}_{k}_{t}_actual_cover.csv"
    )

    print(f"Loading covering C({v},{k},{t}) from local file: {path}")

    if not os.path.exists(path):
        raise FileNotFoundError(f"No local covering file found for C({v},{k},{t}): {path}")

    blocks = []

    with open(path, newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            block_str = row["block"].strip()
            block = [int(x) for x in block_str.split()]
            blocks.append(block)

    if not blocks:
        raise RuntimeError(f"No covering blocks found in file for C({v},{k},{t})")

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


# Unified covering loader
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

def run_for_k_value_monomer(k, all_monomers, n, t, min_total_time, log, fallback_greedy, tolerance):
    """
    Covering design over monomers: each block is a k-subset of monomer indices.
    Runs Normaliz on each subset of monomers directly.
    """
    print(f"\n{'='*70}")
    print(f"[monomer mode] Testing k={k}")
    print(f"{'='*70}\n")

    try:
        if n == 57 and t == 5 and k != 57:
            blocks = list(fetch_covering_local(n, k, t, base_dir="/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/pareto_optimal_set_cascade8"))
        else:
            blocks = list(load_covering_blocks(n, k, t, fallback_greedy=fallback_greedy))
    except (RuntimeError, requests.RequestException) as e:
        print(f"Skipping k={k}: {e}")
        log.write(f"\nk={k}: SKIPPED (covering unavailable: {e})\n")
        log.flush()
        return None, min_total_time

    num_subsets = len(blocks)
    print(f"Using covering design C({n},{k},{t}): {num_subsets} blocks\n")

    all_hilbert_vectors: set[tuple[int, ...]] = set()
    times = []
    probe_size = min(PROBE_LIMIT, num_subsets)
    wall_start = time.time()

    for phase, index_range in [("PROBE", range(probe_size)),
                                ("FULL",  range(probe_size, num_subsets))]:

        if phase == "FULL":
            probe_time = sum(times)
            avg_time = probe_time / probe_size if probe_size else 0
            estimated_total = avg_time * num_subsets
            print(f"Probe complete: avg={avg_time:.3f}s, estimated total={estimated_total:.2f}s")
            # write estimated total to log for reference
            log.write(f"\nk={k}: PROBE complete\n"
                      f"  Probe time: {probe_time:.3f}s for {probe_size} subsets\n"
                      f"  Estimated total time for full run: {estimated_total:.2f}s\n")

            if estimated_total > min_total_time * tolerance:
                print(f"PRUNED k={k}: estimated {estimated_total:.2f}s > best {min_total_time:.2f}s * tolerance {tolerance}")
                log.write(f"\nk={k}: PRUNED\n"
                          f"  Probe time: {probe_time:.3f}s for {probe_size} subsets\n"
                          f"  Estimated: {estimated_total:.2f}s vs best {min_total_time:.2f}s * tolerance {tolerance}\n")
                log.flush()
                return None, min_total_time

        for idx in index_range:
            if check_and_clear_skip():
                print(f"Skipping k={k} by user request.")
                log.write(f"\nk={k}: SKIPPED by user\n")
                log.flush()
                return None, min_total_time

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

    return _finish_run(k, num_subsets, times, all_hilbert_vectors, log, wall_start)


# -------------------------
# Core pipeline: DOMAIN MODE
# -------------------------

def run_for_k_value_domain(k, all_monomers, all_domains, n_domains, n_monomers,
                            t, min_total_time, log, fallback_greedy, tolerance):
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
        return None, min_total_time

    num_subsets = len(blocks)
    print(f"Using covering design C({n_domains},{k},{t}): {num_subsets} blocks\n")

    all_hilbert_vectors: set[tuple[int, ...]] = set()
    times = []
    probe_size = min(PROBE_LIMIT, num_subsets)
    wall_start = time.time()

    for phase, index_range in [("PROBE", range(probe_size)),
                                ("FULL",  range(probe_size, num_subsets))]:

        if phase == "FULL":
            probe_time = sum(times)
            avg_time = probe_time / probe_size if probe_size else 0
            estimated_total = avg_time * num_subsets
            print(f"Probe complete: avg={avg_time:.3f}s, estimated total={estimated_total:.2f}s")

            if estimated_total > min_total_time * tolerance:
                print(f"PRUNED k={k}: estimated {estimated_total:.2f}s > best {min_total_time:.2f}s * tolerance {tolerance}")
                log.write(f"\nk={k}: PRUNED\n"
                          f"  Probe time: {probe_time:.3f}s for {probe_size} subsets\n"
                          f"  Estimated: {estimated_total:.2f}s vs best {min_total_time:.2f}s * tolerance {tolerance}\n")
                log.flush()
                return None, min_total_time

        for idx in index_range:
            if check_and_clear_skip():
                print(f"Skipping k={k} by user request.")
                log.write(f"\nk={k}: SKIPPED by user\n")
                log.flush()
                return None, min_total_time

            block = blocks[idx]
            cleanup_normaliz_files()

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

    return _finish_run(k, num_subsets, times, all_hilbert_vectors, log, wall_start)


# -------------------------
# Core pipeline: NAIVE MONOMER MODE
# -------------------------

def run_naive_monomer(k, all_monomers, n, log):
    """
    Naive exhaustive enumeration over monomers: runs Normaliz on every
    k-subset of the n monomers. Recovers P*_k exactly (equivalent to the
    covering strategy with t=k, which requires every k-subset to be its
    own block). No probe-and-prune; all C(n,k) subsets are always run.
    """
    from math import comb
    num_subsets = comb(n, k)

    print(f"\n{'='*70}")
    print(f"[naive / monomer mode] k={k}, running all {num_subsets} subsets")
    print(f"{'='*70}\n")

    log.write(f"\nnaive monomer k={k}: {num_subsets} subsets\n")
    log.flush()

    all_hilbert_vectors: set[tuple[int, ...]] = set()
    times = []
    wall_start = time.time()

    for idx, block in enumerate(itertools.combinations(range(n), k)):
        if check_and_clear_skip():
            print("Naive run aborted by user.")
            log.write("  ABORTED by user.\n")
            log.flush()
            return None, float("inf")

        cleanup_normaliz_files()
        selected_indices = list(block)
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

        if (idx + 1) % 50 == 0 or (idx + 1) == num_subsets:
            print(f"  Progress: {idx + 1}/{num_subsets} subsets done "
                  f"({sum(times):.1f}s elapsed, {len(all_hilbert_vectors)} vectors so far)")

    return _finish_run(k, num_subsets, times, all_hilbert_vectors, log, wall_start)


# -------------------------
# Core pipeline: NAIVE DOMAIN MODE
# -------------------------

def run_naive_domain(k, all_monomers, all_domains, n_domains, n_monomers, log):
    """
    Naive exhaustive enumeration over domain types: runs Normaliz on every
    k-subset of the n_domains domain types. Monomers are filtered per block
    to those whose domains all lie within the selected subset.
    Recovers P*_k (domain-parameterised) exactly with no covering approximation.
    """
    from math import comb
    num_subsets = comb(n_domains, k)

    print(f"\n{'='*70}")
    print(f"[naive / domain mode] k={k}, running all {num_subsets} domain subsets")
    print(f"{'='*70}\n")

    log.write(f"\nnaive domain k={k}: {num_subsets} subsets\n")
    log.flush()

    all_hilbert_vectors: set[tuple[int, ...]] = set()
    times = []
    wall_start = time.time()

    for idx, block in enumerate(itertools.combinations(range(n_domains), k)):
        if check_and_clear_skip():
            print("Naive run aborted by user.")
            log.write("  ABORTED by user.\n")
            log.flush()
            return None, float("inf")

        cleanup_normaliz_files()
        selected_domains = [all_domains[i] for i in block]
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

        if (idx + 1) % 50 == 0 or (idx + 1) == num_subsets:
            print(f"  Progress: {idx + 1}/{num_subsets} subsets done "
                  f"({sum(times):.1f}s elapsed, {len(all_hilbert_vectors)} vectors so far)")

    return _finish_run(k, num_subsets, times, all_hilbert_vectors, log, wall_start)



def _finish_run(k, num_subsets, times, all_hilbert_vectors, log, wall_start):
    wall_time = time.time() - wall_start
    normaliz_time = sum(times)
    avg_normaliz = (normaliz_time / num_subsets) if num_subsets else 0.0
    overhead = wall_time - normaliz_time  # file I/O, subprocess launch, result parsing

    log.write(f"\nk={k}: FULL RUN COMPLETE\n"
              f"  Number of subsets: {num_subsets}\n"
              f"  Total wall time: {wall_time:.3f}s\n"
              f"  Total Normaliz time: {normaliz_time:.3f}s\n"
              f"  Overhead (I/O + subprocess): {overhead:.3f}s\n"
              f"  Avg Normaliz time per subset: {avg_normaliz:.3f}s\n"
              f"  Unique Hilbert basis vectors: {len(all_hilbert_vectors)}\n")
    log.flush()

    print(f"\nk = {k}: wall={wall_time:.3f}s  normaliz={normaliz_time:.3f}s  "
          f"overhead={overhead:.3f}s  avg/subset={avg_normaliz:.3f}s  "
          f"unique vectors={len(all_hilbert_vectors)}")

    return {
        "k": k,
        "num_subsets": num_subsets,
        "total_wall_time": wall_time,
        "total_normaliz_time": normaliz_time,
        "overhead_time": overhead,
        "avg_normaliz_time_per_subset": avg_normaliz,
        "unique_vectors": len(all_hilbert_vectors),
        "vectors": all_hilbert_vectors
    }, normaliz_time


# -------------------------
# Flag validation
# -------------------------

def validate_args(args, n_monomers, n_domains):
    """
    Validate all flag combinations and emit descriptive error messages.
    n is the size of the covering design universe (monomers or domains).
    """
    errors = []
    warnings = []
    n = n_monomers if args.mode == "monomer" else n_domains
    universe_label = "monomers" if args.mode == "monomer" else "unique domain types"

    # --k-start: required in both strategies
    if args.k_start < 1:
        errors.append(
            f"--k-start must be at least 1 (got {args.k_start})."
        )
    if args.k_start > n:
        errors.append(
            f"--k-start cannot exceed n={n} (got {args.k_start}), the number of\n"
            f"  {universe_label} in the input.\n"
            f"  Each block is a k-subset of these {n} elements, so k <= n is required.\n"
            f"  Use --k-start <= {n}."
        )

    if args.strategy == "naive":
        # Flags that are silently ignored in naive mode — warn the user
        if args.t != 5:  # non-default value was explicitly set
            warnings.append(
                f"--t={args.t} has no effect in naive mode. Naive enumeration runs\n"
                f"  Normaliz on every k-subset (C(n,k) total), which is equivalent\n"
                f"  to t=k. --t is only used by the covering strategy."
            )
        if args.tolerance != 1.0:
            warnings.append(
                f"--tolerance={args.tolerance} has no effect in naive mode.\n"
                f"  Probe-and-prune is only used by the covering strategy."
            )
        if args.include_base:
            warnings.append(
                f"--include-base has no effect in naive mode. To run on all\n"
                f"  {universe_label}, set --k-start {n}."
            )
        if args.fallback_greedy:
            warnings.append(
                f"--fallback-greedy has no effect in naive mode (no covering\n"
                f"  design is fetched or computed)."
            )

    else:  # covering strategy
        # --t
        if args.t < 1:
            errors.append(
                f"--t must be at least 1 (got {args.t}).\n"
                f"  t is the covering strength: every t-element subset of the\n"
                f"  {universe_label} is guaranteed to appear in at least one block."
            )

        # --k-start vs t
        if args.k_start <= args.t:
            errors.append(
                f"--k-start must be strictly greater than t (got k-start={args.k_start}, t={args.t}).\n"
                f"  The algorithm sweeps k from k-start down to t+1; if k-start <= t there\n"
                f"  are no k values to try. Use --k-start >= {args.t + 1}."
            )

        # --tolerance
        if args.tolerance < 1.0:
            errors.append(
                f"--tolerance must be >= 1.0 (got {args.tolerance}).\n"
                f"  A value of 1.0 means a k value is pruned as soon as its estimated\n"
                f"  total time exceeds the current best; values > 1.0 add proportional\n"
                f"  slack (e.g. 1.2 allows 20%% extra time before pruning).\n"
                f"  Values below 1.0 would prune the current best itself."
            )

        # Online repository limits
        if not args.fallback_greedy:
            if n >= 100 or args.k_start > 25 or (args.t > 8 and args.t != 12): # t=12 is a special case stored for n=21 in the local directory
                issues = []
                if n >= 100:
                    issues.append(f"n={n} >= 100")
                if args.k_start > 25:
                    issues.append(f"k-start={args.k_start} > 25")
                if args.t > 8:
                    issues.append(f"t={args.t} > 8")
                errors.append(
                    f"The La Jolla Covering Repository only stores designs for\n"
                    f"  n < 100, k <= 25, t <= 8, but your parameters include: {', '.join(issues)}.\n"
                    f"  Add --fallback-greedy to compute covering designs locally when\n"
                    f"  the online lookup fails."
                )

    if warnings:
        print("\nWarning(s):")
        for msg in warnings:
            print(f"  [!] {msg}\n")

    if errors:
        print("\nError: invalid flag combination(s):\n")
        for i, msg in enumerate(errors, 1):
            print(f"  [{i}] {msg}\n")
        sys.exit(1)



# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Pareto-Optimal Polymer Enumeration via Hilbert Basis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Covering strategy (default): sweep k=5..25, t=3\n"
            "  python covering_pipeline.py --t 3 --k-start 25\n\n"
            "  # Covering strategy, domain mode with 20%% pruning slack\n"
            "  python covering_pipeline.py --t 5 --k-start 20 --mode domain --tolerance 1.2\n\n"
            "  # Covering strategy, include full-system base case\n"
            "  python covering_pipeline.py --t 4 --k-start 30 --fallback-greedy --include-base\n\n"
            "  # Naive enumeration: run Normaliz on all C(n,8) monomer subsets\n"
            "  python covering_pipeline.py --strategy naive --k-start 8\n\n"
            "  # Naive enumeration, domain mode\n"
            "  python covering_pipeline.py --strategy naive --k-start 6 --mode domain\n\n"
            "During a run, type 's' and press Enter to skip the current k (covering)\n"
            "or abort the naive run immediately."
        )
    )

    parser.add_argument(
        "--strategy",
        choices=["covering", "naive"],
        default="covering",
        help=(
            "'covering' (default): use a covering design C(n, k, t) so that only\n"
            "a small number of k-subsets need to be run; sweeps k with probe-and-prune.\n"
            "'naive': run Normaliz on every k-subset of the universe (C(n,k) total)\n"
            "for a single fixed k given by --k-start. Recovers P*_k exactly with\n"
            "no approximation from the covering design. --t, --tolerance,\n"
            "--include-base, and --fallback-greedy are ignored in naive mode."
        )
    )

    parser.add_argument(
        "--t",
        type=int,
        default=5,
        metavar="INT",
        help=(
            "Covering strength: every t-element subset of monomers (monomer mode) "
            "or domain types (domain mode) is covered by at least one block. "
            "Must satisfy 1 <= t < k-start. Default: 5."
        )
    )
    parser.add_argument(
        "--k-start",
        type=int,
        default=25,
        metavar="INT",
        dest="k_start",
        help=(
            "Largest block size k to try. The algorithm sweeps k from k-start "
            "down to t+1, pruning values estimated to be slower than the current best. "
            "Must satisfy t < k-start <= n, where n is the number of monomers "
            "(monomer mode) or unique domain types (domain mode). Default: 25."
        )
    )
    parser.add_argument(
        "--include-base",
        action="store_true",
        dest="include_base",
        help=(
            "Also run the base case k=n (all monomers or all domain types), "
            "which recovers the exact full Hilbert basis P* at the cost of "
            "the longest possible runtime."
        )
    )
    parser.add_argument(
        "--fallback-greedy",
        action="store_true",
        dest="fallback_greedy",
        help=(
            "If a covering design C(n,k,t) is not found in the La Jolla Covering "
            "Repository (which only stores n < 100, k <= 25, t <= 8), compute one "
            "locally using a greedy set-cover algorithm. Can be slow for large n and k."
        )
    )
    parser.add_argument(
        "--mode",
        choices=["monomer", "domain"],
        default="monomer",
        help=(
            "'monomer' (default): covering design over monomer subsets; each block "
            "is a k-subset of monomers and Normaliz is run on those monomers directly. "
            "'domain': covering design over binding-site types; each block selects k "
            "domain types and monomers are filtered to those whose domains all lie "
            "within the selected set. Domain mode may be faster when many monomers "
            "share the same domain types (e.g. cascade or binary-tree TBNs)."
        )
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0,
        metavar="FLOAT",
        help=(
            "Pruning tolerance multiplier (>= 1.0). A candidate k is pruned only if "
            "its estimated total runtime exceeds best_so_far * tolerance. "
            "1.0 (default) prunes as soon as the estimate exceeds the best; "
            "1.2 allows 20%% slack, reducing the risk of premature pruning due to "
            "timing variance at the cost of exploring more k values."
        )
    )

    args = parser.parse_args()

    cleanup_normaliz_files()

    # Load monomers
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

    # Validate after loading input so n is known
    validate_args(args, n_monomers, n_domains)

    n_for_covering = n_monomers if args.mode == "monomer" else n_domains

    os.makedirs("logs", exist_ok=True)

    start_input_listener()

    # -------------------------
    # NAIVE STRATEGY
    # -------------------------
    if args.strategy == "naive":
        from math import comb
        k = args.k_start
        num_subsets = comb(n_for_covering, k)
        log_file = f"logs/log_naive_{n_for_covering}_k{k}.txt"

        print(f"Strategy: naive  |  k={k}  |  {num_subsets} subsets to run")
        print(f"Mode: {args.mode}\n")

        with open(log_file, "a") as log:
            log.write(
                f"Naive Enumeration — Pareto-Optimal Polymer Enumeration\n"
                f"Started: {datetime.now()}\n"
                f"Mode: {args.mode}\n"
                f"n_monomers={n_monomers}, n_domains={n_domains}\n"
                f"k={k}, total subsets=C({n_for_covering},{k})={num_subsets}\n"
                + "=" * 70 + "\n"
            )

            try:
                if args.mode == "monomer":
                    result, total_time = run_naive_monomer(
                        k, all_monomers, n_monomers, log
                    )
                else:
                    result, total_time = run_naive_domain(
                        k, all_monomers, all_domains, n_domains, n_monomers, log
                    )
            except KeyboardInterrupt:
                print("\nInterrupted by user.")
                log.write("\nRun interrupted by user (KeyboardInterrupt).\n")
                log.flush()
                cleanup_normaliz_files()
                return

            if result is not None:
                print(f"\nNaive k={k}: {result['total_normaliz_time']:.2f}s total, "
                    f"{result['unique_vectors']} unique Pareto-optimal polymers found")

                if save:
                    output_path = os.path.join(
                        save_dir,
                        f"hilbert_{n_for_covering}_k{k}_{args.t}_naive.txt"
                    )

                    save_polymer_vectors(
                        result["vectors"],
                        output_path,
                        n_monomers=n_monomers,
                        comment=f"naive mode, k={k}"
                    )

        cleanup_normaliz_files()
        return

    # -------------------------
    # COVERING STRATEGY
    # -------------------------
    k_values = list(range(args.k_start, args.t, -1))

    # also add 26, 27, 28, 29, 30, 35, 40, 45, 50 for t=5 since those are available in the local directory for n=57
    if args.t == 5 and n_for_covering == 57:
        extra_k = [26, 27, 28, 29, 30, 35, 40, 45, 50]
        k_values = sorted(set(k_values + extra_k), reverse=True)

    if args.include_base and n_for_covering not in k_values:
        k_values = [n_for_covering] + k_values

    log_file = f"logs/log_{n_for_covering}_{args.k_start}_{args.t}.txt"

    print(f"Strategy: covering  |  sweeping k in {k_values}")
    print(f"t={args.t}, tolerance={args.tolerance}, fallback-greedy={args.fallback_greedy}\n")

    min_total_time = float("inf")
    results = []

    with open(log_file, "a") as log:
        log.write(
            f"Covering Design Strategy — Pareto-Optimal Polymer Enumeration\n"
            f"Started: {datetime.now()}\n"
            f"Strategy: covering\n"
            f"Mode: {args.mode}\n"
            f"n_monomers={n_monomers}, n_domains={n_domains}\n"
            f"t={args.t}, k_start={args.k_start}, include_base={args.include_base}\n"
            f"tolerance={args.tolerance}, fallback_greedy={args.fallback_greedy}\n"
            f"k values: {k_values}\n"
            + "=" * 70 + "\n"
        )

        try:
            for k in k_values:
                if args.mode == "monomer":
                    result, new_time = run_for_k_value_monomer(
                        k, all_monomers, n_monomers, args.t,
                        min_total_time, log, args.fallback_greedy, args.tolerance
                    )
                else:
                    result, new_time = run_for_k_value_domain(
                        k, all_monomers, all_domains, n_domains, n_monomers,
                        args.t, min_total_time, log, args.fallback_greedy, args.tolerance
                    )

                if new_time < min_total_time:
                    min_total_time = new_time

                if result is None:
                    continue

                results.append(result)

                if save:
                    output_path = os.path.join(
                        save_dir,
                        f"hilbert_{n_for_covering}_k{result['k']}_t{args.t}_covering.txt"
                    )

                    save_polymer_vectors(
                        result["vectors"],
                        output_path,
                        n_monomers=n_monomers,
                        comment=f"covering mode, k={result['k']}, t={args.t}"
                    )

                if new_time < min_total_time:
                    min_total_time = new_time

        except KeyboardInterrupt:
            print("\nInterrupted by user.")
            log.write("\nRun interrupted by user (KeyboardInterrupt).\n")
            log.flush()

    if results:
        best = min(results, key=lambda r: r["total_normaliz_time"])
        print(f"\nBest k = {best['k']}  |  "
              f"wall={best['total_wall_time']:.2f}s  "
              f"normaliz={best['total_normaliz_time']:.2f}s  "
              f"overhead={best['overhead_time']:.2f}s  "
              f"({best['unique_vectors']} unique Pareto-optimal polymers found)")
    else:
        print("\nNo k value completed fully. "
              "Try increasing --tolerance, lowering --k-start, or adding --fallback-greedy.")

    cleanup_normaliz_files()

if __name__ == "__main__":
    main()