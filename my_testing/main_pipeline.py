import subprocess
import time
from datetime import datetime
import os
import sys
import requests
from bs4 import BeautifulSoup


# -------------------------
# Config
# -------------------------
monomer_file = "monomers.txt"
python_script = "monomers_to_normaliz.py"
normaliz_exe = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz"
log_file = "vary_k_ljcr_log.txt"

tmp_monomers = "tmp_monomers.txt"
tmp_eqs = "eqs.in"

PROBE_LIMIT = 100  # number of subsets to probe


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
    full_vector = [0] * n
    for i, idx in enumerate(selected_indices):
        if i < len(reduced_vector):
            full_vector[idx] = reduced_vector[i]
    return tuple(full_vector)


# -------------------------
# Online covering fetch
# -------------------------

def fetch_covering_online(v: int, k: int, t: int):
    """
    Fetches a C(v, k, t) covering design from ljcr.dmgordon.org.
    Returns a list of blocks (each block is a list of 1-indexed integers).
    Raises RuntimeError if no covering is available for the given parameters.
    """
    url = f"https://ljcr.dmgordon.org/show_cover.php?v={v}&k={k}&t={t}"
    print(f"  Fetching covering C({v},{k},{t}) from {url} ...")

    r = requests.get(
        url,
        headers={"User-Agent": "Mozilla/5.0"},
        timeout=30
    )
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

    print(f"  Retrieved {len(blocks)} blocks.")
    return blocks


def load_covering_blocks(v: int, k: int, t: int):
    """
    Returns blocks for a covering design. If k == v, yields the trivial
    full block. Otherwise fetches from the online source.
    Raises RuntimeError if unavailable.
    """
    if k == v:
        yield list(range(1, v + 1))
        return

    blocks = fetch_covering_online(v, k, t)
    for block in blocks:
        yield block


# -------------------------
# Core pipeline with pruning
# -------------------------

def run_for_k_value(k, all_monomers, n, t, min_normaliz_time, log):
    print(f"\n{'='*70}")
    print(f"Testing k={k}")
    print(f"{'='*70}\n")

    # Fetch covering blocks upfront — if unavailable, skip this k
    try:
        blocks = list(load_covering_blocks(n, k, t))
    except (RuntimeError, requests.RequestException) as e:
        print(f"Skipping k={k}: {e}")
        log.write(f"\nk={k}: SKIPPED (covering unavailable: {e})\n")
        log.flush()
        return None, min_normaliz_time

    num_subsets = len(blocks)
    print(f"Using covering design: {num_subsets} subsets\n")
    print("Commands during run: 's' = skip current k, 'q' = quit program")

    all_hilbert_vectors = set()
    times = []

    probe_size = min(PROBE_LIMIT, num_subsets)

    # -------------------------
    # PROBE PHASE
    # -------------------------
    for idx in range(probe_size):
        block = blocks[idx]
        cleanup_normaliz_files()

        selected_indices = [x - 1 for x in block]
        subset = [all_monomers[i] for i in selected_indices]

        with open(tmp_monomers, "w") as f:
            for m in subset:
                f.write(m + "\n")

        subprocess.run(
            ["python", python_script],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )

        start = time.time()
        subprocess.run(
            [normaliz_exe, "-d", "-N", tmp_eqs],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        elapsed = time.time() - start
        times.append(elapsed)

        reduced_vectors = read_hilbert_basis()
        for rv in reduced_vectors:
            fv = expand_vector_to_full_space(rv, selected_indices, n)
            all_hilbert_vectors.add(fv)

    probe_time = sum(times)
    avg_time = probe_time / probe_size if probe_size else 0
    estimated_total = avg_time * num_subsets

    print(f"Probe complete: avg={avg_time:.3f}s, estimated total={estimated_total:.2f}s")

    if estimated_total > min_normaliz_time:
        print(f"PRUNED k={k}: estimated {estimated_total:.2f}s > current best {min_normaliz_time:.2f}s")
        log.write(f"\nk={k}: PRUNED\n")
        log.write(f"  Probe time: {probe_time:.3f}s for {probe_size} subsets\n")
        log.write(f"  Estimated total time: {estimated_total:.2f}s vs current best {min_normaliz_time:.2f}s\n")
        log.flush()
        return None, min_normaliz_time

    # -------------------------
    # FULL RUN
    # -------------------------
    for idx in range(probe_size, num_subsets):
        block = blocks[idx]
        cleanup_normaliz_files()

        selected_indices = [x - 1 for x in block]
        subset = [all_monomers[i] for i in selected_indices]

        with open(tmp_monomers, "w") as f:
            for m in subset:
                f.write(m + "\n")

        subprocess.run(
            ["python", python_script],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )

        start = time.time()
        subprocess.run(
            [normaliz_exe, "-d", "-N", tmp_eqs],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        elapsed = time.time() - start
        times.append(elapsed)

        reduced_vectors = read_hilbert_basis()
        for rv in reduced_vectors:
            fv = expand_vector_to_full_space(rv, selected_indices, n)
            all_hilbert_vectors.add(fv)

    total_time = sum(times)

    log.write(f"\nk={k}: FULL RUN COMPLETE\n")
    log.write(f"  Number of subsets: {num_subsets}\n")
    log.write(f"  Total Normaliz time: {total_time:.3f}s\n")
    log.write(f"  Avg Normaliz time per subset: {total_time / num_subsets:.3f}s\n")
    log.write(f"  Unique Hilbert basis vectors: {len(all_hilbert_vectors)}\n")
    log.flush()

    return {
        "k": k,
        "num_subsets": num_subsets,
        "total_normaliz_time": total_time,
        "avg_normaliz_time_per_subset": total_time / num_subsets,
        "unique_vectors": len(all_hilbert_vectors)
    }, total_time


# -------------------------
# Main
# -------------------------

def main():
    cleanup_normaliz_files()

    # ---- Read monomers ----
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

    n = len(all_monomers)
    print(f"\nDetected number of monomers: n = {n}")
    print(f"Note: Online covering lookup is limited to n < 100, k ≤ 25, t ≤ 8\n")

    t = int(input("Enter value of t (default=5): ") or 5)
    k_start = int(input(f"Enter starting k (≤ {n}): ").strip())
    include_base = input(f"Include base case k={n}? (y/n, default=n): ").strip().lower() or "n"

    k_values = list(range(k_start, t, -1))

    if include_base == "y" and n not in k_values:
        k_values = [n] + k_values

    min_total_time = float("inf")
    results = []

    with open(log_file, "w") as log:
        log.write(f"Vary k Benchmark\n")
        log.write(f"Started: {datetime.now()}\n")
        log.write(f"n={n}, t={t}\n")
        log.write(f"k values: {k_values}\n")
        log.write("="*70 + "\n")

        try:
            for k in k_values:
                result, min_total_time = run_for_k_value(k, all_monomers, n, t, min_total_time, log)

                if result is None:
                    continue

                results.append(result)
                if result["total_normaliz_time"] < min_total_time:
                    min_total_time = result["total_normaliz_time"]

        except KeyboardInterrupt:
            log.flush()

    if results:
        best = min(results, key=lambda r: r["total_normaliz_time"])
        print(f"\nOptimal k = {best['k']} ({best['total_normaliz_time']:.2f}s)")
    else:
        print("\nNo k completed fully.")


if __name__ == "__main__":
    main()