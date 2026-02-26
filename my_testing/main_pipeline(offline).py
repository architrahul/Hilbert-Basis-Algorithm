import subprocess
import time
from datetime import datetime
import os
import sqlite3
import sys
import select


# -------------------------
# Config
# -------------------------
monomer_file = "monomers.txt"
python_script = "monomers_to_normaliz.py"
normaliz_exe = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz"
log_file = "vary_k_ljcr_log.txt"
db_file = "coverings.db"

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
# Covering DB access
# -------------------------

def covering_exists(cur, v, k, t):
    cur.execute(
        "SELECT 1 FROM coverings WHERE v=? AND k=? AND t=? LIMIT 1",
        (v, k, t)
    )
    return cur.fetchone() is not None


def load_covering_blocks_with_base(cur, v, k, t):
    if k == v:
        yield list(range(1, v + 1))
        return

    cur.execute(
        "SELECT block FROM coverings WHERE v=? AND k=? AND t=?",
        (v, k, t)
    )
    for (block_str,) in cur:
        yield [int(x) for x in block_str.split(",")]

# -------------------------
# Core pipeline with pruning
# -------------------------

def run_for_k_value(k, all_monomers, cur, n, t, min_normaliz_time, log):
    print(f"\n{'='*70}")
    print(f"Testing k={k}")
    print(f"{'='*70}\n")

    blocks = list(load_covering_blocks_with_base(cur, n, k, t))
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
        return None, min_normaliz_time  # keep min time unchanged

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

    # -------------------------
    # Logging full run stats
    # -------------------------
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

    # ---- Read monomers with parsing rules ----
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

    t = int(input("Enter value of t (default=5): ") or 5)
    k_start = int(input(f"Enter starting k (≤ {n}): ").strip())

    include_base = input(f"Include base case k={n}? (y/n, default=n): ").strip().lower() or "n"

    k_values = list(range(k_start, t, -1))

    if include_base == "y" and n not in k_values:
        k_values = [n] + k_values


    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

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
                if k != n and not covering_exists(cur, n, k, t):
                    print(f"Skipping k={k}: no covering in DB")
                    continue

                result, min_total_time = run_for_k_value(k, all_monomers, cur, n, t, min_total_time, log)

                if result is None:
                    continue

                results.append(result)

        except KeyboardInterrupt:
            log.flush()


    conn.close()

    if results:
        best = min(results, key=lambda r: r["total_normaliz_time"])
        print(f"\nOptimal k = {best['k']} ({best['total_normaliz_time']:.2f}s)")
    else:
        print("\nNo k completed fully.")


if __name__ == "__main__":
    main()
