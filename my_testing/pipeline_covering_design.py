import subprocess
import time
from itertools import combinations
from datetime import datetime
import os
import math
import random

# Config

monomer_file = "monomers.txt"
n = 41  # Number of monomers
t = 5   # Fixed: guarantee every 5-subset is covered
k_values = list(range(n, n-10, -1))  # Test k from 27 down to 18

python_script = "monomers_to_normaliz.py"
normaliz_exe = '/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz'
log_file = "vary_k_log.txt"

# Temporary files
tmp_monomers = "tmp_monomers.txt"
tmp_eqs = "eqs.in"

# For large computations, limit the number of subsets per k
MAX_SUBSETS_TO_TEST = 500  # Adjust based on your patience


def cleanup_normaliz_files():
    """Remove all Normaliz output files from previous runs"""
    patterns = ["eqs.out", "eqs.gen", "eqs.inv", "eqs.cst", "eqs.typ", "eqs.egn"]
    for pattern in patterns:
        if os.path.exists(pattern):
            os.remove(pattern)


def greedy_covering_fast(n, k, t, max_subsets=None):
    """
    Fast greedy algorithm for t-covering design.
    Uses sampling and early stopping for large problems.

    Returns a list of k-subsets that approximately cover all t-subsets.
    """
    all_t_subsets = set(combinations(range(n), t))
    uncovered = all_t_subsets.copy()
    covering_design = []

    # Pre-generate all possible k-subsets
    all_k_subsets = list(combinations(range(n), k))

    iteration = 0
    while uncovered and (max_subsets is None or len(covering_design) < max_subsets):
        iteration += 1

        if len(all_k_subsets) > 1000:
            candidates = random.sample(all_k_subsets, min(1000, len(all_k_subsets)))
        else:
            candidates = all_k_subsets

        best_k_subset = None
        best_coverage_count = 0

        for k_subset in candidates:
            k_set = set(k_subset)
            # Count how many uncovered t-subsets this k-subset covers
            covered_by_this = sum(1 for t_subset in uncovered if set(t_subset).issubset(k_set))

            if covered_by_this > best_coverage_count:
                best_coverage_count = covered_by_this
                best_k_subset = k_subset

        if best_k_subset is None or best_coverage_count == 0:
            break

        covering_design.append(best_k_subset)

        # Remove covered t-subsets
        k_set = set(best_k_subset)
        uncovered = {t_sub for t_sub in uncovered if not set(t_sub).issubset(k_set)}

        # Print progress for long computations
        if iteration % 10 == 0:
            coverage = 100 * (1 - len(uncovered) / len(all_t_subsets))
            print(f"    Progress: {len(covering_design)} subsets, {coverage:.1f}% coverage")

    return covering_design


def estimate_one_run_time(k, all_monomers, n):
    cleanup_normaliz_files()

    # Pick one random k-subset
    selected_indices = tuple(range(k))  # Just use first k for simplicity
    subset = [all_monomers[i] for i in selected_indices]

    with open(tmp_monomers, 'w') as f:
        for m in subset:
            f.write(m + "\n")

    subprocess.run(["python", python_script], check=True,
                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Run Normaliz and time it
    start_time = time.time()
    subprocess.run([normaliz_exe, "-d", "-N", tmp_eqs],
                  capture_output=True, text=True)
    elapsed = time.time() - start_time

    return elapsed


def read_hilbert_basis(base_filename="eqs"):
    """Read the Hilbert basis from Normaliz output file"""
    output_file = f"{base_filename}.out"

    if not os.path.exists(output_file):
        return []

    vectors = []

    try:
        with open(output_file, 'r') as f:
            lines = f.readlines()

        hilbert_start = None
        for i, line in enumerate(lines):
            if "Hilbert basis elements:" in line:
                hilbert_start = i + 1
                break

        if hilbert_start is None:
            return []

        for i in range(hilbert_start, len(lines)):
            line = lines[i].strip()

            if not line or line.startswith("***") or "extreme rays" in line.lower():
                break

            if line.startswith("#"):
                continue

            try:
                vector = [int(x) for x in line.split()]
                if vector:
                    vectors.append(vector)
            except ValueError:
                continue

    except Exception as e:
        print(f"Error reading {output_file}: {e}")

    return vectors


def expand_vector_to_full_space(reduced_vector, selected_indices, k_current, n):
    """Expand a k-dimensional vector to n-dimensional space"""
    monomer_part = reduced_vector[:k_current]
    full_vector = [0] * n

    for i, idx in enumerate(selected_indices):
        if i < len(monomer_part):
            full_vector[idx] = monomer_part[i]

    return tuple(full_vector)


def run_for_k_value(k, all_monomers, n, t):
    """Run the pipeline for a specific k value"""
    print(f"\n{'='*70}")
    print(f"Testing k={k}")
    print(f"{'='*70}\n")

    # First, estimate time for one run
    print(f"Estimating time per subset with k={k}...")
    est_time = estimate_one_run_time(k, all_monomers, n)
    print(f"  Estimated time per subset: {est_time:.3f} sec")

    # Generate covering design
    print(f"Generating {t}-covering design for n={n}, k={k}, t={t}...")
    start_gen = time.time()
    covering_subsets = greedy_covering_fast(n, k, t, max_subsets=MAX_SUBSETS_TO_TEST)
    gen_time = time.time() - start_gen

    num_subsets = len(covering_subsets)
    theoretical_min = math.comb(n, t) / math.comb(k, t) if k >= t else float('inf')

    print(f"Generated {num_subsets} subsets in {gen_time:.3f} sec")
    print(f"  (Theoretical minimum: ~{theoretical_min:.0f} subsets)")
    print(f"  Estimated total Normaliz time: {num_subsets * est_time:.1f} sec\n")

    # Decide whether to run all or sample
    if num_subsets > MAX_SUBSETS_TO_TEST:
        print(f"WARNING: {num_subsets} subsets is too many. Sampling {MAX_SUBSETS_TO_TEST} subsets.")
        covering_subsets = random.sample(covering_subsets, MAX_SUBSETS_TO_TEST)
        num_subsets = len(covering_subsets)

    # Run Normaliz on all subsets
    all_hilbert_vectors = set()
    times = []

    for combo_idx, selected_indices in enumerate(covering_subsets):
        cleanup_normaliz_files()

        subset = [all_monomers[i] for i in selected_indices]
        k_current = len(selected_indices)

        with open(tmp_monomers, 'w') as f:
            for m in subset:
                f.write(m + "\n")

        subprocess.run(["python", python_script], check=True,
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        start_time = time.time()
        subprocess.run([normaliz_exe, "-d", "-N", tmp_eqs],
                      capture_output=True, text=True)
        elapsed = time.time() - start_time

        times.append(elapsed)

        reduced_vectors = read_hilbert_basis("eqs")

        num_new_vectors = 0
        for reduced_vec in reduced_vectors:
            full_vec = expand_vector_to_full_space(reduced_vec, selected_indices, k_current, n)
            if full_vec not in all_hilbert_vectors:
                all_hilbert_vectors.add(full_vec)
                num_new_vectors += 1

        if (combo_idx + 1) % 10 == 0 or combo_idx == num_subsets - 1:
            print(f"  Subset {combo_idx+1}/{num_subsets}: {elapsed:.3f} sec, "
                  f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, "
                  f"total unique: {len(all_hilbert_vectors)}")

    total_time = sum(times)
    avg_time = total_time / len(times) if times else 0

    result_summary = {
        'k': k,
        'num_subsets': num_subsets,
        'generation_time': gen_time,
        'total_normaliz_time': total_time,
        'total_time': gen_time + total_time,
        'avg_time_per_subset': avg_time,
        'estimated_time_per_subset': est_time,
        'unique_vectors': len(all_hilbert_vectors),
        'theoretical_min_subsets': theoretical_min
    }

    return result_summary


def main():
    cleanup_normaliz_files()

    # Read all monomers
    with open(monomer_file, 'r') as f:
        all_monomers = [line.strip() for line in f if line.strip()]

    actual_n = len(all_monomers)

    global n
    if actual_n < n:
        print(f"WARNING: Found only {actual_n} monomers, expected {n}. Using n={actual_n}")
        n = actual_n
    else:
        # Use only first n monomers
        all_monomers = all_monomers[:n]

    print(f"Total monomers: {n}")
    print(f"Fixed t: {t}")
    print(f"Testing k values: {k_values}")
    print(f"Max subsets to test per k: {MAX_SUBSETS_TO_TEST}")

    results = []

    with open(log_file, 'w') as log:
        log.write(f"Vary k Benchmark (Fixed t={t})\n")
        log.write(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write(f"n={n}, t={t}\n")
        log.write(f"Testing k values: {k_values}\n")
        log.write("="*70 + "\n\n")

        for k in k_values:
            print(f"\n{'#'*70}")
            print(f"# TESTING k={k}")
            print(f"{'#'*70}")

            result = run_for_k_value(k, all_monomers, n, t)
            results.append(result)

            log.write(f"\nk={k} Results:\n")
            log.write(f"  Subsets generated: {result['num_subsets']}\n")
            log.write(f"  Theoretical minimum: {result['theoretical_min_subsets']:.0f}\n")
            log.write(f"  Generation time: {result['generation_time']:.3f} sec\n")
            log.write(f"  Total Normaliz time: {result['total_normaliz_time']:.3f} sec\n")
            log.write(f"  Total time: {result['total_time']:.3f} sec\n")
            log.write(f"  Avg time per subset: {result['avg_time_per_subset']:.3f} sec\n")
            log.write(f"  Unique Hilbert basis vectors: {result['unique_vectors']}\n")
            log.write("-"*70 + "\n")
            log.flush()

        # Find optimal k
        optimal = min(results, key=lambda x: x['total_time'])

        summary = f"""
{'='*70}
SUMMARY - Optimal k Value (with t={t} fixed)
{'='*70}
"""
        for r in results:
            marker = " <-- OPTIMAL" if r['k'] == optimal['k'] else ""
            summary += f"k={r['k']:2d}: {r['total_time']:7.2f} sec total ({r['num_subsets']:4d} subsets, {r['unique_vectors']:4d} vectors){marker}\n"

        summary += f"""
{'='*70}
Best k={optimal['k']} with total time {optimal['total_time']:.2f} sec
  - {optimal['num_subsets']} subsets required
  - {optimal['unique_vectors']} unique Hilbert basis vectors found
  - Avg {optimal['avg_time_per_subset']:.3f} sec per subset
{'='*70}
Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

        log.write(summary)
        print(summary)


if __name__ == "__main__":
    main()
