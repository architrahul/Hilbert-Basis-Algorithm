import subprocess
import time
from itertools import combinations
from datetime import datetime
import os
import math

# -------------------------
# Config
# -------------------------
monomer_file = "monomers.txt"
k = 5  # Fixed: size of each subset
n = 21  # Number of monomers (from damien.tbn)
t_values = [1, 2, 3, 4, 5]  # Different t values to test

python_script = "monomers_to_normaliz.py"
normaliz_exe = '/Users/archit/Desktop/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz'
log_file = "covering_design_log.txt"  # Log file path
hilbert_basis_file = "all_hilbert_basis.txt"

# Temporary files for each loop
tmp_monomers = "tmp_monomers.txt"
tmp_vectors = "vectors.txt"
tmp_eqs = "eqs.in"


def cleanup_normaliz_files():
    """Remove all Normaliz output files from previous runs"""
    patterns = ["eqs.out", "eqs.gen", "eqs.inv", "eqs.cst", "eqs.typ", "eqs.egn"]
    for pattern in patterns:
        if os.path.exists(pattern):
            os.remove(pattern)


def greedy_t_covering(n, k, t):
    """
    Generate a t-covering design C(n, k, t) using a greedy algorithm.

    Returns a list of k-subsets such that every t-subset is covered by at least one k-subset.

    Args:
        n: Total number of elements (monomers)
        k: Size of each subset
        t: Coverage parameter (every t-subset must be covered)

    Returns:
        List of k-subsets (each subset is a tuple of indices)
    """
    # Generate all t-subsets that need to be covered
    all_t_subsets = set(combinations(range(n), t))
    uncovered = all_t_subsets.copy()

    # Selected k-subsets
    covering_design = []

    # Greedy algorithm: repeatedly select the k-subset that covers the most uncovered t-subsets
    while uncovered:
        best_k_subset = None
        best_coverage_count = 0

        # Try all possible k-subsets
        for k_subset in combinations(range(n), k):
            # Count how many uncovered t-subsets this k-subset covers
            covered_by_this = 0
            for t_subset in uncovered:
                # Check if t_subset is contained in k_subset
                if set(t_subset).issubset(set(k_subset)):
                    covered_by_this += 1

            if covered_by_this > best_coverage_count:
                best_coverage_count = covered_by_this
                best_k_subset = k_subset

        # Add the best k-subset to our design
        if best_k_subset is None:
            break

        covering_design.append(best_k_subset)

        # Remove all t-subsets covered by this k-subset
        newly_covered = set()
        for t_subset in uncovered:
            if set(t_subset).issubset(set(best_k_subset)):
                newly_covered.add(t_subset)

        uncovered -= newly_covered

    return covering_design


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


def run_covering_design_for_t(t, all_monomers, n):
    """Run the covering design algorithm for a specific value of t"""
    print(f"\n{'='*70}")
    print(f"Testing t={t}")
    print(f"{'='*70}\n")

    # Generate the t-covering design
    print(f"Generating {t}-covering design for n={n}, k={k}...")
    start_gen = time.time()
    covering_subsets = greedy_t_covering(n, k, t)
    gen_time = time.time() - start_gen

    num_subsets = len(covering_subsets)
    all_possible = math.comb(n, k)

    print(f"Generated {num_subsets} subsets in {gen_time:.3f} sec")
    print(f"(vs. {all_possible} total C({n},{k}) combinations)")
    print(f"Reduction: {100 * (1 - num_subsets/all_possible):.1f}%\n")

    # Store all unique Hilbert basis vectors
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

        # Run Normaliz dual algorithm
        start_time = time.time()
        result = subprocess.run([normaliz_exe, "-d", "-N", tmp_eqs],
                              capture_output=True, text=True)
        elapsed = time.time() - start_time

        times.append(elapsed)

        # Read the Hilbert basis output from Normaliz
        reduced_vectors = read_hilbert_basis("eqs")

        num_new_vectors = 0
        for reduced_vec in reduced_vectors:
            full_vec = expand_vector_to_full_space(reduced_vec, selected_indices, k_current, n)
            if full_vec not in all_hilbert_vectors:
                all_hilbert_vectors.add(full_vec)
                num_new_vectors += 1

        print(f"  Subset {combo_idx+1}/{num_subsets}: {elapsed:.3f} sec, "
              f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, "
              f"total unique: {len(all_hilbert_vectors)}")

    # Calculate statistics
    total_time = sum(times)
    avg_time = total_time / len(times) if times else 0
    min_time = min(times) if times else 0
    max_time = max(times) if times else 0

    result_summary = {
        't': t,
        'num_subsets': num_subsets,
        'generation_time': gen_time,
        'total_normaliz_time': total_time,
        'total_time': gen_time + total_time,
        'avg_time_per_subset': avg_time,
        'min_time': min_time,
        'max_time': max_time,
        'unique_vectors': len(all_hilbert_vectors),
        'reduction_percent': 100 * (1 - num_subsets/all_possible)
    }

    return result_summary, all_hilbert_vectors


def main():
    cleanup_normaliz_files()

    # Read all monomers
    with open(monomer_file, 'r') as f:
        all_monomers = [line.strip() for line in f if line.strip()]

    actual_n = len(all_monomers)

    if actual_n != n:
        print(f"WARNING: Expected {n} monomers but found {actual_n}.")
        print(f"Using actual n={actual_n}")
        global n
        n = actual_n

    print(f"Total monomers: {n}")
    print(f"Fixed k: {k}")
    print(f"Testing t values: {t_values}")
    print(f"Total possible combinations: {math.comb(n, k)}")

    # Run for each t value
    results = []

    with open(log_file, 'w') as log:
        log.write(f"Covering Design Benchmark\n")
        log.write(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write(f"n={n}, k={k}\n")
        log.write(f"Testing t values: {t_values}\n")
        log.write("="*70 + "\n\n")

        for t in t_values:
            print(f"\n{'#'*70}")
            print(f"# TESTING t={t}")
            print(f"{'#'*70}")

            result, hilbert_vectors = run_covering_design_for_t(t, all_monomers, n)
            results.append(result)

            # Log results
            log.write(f"\nt={t} Results:\n")
            log.write(f"  Subsets generated: {result['num_subsets']}\n")
            log.write(f"  Generation time: {result['generation_time']:.3f} sec\n")
            log.write(f"  Total Normaliz time: {result['total_normaliz_time']:.3f} sec\n")
            log.write(f"  Total time: {result['total_time']:.3f} sec\n")
            log.write(f"  Avg time per subset: {result['avg_time_per_subset']:.3f} sec\n")
            log.write(f"  Unique Hilbert basis vectors: {result['unique_vectors']}\n")
            log.write(f"  Reduction: {result['reduction_percent']:.1f}%\n")
            log.write("-"*70 + "\n")
            log.flush()

        # Find optimal t
        optimal = min(results, key=lambda x: x['total_time'])

        summary = f"""
{'='*70}
SUMMARY - Optimal t Value
{'='*70}
"""
        for r in results:
            marker = " <-- OPTIMAL" if r['t'] == optimal['t'] else ""
            summary += f"t={r['t']}: {r['total_time']:.3f} sec total ({r['num_subsets']} subsets, {r['unique_vectors']} vectors){marker}\n"

        summary += f"""
{'='*70}
Best t={optimal['t']} with total time {optimal['total_time']:.3f} sec
  - {optimal['num_subsets']} subsets (vs. {math.comb(n, k)} total)
  - {optimal['reduction_percent']:.1f}% reduction in subsets
  - {optimal['unique_vectors']} unique Hilbert basis vectors found
{'='*70}
Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

        log.write(summary)
        print(summary)


if __name__ == "__main__":
    main()
