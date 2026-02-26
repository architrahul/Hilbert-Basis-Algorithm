import subprocess
import time
from itertools import combinations
from datetime import datetime
import os
from itertools import combinations, islice
import math

# -------------------------
# Config
# -------------------------
monomer_file = "monomers.txt" 
k = 7
max_loops = 100
python_script = "monomers_to_normaliz.py"
normaliz_exe = '/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz'
log_file = "benchmark_log.txt"  # Log file path
hilbert_basis_file = "all_hilbert_basis.txt" 

# Temporary files for each loop
tmp_monomers = "tmp_monomers.txt"
tmp_vectors = "vectors.txt"
tmp_eqs = "eqs.in"

# Clean up any previous Normaliz output files
def cleanup_normaliz_files():
    """Remove all Normaliz output files from previous runs"""
    patterns = ["eqs.out", "eqs.gen", "eqs.inv", "eqs.cst", "eqs.typ", "eqs.egn"]
    for pattern in patterns:
        if os.path.exists(pattern):
            os.remove(pattern)
            print(f"Cleaned up: {pattern}")

# Clean up at the start
cleanup_normaliz_files()

# Read all monomers
with open(monomer_file, 'r') as f:
    all_monomers = [line.strip() for line in f if line.strip()]

n = len(all_monomers)
if k > n:
    raise ValueError(f"k={k} is larger than total number of monomers={n}")

# Generate all combinations
total_combinations = math.comb(n, k)
loops_to_run = min(total_combinations, max_loops)
combinations_iter = combinations(range(n), k)

# Store all unique Hilbert basis vectors (in original n-dimensional space)
all_hilbert_vectors = set()

def read_hilbert_basis(base_filename="eqs"):
    """
    Read the Hilbert basis from Normaliz output file (eqs.out).
    Returns a list of vectors.
    """
    output_file = f"{base_filename}.out"
    
    if not os.path.exists(output_file):
        print(f"Warning: {output_file} not found")
        return []
    
    vectors = []
    
    try:
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        # Find the line that says "X Hilbert basis elements:"
        hilbert_start = None
        for i, line in enumerate(lines):
            if "Hilbert basis elements:" in line:
                hilbert_start = i + 1
                break
        
        if hilbert_start is None:
            print(f"Warning: Could not find 'Hilbert basis elements:' in {output_file}")
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
    """
    Expand a k-dimensional vector to n-dimensional space.
    The reduced vector includes both monomers and singleton vectors.
    We only care about the monomer components (first k_current elements).
    
    Args:
        reduced_vector: Vector that includes k_current monomers + singleton vectors
        selected_indices: List of k_current indices that were selected (0-indexed)
        k_current: Number of monomers in this combination
        n: Total number of monomers in original space
    Returns:
        Vector in n-dimensional space with zeros in unselected positions
    """
    # Only take the first k_current elements (corresponding to the actual monomers)
    monomer_part = reduced_vector[:k_current]
    
    full_vector = [0] * n
    
    for i, idx in enumerate(selected_indices):
        if i < len(monomer_part):
            full_vector[idx] = monomer_part[i]
    
    return tuple(full_vector)


with open(log_file, 'w') as log:
    log.write(f"Benchmark started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log.write(f"Total monomers: {n}\n")
    log.write(f"k value: {k}\n")
    log.write(f"Total combinations to test: {total_combinations} (C({n},{k}))\n")
    log.write("="*70 + "\n\n")
    log.flush()
    
    print(f"Total combinations to test: {total_combinations} (C({n},{k}))")
    print(f"Logging to: {log_file}")
    
    times = []
    
    for combo_idx, selected_indices in enumerate(islice(combinations_iter, loops_to_run)):
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
        
        if combo_idx == 0:
            print(f"\nDEBUG - First combination:")
            print(f"  Selected indices: {selected_indices}")
            print(f"  k_current: {k_current}")
            print(f"  Found {len(reduced_vectors)} Hilbert basis vectors")
            if reduced_vectors:
                print(f"  First vector length: {len(reduced_vectors[0])}")
                print(f"  First 5 vectors (full):")
                for i in range(min(5, len(reduced_vectors))):
                    vec = reduced_vectors[i]
                    print(f"    Vector {i}: {vec}")
                    expanded = expand_vector_to_full_space(vec, selected_indices, k_current, n)
                    print(f"    Expanded: {expanded}")
                    print()
        
        num_new_vectors = 0
        for reduced_vec in reduced_vectors:
            full_vec = expand_vector_to_full_space(reduced_vec, selected_indices, k_current, n)
            if full_vec not in all_hilbert_vectors:
                all_hilbert_vectors.add(full_vec)
                num_new_vectors += 1
        
        log_line = f"Combination {combo_idx+1}/{total_combinations}: {elapsed:.3f} sec, " \
                   f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, " \
                   f"total unique: {len(all_hilbert_vectors)}\n"
        log.write(log_line)
        log.flush()
        
        print(f"Combination {combo_idx+1}/{total_combinations}: {elapsed:.3f} sec, "
              f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, "
              f"total unique: {len(all_hilbert_vectors)}")
    
    if times:
        avg_time = sum(times)/len(times)
        min_time = min(times)
        max_time = max(times)
    else:
        avg_time = min_time = max_time = 0.0
    
    summary = f"""
{'='*70}
Statistics for k={k} (all {total_combinations} combinations):
Average time: {avg_time:.3f} sec
Min time: {min_time:.3f} sec
Max time: {max_time:.3f} sec
Total time: {sum(times):.3f} sec
Total unique Hilbert basis vectors found: {len(all_hilbert_vectors)}
{'='*70}
Benchmark completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""
    
    log.write(summary)
    print(summary)

with open(hilbert_basis_file, 'w') as f:
    f.write(f"# Total unique Hilbert basis vectors: {len(all_hilbert_vectors)}\n")
    f.write(f"# Format: space-separated vector in {n}-dimensional space\n")
    f.write(f"# Computed using k={k} combinations\n\n")
    
    for vec in sorted(all_hilbert_vectors):
        f.write(' '.join(map(str, vec)) + '\n')

print(f"\nAll unique Hilbert basis vectors written to: {hilbert_basis_file}")