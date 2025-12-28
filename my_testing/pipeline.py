import subprocess
import time
from itertools import combinations
from datetime import datetime
import os

# -------------------------
# Config
# -------------------------
monomer_file = "monomers.txt"   # Your full monomers file
k = 16                         # Number of monomers to pick per combination
python_script = "monomers_to_normaliz.py"  # Path to your conversion script
normaliz_exe = '/Users/archit/Desktop/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz'
log_file = "benchmark_log.txt"  # Log file path
hilbert_basis_file = "all_hilbert_basis.txt"  # File to store all unique Hilbert basis vectors

# Temporary files for each loop
tmp_monomers = "tmp_monomers.txt"
tmp_vectors = "vectors.txt"
tmp_eqs = "eqs.in"

# -------------------------
# Read all monomers
# -------------------------
with open(monomer_file, 'r') as f:
    all_monomers = [line.strip() for line in f if line.strip()]

n = len(all_monomers)
if k > n:
    raise ValueError(f"k={k} is larger than total number of monomers={n}")

# -------------------------
# Generate all combinations
# -------------------------
all_combinations = list(combinations(range(n), k))  # Store indices instead of monomers
total_combinations = len(all_combinations)

# Store all unique Hilbert basis vectors (in original n-dimensional space)
all_hilbert_vectors = set()

# -------------------------
# Helper function to read Normaliz output
# -------------------------
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
        
        # Read vectors until we hit a blank line or another section
        for i in range(hilbert_start, len(lines)):
            line = lines[i].strip()
            
            # Stop at empty line or next section
            if not line or line.startswith("***") or "extreme rays" in line.lower():
                break
            
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Parse the vector
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
    
    # Create full vector with zeros
    full_vector = [0] * n
    
    # Place the monomer values at their correct positions
    for i, idx in enumerate(selected_indices):
        if i < len(monomer_part):
            full_vector[idx] = monomer_part[i]
    
    return tuple(full_vector)

# -------------------------
# Open log file
# -------------------------
with open(log_file, 'w') as log:
    # Write header
    log.write(f"Benchmark started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log.write(f"Total monomers: {n}\n")
    log.write(f"k value: {k}\n")
    log.write(f"Total combinations to test: {total_combinations} (C({n},{k}))\n")
    log.write("="*70 + "\n\n")
    log.flush()
    
    print(f"Total combinations to test: {total_combinations} (C({n},{k}))")
    print(f"Logging to: {log_file}")
    
    # -------------------------
    # Run benchmarking loop
    # -------------------------
    times = []
    
    for combo_idx, selected_indices in enumerate(all_combinations):
        # Get the selected monomers
        subset = [all_monomers[i] for i in selected_indices]
        k_current = len(selected_indices)  # This should always equal k, but be explicit
        
        # Write temporary monomers file
        with open(tmp_monomers, 'w') as f:
            for m in subset:
                f.write(m + "\n")
        
        # Run the conversion script to generate eqs.in
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
        
        # Debug: Print first combination details
if combo_idx == 0:
    print(f"\nDEBUG - First combination:")
    print(f"  Selected indices: {selected_indices}")
    print(f"  k_current: {k_current}")
    print(f"  Found {len(reduced_vectors)} Hilbert basis vectors")
    if reduced_vectors:
        print(f"  First vector length: {len(reduced_vectors[0])}")
        print(f"  First 5 vectors:")
        for i in range(min(5, len(reduced_vectors))):
            vec = reduced_vectors[i]
            print(f"    Vector {i}: {vec}")
            expanded = expand_vector_to_full_space(vec, selected_indices, k_current, n)
            print(f"    Expanded: {expanded}")
        
        # Expand each vector to full space and add to set
        num_new_vectors = 0
        for reduced_vec in reduced_vectors:
            full_vec = expand_vector_to_full_space(reduced_vec, selected_indices, k_current, n)
            if full_vec not in all_hilbert_vectors:
                all_hilbert_vectors.add(full_vec)
                num_new_vectors += 1
        
        # Log to file
        log_line = f"Combination {combo_idx+1}/{total_combinations}: {elapsed:.3f} sec, " \
                   f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, " \
                   f"total unique: {len(all_hilbert_vectors)}\n"
        log.write(log_line)
        log.flush()
        
        # Print to console
        print(f"Combination {combo_idx+1}/{total_combinations}: {elapsed:.3f} sec, "
              f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, "
              f"total unique: {len(all_hilbert_vectors)}")
    
    # -------------------------
    # Summary
    # -------------------------
    avg_time = sum(times)/len(times)
    min_time = min(times)
    max_time = max(times)
    
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

# -------------------------
# Write all Hilbert basis vectors to file
# -------------------------
with open(hilbert_basis_file, 'w') as f:
    f.write(f"# Total unique Hilbert basis vectors: {len(all_hilbert_vectors)}\n")
    f.write(f"# Format: space-separated vector in {n}-dimensional space\n")
    f.write(f"# Computed using k={k} combinations\n\n")
    
    for vec in sorted(all_hilbert_vectors):  # Sort for consistent output
        f.write(' '.join(map(str, vec)) + '\n')

print(f"\nAll unique Hilbert basis vectors written to: {hilbert_basis_file}")