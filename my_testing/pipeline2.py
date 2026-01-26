import subprocess
import time
from itertools import combinations
from datetime import datetime
import os
from collections import OrderedDict

# -------------------------
# Config
# -------------------------
monomer_file = "monomers.txt"   # Your full monomers file
k = 9                          # Number of DOMAINS (binding sites) to pick per combination
python_script = "monomers_to_normaliz.py"  # Path to your conversion script
normaliz_exe = '/Users/archit/Desktop/Hilbert Basis Algorithm/my_testing/Normaliz/source/normaliz'
log_file = "benchmark_log.txt"  # Log file path
hilbert_basis_file = "all_hilbert_basis.txt"  # File to store all unique Hilbert basis vectors

# Temporary files for each loop
tmp_monomers = "tmp_monomers.txt"
tmp_vectors = "vectors.txt"
tmp_eqs = "eqs.in"

# -------------------------
# Clean up any previous Normaliz output files
# -------------------------
def cleanup_normaliz_files():
    """Remove all Normaliz output files from previous runs"""
    patterns = ["eqs.out", "eqs.gen", "eqs.inv", "eqs.cst", "eqs.typ", "eqs.egn"]
    for pattern in patterns:
        if os.path.exists(pattern):
            os.remove(pattern)

# Clean up at the start
cleanup_normaliz_files()

# -------------------------
# Read all monomers and extract unique domains
# -------------------------
with open(monomer_file, 'r') as f:
    all_monomers = [line.strip() for line in f if line.strip()]

def get_domains_from_monomer(monomer):
    """Extract unique domains from a monomer (without stars)"""
    domains = set()
    for domain in monomer.split():
        domains.add(domain.rstrip('*'))
    return domains

def get_all_unique_domains(monomers):
    """Get all unique domains across all monomers"""
    all_domains = OrderedDict()
    for monomer in monomers:
        for domain in get_domains_from_monomer(monomer):
            all_domains[domain] = None
    return list(all_domains.keys())

# Get all unique domains
all_domains = get_all_unique_domains(all_monomers)
n_domains = len(all_domains)
n_monomers = len(all_monomers)

print(f"Total monomers: {n_monomers}")
print(f"Total unique domains: {n_domains}")
print(f"Domains: {all_domains}")

if k > n_domains:
    raise ValueError(f"k={k} is larger than total number of domains={n_domains}")

# -------------------------
# Generate all domain combinations and filter monomers
# -------------------------
def filter_monomers_by_domains(monomers, selected_domains):
    """
    Return only monomers that contain ONLY the selected domains.
    A monomer is included if all its domains are in the selected set.
    """
    selected_set = set(selected_domains)
    filtered = []
    for monomer in monomers:
        monomer_domains = get_domains_from_monomer(monomer)
        # Include monomer if all its domains are in the selected set
        if monomer_domains.issubset(selected_set):
            filtered.append(monomer)
    return filtered

all_domain_combinations = list(combinations(range(n_domains), k))
total_combinations = len(all_domain_combinations)

# Store all unique Hilbert basis vectors (in original n-dimensional monomer space)
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

def expand_vector_to_full_monomer_space(reduced_vector, filtered_monomers_indices, n_monomers):
    """
    Expand a vector from the filtered monomer space to the full monomer space.
    
    Args:
        reduced_vector: Vector in the space of filtered monomers (includes monomers + singletons)
        filtered_monomers_indices: Indices of filtered monomers in the original monomer list
        n_monomers: Total number of monomers in original space
    Returns:
        Vector in full n_monomers-dimensional space with zeros for non-selected monomers
    """
    k_filtered = len(filtered_monomers_indices)
    
    # Only take the first k_filtered elements (corresponding to the actual monomers)
    monomer_part = reduced_vector[:k_filtered]
    
    # Create full vector with zeros
    full_vector = [0] * n_monomers
    
    # Place the monomer values at their correct positions
    for i, idx in enumerate(filtered_monomers_indices):
        if i < len(monomer_part):
            full_vector[idx] = monomer_part[i]
    
    return tuple(full_vector)

# -------------------------
# Open log file
# -------------------------
with open(log_file, 'w') as log:
    # Write header
    log.write(f"Benchmark started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log.write(f"Total monomers: {n_monomers}\n")
    log.write(f"Total domains: {n_domains}\n")
    log.write(f"k value (domains): {k}\n")
    log.write(f"Total combinations to test: {total_combinations} (C({n_domains},{k}))\n")
    log.write("="*70 + "\n\n")
    log.flush()
    
    print(f"\nTotal domain combinations to test: {total_combinations} (C({n_domains},{k}))")
    print(f"Logging to: {log_file}\n")
    
    # -------------------------
    # Run benchmarking loop
    # -------------------------
    times = []
    
    for combo_idx, selected_domain_indices in enumerate(all_domain_combinations):
        # Clean up Normaliz output files before each run
        cleanup_normaliz_files()
        
        # Get the selected domains
        selected_domains = [all_domains[i] for i in selected_domain_indices]
        
        # Filter monomers that only contain these domains
        filtered_monomers = filter_monomers_by_domains(all_monomers, selected_domains)
        
        # Get indices of filtered monomers in original list
        filtered_monomers_indices = [
            i for i, monomer in enumerate(all_monomers) 
            if monomer in filtered_monomers
        ]
        
        if len(filtered_monomers) == 0:
            print(f"Combination {combo_idx+1}/{total_combinations}: Skipped (no monomers with only these domains)")
            log.write(f"Combination {combo_idx+1}/{total_combinations}: Skipped (no monomers)\n")
            log.flush()
            continue
        
        # Write temporary monomers file
        with open(tmp_monomers, 'w') as f:
            for m in filtered_monomers:
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
            print(f"DEBUG - First combination:")
            print(f"  Selected domains: {selected_domains}")
            print(f"  Filtered to {len(filtered_monomers)} monomers: {filtered_monomers}")
            print(f"  Filtered monomer indices: {filtered_monomers_indices}")
            print(f"  Found {len(reduced_vectors)} Hilbert basis vectors")
            if reduced_vectors and len(reduced_vectors) > 0:
                print(f"  First vector length: {len(reduced_vectors[0])}")
                print(f"  First 3 vectors (showing first {min(len(filtered_monomers)+2, len(reduced_vectors[0]))} elements):")
                for i in range(min(3, len(reduced_vectors))):
                    vec = reduced_vectors[i]
                    print(f"    Vector {i}: {vec[:min(len(filtered_monomers)+2, len(vec))]}")
                    expanded = expand_vector_to_full_monomer_space(vec, filtered_monomers_indices, n_monomers)
                    print(f"    Expanded: {expanded}")
                print()
        
        # Expand each vector to full space and add to set
        num_new_vectors = 0
        for reduced_vec in reduced_vectors:
            full_vec = expand_vector_to_full_monomer_space(reduced_vec, filtered_monomers_indices, n_monomers)
            if full_vec not in all_hilbert_vectors:
                all_hilbert_vectors.add(full_vec)
                num_new_vectors += 1
        
        # Log to file
        log_line = f"Combination {combo_idx+1}/{total_combinations} (domains: {selected_domains}): " \
                   f"{elapsed:.3f} sec, {len(filtered_monomers)} monomers, " \
                   f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, " \
                   f"total unique: {len(all_hilbert_vectors)}\n"
        log.write(log_line)
        log.flush()
        
        # Print to console
        print(f"Combination {combo_idx+1}/{total_combinations}: {elapsed:.3f} sec, "
              f"{len(filtered_monomers)} monomers, "
              f"found {len(reduced_vectors)} vectors, {num_new_vectors} new, "
              f"total unique: {len(all_hilbert_vectors)}")
    
    # -------------------------
    # Summary
    # -------------------------
    if times:
        avg_time = sum(times)/len(times)
        min_time = min(times)
        max_time = max(times)
        total_time = sum(times)
    else:
        avg_time = min_time = max_time = total_time = 0
    
    summary = f"""
{'='*70}
Statistics for k={k} domains (all {total_combinations} combinations):
Average time: {avg_time:.3f} sec
Min time: {min_time:.3f} sec
Max time: {max_time:.3f} sec
Total time: {total_time:.3f} sec
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
    f.write(f"# Format: space-separated vector in {n_monomers}-dimensional monomer space\n")
    f.write(f"# Computed using k={k} domain combinations\n\n")
    
    for vec in sorted(all_hilbert_vectors):  # Sort for consistent output
        f.write(' '.join(map(str, vec)) + '\n')

print(f"\nAll unique Hilbert basis vectors written to: {hilbert_basis_file}")