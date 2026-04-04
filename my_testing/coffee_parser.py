import random
import re
import os


# -------------------------
# Parse monomer file
# -------------------------
def parse_monomers(path):
    monomers = []
    domain_pattern = re.compile(r"[a-zA-Z]+\d+_\d+\*?")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if ":" in line:
                line = line.split(":", 1)[1].strip()

            parts = line.split(",")
            domains = domain_pattern.findall(parts[0])
            monomers.append(domains)

    return monomers


# -------------------------
# Assign bond energies ONCE
# -------------------------
def assign_domain_energies(monomers, seed=42):
    random.seed(seed)   # ensures reproducibility

    domain_types = set()

    for m in monomers:
        for d in m:
            base = d.replace("*", "")
            domain_types.add(base)

    energies = {}
    for d in sorted(domain_types):
        # instead of random all of them have -20
        energies[d] = -20.0
    return energies


# -------------------------
# Count bonds in polymer
# -------------------------
def compute_polymer_energy(polymer, monomers, domain_energy):
    domain_count = {}

    for i, count in enumerate(polymer):
        if count == 0:
            continue

        for _ in range(count):
            for d in monomers[i]:
                domain_count[d] = domain_count.get(d, 0) + 1

    energy = 0.0

    # Count complementary matches s and s*
    for d in list(domain_count.keys()):
        if d.endswith("*"):
            base = d[:-1]
            if base in domain_count:
                bonds = min(domain_count[d], domain_count[base])
                energy += bonds * domain_energy[base]

    return energy


# -------------------------
# Read polymer file
# -------------------------
def read_polymers(path):
    polymers = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            vec = list(map(int, line.split()))
            polymers.append(vec)

    return polymers


# -------------------------
# Write input.ocx
# -------------------------
def write_ocx(polymers, energies, output_path):
    with open(output_path, "w") as f:
        for i, (vec, energy) in enumerate(zip(polymers, energies), start=1):
            row = [str(i), "1"] + [str(x) for x in vec] + [f"{energy:.6e}"]
            f.write("\t".join(row) + "\n")


# -------------------------
# Write input.con
# -------------------------
# def write_con(n_monomers, output_path, seed=123):
#     random.seed(seed)   # optional reproducibility for concentrations too
#     with open(output_path, "w") as f:
#         for _ in range(n_monomers):
#             val = random.uniform(1e-7, 1e-6)
#             f.write(f"{val:.6e}\n")
def write_con(n_monomers, output_path):
    """
    Writes a .con file with all monomer concentrations set to 1 μM.
    """
    with open(output_path, "w") as f:
        for _ in range(n_monomers):
            f.write("1.000000e-6\n")


# -------------------------
# Optional: save energies for reference
# -------------------------
def write_domain_energies(domain_energy, output_path):
    with open(output_path, "w") as f:
        for d in sorted(domain_energy):
            f.write(f"{d}\t{domain_energy[d]:.6f}\n")


# -------------------------
# Main pipeline
# -------------------------
def generate_coffee_inputs(monomers, polymer_file, out_dir, domain_energy):
    os.makedirs(out_dir, exist_ok=True)

    polymers = read_polymers(polymer_file)

    energies = []
    for p in polymers:
        e = compute_polymer_energy(p, monomers, domain_energy)
        energies.append(e)

    write_ocx(polymers, energies, os.path.join(out_dir, "input.ocx"))
    write_con(len(monomers), os.path.join(out_dir, "input.con"))
    write_domain_energies(domain_energy, os.path.join(out_dir, "domain_energies.txt"))

    print(f"Generated COFFEE inputs in {out_dir}")


# -------------------------
# Usage
# -------------------------
if __name__ == "__main__":

    monomer_file = "/Users/archit/Projects/Hilbert Basis Algorithm/example-tbns/cascade_n10_incomplete.txt"
    polymer_dir = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/pareto_optimal_set_cascade10"
    base_out_dir = "/Users/archit/Projects/Hilbert Basis Algorithm/my_testing/pareto_optimal_set_cascade10/coffee_input"

    # Parse monomers ONCE
    monomers = parse_monomers(monomer_file)

    # Assign energies ONCE and reuse everywhere
    domain_energy = assign_domain_energies(monomers, seed=42)

    # for (n) in [81, 80]:
    for (n) in [80]:
        polymer_file = os.path.join(polymer_dir, f"hilbert_{n}_k25_t5_covering.txt")
        out_dir = os.path.join(base_out_dir, f"n{n}_k25_t5")

        generate_coffee_inputs(
            monomers=monomers,
            polymer_file=polymer_file,
            out_dir=out_dir,
            domain_energy=domain_energy
        )