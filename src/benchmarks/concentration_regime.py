"""Shared helpers for the 10 nM / 100 nM leakage-concentration regime."""
from __future__ import annotations

import csv
import importlib.util
import os
import subprocess
import time
from pathlib import Path
from typing import Sequence

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from coffee_parser import assign_domain_energies, find_input_signal_indices, generate_coffee_inputs, modify_output, parse_monomers, read_polymers
from paths import EXAMPLE_TBNS_DIR, REPO_ROOT

RESULTS_NEW = Path(REPO_ROOT) / "results"
OLD_HB_DIR = RESULTS_NEW / "common" / "hilbert_basis"
COMMON_COFFEE = RESULTS_NEW / "common" / "coffee"
COMMON_CONCENTRATIONS = RESULTS_NEW / "common" / "concentrations"
COMMON_MONOMERS = RESULTS_NEW / "common" / "monomer_files"
COFFEE_CLI = Path(REPO_ROOT) / "coffee" / "crates" / "coffee-cli" / "target" / "release" / "coffee-cli"
CONCENTRATIONS_EXE = Path(REPO_ROOT) / "concentrations"

INPUT_CONC = 1e-8       # 10 nM
NON_INPUT_CONC = 1e-7   # 100 nM
EXTRA_CONC = 5e-8       # 50 nM in regime B
DEFAULT_BOND_ENERGY = -10.0


def energy_tag(energy: float) -> str:
    sign = "m" if energy < 0 else "p"
    mag = str(abs(float(energy))).rstrip("0").rstrip(".").replace(".", "p")
    return f"e_{sign}{mag}"


def _load_tbn_builder():
    path = Path(REPO_ROOT) / "example-tbns" / "tbn_builder.py"
    spec = importlib.util.spec_from_file_location("tbn_builder", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ensure_incomplete_cascade_m(m: int) -> Path:
    COMMON_MONOMERS.mkdir(parents=True, exist_ok=True)
    out = COMMON_MONOMERS / f"monomers_cascade_m{m}_incomplete_x1.txt"
    if out.exists():
        return out
    builder = _load_tbn_builder()
    monomers = builder.generate_cascade(m)[1:]
    if len(monomers) != 8 * m:
        raise RuntimeError(f"expected 8m={8*m} monomers for incomplete m={m}, got {len(monomers)}")
    with out.open("w") as f:
        f.write(f"# Linear cascade m={m}, first input x1 removed\n")
        for mon in monomers:
            f.write(mon + "\n")
    return out


def monomer_file_n7_removed_input(removed_input: int) -> Path:
    if removed_input == 1:
        return Path(EXAMPLE_TBNS_DIR) / "monomers_cascade_n7_incomplete.txt"
    return Path(EXAMPLE_TBNS_DIR) / f"monomers_cascade_n7_missing_input{removed_input}.txt"


def concentration_profile(monomers: list[list[str]], regime: str = "A") -> list[float]:
    """Initial concentrations for regime A/B.

    A: input signal monomers at 10 nM; all non-input monomers at 100 nM.
    B: same as A plus 50 nM of every monomer except the final output monomer.
    """
    if regime not in {"A", "B"}:
        raise ValueError("regime must be A or B")
    input_idxs = set(find_input_signal_indices(monomers))
    conc = [INPUT_CONC if i in input_idxs else NON_INPUT_CONC for i in range(len(monomers))]
    if regime == "B":
        final_idx = len(monomers) - 1
        for i in range(len(conc)):
            if i != final_idx:
                conc[i] += EXTRA_CONC
    return conc


def domain_energy_for(monomers: list[list[str]], bond_energy: float) -> dict[str, float]:
    energies = assign_domain_energies(monomers, seed=42)
    for key in list(energies):
        energies[key] = float(bond_energy)
    return energies


def run_coffee(monomer_file: Path, hb_path: Path, out_dir: Path, *, regime: str, bond_energy: float = DEFAULT_BOND_ENERGY, force: bool = False) -> tuple[Path, Path, float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    coffee_output = out_dir / "coffee_output.txt"
    sorted_csv = out_dir / "coffee_output_sorted.csv"
    if coffee_output.exists() and sorted_csv.exists() and not force:
        return coffee_output, sorted_csv, 0.0
    if not COFFEE_CLI.is_file():
        raise RuntimeError(f"coffee-cli not built at {COFFEE_CLI}")
    monomers = parse_monomers(str(monomer_file))
    generate_coffee_inputs(
        monomers=monomers,
        polymer_file=str(hb_path),
        out_dir=str(out_dir),
        domain_energy=domain_energy_for(monomers, bond_energy),
        concentrations=concentration_profile(monomers, regime),
        label=f"{out_dir.name}_{regime}_{energy_tag(bond_energy)}",
    )
    t0 = time.time()
    subprocess.run([str(COFFEE_CLI), str(out_dir / "input.ocx"), str(out_dir / "input.con"), "-o", str(coffee_output)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
    dt = time.time() - t0
    modify_output(str(coffee_output), str(hb_path), str(sorted_csv))
    return coffee_output, sorted_csv, dt


def run_concentrations(monomer_file: Path, hb_path: Path, out_dir: Path, *, regime: str, bond_energy: float = DEFAULT_BOND_ENERGY, force_inputs: bool = False, force_solver: bool = False, exe: Path = CONCENTRATIONS_EXE) -> tuple[Path, Path, float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    ocx = out_dir / "input.ocx"
    con = out_dir / "input.con"
    if force_inputs or not (ocx.exists() and con.exists()):
        monomers = parse_monomers(str(monomer_file))
        generate_coffee_inputs(
            monomers=monomers,
            polymer_file=str(hb_path),
            out_dir=str(out_dir),
            domain_energy=domain_energy_for(monomers, bond_energy),
            concentrations=concentration_profile(monomers, regime),
            label=f"conc_{out_dir.name}_{regime}_{energy_tag(bond_energy)}",
        )
    if not exe.is_file():
        raise RuntimeError(f"concentrations executable not found at {exe}")
    eq = out_dir / "input.eq"
    if eq.exists() and not force_solver:
        dt = 0.0
    else:
        t0 = time.time()
        with (out_dir / "concentrations_stdout.txt").open("w") as so, (out_dir / "concentrations_stderr.txt").open("w") as se:
            subprocess.run([str(exe), "input"], cwd=str(out_dir), stdout=so, stderr=se, check=True)
        dt = time.time() - t0
    rows = parse_eq(eq, len(parse_monomers(str(monomer_file))))
    sorted_csv = out_dir / "concentrations_sorted.csv"
    write_sorted_rows(sorted_csv, rows)
    return eq, sorted_csv, dt


def parse_eq(eq_path: Path, n_monomers: int) -> list[dict]:
    rows = []
    with eq_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) < n_monomers + 4:
                continue
            try:
                vec = tuple(int(x) for x in parts[2:2+n_monomers])
                energy = float(parts[2+n_monomers])
                conc = float(parts[3+n_monomers])
            except ValueError:
                continue
            rows.append({"concentration_M": conc, "polymer_vector": vec, "complex_energy": energy})
    if not rows:
        raise RuntimeError(f"could not parse concentrations from {eq_path}")
    rows.sort(key=lambda r: -r["concentration_M"])
    return rows


def write_sorted_rows(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["concentration_M", "polymer_vector", "complex_energy"])
        for row in rows:
            w.writerow([f"{row['concentration_M']:.12e}", " ".join(str(x) for x in row["polymer_vector"]), f"{row.get('complex_energy', 0.0):.12e}"])


def read_sorted_conc_csv(path: Path) -> list[dict]:
    rows=[]
    with path.open() as f:
        r=csv.DictReader(f)
        for row in r:
            rows.append({"concentration_M": float(row["concentration_M"]), "polymer_vector": tuple(int(x) for x in row["polymer_vector"].split())})
    return rows


def final_output_metrics_from_csv(sorted_csv: Path) -> tuple[float, float]:
    rows = read_sorted_conc_csv(sorted_csv)
    if not rows:
        return 0.0, 0.0
    n = len(rows[0]["polymer_vector"])
    final_idx = n - 1
    singleton = tuple(1 if i == final_idx else 0 for i in range(n))
    free = 0.0
    total = 0.0
    for row in rows:
        vec = row["polymer_vector"]
        c = row["concentration_M"]
        if vec == singleton:
            free = c
        if vec[final_idx] > 0:
            total += c * vec[final_idx]
    return free, total


def hb_exp5_path(m: int, method: str) -> Path:
    label = f"linear_cascade_m{m}_incomplete_x1"
    if method == "full_hb":
        return OLD_HB_DIR / f"exp5_full_hb_{label}.txt"
    return OLD_HB_DIR / f"exp5_{label}_k25_t{method.split('t')[-1]}.txt"


def hb_exp7_path(removed_input: int) -> Path:
    if removed_input == 1:
        canonical = OLD_HB_DIR / "hilbert_k25_t5_monomer_n7_incomplete.txt"
        exp7 = OLD_HB_DIR / "exp7_hilbert_k25_t5_monomer_n7_incomplete.txt"
        return canonical if canonical.exists() else exp7
    return OLD_HB_DIR / f"exp7_hilbert_k25_t5_monomer_n7_missing_input{removed_input}.txt"

# ---------------------------------------------------------------------------
# Correct-output candidate sets (full cascade, X1 present)
# ---------------------------------------------------------------------------

def ensure_full_cascade_m(m: int) -> Path:
    COMMON_MONOMERS.mkdir(parents=True, exist_ok=True)
    if m == 7:
        p = Path(EXAMPLE_TBNS_DIR) / "monomers_cascade_n7.txt"
        if p.exists():
            return p
    out = COMMON_MONOMERS / f"monomers_cascade_m{m}_full.txt"
    if out.exists():
        return out
    builder = _load_tbn_builder()
    monomers = builder.generate_cascade(m)
    if len(monomers) != 8 * m + 1:
        raise RuntimeError(f"expected 8m+1={8*m+1} monomers for full m={m}, got {len(monomers)}")
    with out.open("w") as f:
        f.write(f"# Linear cascade m={m}, full/correct-output system with x1 present\n")
        for mon in monomers:
            f.write(mon + "\n")
    return out


def _saved_full_covering_hb(m: int, t: int) -> Path | None:
    candidates = [
        OLD_HB_DIR / f"hilbert_k25_t{t}_monomer_n{m}_full.txt",
        OLD_HB_DIR / f"hilbert_k25_t{t}_monomer_linear_cascade_m{m}_full.txt",
        OLD_HB_DIR / f"exp5_linear_cascade_m{m}_full_k25_t{t}.txt",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def ensure_full_covering_hb(m: int, *, preferred_t: int = 5, fallback_t: int = 4, cap_s: float = 1800.0, force: bool = False, logs: bool = False) -> tuple[Path, dict]:
    """Return/generate a covering HB for the full cascade (X1 present).

    Try t=5 first. If the probe-projected runtime is over cap_s, use t=4.
    This intentionally computes a covering candidate set, not Full HB.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from hilbert_pipeline import cleanup_normaliz_files, full_run_k, get_all_unique_domains, load_covering_blocks, load_monomers, probe_k, run_normaliz_on_subset, save_polymer_vectors

    monomer_file = ensure_full_cascade_m(m)
    all_monomers = load_monomers(str(monomer_file))
    n_monomers = len(all_monomers)
    domains = get_all_unique_domains(all_monomers)
    meta_dir = RESULTS_NEW / "common" / "hilbert_basis" / "metadata"
    out_dir = RESULTS_NEW / "common" / "hilbert_basis"
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    # For small full cascades, k=25 is at least the full monomer set size.
    # The covering consists of a single full-system block, so do not query LJCR
    # for impossible objects such as C(9,25,5).
    if 25 >= n_monomers:
        t = min(preferred_t, n_monomers)
        out = out_dir / f"linear_cascade_m{m}_full_single_block.txt"
        meta_path = meta_dir / f"linear_cascade_m{m}_full_single_block.json"
        if out.exists() and meta_path.exists() and not force:
            import json
            with meta_path.open() as f:
                return out, json.load(f)
        print(f"[full/correct HB] m={m}, |M|={n_monomers} <= k=25: running single full-system block")
        wall0 = time.time()
        elapsed, raw = run_normaliz_on_subset(all_monomers)
        wall = time.time() - wall0
        if not raw:
            raise RuntimeError(f"Normaliz returned no vectors for full cascade m={m}")
        vecs = {tuple(v[:n_monomers]) for v in raw}
        vecs.discard(tuple([0] * n_monomers))
        save_polymer_vectors(vecs, str(out), n_monomers=n_monomers, comment=f"new regime correct-output full cascade m={m}, single full-system block")
        meta = {"m": m, "t": t, "k": n_monomers, "source": "single_full_block", "hb_path": str(out), "hb_runtime_s": wall, "hb_normaliz_s": elapsed, "num_vectors": len(vecs)}
        import json
        with meta_path.open("w") as f:
            json.dump(meta, f, indent=2)
        cleanup_normaliz_files()
        return out, meta

    for t in [preferred_t, fallback_t]:
        saved = _saved_full_covering_hb(m, t)
        if saved is not None and not force:
            return saved, {"m": m, "t": t, "source": "saved", "hb_path": str(saved)}

        out = out_dir / f"linear_cascade_m{m}_full_k25_t{t}.txt"
        meta_path = meta_dir / f"linear_cascade_m{m}_full_k25_t{t}.json"
        if out.exists() and meta_path.exists() and not force:
            import json
            with meta_path.open() as f:
                return out, json.load(f)

        print(f"[full/correct HB] m={m}, k=25, t={t}: building/probing covering design")
        cover0 = time.time()
        blocks = load_covering_blocks(n_monomers, 25, t, fallback_dp=True)
        covering_s = time.time() - cover0
        with open(os.devnull, "w") as log:
            probe0 = time.time()
            projected, probe_times, num_blocks = probe_k(25, t, blocks, n_monomers, all_monomers, "monomer", domains, n_monomers, best_projected=None, log=log)
            probe_s = time.time() - probe0
        if projected is None:
            raise RuntimeError(f"probe failed for full cascade m={m}, t={t}")
        if t == preferred_t and projected > cap_s:
            print(f"  t={t} projected {projected:.1f}s > {cap_s:.1f}s; switching to t={fallback_t}")
            continue
        print(f"  t={t} projected {projected:.1f}s; running covering enumeration")
        with open(os.devnull, "w") as log:
            result, _ = full_run_k(25, blocks, n_monomers, all_monomers, "monomer", domains, n_monomers, log)
        if result is None:
            raise RuntimeError(f"covering enumeration failed for full cascade m={m}, t={t}")
        save_polymer_vectors(result["vectors"], str(out), n_monomers=n_monomers, comment=f"new regime correct-output full cascade m={m}, k=25, t={t}")
        meta = {
            "m": m, "t": t, "k": 25, "source": "computed_covering", "hb_path": str(out),
            "projected_runtime_s": projected, "probe_wall_s": probe_s, "probe_normaliz_s": sum(probe_times),
            "covering_construction_s": covering_s, "num_blocks": num_blocks,
            "hb_runtime_s": result["total_wall_time"], "hb_normaliz_s": result["total_normaliz_time"], "num_vectors": result["unique_vectors"],
        }
        import json
        with meta_path.open("w") as f:
            json.dump(meta, f, indent=2)
        cleanup_normaliz_files()
        return out, meta
    raise RuntimeError(f"could not produce full/correct-output covering HB for m={m}")


def ensure_complete_full_hb(m: int, *, force: bool = False) -> tuple[Path, dict]:
    """Return/generate Full-HB P* for the complete cascade with X1 present."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from hilbert_pipeline import cleanup_normaliz_files, get_all_unique_domains, load_monomers, run_normaliz_on_subset, save_polymer_vectors

    # Experiment 2 stores complete linear-cascade Full-HB baselines as size{m}.
    saved = OLD_HB_DIR / f"exp2_full_hb_linear_cascade_size{m}.txt"
    if saved.exists() and not force:
        return saved, {"m": m, "source": "saved_exp2_full_hb", "hb_path": str(saved)}

    monomer_file = ensure_full_cascade_m(m)
    all_monomers = load_monomers(str(monomer_file))
    n_monomers = len(all_monomers)
    out_dir = RESULTS_NEW / "common" / "hilbert_basis"
    meta_dir = out_dir / "metadata"
    out_dir.mkdir(parents=True, exist_ok=True); meta_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"linear_cascade_m{m}_full_full_hb.txt"
    meta_path = meta_dir / f"linear_cascade_m{m}_full_full_hb.json"
    if out.exists() and meta_path.exists() and not force:
        import json
        with meta_path.open() as f:
            return out, json.load(f)
    print(f"[full/correct Full HB] m={m}, |M|={n_monomers}: running Normaliz on complete system")
    wall0 = time.time()
    elapsed, raw = run_normaliz_on_subset(all_monomers)
    wall = time.time() - wall0
    if not raw:
        raise RuntimeError(f"Normaliz returned no Full-HB vectors for complete cascade m={m}")
    vecs = {tuple(v[:n_monomers]) for v in raw}
    vecs.discard(tuple([0] * n_monomers))
    save_polymer_vectors(vecs, str(out), n_monomers=n_monomers, comment=f"new regime correct-output complete cascade m={m}, Full HB")
    meta = {"m": m, "source": "computed_full_hb", "hb_path": str(out), "hb_runtime_s": wall, "hb_normaliz_s": elapsed, "num_vectors": len(vecs)}
    import json
    with meta_path.open("w") as f:
        json.dump(meta, f, indent=2)
    cleanup_normaliz_files()
    return out, meta


def ensure_full_covering_hb_t(m: int, t: int, *, force: bool = False) -> tuple[Path, dict]:
    """Return/generate complete-cascade covering HB at exact t."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from hilbert_pipeline import cleanup_normaliz_files, full_run_k, get_all_unique_domains, load_covering_blocks, load_monomers, save_polymer_vectors

    monomer_file = ensure_full_cascade_m(m)
    all_monomers = load_monomers(str(monomer_file))
    n_monomers = len(all_monomers)
    domains = get_all_unique_domains(all_monomers)
    out_dir = RESULTS_NEW / "common" / "hilbert_basis"
    meta_dir = out_dir / "metadata"
    out_dir.mkdir(parents=True, exist_ok=True); meta_dir.mkdir(parents=True, exist_ok=True)

    saved = _saved_full_covering_hb(m, t)
    if saved is not None and not force:
        return saved, {"m": m, "t": t, "source": "saved", "hb_path": str(saved)}

    if 25 >= n_monomers:
        # Exact t is irrelevant when the single full-system block is used.
        return ensure_full_covering_hb(m, preferred_t=t, fallback_t=t, force=force)

    out = out_dir / f"linear_cascade_m{m}_full_k25_t{t}.txt"
    meta_path = meta_dir / f"linear_cascade_m{m}_full_k25_t{t}.json"
    if out.exists() and meta_path.exists() and not force:
        import json
        with meta_path.open() as f:
            return out, json.load(f)

    print(f"[full/correct covering HB] m={m}, k=25, t={t}: running exact-t covering enumeration")
    cover0 = time.time()
    blocks = load_covering_blocks(n_monomers, 25, t, fallback_dp=True)
    covering_s = time.time() - cover0
    with open(os.devnull, "w") as log:
        result, _ = full_run_k(25, blocks, n_monomers, all_monomers, "monomer", domains, n_monomers, log)
    if result is None:
        raise RuntimeError(f"covering enumeration failed for complete cascade m={m}, t={t}")
    save_polymer_vectors(result["vectors"], str(out), n_monomers=n_monomers, comment=f"new regime correct-output complete cascade m={m}, k=25, t={t}")
    meta = {"m": m, "t": t, "k": 25, "source": "computed_covering_exact_t", "hb_path": str(out), "covering_construction_s": covering_s, "hb_runtime_s": result["total_wall_time"], "hb_normaliz_s": result["total_normaliz_time"], "num_vectors": result["unique_vectors"]}
    import json
    with meta_path.open("w") as f:
        json.dump(meta, f, indent=2)
    cleanup_normaliz_files()
    return out, meta
