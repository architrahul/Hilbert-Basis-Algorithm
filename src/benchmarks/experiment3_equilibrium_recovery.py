#!/usr/bin/env python3
"""Experiment 3 under the new leakage concentration regime.

Uses saved Hilbert bases for the X1=0 leak case. For the X1=10 nM
correct-output case, it uses a saved full-system covering HB if available,
or computes a covering HB at t=5; if the projected t=5 runtime exceeds
1800s, it uses t=4.
  A: inputs X2..X{k+1}=10 nM, non-input monomers=100 nM
  B: A plus 50 nM of all monomers except the final output monomer
Bond energy defaults to -10 kcal/mol.
"""
from __future__ import annotations
import argparse, csv, json, os, sys
from datetime import datetime
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "results" / ".matplotlib"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import concentration_regime as nr
from paths import EXAMPLE_TBNS_DIR

BENCH_DIR = nr.RESULTS_NEW / "benchmarks" / "03_equilibrium_recovery"
CSV_DIR = BENCH_DIR / "csv"; FIG_DIR = BENCH_DIR / "figures"
SETS = {
    "full_pstar": "hilbert_full_p_star_n7_incomplete.txt",
    "k25_t3": "hilbert_k25_t3_monomer_n7_incomplete.txt",
    "k25_t5": "hilbert_k25_t5_monomer_n7_incomplete.txt",
}

CUTOFF_M = 1e-10

def energy_label(energy: float) -> str:
    mag = str(abs(float(energy))).rstrip("0").rstrip(".").replace(".", "p")
    return f"E{mag}"

def scenario_name(regime: str, case: str, energy: float) -> str:
    kind = "leak" if case == "x1_0_leak" else "correct"
    return f"regime{regime}_{kind}_{energy_label(energy)}"

def concentration_rank_plot(sorted_csv: Path, out_png: Path, title: str) -> None:
    rows = load_rows(sorted_csv)
    concs = sorted([c for c, _v in rows if c > 0], reverse=True)
    x = np.arange(len(concs))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(x, concs, color="#2ca02c", lw=1.6)
    ax.set_yscale("log")
    ax.set_xlabel("polymer index (decreasing concentration)")
    ax.set_ylabel("concentration (M)")
    ax.set_title(title)
    ax.grid(True, which="both", ls=":", alpha=.4)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(out_png, dpi=160); plt.close(fig)

def load_rows(path: Path):
    rows=[]
    with path.open() as f:
        r=csv.DictReader(f)
        for row in r:
            rows.append((float(row["concentration_M"]), tuple(int(x) for x in row["polymer_vector"].split())))
    return rows

def relative_error(full_csv: Path, red_csv: Path, out_csv: Path, out_png: Path, title: str):
    full=load_rows(full_csv); red=dict((v,c) for c,v in load_rows(red_csv))
    sig=[(c,v) for c,v in full if c>=CUTOFF_M]
    sig.sort(reverse=True, key=lambda x:x[0])
    xs=[]; errs=[]; present=[]; full_concs=[]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w=csv.writer(f); w.writerow(["rank","full_pstar_concentration_M","reduced_concentration_M","relative_error","present","polymer_vector"])
        for i,(c,v) in enumerate(sig):
            rc=red.get(v)
            ok=rc is not None
            err=float("inf") if rc is None else abs(rc-c)/c
            w.writerow([i,f"{c:.12e}","" if rc is None else f"{rc:.12e}","Inf" if not np.isfinite(err) else f"{err:.12e}",ok," ".join(map(str,v))])
            xs.append(i); errs.append(err); present.append(ok); full_concs.append(c)
    x=np.arange(len(sig)); errs=np.array(errs); present=np.array(present, dtype=bool); full_concs=np.array(full_concs)
    fig, ax=plt.subplots(figsize=(7,4.5)); ax2=ax.twinx()
    ax2.plot(x, full_concs, color="#2ca02c", alpha=.32, lw=1.8, label="full P* concentration"); ax2.set_yscale("log"); ax2.set_ylabel("full P* concentration (M)", color="#2ca02c"); ax2.tick_params(axis="y", labelcolor="#2ca02c")
    finite=present & np.isfinite(errs)
    ax.scatter(x[finite], np.clip(errs[finite], 1e-16, None), s=14, label="present in reduced set")
    absent=~present
    if absent.any():
        ymax=max(np.nanmax(errs[finite]) if finite.any() else 1.0,1.0)*10
        for xi in x[absent]: ax.plot([xi,xi],[CUTOFF_M,ymax], color="#d62728", lw=.8)
        ax.plot([],[], color="#d62728", lw=2, label=f"absent ({absent.sum()})")
    ax.set_yscale("log"); ax.set_xlabel("polymer index (decreasing full-P* concentration)"); ax.set_ylabel("relative concentration error"); ax.set_title(title); ax.grid(True, which="both", ls=":", alpha=.4)
    h,l=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels(); ax.legend(h+h2,l+l2, fontsize=8)
    fig.tight_layout(); out_png.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out_png, dpi=160); plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bond-energy", type=float, default=nr.DEFAULT_BOND_ENERGY)
    ap.add_argument("--regimes", nargs="+", default=["A","B"], choices=["A","B"])
    ap.add_argument("--force", action="store_true", help="rerun COFFEE even when cached")
    args=ap.parse_args()
    CSV_DIR.mkdir(parents=True, exist_ok=True); FIG_DIR.mkdir(parents=True, exist_ok=True)
    leak_monomer_file=Path(EXAMPLE_TBNS_DIR)/"monomers_cascade_n7_incomplete.txt"
    correct_monomer_file=nr.ensure_full_cascade_m(7)
    rows=[]
    for regime in args.regimes:
        # X1 = 0 leak case: compare full P* baseline to t=3/t=5 reduced sets.
        for set_name, hb_name in SETS.items():
            hb=nr.OLD_HB_DIR/hb_name
            if not hb.exists(): raise FileNotFoundError(f"missing saved HB: {hb}")
            out_dir=nr.COMMON_COFFEE/"03_equilibrium_recovery"/f"n7_x1_0"/f"{set_name}_regime_{regime}_{nr.energy_tag(args.bond_energy)}"
            cof, sorted_csv, coffee_s=nr.run_coffee(leak_monomer_file, hb, out_dir, regime=regime, bond_energy=args.bond_energy, force=args.force)
            free,total=nr.final_output_metrics_from_csv(sorted_csv)
            rows.append({"case":"x1_0_leak","regime":regime,"set":set_name,"bond_energy":args.bond_energy,"free_final_output_M":free,"total_final_output_containing_M":total,"hilbert_basis":str(hb),"coffee_output":str(cof),"sorted_concentration_csv":str(sorted_csv),"coffee_runtime_s":coffee_s})
        for t in [3,5]:
            relative_error(
                nr.COMMON_COFFEE/"03_equilibrium_recovery"/"n7_x1_0"/f"full_pstar_regime_{regime}_{nr.energy_tag(args.bond_energy)}"/"coffee_output_sorted.csv",
                nr.COMMON_COFFEE/"03_equilibrium_recovery"/"n7_x1_0"/f"k25_t{t}_regime_{regime}_{nr.energy_tag(args.bond_energy)}"/"coffee_output_sorted.csv",
                CSV_DIR/scenario_name(regime, "x1_0_leak", args.bond_energy)/f"relative_error_t{t}.csv",
                FIG_DIR/scenario_name(regime, "x1_0_leak", args.bond_energy)/f"relative_error_t{t}.png",
                f"Experiment 3: X1=0, regime {regime}, t={t}, E={args.bond_energy:g}")

        # X1 = 10 nM correct-output case: same treatment as leak.
        correct_sources = {
            "full_pstar": nr.ensure_complete_full_hb(7),
            "k25_t3": nr.ensure_full_covering_hb_t(7, 3),
            "k25_t5": nr.ensure_full_covering_hb_t(7, 5),
        }
        correct_csvs = {}
        scenario = scenario_name(regime, "x1_10nM_correct", args.bond_energy)
        scen_csv = CSV_DIR/scenario/"outputs.csv"
        scen_csv.parent.mkdir(parents=True, exist_ok=True)
        with scen_csv.open("w", newline="") as sf:
            sw = csv.writer(sf); sw.writerow(["set", "free_final_output_M", "total_final_output_containing_M", "sorted_concentration_csv", "hilbert_basis"])
            for correct_set, (correct_hb, correct_meta) in correct_sources.items():
                out_dir=nr.COMMON_COFFEE/"03_equilibrium_recovery"/"n7_x1_10nM"/f"{correct_set}_regime_{regime}_{nr.energy_tag(args.bond_energy)}"
                cof, sorted_csv, coffee_s=nr.run_coffee(correct_monomer_file, correct_hb, out_dir, regime=regime, bond_energy=args.bond_energy, force=args.force)
                correct_csvs[correct_set] = sorted_csv
                free,total=nr.final_output_metrics_from_csv(sorted_csv)
                rows.append({"case":"x1_10nM_correct","regime":regime,"set":correct_set,"bond_energy":args.bond_energy,"free_final_output_M":free,"total_final_output_containing_M":total,"hilbert_basis":str(correct_hb),"coffee_output":str(cof),"sorted_concentration_csv":str(sorted_csv),"coffee_runtime_s":coffee_s,"hb_t":correct_meta.get('t',''),"hb_projected_runtime_s":correct_meta.get('projected_runtime_s','')})
                sw.writerow([correct_set, f"{free:.12e}", f"{total:.12e}", str(sorted_csv), str(correct_hb)])
        for t in [3,5]:
            relative_error(
                correct_csvs["full_pstar"],
                correct_csvs[f"k25_t{t}"],
                CSV_DIR/scenario/f"relative_error_t{t}.csv",
                FIG_DIR/scenario/f"relative_error_t{t}.png",
                f"Experiment 3: X1=10 nM correct, regime {regime}, t={t}, E={args.bond_energy:g}",
            )
    with (CSV_DIR/f"summary_{nr.energy_tag(args.bond_energy)}.csv").open("w", newline="") as f:
        preferred = [
            "case", "regime", "set", "bond_energy", "free_final_output_M",
            "total_final_output_containing_M", "hb_t", "hb_projected_runtime_s",
            "hilbert_basis", "coffee_output", "sorted_concentration_csv", "coffee_runtime_s",
        ]
        extras = sorted({k for row in rows for k in row} - set(preferred))
        fields = [k for k in preferred if any(k in row for row in rows)] + extras
        w=csv.DictWriter(f, fieldnames=fields, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with (BENCH_DIR/f"summary_{nr.energy_tag(args.bond_energy)}.json").open("w") as f: json.dump({"generated":datetime.now().isoformat(),"rows":rows}, f, indent=2)
    print(f"Experiment 3 new-regime outputs: {BENCH_DIR}")
    print("Includes X1=0 leak and X1=10 nM correct-output cases. Correct-output uses covering HB t=5 unless projected over 1800s, then t=4.")
if __name__=="__main__": main()
