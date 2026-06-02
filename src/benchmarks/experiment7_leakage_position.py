#!/usr/bin/env python3
"""Experiment 7 under the new 10 nM/100 nM regime.

Seven-module cascade, one input absent at a time, using saved covering k=25,t=5 Hilbert bases and COFFEE. Also includes the X1=10 nM correct-output full system via covering HB t=5 unless projected over 1800s, then t=4. Bond energy defaults to -10.
"""
from __future__ import annotations
import argparse, csv, json, os, sys
from datetime import datetime
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2]/"results"/".matplotlib")); os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import concentration_regime as nr
BENCH_DIR=nr.RESULTS_NEW/"benchmarks"/"07_leakage_position"; CSV_DIR=BENCH_DIR/"csv"; FIG_DIR=BENCH_DIR/"figures"

def module_to_removed_input(module:int)->int: return 1 if module==1 else module+1

def run_cell(module, regime, energy, force):
    removed=module_to_removed_input(module)
    monomer=nr.monomer_file_n7_removed_input(removed)
    hb=nr.hb_exp7_path(removed)
    if not hb.exists():
        # For x1 the older canonical name is expected; for x3..x8 exp7 files exist.
        raise FileNotFoundError(f"missing saved HB for removed x{removed}: {hb}")
    label="n7_incomplete" if removed==1 else f"n7_missing_input{removed}"
    out_dir=nr.COMMON_COFFEE/"07_leakage_position"/label/f"k25_t5_regime_{regime}_{nr.energy_tag(energy)}"
    cof, sorted_csv, coffee_s=nr.run_coffee(monomer,hb,out_dir,regime=regime,bond_energy=energy,force=force)
    free,total=nr.final_output_metrics_from_csv(sorted_csv)
    return {"missing_input_module":module,"removed_input":removed,"case":f"x{removed}_0_leak","regime":regime,"bond_energy":energy,"final_output_concentration_M":free,"total_final_output_containing_M":total,"coffee_runtime_s":coffee_s,"hilbert_basis":str(hb),"monomer_file":str(monomer),"coffee_output":str(cof),"sorted_concentration_csv":str(sorted_csv)}

def run_correct(regime, energy, force):
    hb, meta = nr.ensure_full_covering_hb(7)
    monomer=nr.ensure_full_cascade_m(7)
    t=meta.get("t", "")
    out_dir=nr.COMMON_COFFEE/"07_leakage_position"/"n7_x1_10nM"/f"correct_k25_t{t}_regime_{regime}_{nr.energy_tag(energy)}"
    cof, sorted_csv, coffee_s=nr.run_coffee(monomer,hb,out_dir,regime=regime,bond_energy=energy,force=force)
    free,total=nr.final_output_metrics_from_csv(sorted_csv)
    return {"missing_input_module":"correct", "removed_input":"", "case":"x1_10nM_correct", "regime":regime, "bond_energy":energy, "final_output_concentration_M":free, "total_final_output_containing_M":total, "coffee_runtime_s":coffee_s, "hilbert_basis":str(hb), "monomer_file":str(monomer), "coffee_output":str(cof), "sorted_concentration_csv":str(sorted_csv), "hb_t":t, "hb_projected_runtime_s":meta.get("projected_runtime_s", "")}

def write_outputs(rows):
    CSV_DIR.mkdir(parents=True, exist_ok=True); FIG_DIR.mkdir(parents=True, exist_ok=True)
    energy_tag = nr.energy_tag(float(rows[0]["bond_energy"])) if rows else "energy"
    fields=["missing_input_module","removed_input","case","regime","bond_energy","final_output_concentration_M","total_final_output_containing_M","coffee_runtime_s","hilbert_basis","monomer_file","coffee_output","sorted_concentration_csv"]
    with (CSV_DIR/f"leakage_position_final_output_{energy_tag}.csv").open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=fields, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with (BENCH_DIR/f"summary_{energy_tag}.json").open("w") as f: json.dump({"generated":datetime.now().isoformat(),"rows":rows}, f, indent=2)
    for regime in sorted({r["regime"] for r in rows}):
        leak_sub=sorted([r for r in rows if r["regime"]==regime and isinstance(r["missing_input_module"], int)], key=lambda r:r["missing_input_module"])
        correct=[r for r in rows if r["regime"]==regime and r["case"]=="x1_10nM_correct"]
        fig,ax=plt.subplots(figsize=(7,4.2)); ax.plot([r["missing_input_module"] for r in leak_sub],[r["final_output_concentration_M"] for r in leak_sub], marker="o", lw=1.8, label="Xj=0 leak")
        if correct:
            ax.axhline(correct[-1]["final_output_concentration_M"], color="#ff7f0e", lw=1.8, ls="-", label="X1=10 nM correct")
        ax.set_xlabel("missing-input module"); ax.set_ylabel("free final output concentration (M)"); ax.set_title(f"Experiment 7: leakage position, regime {regime}, E={leak_sub[0]['bond_energy'] if leak_sub else correct[-1]['bond_energy']:g}"); ax.legend(); ax.grid(True, ls=":", alpha=.4); fig.tight_layout(); fig.savefig(FIG_DIR/f"leakage_position_regime_{regime}_{energy_tag}.png", dpi=160); plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bond-energy", type=float, default=nr.DEFAULT_BOND_ENERGY)
    ap.add_argument("--regimes", nargs="+", default=["A","B"], choices=["A","B"])
    ap.add_argument("--modules", nargs="+", type=int, default=list(range(1,8)))
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-correct", action="store_true", help="skip X1=10 nM correct-output full-system run")
    args=ap.parse_args(); rows=[]
    for regime in args.regimes:
        if not args.no_correct:
            print(f"[exp7] correct-output full system regime={regime} E={args.bond_energy:g}")
            rows.append(run_correct(regime,args.bond_energy,args.force)); write_outputs(rows)
        for module in args.modules:
            print(f"[exp7] module={module} regime={regime} E={args.bond_energy:g}")
            row=run_cell(module,regime,args.bond_energy,args.force)
            rows.append(row); write_outputs(rows)
    print(f"Experiment 7 new-regime outputs: {BENCH_DIR}")
if __name__=="__main__": main()
