#!/usr/bin/env python3
"""Experiment 6 under the new regime using the `concentrations` executable.

Uses saved Exp5 Hilbert bases only. Defaults to bond energy -10 and regimes A/B;
pass --bond-energies to sweep more values.
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
BENCH_DIR=nr.RESULTS_NEW/"benchmarks"/"06_bond_energy_concentrations"; CSV_DIR=BENCH_DIR/"csv"; FIG_DIR=BENCH_DIR/"figures"

def run_cell(m, energy, regime, exe, force_inputs, force_solver, correct=False):
    if correct:
        hb, meta = nr.ensure_full_covering_hb(m)
        method=f"correct_k25_t{meta.get('t','')}"
        monomer=nr.ensure_full_cascade_m(m)
        out_dir=nr.COMMON_CONCENTRATIONS/"06_bond_energy"/f"m{m}_x1_10nM"/f"{method}_regime_{regime}_{nr.energy_tag(energy)}"
    else:
        method="full_hb" if m<=3 else "k25_t5"
        hb=nr.hb_exp5_path(m, method)
        if not hb.exists(): print(f"missing saved HB; skip m={m}: {hb}"); return None
        monomer=nr.ensure_incomplete_cascade_m(m)
        out_dir=nr.COMMON_CONCENTRATIONS/"06_bond_energy"/f"m{m}_x1_0"/f"{method}_regime_{regime}_{nr.energy_tag(energy)}"
    eq, sorted_csv, runtime=nr.run_concentrations(monomer,hb,out_dir,regime=regime,bond_energy=energy,force_inputs=force_inputs,force_solver=force_solver,exe=exe)
    free,total=nr.final_output_metrics_from_csv(sorted_csv)
    return {"m":m,"case":"x1_10nM_correct" if correct else "x1_0_leak","regime":regime,"method":method,"bond_energy":energy,"free_final_output_M":free,"total_final_output_containing_M":total,"solver_runtime_s":runtime,"hilbert_basis":str(hb),"input_dir":str(out_dir),"eq_output":str(eq),"sorted_concentration_csv":str(sorted_csv)}

def write_outputs(rows):
    CSV_DIR.mkdir(parents=True, exist_ok=True); FIG_DIR.mkdir(parents=True, exist_ok=True)
    fields=["m","case","regime","method","bond_energy","free_final_output_M","total_final_output_containing_M","solver_runtime_s","hilbert_basis","input_dir","eq_output","sorted_concentration_csv"]
    with (CSV_DIR/"final_output_by_bond_energy.csv").open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=fields, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with (BENCH_DIR/"summary.json").open("w") as f: json.dump({"generated":datetime.now().isoformat(),"rows":rows}, f, indent=2)
    colors={-10.0:"#9467bd",-20.0:"#000000",-50.0:"#d62728",-60.0:"#ff7f0e",-65.0:"#2ca02c",-70.0:"#17becf",-75.0:"#8c564b",-100.0:"#e377c2"}
    markers={"x1_0_leak":"o","x1_10nM_correct":"s"}
    linestyles={"x1_0_leak":"--","x1_10nM_correct":"-"}
    labels={"x1_0_leak":"leak, X1=0","x1_10nM_correct":"correct, X1=10 nM"}
    for regime in sorted({r["regime"] for r in rows}):
        fig, ax=plt.subplots(figsize=(7.5,4.8))
        for e in sorted({float(r["bond_energy"]) for r in rows if r["regime"]==regime}):
            for case in ["x1_0_leak", "x1_10nM_correct"]:
                sub=sorted([r for r in rows if r["regime"]==regime and float(r["bond_energy"])==e and r.get("case")==case], key=lambda r:int(r["m"]))
                if not sub:
                    continue
                ax.plot([int(r["m"]) for r in sub],[float(r["free_final_output_M"]) for r in sub], marker=markers[case], lw=1.9, linestyle=linestyles[case], label=f"E={e:g}, {labels[case]}", color=colors.get(e))
        ax.set_xlabel("cascade length m")
        ax.set_ylabel("free final output concentration (M)")
        ax.set_title(f"Experiment 6: concentrations executable, regime {regime}")
        ax.grid(True, ls=":", alpha=.4)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR/f"final_output_leak_vs_correct_by_bond_energy_regime_{regime}.png", dpi=160)
        plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--concentrations-exe", type=Path, default=nr.CONCENTRATIONS_EXE)
    ap.add_argument("--m-values", nargs="+", type=int, default=list(range(1,17)))
    ap.add_argument("--bond-energies", nargs="+", type=float, default=[nr.DEFAULT_BOND_ENERGY])
    ap.add_argument("--regimes", nargs="+", default=["A","B"], choices=["A","B"])
    ap.add_argument("--force-inputs", action="store_true"); ap.add_argument("--force-solver", action="store_true")
    ap.add_argument("--include-correct", action="store_true", help="also run X1=10 nM correct-output full-system covering candidate sets")
    args=ap.parse_args(); rows=[]
    for regime in args.regimes:
        for e in args.bond_energies:
            for m in args.m_values:
                print(f"[exp6] m={m} regime={regime} E={e:g}")
                row=run_cell(m,e,regime,args.concentrations_exe,args.force_inputs,args.force_solver, correct=False)
                if row: rows.append(row); write_outputs(rows)
                if args.include_correct:
                    print(f"[exp6] m={m} correct-output regime={regime} E={e:g}")
                    row=run_cell(m,e,regime,args.concentrations_exe,args.force_inputs,args.force_solver, correct=True)
                    if row: rows.append(row); write_outputs(rows)
    print(f"Experiment 6 new-regime outputs: {BENCH_DIR}")
if __name__=="__main__": main()
