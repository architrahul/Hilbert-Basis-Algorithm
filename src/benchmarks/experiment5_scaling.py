#!/usr/bin/env python3
"""Experiment 5 under the new 10 nM/100 nM leakage regime.

Uses saved Hilbert bases from results/common/hilbert_basis; does not compute HB.
Runs COFFEE for incomplete cascades (X1=0 leak case) using saved HBs, and for the X1=10 nM correct-output case using saved/computed full-system covering HBs. Full-system covering uses t=5 unless the probe projection is over 1800s, then t=4.
"""
from __future__ import annotations
import argparse, csv, json, os, sys
from datetime import datetime
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "results" / ".matplotlib")); os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import concentration_regime as nr

BENCH_DIR=nr.RESULTS_NEW/"benchmarks"/"05_scaling"; CSV_DIR=BENCH_DIR/"csv"; FIG_DIR=BENCH_DIR/"figures"
METHODS=[("full_hb","Full HB",None), ("k25_t5","t=5",5), ("k25_t4","t=4",4), ("k25_t3","t=3",3)]

def method_m_values(method):
    return range(1,8) if method=="full_hb" else range(4,19)

def run_cell(m, method, label, t, regime, energy, force):
    hb=nr.hb_exp5_path(m, method)
    if not hb.exists():
        print(f"  missing saved HB; skip {method} m={m}: {hb}")
        return None
    monomer_file=nr.ensure_incomplete_cascade_m(m)
    out_dir=nr.COMMON_COFFEE/"05_scaling"/f"m{m}_x1_0"/f"{method}_regime_{regime}_{nr.energy_tag(energy)}"
    cof, sorted_csv, coffee_s=nr.run_coffee(monomer_file, hb, out_dir, regime=regime, bond_energy=energy, force=force)
    free,total=nr.final_output_metrics_from_csv(sorted_csv)
    return {"m":m,"case":"x1_0_leak","regime":regime,"method":method,"label":label,"t":t or "","bond_energy":energy,"free_final_output_M":free,"total_final_output_containing_M":total,"coffee_runtime_s":coffee_s,"hilbert_basis":str(hb),"monomer_file":str(monomer_file),"coffee_output":str(cof),"sorted_concentration_csv":str(sorted_csv)}

def run_correct_cell(m, regime, energy, force):
    hb, meta = nr.ensure_full_covering_hb(m)
    monomer_file=nr.ensure_full_cascade_m(m)
    t=meta.get("t", "")
    method=f"correct_k25_t{t}"
    out_dir=nr.COMMON_COFFEE/"05_scaling"/f"m{m}_x1_10nM"/f"{method}_regime_{regime}_{nr.energy_tag(energy)}"
    cof, sorted_csv, coffee_s=nr.run_coffee(monomer_file, hb, out_dir, regime=regime, bond_energy=energy, force=force)
    free,total=nr.final_output_metrics_from_csv(sorted_csv)
    return {"m":m,"case":"x1_10nM_correct","regime":regime,"method":method,"label":f"correct t={t}","t":t,"bond_energy":energy,"free_final_output_M":free,"total_final_output_containing_M":total,"coffee_runtime_s":coffee_s,"hilbert_basis":str(hb),"monomer_file":str(monomer_file),"coffee_output":str(cof),"sorted_concentration_csv":str(sorted_csv),"hb_projected_runtime_s":meta.get("projected_runtime_s", ""),"hb_runtime_s":meta.get("hb_runtime_s", "") }

def write_outputs(rows):
    CSV_DIR.mkdir(parents=True, exist_ok=True); FIG_DIR.mkdir(parents=True, exist_ok=True)
    energy_tag = nr.energy_tag(float(rows[0]["bond_energy"])) if rows else "energy"
    fields=["m","case","regime","method","label","t","bond_energy","free_final_output_M","total_final_output_containing_M","coffee_runtime_s","hilbert_basis","monomer_file","coffee_output","sorted_concentration_csv"]
    with (CSV_DIR/f"summary_{energy_tag}.csv").open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=fields, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with (CSV_DIR/f"final_output_by_m_{energy_tag}.csv").open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=fields[:9]+["sorted_concentration_csv"], extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with (BENCH_DIR/f"summary_{energy_tag}.json").open("w") as f: json.dump({"generated":datetime.now().isoformat(),"rows":rows}, f, indent=2)
    colors={"full_hb":"#000000","k25_t5":"#1f77b4","k25_t4":"#ff7f0e","k25_t3":"#2ca02c","correct":"#d62728"}
    for regime in sorted({r["regime"] for r in rows}):
        fig, ax=plt.subplots(figsize=(7.5,4.8))
        # Leak case: one line per candidate-set method.
        for method,label,_ in METHODS:
            sub=sorted([r for r in rows if r["regime"]==regime and r["case"]=="x1_0_leak" and r["method"]==method], key=lambda r:int(r["m"]))
            if sub:
                ax.plot([int(r["m"]) for r in sub], [float(r["free_final_output_M"]) for r in sub], marker="o", lw=1.8, color=colors[method], linestyle="--", label=f"leak, {label}")
        # Correct-output case: full system with X1=10 nM. Its covering t may switch
        # from 5 to 4 per m, so keep it as one physical reference line.
        correct=sorted([r for r in rows if r["regime"]==regime and r["case"]=="x1_10nM_correct"], key=lambda r:int(r["m"]))
        if correct:
            ax.plot([int(r["m"]) for r in correct], [float(r["free_final_output_M"]) for r in correct], marker="s", lw=2.2, color=colors["correct"], linestyle="-", label="correct output, X1=10 nM")
        ax.set_xlabel("cascade length m")
        ax.set_ylabel("free final output concentration (M)")
        ax.set_title(f"Experiment 5: leak vs correct output, regime {regime}, E={float(rows[0]['bond_energy']):g}")
        ax.grid(True, ls=":", alpha=.4)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR/f"final_output_leak_vs_correct_regime_{regime}_{energy_tag}.png", dpi=160)
        plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bond-energy", type=float, default=nr.DEFAULT_BOND_ENERGY)
    ap.add_argument("--regimes", nargs="+", default=["A","B"], choices=["A","B"])
    ap.add_argument("--methods", nargs="+", default=[m[0] for m in METHODS], choices=[m[0] for m in METHODS])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--include-correct", action="store_true", default=True, help="include X1=10 nM correct-output case using full-system covering HB")
    ap.add_argument("--leak-only", action="store_true", help="only run X1=0 leak case")
    args=ap.parse_args()
    rows=[]
    for regime in args.regimes:
        for method,label,t in METHODS:
            if method not in args.methods: continue
            for m in method_m_values(method):
                print(f"[exp5] m={m} {method} regime={regime} E={args.bond_energy:g}")
                row=run_cell(m,method,label,t,regime,args.bond_energy,args.force)
                if row: rows.append(row); write_outputs(rows)
        if args.include_correct and not args.leak_only:
            for m in range(1,19):
                print(f"[exp5] m={m} correct-output full system regime={regime} E={args.bond_energy:g}")
                row=run_correct_cell(m,regime,args.bond_energy,args.force)
                rows.append(row); write_outputs(rows)
    print(f"Experiment 5 new-regime outputs: {BENCH_DIR}")
    print("Includes X1=0 leak plus X1=10 nM correct-output where requested; correct-output uses covering HB t=5 unless projected over 1800s, then t=4.")
if __name__=="__main__": main()
