"""Run baseline, repeated communities, and extinction-cutoff sensitivity."""
import argparse
from dataclasses import replace
from pathlib import Path
import pandas as pd
from model import Config,generate,simulate,export,summarize
from plots import make_plots,ensemble,ensemble_plot

def main():
    a=argparse.ArgumentParser();a.add_argument('--seed',type=int,default=42);a.add_argument('--replicates',type=int,default=20);a.add_argument('--output',default='outputs');args=a.parse_args()
    c=Config(seed=args.seed);p=generate(c);r=simulate(p);out=Path(args.output)
    export(p,r,out);make_plots(p,r,out)
    rows=[]
    for f in [1e-4,1e-6,1e-8]:
        pp=generate(replace(c,extinction_fraction=f));rr=simulate(pp)
        rows.append({'cutoff_fraction':f,'survivors':int((rr['N'][-1]>0).sum()),'final_total_biomass':rr['N'][-1].sum(),'max_budget_error':abs(rr['budget_error']).max()})
    pd.DataFrame(rows).to_csv(out/'cutoff_sensitivity.csv',index=False)
    if args.replicates:
        df=ensemble(c,args.replicates);df.to_csv(out/'ensemble_species_summary.csv',index=False);ensemble_plot(df,out)
    print(f"Survivors at t={c.t_end}: {(r['N'][-1]>0).sum()}/{c.n_species}; maximum budget error {abs(r['budget_error']).max():.2g}")
    print(summarize(p,r)[['growth_score','death_rate','final_biomass']].describe().to_string())
if __name__=='__main__':main()
