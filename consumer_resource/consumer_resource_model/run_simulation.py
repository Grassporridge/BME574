import argparse
from model import Config
from analysis_tools import run_scenarios,make_figures,repeated_scenarios,repeated_plot

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--seed',type=int,default=42);parser.add_argument('--replicates',type=int,default=10);parser.add_argument('--output',default='outputs');args=parser.parse_args()
    c=Config(seed=args.seed);runs,summary=run_scenarios(c,args.output);figs=make_figures(runs,args.output)
    if args.replicates:
        df=repeated_scenarios(c,args.replicates);df.to_csv(f'{args.output}/repeated_species_summary.csv',index=False);figs.append(repeated_plot(df,args.output))
    from pathlib import Path
    Path(args.output,'figure_captions.txt').write_text('\n\n'.join(caption for _,caption in figs))
    print(summary.to_string(index=False))
if __name__=='__main__':main()
