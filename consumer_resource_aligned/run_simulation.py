import argparse
from model import Config
from analysis_tools import run_scenarios,make_figures,repeated_scenarios,repeated_plot

def main():
    print("Aligned control: zero mortality; balanced feeding categories; cost/affinity scenarios.")
    parser=argparse.ArgumentParser();parser.add_argument('--seed',type=int,default=42);parser.add_argument('--replicates',type=int,default=10);parser.add_argument('--output',default='outputs');parser.add_argument('--composition',action='store_true',help='Also run the 13 starting-composition conditions');args=parser.parse_args()
    c=Config(seed=args.seed);runs,summary=run_scenarios(c,args.output);figs=make_figures(runs,args.output)
    if args.replicates:
        df=repeated_scenarios(c,args.replicates);df.to_csv(f'{args.output}/repeated_species_summary.csv',index=False);figs.append(repeated_plot(df,args.output))
    from pathlib import Path
    Path(args.output,'figure_captions.txt').write_text('\n\n'.join(caption for _,caption in figs))
    if args.composition:
        from run_initial_composition import run
        run(Config(seed=args.seed),Path(args.output)/'initial_composition_experiment')
    print(summary.to_string(index=False))
if __name__=='__main__':main()
