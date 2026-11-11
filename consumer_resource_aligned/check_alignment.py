"""Compare this control with the antibiotic model at initial_antibiotic=0."""
import argparse,importlib.util,sys
from pathlib import Path
import numpy as np
import model

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--antibiotic-model',required=True,type=Path)
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('alignment_antibiotic',args.antibiotic_model)
    antibiotic=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=antibiotic;spec.loader.exec_module(antibiotic)
    for seed in [42,43]:
        for settings in [{},{'variable_affinity':False,'resistance_cost_max':0.},
                         {'structured_community':True}]:
            common=dict(seed=seed,t_end=300.,n_times=301,**settings)
            p=model.generate(model.Config(**common))
            q=antibiotic.generate(antibiotic.Config(initial_antibiotic=0.,**common))
            for key in ['initial','K','uptake','split','retain','leak','resistance','IC50','growth_multiplier']:
                np.testing.assert_array_equal(p[key],q[key],err_msg=key)
            np.testing.assert_array_equal(p['table'].to_numpy(),q['table'].to_numpy())
            r=model.simulate(p);s=antibiotic.simulate(q)
            for key in ['N','R','gross','net','death_rates','death_loss','cutoff_loss']:
                # Removing the constant-zero antibiotic state changes adaptive
                # solver steps slightly. Traits match exactly; trajectories are
                # compared within numerical integration tolerance.
                np.testing.assert_allclose(r[key],s[key],rtol=3e-6,atol=2e-8,err_msg=key)
            print('Matched:',seed,settings,'max biomass difference',np.max(abs(r['N']-s['N'])))
    print('All six zero-antibiotic alignment checks passed.')

if __name__=='__main__':main()
