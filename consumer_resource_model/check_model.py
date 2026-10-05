"""Capacity, matched-draw, accounting, mortality and extinction checks."""
import numpy as np
from dataclasses import replace
from model import Config,generate,simulate,death_rates
from analysis_tools import SCENARIOS,diversity
c=Config();p=generate(c)
assert np.allclose(np.sum(p['uptake']*p['retain'],axis=1),p['table'].mu_max)
assert np.allclose(death_rates(p['table'].mu_max.to_numpy(),p),p['table'].death_rate)
assert np.allclose(death_rates(np.zeros(25),p),p['table'].death_rate+p['starve'])
base=generate(replace(c,**SCENARIOS['baseline']))
for name in ['affinity','starvation','both']:
    pp=generate(replace(c,**SCENARIOS[name]))
    for key in ['uptake','split','initial','leak']:assert np.array_equal(base[key],pp[key])
    assert np.array_equal(base['table'].death_rate,pp['table'].death_rate)
for name,settings in SCENARIOS.items():
    pp=generate(replace(c,**settings));r=simulate(pp)
    assert abs(r['budget_error']).max()<1e-6
    for i,te in enumerate(r['extinction']):
        if np.isfinite(te):assert (r['N'][r['t']>te,i]==0).all()
p0=generate(replace(c,initial_resource=0,t_end=10,n_times=41));r0=simulate(p0)
expected=p0['initial']*np.exp(-r0['t'][:,None]*(p0['table'].death_rate.to_numpy()+p0['starve']))
assert np.allclose(r0['N'],expected,rtol=1e-6,atol=1e-12)
assert np.isclose(diversity(r0).shannon.iloc[0],np.log(25))
fake={'N':np.zeros((2,25)),'t':np.array([0,1])};assert diversity(fake).shannon.isna().all()
print('Passed capacity, matched draws, death bounds, all-scenario budgets, extinction, analytic starvation and Shannon checks.')
