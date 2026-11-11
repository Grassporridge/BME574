"""Meaningful accounting, capacity, no-food and extinction checks."""
import numpy as np
from model import Config,generate,simulate
p=generate();assert np.allclose(p['uptake']@p['retain'],p['table'].mu_max)
r=simulate(p);assert abs(r['budget_error']).max()<1e-6
p0=generate(Config(initial_resource=0,n_times=101));r0=simulate(p0)
expected=p0['config'].initial_biomass/25*np.exp(-r0['t'][:,None]*p0['table'].death_rate.to_numpy())
assert np.allclose(r0['N'],expected,rtol=1e-6,atol=1e-12)
p1=generate(Config(death_mean=.2,death_sd=.02));r1=simulate(p1)
assert np.isfinite(r1['extinction']).all()
for i,te in enumerate(r1['extinction']): assert (r1['N'][r1['t']>te,i]==0).all()
print('Passed: saturation capacity, finite budget, analytic no-food mortality, irreversible extinction.')
