"""Checks for the antibiotic-aligned, mortality-free consumer-resource model."""
from dataclasses import replace
import numpy as np
from model import Config,generate,simulate,death_rates
from analysis_tools import SCENARIOS,diversity

c=Config();p=generate(c)
assert np.allclose(np.sum(p['uptake']*p['retain'],axis=1),p['table'].mu_max)
assert np.all(death_rates(np.zeros(c.n_species),p)==0)
counts=p['table'].food_category_code.value_counts()
assert len(counts)==7 and counts.max()-counts.min()<=1
base=generate(replace(c,**SCENARIOS['baseline']))
for name in ['affinity','cost','both']:
    pp=generate(replace(c,**SCENARIOS[name]))
    for key in ['split','initial','leak','resistance','food_category_code']:
        np.testing.assert_array_equal(base[key],pp[key])
    np.testing.assert_array_equal(base['table'].intrinsic_mu_max,pp['table'].intrinsic_mu_max)
for name,settings in SCENARIOS.items():
    pp=generate(replace(c,**settings));r=simulate(pp)
    assert abs(r['budget_error']).max()<1e-6
    assert np.all(np.diff(r['N'],axis=0)>=-1e-9)
    assert not np.isfinite(r['extinction']).any()
    assert np.max(abs(r['death_rates']))==0
p0=generate(replace(c,initial_resource=0,t_end=10,n_times=41));r0=simulate(p0)
np.testing.assert_allclose(r0['N'],np.broadcast_to(p0['initial'],r0['N'].shape),atol=1e-12)
assert np.isclose(diversity(r0).shannon.iloc[0],np.log(25))
fake={'N':np.zeros((2,25)),'t':np.array([0,1])}
assert diversity(fake).shannon.isna().all()
# Legacy starvation flag is deliberately inert.
legacy=generate(replace(c,starvation_mortality=True))
np.testing.assert_array_equal(p['uptake'],legacy['uptake'])
assert np.all(death_rates(np.zeros(c.n_species),legacy)==0)
print('Passed capacity, balanced categories, matched traits, all-scenario budgets, no-food/no-death and diversity checks.')
