"""Energy-accounted, well-mixed batch consumer-resource model."""
from dataclasses import dataclass, asdict
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import truncnorm

@dataclass
class Config:
    seed: int = 42
    n_species: int = 25
    t_end: float = 300.
    n_times: int = 1201
    growth_mean: float = .5
    growth_sd: float = .15
    growth_reference: float = 1.
    death_mean: float = .01
    death_sd: float = .003
    initial_biomass: float = .01
    initial_resource: float = 1.
    half_saturation: float = .1
    leakage: float = .2
    extinction_fraction: float = 1e-6

FOODS = ['R0', 'M1', 'M2']
def truncated(rng, mean, sd, low, high, size):
    return truncnorm.rvs((low-mean)/sd, (high-mean)/sd, loc=mean, scale=sd,
                         size=size, random_state=rng)

def generate(c=Config()):
    if c.n_species < 3 or c.n_times < 2 or c.t_end <= 0:
        raise ValueError('Require >=3 species, >=2 output times, and positive duration')
    if not (0 <= c.leakage < 1 and 0 < c.extinction_fraction < 1):
        raise ValueError('Require 0<=leakage<1 and 0<extinction_fraction<1')
    if min(c.growth_sd,c.death_sd,c.growth_reference,c.initial_biomass,c.half_saturation) <= 0 or c.initial_resource < 0:
        raise ValueError('Invalid positive scale or initial resource')
    rng = np.random.default_rng(c.seed)
    score = truncated(rng,c.growth_mean,c.growth_sd,0,1,c.n_species)
    death = truncated(rng,c.death_mean,c.death_sd,0,np.inf,c.n_species)
    # Sample uniformly among the seven nonempty food subsets. Guarantee coverage.
    codes = rng.integers(1,8,c.n_species)
    if c.n_species >= 3: codes[:3] = [1,2,4]
    mask = ((codes[:,None] >> np.arange(3)) & 1).astype(bool)
    weights = rng.dirichlet(np.ones(3),c.n_species)*mask
    weights /= weights.sum(axis=1,keepdims=True)
    # Retained energy fractions: R0 leaks, metabolites do not.
    retain = np.array([1-c.leakage,1.,1.])
    # Normalize capacity so growth at saturation of all permitted foods is mu_max.
    capacity = score*c.growth_reference/(weights@retain)
    uptake = capacity[:,None]*weights
    split = rng.dirichlet([1.,1.],c.n_species)
    table = pd.DataFrame({'species':[f'S{i+1:02d}' for i in range(c.n_species)],
        'growth_score':score,'mu_max':score*c.growth_reference,'death_rate':death,
        'total_uptake_capacity':capacity,'M1_secretion_share':split[:,0],
        'M2_secretion_share':split[:,1], 'initial_biomass':c.initial_biomass/c.n_species})
    for j,f in enumerate(FOODS):
        table[f'can_consume_{f}']=mask[:,j]
        table[f'preference_{f}']=weights[:,j]
        table[f'uptake_max_{f}']=uptake[:,j]
    return {'config':c,'table':table,'uptake':uptake,'split':split,'retain':retain}

def fluxes(N,R,p):
    c=p['config']
    u=p['uptake']*np.maximum(R,0)/(c.half_saturation+np.maximum(R,0))
    gross=u@p['retain']
    consumed=np.maximum(N,0)[:,None]*u
    secretion=np.sum(consumed[:,0,None]*c.leakage*p['split'],axis=0)
    return gross,consumed,secretion

def simulate(p):
    c=p['config']; n=c.n_species; death=p['table'].death_rate.to_numpy()
    N0=np.full(n,c.initial_biomass/n); cutoff=N0[0]*c.extinction_fraction
    # State: N, R0/M1/M2, cumulative mortality, biomass removed at cutoffs.
    y=np.r_[N0,c.initial_resource,0.,0.,0.,0.]
    active=np.ones(n,dtype=bool); extinction=np.full(n,np.nan)
    times=np.linspace(0,c.t_end,c.n_times); states=np.empty((len(times),n+5))
    t=0.; filled=np.zeros(len(times),dtype=bool)
    while t < c.t_end:
        def rhs(t,y):
            N=np.maximum(y[:n],0)*active; R=y[n:n+3]
            gross,consumed,secreted=fluxes(N,R,p)
            return np.r_[N*(gross-death),-consumed.sum(axis=0)+np.r_[0.,secreted],
                         N@death,0.]
        events=[]; ids=np.flatnonzero(active)
        for i in ids:
            def event(t,y,i=i): return y[i]-cutoff
            event.terminal=True; event.direction=-1; events.append(event)
        sol=solve_ivp(rhs,(t,c.t_end),y,method='DOP853',rtol=1e-8,atol=1e-13,
                      dense_output=True,events=events or None,max_step=2.)
        if not sol.success: raise RuntimeError(sol.message)
        end=sol.t[-1]; sel=(times>=t)&(times<=end)&(~filled)
        if sel.any():
            states[sel]=sol.sol(times[sel]).T; filled[sel]=True
        y=sol.y[:,-1].copy(); t=end
        if t>=c.t_end: break
        hits=[ids[k] for k,e in enumerate(sol.t_events) if len(e)]
        hits=list(set(hits+list(np.flatnonzero(active & (y[:n]<=cutoff*(1+1e-8))))))
        if not hits: raise RuntimeError('Extinction event without a crossing')
        y[-1]+=y[hits].sum(); y[hits]=0; active[hits]=False; extinction[hits]=t
    assert filled.all()
    N=np.maximum(states[:,:n],0); R=np.maximum(states[:,n:n+3],0)
    gross=np.array([fluxes(N[k],R[k],p)[0] for k in range(len(times))])
    result={'t':times,'N':N,'R':R,'gross':gross,'net':gross-death,
            'mortality':N*death,'extinction':extinction,'death_loss':states[:,-2],
            'cutoff_loss':states[:,-1],'cutoff':cutoff}
    budget=N.sum(axis=1)+R.sum(axis=1)+states[:,-2]+states[:,-1]
    result['budget_error']=budget-(c.initial_biomass+c.initial_resource)
    assert np.max(np.abs(result['budget_error']))<1e-6
    assert states[:,:n+3].min()>-1e-9
    return result

def summarize(p,r):
    df=p['table'].copy(); df['final_biomass']=r['N'][-1]
    df['peak_biomass']=r['N'].max(axis=0)
    df['biomass_time_integral']=np.trapz(r['N'],r['t'],axis=0)
    df['extinction_time']=r['extinction']; df['survived']=np.isnan(r['extinction'])
    df['final_percent']=100*r['N'][-1]/max(r['N'][-1].sum(),1e-300)
    return df

def export(p,r,out):
    from pathlib import Path
    import json
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    summarize(p,r).to_csv(out/'species_summary.csv',index=False)
    frames=[]
    for i,name in enumerate(p['table'].species):
        frames.append(pd.DataFrame({'time':r['t'],'species':name,'biomass':r['N'][:,i],
          'percent':100*r['N'][:,i]/np.maximum(r['N'].sum(axis=1),1e-300),
          'gross_growth_rate':r['gross'][:,i],'death_rate':p['table'].death_rate.iloc[i],
          'net_growth_rate':r['net'][:,i], 'biomass_death_flux':r['mortality'][:,i]}))
    pd.concat(frames).to_csv(out/'species_timeseries.csv',index=False)
    pd.DataFrame({'time':r['t'],**{f:r['R'][:,j] for j,f in enumerate(FOODS)},
      'total_biomass':r['N'].sum(axis=1),'cumulative_death_loss':r['death_loss'],
      'cumulative_cutoff_loss':r['cutoff_loss'],'budget_error':r['budget_error']}).to_csv(out/'community_timeseries.csv',index=False)
    (out/'config.json').write_text(json.dumps(asdict(p['config']),indent=2))
