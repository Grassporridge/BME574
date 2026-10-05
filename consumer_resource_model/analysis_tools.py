"""Matched scenario comparisons and captions."""
from pathlib import Path
from dataclasses import replace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from model import Config,generate,simulate,summarize,export,FOODS

SCENARIOS={
 'baseline':dict(variable_affinity=False,starvation_mortality=False,structured_community=False),
 'affinity':dict(variable_affinity=True,starvation_mortality=False,structured_community=False),
 'starvation':dict(variable_affinity=False,starvation_mortality=True,structured_community=False),
 'both':dict(variable_affinity=True,starvation_mortality=True,structured_community=False),
 'structured':dict(variable_affinity=True,starvation_mortality=True,structured_community=True)}

def diversity(r):
    total=r['N'].sum(axis=1);alive=total>0
    shares=np.divide(r['N'],total[:,None],out=np.zeros_like(r['N']),where=total[:,None]>0)
    H=-np.sum(shares*np.log(np.maximum(shares,1e-300)),axis=1);H[~alive]=np.nan
    effective=np.exp(H);richness=(r['N']>0).sum(axis=1)
    return pd.DataFrame({'time':r['t'],'total_biomass':total,'shannon':H,'effective_species':effective,'richness':richness})

def run_scenarios(config=Config(),out='outputs'):
    root=Path(out);root.mkdir(exist_ok=True,parents=True);runs={};rows=[]
    for name,settings in SCENARIOS.items():
        p=generate(replace(config,**settings));r=simulate(p);runs[name]=(p,r)
        folder=root/name;export(p,r,folder);d=diversity(r);d.to_csv(folder/'diversity_timeseries.csv',index=False)
        p['table'].select_dtypes(include='number').describe().T.to_csv(folder/'realized_parameter_spread.csv')
        df=summarize(p,r);rows.append({'scenario':name,'survivors':int(df.survived.sum()),'final_biomass':r['N'][-1].sum(),'final_shannon':d.shannon.iloc[-1],'max_budget_error':abs(r['budget_error']).max()})
    summary=pd.DataFrame(rows);summary.to_csv(root/'scenario_summary.csv',index=False)
    return runs,summary

def make_figures(runs,out='outputs'):
    out=Path(out);figures=[];plt.rcParams.update({'figure.dpi':120,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name,caption):
        fig.tight_layout();path=out/(name+'.png');fig.savefig(path,bbox_inches='tight');plt.close(fig);figures.append((path,caption))
    fig,axs=plt.subplots(2,2,figsize=(11,7))
    for name,(p,r) in runs.items():
        d=diversity(r)
        for ax,col in zip(axs.flat,['shannon','effective_species','richness','total_biomass']):ax.plot(d.time,d[col],label=name);ax.set(xlabel='Time',ylabel=col.replace('_',' '))
    axs[0,0].legend(fontsize=8);axs[1,1].set_yscale('log')
    save(fig,'01_scenario_diversity','Figure 1. Shannon H (natural logarithm), effective diversity exp(H), species remaining above the cutoff, and total living biomass. The first four scenarios share growth, baseline death, permissions, secretion and initial abundance draws. The structured scenario changes several traits jointly, so its difference cannot be attributed to one mechanism. Shannon is undefined when total biomass is zero.')
    p,r=runs['both'];df=p['table'];colors=plt.cm.turbo(np.linspace(0,1,len(df)))
    fig,axs=plt.subplots(2,3,figsize=(13,8))
    axs[0,0].hist(df.mu_max,bins=8);axs[0,0].set(xlabel='Maximum gross growth',ylabel='Species count')
    axs[0,1].scatter(df.death_rate,df.starvation_death_increment);axs[0,1].set(xlabel='Baseline mortality',ylabel='Additional starvation mortality')
    im=axs[0,2].imshow(p['K'],aspect='auto',norm=plt.matplotlib.colors.LogNorm());fig.colorbar(im,ax=axs[0,2],label='Half-saturation K');axs[0,2].set(xticks=range(3),xticklabels=FOODS,ylabel='Species index')
    im=axs[1,0].imshow(p['uptake'],aspect='auto');fig.colorbar(im,ax=axs[1,0],label='Maximum uptake');axs[1,0].set(xticks=range(3),xticklabels=FOODS,ylabel='Species index')
    axs[1,1].bar(range(len(df)),p['split'][:,0],label='M1');axs[1,1].bar(range(len(df)),p['split'][:,1],bottom=p['split'][:,0],label='M2');axs[1,1].legend();axs[1,1].set(xlabel='Species index',ylabel='Leaked-energy allocation')
    groups=df[[f'can_consume_{f}' for f in FOODS]].apply(lambda row:'+'.join(f for f,v in zip(FOODS,row) if v),axis=1).value_counts()
    axs[1,2].barh(groups.index,groups.values);axs[1,2].set(xlabel='Species count',ylabel='Permitted foods')
    save(fig,'02_species_traits','Figure 2. Realized traits in the both scenario. Lower K means uptake remains more efficient at low chemical concentration. Maximum uptake is zero for forbidden foods. Growth capacity, death parameters and K were sampled independently; no growth-survival or growth-affinity tradeoff was imposed. Secretion allocations describe potential output when R0 is actually consumed.')
    for name in ['both','structured']:
        p,r=runs[name];N=r['N'];t=r['t'];total=N.sum(axis=1)
        fig,axs=plt.subplots(2,2,figsize=(12,8))
        for i in range(N.shape[1]):axs[0,0].plot(t,np.maximum(N[:,i],r['cutoff']/10),color=colors[i],lw=.9)
        axs[0,0].axhline(r['cutoff'],color='black',ls='--');axs[0,0].set(yscale='log',xlabel='Time',ylabel='Absolute biomass')
        perc=np.divide(100*N,total[:,None],out=np.zeros_like(N),where=total[:,None]>0)
        axs[0,1].stackplot(t,perc.T,colors=colors);axs[0,1].set(xlabel='Time',ylabel='Community percentage',ylim=(0,100))
        for j,f in enumerate(FOODS):axs[1,0].plot(t,r['R'][:,j],label=f)
        axs[1,0].legend();axs[1,0].set(xlabel='Time',ylabel='Chemical energy concentration')
        for j,f in enumerate(FOODS):axs[1,1].plot(t,r['R'][:,j],label=f)
        axs[1,1].set(xlabel='Time',ylabel='Chemical energy concentration',xlim=(0,50));axs[1,1].legend()
        fig.legend(handles=[plt.Line2D([0],[0],color=colors[i],label=f'S{i+1:02d}') for i in range(N.shape[1])],loc='lower center',ncol=13,bbox_to_anchor=(.5,-.06),fontsize=7)
        save(fig,'03_dynamics_'+name,f'Figure 3 ({name}). Species biomass, composition, chemical depletion and an early-time chemical zoom. M1/M2 begin at zero and accumulate only through R0 consumption. Log plots show extinct species at a display floor, not positive simulated biomass. Percentages can rise while absolute biomass falls. Species colors are consistent across these panels.')
    p,r=runs['both'];t=r['t'];N=r['N'];fig,axs=plt.subplots(2,2,figsize=(11,7))
    for i in range(N.shape[1]):
        for ax,key in zip([axs[0,0],axs[0,1]],['gross','death_rates']):ax.plot(t,np.where(N[:,i]>0,r[key][:,i],np.nan),color=colors[i],lw=.8)
    axs[0,0].set(xlabel='Time',ylabel='Gross per-capita growth');axs[0,1].set(xlabel='Time',ylabel='Per-capita death')
    for i in range(N.shape[1]):axs[1,0].plot(t,np.where(N[:,i]>0,r['net'][:,i],np.nan),color=colors[i],lw=.8)
    axs[1,0].axhline(0,color='black',ls='--');axs[1,0].set(xlabel='Time',ylabel='Net per-capita growth')
    axs[1,1].plot(t,(N*r['gross']).sum(axis=1),label='Growth flux');axs[1,1].plot(t,r['mortality'].sum(axis=1),label='Death flux');axs[1,1].legend();axs[1,1].set(xlabel='Time',ylabel='Biomass / time')
    save(fig,'04_growth_death','Figure 4. Resource-limited gross growth, starvation-dependent mortality, net growth and community biomass fluxes in the both scenario. Death rises toward baseline plus starvation increment as food-supported growth approaches zero. Negative net growth means biomass decreases; it does not mean immediate extinction. Rates after irreversible removal are hidden.')
    df=summarize(p,r);fig,axs=plt.subplots(2,2,figsize=(11,7))
    for ax,x in zip(axs.flat,['mu_max','death_rate','starvation_death_increment','K_R0']):
        ax.scatter(df[x],np.log10(np.maximum(df.final_biomass,r['cutoff']/10)),c=df.survived.map({True:'tab:blue',False:'tab:red'}));ax.set(xlabel=x,ylabel='log10 final biomass (floor for extinct)')
    save(fig,'05_trait_outcomes','Figure 5. Single-community trait associations: blue species persist above the cutoff, red species were removed. Extinct final biomass is displayed at a floor. K_R0 has no uptake effect for species forbidden from using R0. These plots are descriptive; food permissions and interactions can confound apparent trait effects.')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for name,(p,r) in runs.items():
        axs[0].plot(r['t'],r['budget_error'],label=name)
        e=np.sort(r['extinction'][np.isfinite(r['extinction'])]);axs[1].step(np.r_[0,e,300],np.r_[0,np.arange(1,len(e)+1),len(e)],where='post',label=name)
    axs[0].set(xlabel='Time',ylabel='Energy-accounting error');axs[1].set(xlabel='Time',ylabel='Cumulative threshold extinctions');axs[1].legend(fontsize=8)
    save(fig,'06_accounting_extinctions','Figure 6. Accounting residuals and cumulative threshold extinctions. Energy-equivalent resources plus living biomass, accumulated mortality and cutoff losses remain equal to the initial total. The extinction curve records numerical cutoff crossings, not experimental viability measurements.')
    return figures

def repeated_scenarios(config=Config(),replicates=10):
    rows=[]
    for seed in range(config.seed,config.seed+replicates):
        for name,settings in SCENARIOS.items():
            p=generate(replace(config,seed=seed,**settings));r=simulate(p);df=summarize(p,r);df['seed']=seed;df['scenario']=name;df['shannon_final']=diversity(r).shannon.iloc[-1];rows.append(df)
    return pd.concat(rows,ignore_index=True)

def repeated_plot(df,out):
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    grouped=df.groupby(['seed','scenario']).agg(survivors=('survived','sum'),shannon=('shannon_final','first')).reset_index()
    names=list(SCENARIOS)
    for ax,key in zip(axs,['survivors','shannon']):
        for _,g in grouped.groupby('seed'):
            ax.plot(range(len(names)),g.set_index('scenario').reindex(names)[key],color='gray',alpha=.4,marker='o',lw=.7)
        ax.set(xticks=range(len(names)),xticklabels=names,ylabel=key);ax.tick_params(axis='x',rotation=25)
    fig.tight_layout();path=Path(out)/'07_repeated_comparisons.png';fig.savefig(path,dpi=120);plt.close(fig)
    return path,'Figure 7. Matched scenario comparisons over ten random communities. Each gray line connects conditions using one seed. Species counts and final Shannon need not move together: losing rare species may change richness more than Shannon. The structured condition changes several assumptions together. Lines summarize simulations, not biological confidence intervals.'
