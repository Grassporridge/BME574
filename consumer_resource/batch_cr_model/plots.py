"""Visualization and repeated-community diagnostics."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from model import FOODS, Config, generate, simulate, summarize
plt.rcParams.update({'figure.dpi':120,'axes.spines.top':False,'axes.spines.right':False})

def make_plots(p,r,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    df=summarize(p,r); t=r['t'];N=r['N']; colors=plt.cm.turbo(np.linspace(0,1,N.shape[1]))
    paths=[]
    def save(fig,name):
        fig.tight_layout();path=out/(name+'.png');fig.savefig(path,bbox_inches='tight');plt.close(fig);paths.append(path)
    fig,ax=plt.subplots(2,2,figsize=(11,8))
    ax[0,0].hist(df.growth_score,bins=np.linspace(0,1,11));ax[0,0].set(xlabel='Growth score',ylabel='Species count')
    ax[0,1].hist(df.death_rate,bins=8);ax[0,1].set(xlabel='Death rate / time unit',ylabel='Species count')
    im=ax[1,0].imshow(p['uptake'],aspect='auto');fig.colorbar(im,ax=ax[1,0],label='Maximum uptake / biomass / time')
    ax[1,0].set(xticks=range(3),xticklabels=FOODS,yticks=range(len(df)),yticklabels=df.species,title='Uptake profiles')
    ax[1,1].barh(df.species,df.M1_secretion_share,label='M1');ax[1,1].barh(df.species,df.M2_secretion_share,left=df.M1_secretion_share,label='M2')
    ax[1,1].set(xlabel='Allocation of leaked R0 energy');ax[1,1].legend()
    save(fig,'01_species_characterization')
    fig,ax=plt.subplots(2,2,figsize=(12,8))
    for i in range(N.shape[1]): ax[0,0].plot(t,np.maximum(N[:,i],r['cutoff']/10),color=colors[i],lw=1)
    ax[0,0].axhline(r['cutoff'],color='black',ls='--',label='Extinction cutoff');ax[0,0].set(yscale='log',ylabel='Absolute biomass',xlabel='Time');ax[0,0].legend()
    ax[0,1].stackplot(t,(100*N/np.maximum(N.sum(axis=1),1e-300)[:,None]).T,colors=colors)
    ax[0,1].set(ylabel='Community composition (%)',xlabel='Time',ylim=(0,100))
    for j,f in enumerate(FOODS): ax[1,0].plot(t,r['R'][:,j],label=f)
    ax[1,0].legend();ax[1,0].set(xlabel='Time',ylabel='Resource energy-equivalent concentration')
    ax[1,1].plot(t,N.sum(axis=1),label='Total biomass');ax[1,1].set(xlabel='Time',ylabel='Total biomass')
    a=ax[1,1].twinx();a.plot(t,(N>0).sum(axis=1),color='black',ls='--');a.set(ylabel='Species above extinction cutoff',ylim=(0,26))
    fig.legend(handles=[plt.Line2D([0],[0],color=colors[i],label=df.species.iloc[i]) for i in range(len(df))],loc='lower center',ncol=13,bbox_to_anchor=(.5,-.055),fontsize=7)
    save(fig,'02_batch_dynamics')
    fig,ax=plt.subplots(2,2,figsize=(12,7))
    for i in range(N.shape[1]):
        alive=N[:,i]>0
        ax[0,0].plot(t,np.where(alive,r['gross'][:,i],np.nan),color=colors[i],lw=.8)
        ax[0,1].plot(t,np.where(alive,r['net'][:,i],np.nan),color=colors[i],lw=.8)
    ax[0,0].set(xlabel='Time',ylabel='Gross per-capita growth rate')
    ax[0,1].axhline(0,color='black',ls='--');ax[0,1].set(xlabel='Time',ylabel='Net per-capita growth rate')
    ax[1,0].plot(t,(N*r['gross']).sum(axis=1),label='Growth flux');ax[1,0].plot(t,r['mortality'].sum(axis=1),label='Death flux')
    ax[1,0].legend();ax[1,0].set(xlabel='Time',ylabel='Biomass / time')
    ax[1,1].plot(t,r['death_loss'],label='Cumulative mortality');ax[1,1].plot(t,r['cutoff_loss'],label='Cutoff removals');ax[1,1].legend();ax[1,1].set(xlabel='Time',ylabel='Cumulative lost biomass')
    save(fig,'03_growth_and_death')
    fig,ax=plt.subplots(2,2,figsize=(11,8))
    for a,x in zip(ax.flat,['mu_max','death_rate','preference_R0','M1_secretion_share']):
        a.scatter(df[x],np.log10(np.maximum(df.final_biomass,r['cutoff']/10)),c=df.survived.map({True:'tab:blue',False:'tab:red'}))
        a.set(xlabel=x,ylabel='log10(final biomass; extinct at floor)')
    fig.suptitle('Single-community associations: blue survives, red extinct; not causal')
    save(fig,'04_trait_outcomes')
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    shares=N/np.maximum(N.sum(axis=1),1e-300)[:,None]
    shannon=-np.sum(shares*np.log(np.maximum(shares,1e-300)),axis=1)
    ax[0].plot(t,np.exp(shannon));ax[0].set(xlabel='Time',ylabel='Effective species number exp(Shannon)')
    ax[1].plot(t,r['budget_error']);ax[1].set(xlabel='Time',ylabel='Energy accounting error')
    save(fig,'05_diversity_and_budget')
    return paths

def ensemble(c,replicates=20):
    from dataclasses import replace
    frames=[]
    for k in range(replicates):
        p=generate(replace(c,seed=c.seed+k));r=simulate(p);df=summarize(p,r);df['community_seed']=c.seed+k;frames.append(df)
    return pd.concat(frames,ignore_index=True)

def ensemble_plot(df,out):
    fig,ax=plt.subplots(1,3,figsize=(13,4))
    # Continuous biomass-time outcome avoids treating zero biomass as a precise log value.
    for a,x in zip(ax,['mu_max','death_rate','preference_R0']):
        y=np.log10(np.maximum(df.biomass_time_integral,1e-300));a.scatter(df[x],y,s=8,alpha=.3)
        rho=spearmanr(df[x],y).statistic;a.set(xlabel=x,ylabel='log10 biomass-time integral',title=f'Pooled descriptive Spearman rho={rho:.2f}')
    fig.tight_layout();path=Path(out)/'06_repeated_community_associations.png';fig.savefig(path);plt.close(fig);return path
