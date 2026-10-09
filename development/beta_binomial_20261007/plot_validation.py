"""Standalone scientific figures from the completed, untouched-seed evaluation."""
from pathlib import Path
import json
import os
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
OUT = Path(__file__).parent/'holdout_v3'
sys.path[:0] = [str(ROOT/'vendor310')]
os.environ.setdefault('MPLCONFIGDIR','/tmp/epykit_bb_validation_mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


def main():
    context = json.loads((OUT/'context.json').read_text())
    assert context['status']=='complete'
    data = pd.read_csv(OUT/'metrics.csv')
    signal = data[data.signal].reset_index(drop=True)
    names = {
        'balanced_signalTrue':'Balanced coverage · 5 per group',
        'high_depth_signalTrue':'Higher depth · 5 per group',
        'many_replicates_signalTrue':'Unequal coverage · 10 per group',
        'mean_rho_association_signalTrue':'Mean–dispersion association · 5 per group',
        'unequal_group_rho_signalTrue':'Unequal group dispersions · 5 per group',
        'missing_coverage_signalTrue':'15% missing coverage · 5 per group',
        'strong_effect_signalTrue':'Strong effects (0.6) · 5 per group',
    }
    assert set(signal.label)==set(names)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(1,3,figsize=(13.5,5.6),sharey=True)
    y = np.arange(len(signal))
    for ax,column,title in zip(axes,['recall','average_precision','fdp'],
                              ['Recall at BH q ≤ 0.05','Average precision (ranking)','Observed false-discovery proportion']):
        values = signal[column].to_numpy().copy()
        if column=='fdp':values[signal.tp+signal.fp==0] = np.nan
        ax.barh(y,values,color='#184e87',height=.6)
        for index,value in enumerate(values):
            if np.isnan(value):
                ax.text(.02,index,'No discoveries',va='center',color='#626e7b')
            else:
                ax.annotate(f'{value:.1%}',(value,index),xytext=(4,0),textcoords='offset points',va='center',fontsize=9)
        ax.set_xlim(0,1.1 if column!='fdp' else .10)
        ax.xaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
        ax.set_title(title,fontsize=10)
        if column=='fdp':ax.axvline(.05,color='#626e7b',linestyle='--',linewidth=1)
    axes[0].set_yticks(y,[names[label] for label in signal.label])
    axes[0].invert_yaxis()
    fig.suptitle('Experimental BB-F engine: held-out CpG power and ranking',fontsize=14)
    fig.text(.01,.02,'75,000 sites per scenario; approximately 1% signals. Five-per-group scenarios use effects of 0.2 unless labeled strong.\n'
                     'The F reference is approximate. One realized false-discovery proportion does not establish FDR control.',fontsize=9)
    fig.tight_layout(rect=(0,.08,1,.95))
    for suffix in ['png','svg']:
        fig.savefig(OUT/f'power_and_ranking.{suffix}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    null = data[~data.signal].reset_index(drop=True)
    fig,ax = plt.subplots(figsize=(10,7))
    y = np.arange(len(null))
    values = null['null_tail_0.05'].to_numpy()
    low,high = null['null_tail_0.05_lo'].to_numpy(),null['null_tail_0.05_hi'].to_numpy()
    ax.errorbar(values,y,xerr=np.vstack([values-low,high-values]),fmt='o',color='#184e87',
                ecolor='#79a9cf',capsize=3,markersize=5)
    ax.axvline(.05,color='#626e7b',linestyle='--',linewidth=1,label='Nominal raw p threshold: 5%')
    ax.set_yticks(y,[label.replace('_signalFalse','').replace('_',' ') for label in null.label])
    ax.invert_yaxis();ax.set_xlim(0,.06);ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlabel('Fraction of true-null sites with raw p ≤ 0.05')
    ax.set_title('Held-out null tails: approximate F reference')
    ax.grid(axis='x',alpha=.2);ax.legend(loc='lower right',fontsize=9)
    fig.text(.01,.015,'Nominal 95% site-binomial intervals omit uncertainty from the shared learned prior.\n'
                     'Raw p tails and genome-wide BH error rates are distinct measurements.',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1))
    for suffix in ['png','svg']:
        fig.savefig(OUT/f'null_tails.{suffix}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    print(f'Held-out figures written to {OUT}')


if __name__=='__main__':
    main()
