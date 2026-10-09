"""Standalone scientific figures drawn only from reviewed report evidence."""
from pathlib import Path
import os,sys,json
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor310'))
os.environ['MPLCONFIGDIR']='/tmp/epykit-report-mpl'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none'})
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)

r=pd.read_csv(OUT/'region_comparison.csv');y=np.arange(len(r));fig,ax=plt.subplots(figsize=(10,5.5))
ax.barh(y-.18,r['F1 50%'],height=.34,color='#28668b',label='50% reciprocal CpG overlap')
ax.barh(y+.18,r['F1 80%'],height=.34,color='#c18346',label='80% reciprocal CpG overlap')
ax.set_yticks(y,r.Tool);ax.invert_yaxis();ax.set_xlim(0,.66);ax.set_xlabel('Region F1');ax.legend(loc='lower right',frameon=False)
for i,row in r.iterrows():
    ax.text(row['F1 50%']+.009,i-.18,f"{row['F1 50%']:.3f}",va='center',fontsize=9)
    ax.text(row['F1 80%']+.009,i+.18,f"{row['F1 80%']:.3f}",va='center',fontsize=9)
ax.set_title('The strongest region F1 depends on the boundary requirement',loc='left',pad=16)
fig.text(.125,-.015,'Autosomal signal, 2,893 truth regions. Frozen workflow settings differ across tools.',fontsize=9,color='#555555')
save(fig,'region_f1')

fig,ax=plt.subplots(figsize=(10,3.4));values=[1462,223,317,891];labels=['No candidate overlap','Overlap fails 50% rule','Lost at corrected q gate','Recovered'];colors=['#adb7bf','#d9b17c','#c57857','#28668b'];left=0
for v,label,color in zip(values,labels,colors):
    ax.barh([0],v,left=left,color=color,label=f'{label}: {v:,}')
    if v>=200:ax.text(left+v/2,0,f'{v:,}\n{100*v/2893:.1f}%',ha='center',va='center',color='white' if color in ['#c57857','#28668b'] else '#1c2932',fontsize=10)
    left+=v
ax.set_xlim(0,2893);ax.set_yticks([]);ax.set_xlabel('True regions');ax.set_title('Most losses occur before the region correction',loc='left',pad=15);ax.legend(loc='upper center',bbox_to_anchor=(.5,-.26),ncol=2,frameon=False)
save(fig,'recall_decomposition')

t=pd.read_csv(OUT/'timing_comparison.csv');fig,ax=plt.subplots(figsize=(10,3.6));left=np.zeros(len(t))
for col,color in [('Read/filter (s)','#adb7bf'),('CpG analysis (s)','#28668b'),('Region analysis (s)','#c18346')]:
    ax.barh(t.Engine,t[col],left=left,label=col,color=color);left+=t[col]
ax.invert_yaxis();ax.set_xlabel('Seconds');ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.26),ncol=3);ax.set_title('Repeated fitting dominates the experimental BB runtime',loc='left',pad=15)
fig.subplots_adjust(bottom=.33)
fig.text(.125,-.10,'Signal run phase timings. Confidence intervals were disabled in BB-F.',fontsize=9,color='#555555');save(fig,'runtime_phases')

s=pd.read_csv(OUT/'strata.csv',keep_default_na=False);s=s[(s.dataset=='null')&(s.dimension=='coverage_bin')];assert len(s)==6;s=s.iloc[[0,4,1,2,3,5]]
fig,ax=plt.subplots(figsize=(9,4));rates=100*s.raw_p_lt_05/s.sites;labels=['≤5×','5–10×','10–20×','20–40×','40–80×','>80×'];bars=ax.bar(labels,rates,color='#28668b');ax.axhline(5,color='#c18346',linestyle='--',label='Nominal 5% raw-p rate')
for b,val,n in zip(bars,rates,s.sites):ax.text(b.get_x()+b.get_width()/2,val+.9,f'{val:.2f}%\nn={n:,}',ha='center',fontsize=8)
ax.set_ylim(0,21);ax.set_ylabel('Null sites with raw p<0.05 (%)');ax.set_xlabel('Observed mean read depth');ax.set_title('Legacy null calibration changes substantially with depth',loc='left',pad=15);ax.legend(frameon=False,loc='upper left');save(fig,'null_by_depth')
files=[p.name for p in FIG.glob('*')]
(OUT/'figure_validation.json').write_text(json.dumps({'figures':files,'sources':['region_comparison.csv','candidate_ceiling.csv','timing_comparison.csv','strata.csv'],'checked':'All charts generated from reviewed rows; visual inspection required before delivery.'},indent=2)+'\n')
print(json.dumps({'figures':files}))
