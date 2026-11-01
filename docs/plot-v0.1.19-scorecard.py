"""Staged scorecard only; publication is owned by the root agent."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np

FIELDS = ('t850','t500','t200','u850','u200','cwv','precip','olr','osr')
CORE = tuple(f for f in FIELDS if f != 'cwv')
REGIONS = ('tropics','subtropics','nh_mid','sh_mid','nh_polar','sh_polar')
WEIGHTS = np.array([np.sin(np.deg2rad(20)),.5-np.sin(np.deg2rad(20)),
                    (np.sin(np.deg2rad(60))-.5)/2,(np.sin(np.deg2rad(60))-.5)/2,
                    (1-np.sin(np.deg2rad(60)))/2,(1-np.sin(np.deg2rad(60)))/2])
NAMES = ('Temperature 850 hPa','Temperature 500 hPa','Temperature 200 hPa',
         'Zonal wind 850 hPa','Zonal wind 200 hPa','Column water','Precipitation',
         'Outgoing longwave','Reflected shortwave')
REGION_NAMES = ('20°S–20°N','20–30° N/S','30–60°N','30–60°S','60–90°N','60–90°S')
COLOURS = {'control':'#89939c','candidate':'#245f91'}


def require(condition,message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def recompute(score,label):
    """Independent calculation from regional errors, not the stored aggregator."""
    require(score['schema']=='scorecard_v4' and len(score['cmip6']['cohort'])==59,'wrong cohort')
    summary=score['cmip6']['summary']
    med=np.array([[summary[f][r]['rmse']['median'] for r in REGIONS] for f in FIELDS])
    require(np.isfinite(med).all() and (med>0).all(),'invalid normalizers')
    for fi,f in enumerate(FIELDS):
        for ri,r in enumerate(REGIONS):
            available=[]
            for model in score['cmip6']['models'].values():
                cell=model['fields'].get(f,{}).get(r)
                if cell and cell['rmse'] is not None and np.isfinite(cell['rmse']):
                    available.append(cell['rmse'])
            require(len(available)==summary[f][r]['rmse']['n'],'wrong median population')
            require(np.isclose(np.median(available),med[fi,ri],rtol=0,atol=1e-12),'median mismatch')
    errors=np.array([[score['runs'][label]['fields'][f][r]['rmse'] for r in REGIONS] for f in FIELDS])
    require(np.isfinite(errors).all() and (errors>=0).all(),'invalid run errors')
    ratios=errors/med
    field_index=ratios@WEIGHTS
    onr=float(np.mean(field_index))
    core_indices=[FIELDS.index(f) for f in CORE]
    onr8=float(np.mean(field_index[core_indices]))
    cohort=[]
    for model in score['cmip6']['models'].values():
        values=np.array([[model['fields'][f][r]['rmse'] for r in REGIONS] for f in CORE])
        cohort.append(float(np.mean((values/med[core_indices])@WEIGHTS)))
    require(len(cohort)==59 and np.isfinite(cohort).all(),'incomplete position cohort')
    position=int(np.count_nonzero(np.array(cohort)>onr8))
    saved=score['aggregate']['runs'][label]
    require(np.isclose(onr,saved['onr9'],rtol=0,atol=1e-12),'ONR9 mismatch')
    require(np.isclose(onr8,saved['onr8_position'],rtol=0,atol=1e-12),'ONR8 mismatch')
    require(position==saved['cmip6_models_worse'],'rank mismatch')
    require(np.allclose(field_index,[saved['field_index'][f] for f in FIELDS],rtol=0,atol=1e-12),'field index mismatch')
    return dict(onr9=onr,position=position,ratios=ratios.tolist(),rmse=errors.tolist(),
                field_index=field_index.tolist(),toa_bias=score['runs'][label]['toa_bias'])


def clean_axis(ax):
    for side in ('top','right'):
        ax.spines[side].set_visible(False)
    ax.spines['left'].set_color('#b3bac0')
    ax.spines['bottom'].set_color('#b3bac0')
    ax.grid(axis='y',color='#e7e9ec',lw=.7)
    ax.set_axisbelow(True)


def draw(data,png,pdf):
    LABELS = data["labels"]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#202b34',
                         'axes.labelcolor':'#202b34','xtick.color':'#33434f','ytick.color':'#33434f'})
    fig=plt.figure(figsize=(14.5,16.9),facecolor='white')
    gs=fig.add_gridspec(4,2,height_ratios=(1.25,2.65,1.1,1.35),hspace=.53,wspace=.36,
                        left=.16,right=.94,bottom=.09,top=.865)
    fig.suptitle(data['title'],x=.5,y=.976,fontsize=20,weight='bold')
    fig.text(.5,.955,'31-day runs from observed starts • January and July 1993 • three matched pairs per month',ha='center',fontsize=11)
    fig.text(.5,.937,'Weather and rainfall compared with 59 climate models',ha='center',fontsize=10,color='#52606b')
    # Preserve the frozen scorer verdicts and the separate neutral decision.
    scored = data['frozen_verdicts_as_scored']
    neutral = data['neutral_release_decision']
    v3, v4 = scored['v3_season_verdicts'], scored['v4_season_verdicts']
    decision = neutral['skill_adjudication']
    banner = ("Seven-field score: January " + str(v3['JAN']) + ", July " + str(v3['JUL']) +
              ", overall " + str(scored['v3_metric_verdict']) +
              "   |   Nine-field score: January " + str(v4['JAN']) + ", July " +
              str(v4['JUL']) + ", overall " + str(scored['v4_metric_verdict']) +
              "\nPlanned release rule: " + decision +
              ". The energy or sea-ice check failed.")
    fig.text(.5,.9,banner,ha='center',va='center',fontsize=10.5,weight='bold',color='#202b34',
             bbox=dict(boxstyle='round,pad=0.55',facecolor='#f4f6f8',edgecolor='#6b7782',linewidth=1.3))
    for col,(metric,title,ylabel) in enumerate((('onr9','Overall error ↓','Error across nine weather fields'),
                                               ('position','CMIP6 position ↑','Models with larger error\n(out of 59)'))):
        ax=fig.add_subplot(gs[0,col])
        for x,season in enumerate(('JAN','JUL')):
            old=[r[metric] for r in data['seasons'][season]['control']]
            new=[r[metric] for r in data['seasons'][season]['candidate']]
            for i in range(3):
                offset=(i-1)*.025
                ax.plot([x-.15+offset,x+.15+offset],[old[i],new[i]],color='#d1d7dc',lw=.9,zorder=1)
            for arm,offset,marker in (('control',-.15,'o'),('candidate',.15,'s')):
                values=[r[metric] for r in data['seasons'][season][arm]]
                ax.scatter(x+offset+np.array([-.025,0,.025]),values,s=25,color=COLOURS[arm],marker=marker,
                           label=LABELS[arm] if x==0 else None,zorder=3)
                mean=float(np.mean(values))
                ax.annotate(f'{mean:.3f}' if metric=='onr9' else f'{mean:.2f}',(x+offset,mean),
                            xytext=(0,12),textcoords='offset points',ha='center',fontsize=10,weight='bold')
        ax.set_xticks([0,1],['January','July'])
        ax.set_xlim(-.45,1.45)
        ax.set_ylabel(ylabel)
        ax.set_title(title,loc='left',fontsize=13,weight='bold',pad=12)
        ax.set_ylim((min(.85, min(r['onr9'] for season in data['seasons'].values() for rows in season.values() for r in rows)-.06), max(1.15, max(r['onr9'] for season in data['seasons'].values() for rows in season.values() for r in rows)+.06)) if metric=='onr9' else (0,66))
        if metric == 'position':
            ax.set_yticks([0,10,20,30,40,50,59])
        if metric=='onr9':
            ax.axhline(1,color='#48525b',ls='--',lw=.8)
            ax.legend(frameon=False,ncol=2,loc='upper left',fontsize=9)
        clean_axis(ax)
    cmap=LinearSegmentedColormap.from_list('skill_reference',['#3573a1','#f9fafb','#b85b25'])
    norm=TwoSlopeNorm(vmin=0,vcenter=1,vmax=2.5)
    heat_axes=[]
    for col,season in enumerate(('JAN','JUL')):
        ax=fig.add_subplot(gs[1,col])
        heat_axes.append(ax)
        values=np.mean([r['ratios'] for r in data['seasons'][season]['candidate']],axis=0)
        im=ax.imshow(values,cmap=cmap,norm=norm,aspect='auto')
        for i in range(9):
            for j in range(6):
                value=values[i,j]
                ax.text(j,i,f'{value:.2f}',ha='center',va='center',fontsize=9,
                        color='white' if value<.35 or value>2 else '#202b34')
        ax.set_yticks(np.arange(9),NAMES if col==0 else ['']*9,fontsize=10)
        ax.set_xticks(np.arange(6),REGION_NAMES,rotation=45,ha='right',fontsize=9)
        ax.tick_params(length=0,pad=7)
        for spine in ax.spines.values(): spine.set_visible(False)
        ax.set_title(('January' if season=='JAN' else 'July')+' regional errors',loc='left',fontsize=13,weight='bold',pad=10)
    cax=fig.add_axes([.955,.425,.012,.27])
    cb=fig.colorbar(im,cax=cax,extend='max',ticks=[0,.5,1,1.5,2,2.5])
    cb.set_label('Error relative to climate-model median',fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.text(.55,.405,data['labels']['candidate']+' member means. Below 1: smaller error than the CMIP6 median. Above 1: larger error.',ha='center',fontsize=9)
    rgrid=gs[2,:].subgridspec(1,3,wspace=.6)
    for col,(key,label,unit) in enumerate((('olr','Outgoing-longwave bias','W m⁻²'),
                                          ('planetary_albedo','Planetary-albedo bias','fraction'),
                                          ('net_down','Energy imbalance at top of atmosphere','W m⁻²'))):
        ax=fig.add_subplot(rgrid[0,col])
        for i in range(3):
            ax.plot([0,1],[data['radiation'][key][arm][i] for arm in LABELS],color='#d1d7dc',lw=.9)
        for x,arm in enumerate(LABELS):
            vals=data['radiation'][key][arm]
            ax.scatter(x+np.array([-.035,0,.035]),vals,color=COLOURS[arm],marker='o' if x==0 else 's',s=22,zorder=3)
            ax.annotate(f'{np.mean(vals):.4f}' if key=='planetary_albedo' else f'{np.mean(vals):.2f}',
                        (x,np.mean(vals)),xytext=(0,10),textcoords='offset points',ha='center',fontsize=10)
        ax.set_xticks([0,1],list(LABELS.values()),fontsize=9)
        ax.set_xlim(-.4,1.4)
        ax.set_ylim(0,max(max(v) for v in data['radiation'][key].values())*1.35)
        ax.set_title(label,fontsize=11,loc='left',weight='bold')
        ax.set_ylabel(unit,fontsize=9)
        clean_axis(ax)
    ax=fig.add_subplot(gs[3,:])
    ax.axis('off')
    table_rows=[]
    for metric in ('Day-31 mean SST (°C)','Ice-area excess (million km²)','850-hPa wind error / CMIP6 (all regions)','Day-31 maximum ice h × A (m)'):
        row=[metric]
        for season in ('JAN','JUL'):
            values=[]
            for arm in LABELS:
                if metric.startswith('Day-31 mean'):
                    number=np.mean([r['safety']['ocean_area_weighted_surface_temperature_C'] for r in data['seasons'][season][arm]])
                    values.append(f'{number:+.2f}')
                elif metric.startswith('Ice-area'):
                    number=data['physical_review']['seasons'][season][arm+'_mean_bias_million_km2']
                    values.append(f'{number:+.2f}')
                elif metric.startswith('850'):
                    values.append(f"{np.mean([r['field_index'][3] for r in data['seasons'][season][arm]]):.2f}")
                else:
                    peaks=[r['thickness'][arm]['maximum_m'] for r in data['physical_review']['pairs'] if r['season']==season]
                    values.append(f'{min(peaks):.2f}–{max(peaks):.2f}')
            row.append(' → '.join(values))
        table_rows.append(row)
    table=ax.table(cellText=table_rows,colLabels=['Physical checks', 'January, '+LABELS['control']+' → '+LABELS['candidate'], 'July, '+LABELS['control']+' → '+LABELS['candidate']],
                   cellLoc='left',colLoc='left',colWidths=[.40,.30,.30],bbox=[0,.12,1,.86])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for (row,col),cell in table.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_facecolor('#e9eef2' if row==0 else ('#f4f6f8' if row%2 else 'white'))
        if row==0: cell.set_text_props(weight='bold')
    fig.text(.16,.065,'Energy panels use the absolute two-season mean of the three matched runs.',fontsize=9)
    fig.text(.16,.051,'Ocean temperature is the day-31 mean. Ice excess is the monthly difference from observations. Ice peaks show the day-31 range.',fontsize=9)
    fig.text(.16,.037,'These runs start from observations. Simple forecasts beat all 59 climate models in a separate test. Ranks do not show long-term skill.',fontsize=9)
    fig.text(.16,.021,data['physical_note'],fontsize=9,weight='bold')
    fig.text(.16,.007,data['scope'],fontsize=9,color='#52606b')
    fig.savefig(png,dpi=160,facecolor='white')
    fig.savefig(pdf,facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    data=load(args.data)
    require(data['schema']=='readiesm_staged_skill_scorecard_v019','wrong staged schema')
    draw(data,args.output/'scorecard.png',args.output/'scorecard.pdf')
    (args.output/'data.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print('STAGED_SCORECARD_RENDER_PASS',args.output)
