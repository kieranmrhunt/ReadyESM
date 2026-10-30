#!/usr/bin/env python3
"""Render the source-backed development inventory to standalone PNG/PDF."""
import argparse
from pathlib import Path
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter

STAGES=['design','component','pilot','scored','release','accepted']
STAGE_LABELS=['Idea','Code','Pilot','Score','Check','In use']
FLAGS=['Off','Resume','Balance','Limits','Repeat','Climate']
PLAIN_RELEASE_LABELS={
    'component and coupled-step gates':'component and coupled-step checks',
    'coupled output smoke':'coupled output check',
    'component gates':'component checks',
    '390 d prototype: identity FAILED':'390 d test run: source mismatch',
    '31 d · optional ON/OFF pair':'31 d · optional feature on/off pair',
    '120 d output; drift review pending':'120 d output. Climate drift check pending',
}
PLAIN_RELEASE_CHANGES={
    'Darcy boundary fix':'Soil-water boundary fix',
    'Sigma-RH convection profile':'Height-varying convection humidity',
    'Momentum-only K profile':'Boundary-layer wind mixing',
    'CATKE energy and scored outputs':'Ocean surface-mixing energy and score files',
}
COLOURS={'accepted':'#285c8d','active':'#247c8b','supported':'#b07b13',
         'design':'#6e62a3','reserve':'#8a929c','closed':'#485463'}
FLAG_COLOURS={'P':'#285c8d','F':'#a83b55','~':'#b07b13','.':'#a4adb8','-':'#8a929c'}

def load_data(path):
    data=tomllib.loads(path.read_text())
    if [r['version'] for r in data['release']]!=[f'v0.1.{i}' for i in range(19)]:
        raise ValueError('Expected every version v0.1.0–v0.1.18 in order')
    meta=data['meta']
    if meta['candidate_version']!='v0.1.18' or meta['release_status'] not in ('candidate_preview','evidence_verified_for_v018'):
        raise ValueError('Expected v0.1.18 preview or evidence-verified inventory')
    if meta['release_status']=='evidence_verified_for_v018' and data['release'][-1]['status']!='qualified':
        raise ValueError('Evidence-verified inventory needs a qualified v0.1.18 row')
    if [r['id'] for r in data['experiment']]!=[f'{p}{i:02d}' for p in ('BS','C') for i in range(1,11)]:
        raise ValueError('Expected all twenty named programme tracks')
    for row in data['release']:
        if row['status'] not in ('released','unreleased','qualified'):raise ValueError('Unknown release status')
        if row.get('direct_days',1)<=0:raise ValueError('Duration must be positive')
        if row.get('provisional_days') and not row.get('provisional_label'):raise ValueError('Missing provisional scope')
    for row in data['experiment']:
        if row['status'] not in COLOURS or row['reached'] not in STAGES:raise ValueError('Unknown track disposition')
        if len(row['flags'])!=6 or any(f not in FLAG_COLOURS for f in row['flags']):raise ValueError('Invalid gate flags')
        if not row['evidence_scope'] or not row['source']:raise ValueError('Evidence scope/source missing')
    return data

def releases(ax,rows):
    ax.set_xscale('log');ax.set_xlim(.055,2400);ax.set_ylim(len(rows)-.4,-.8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{r['version']}   {r['date']:>7}   {PLAIN_RELEASE_CHANGES.get(r['change'],r['change'])}" for r in rows],fontsize=9.2)
    ax.tick_params(axis='y',length=0,pad=10)
    ticks=[1/12,1,10,31,120,390]
    ax.xaxis.set_major_locator(FixedLocator(ticks));ax.xaxis.set_major_formatter(FixedFormatter(['2 h','1 d','10 d','31 d','120 d','390 d']))
    ax.tick_params(axis='x',labelsize=9)
    ax.grid(axis='x',color='#e1e5e9',lw=.7);ax.set_axisbelow(True)
    for i,row in enumerate(rows):
        if i%2==0:ax.axhspan(i-.45,i+.45,color='#f5f7f9',zorder=0)
        candidate=row['status']=='unreleased'
        colour='#b07b13' if candidate else '#285c8d' if row['kind']=='model' else '#687482'
        day=row.get('direct_days')
        if day is not None:
            ax.scatter(day,i,s=48,marker='o',facecolor='white' if candidate else colour,edgecolor=colour,lw=1.5,zorder=4)
            # v0.1.14 also has a 120-day output marker. Put the 31-day label
            # left of its marker so the two labels cannot collide.
            annotation_kwargs={'xytext':(-7,0),'ha':'right'} if row['version']=='v0.1.14' else {'xytext':(7,0),'ha':'left'}
            ax.annotate(PLAIN_RELEASE_LABELS.get(row['direct_label'],row['direct_label']),(day,i),
                        textcoords='offset points',va='center',fontsize=8.6,color=colour,**annotation_kwargs)
        else:ax.text(.063,i,PLAIN_RELEASE_LABELS.get(row['direct_label'],row['direct_label']),va='center',fontsize=8.5,color='#647183')
        inherited=row.get('inherited_days')
        if inherited is not None:
            ax.scatter(inherited,i,s=50,facecolor='white',edgecolor=colour,lw=1.5,zorder=4)
            ax.annotate(f"{inherited:g} d from {row['inherited_from']}",(inherited,i),xytext=(7,0),textcoords='offset points',va='center',fontsize=8.6,color='#647183')
        provisional=row.get('provisional_days')
        if provisional is not None:
            failed=row.get('provisional_status')=='rejected_identity';col='#a83b55' if failed else '#b07b13'
            if failed:ax.scatter(provisional,i,s=48,marker='x',color=col,lw=1.4,zorder=4)
            else:ax.scatter(provisional,i,s=48,marker='D',facecolor='white',edgecolor=col,lw=1.4,zorder=4)
            ax.annotate(PLAIN_RELEASE_LABELS.get(row['provisional_label'],row['provisional_label']),(provisional,i),xytext=(7,0),textcoords='offset points',va='center',fontsize=8.5,color=col)
    ax.set_title('A. Release tests: simulated time reached',loc='left',fontweight='bold',fontsize=13,pad=31)
    ax.text(0,1.035,'Logarithmic duration axis · starting from observations does not establish long-term climate skill',transform=ax.transAxes,fontsize=9,color='#526174')
    for spine in ('top','right','left'):ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color('#bcc5cf')

def gate(ax,x,y,flag):
    col=FLAG_COLOURS[flag]
    if flag in ('P','F','~'):
        ax.scatter(x,y,s=73 if flag!='~' else 83,marker={'P':'o','F':'s','~':'D'}[flag],color=col,zorder=3)
        ax.text(x,y,{'P':'✓','F':'×','~':'~'}[flag],ha='center',va='center',color='white',fontsize=9,zorder=4)
    elif flag=='.':ax.scatter(x,y,s=52,facecolor='white',edgecolor=col,lw=1,zorder=3)
    else:ax.text(x,y,'—',ha='center',va='center',fontsize=10,color=col)

def programme(ax,rows,title):
    ax.set_xlim(0,100);ax.set_ylim(9.95,-1.2);ax.axis('off')
    ax.text(0,-.95,title,fontsize=13,fontweight='bold',va='bottom')
    stage_x=[49+i*3.6 for i in range(6)];flag_x=[76+i*4.5 for i in range(6)]
    ax.text(0,-.48,'Current work and test results',fontsize=9,color='#526174')
    ax.text(58,-.48,'Furthest test completed',fontsize=9,color='#526174',ha='center')
    ax.text(87.25,-.48,'Checks completed for each item',fontsize=9,color='#526174',ha='center')
    for x,name in zip(stage_x,STAGE_LABELS):ax.text(x,-.05,name,fontsize=7.6,ha='center',color='#526174')
    for x,name in zip(flag_x,FLAGS):ax.text(x,-.05,name,fontsize=8,ha='center',color='#526174')
    for i,row in enumerate(rows):
        y=i+.42
        if i%2==0:ax.add_patch(Rectangle((0,y-.46),100,.92,facecolor='#f5f7f9',edgecolor='none',zorder=0))
        col=COLOURS[row['status']]
        ax.add_patch(Rectangle((.2,y-.3),.28,.6,facecolor=col,edgecolor='none',zorder=2))
        label=row['label'].replace('EDMF','mixing')
        scope=row['evidence_scope'].replace('gate','test').replace('UNDECIDED','undecided').replace('; ','. ')
        if row['id']=='C04':scope='Slanted convection appears in less than 1% of inspected cases.'
        if row['id']=='C10':scope='Only wind mixing was accepted in v0.1.15.'
        if row['id']=='BS05':scope='Only one-column checks of which settings can be told apart.'
        ax.text(1,y-.17,f"{row['id']}  {label}",fontsize=9.5,fontweight='bold',va='center')
        ax.text(1,y+.14,scope,fontsize=8.7,color='#526174',va='center')
        ax.text(45,y-.17,row['status'],fontsize=7.6,ha='right',va='center',color=col)
        reached=STAGES.index(row['reached'])
        ax.plot([stage_x[0],stage_x[-1]],[y,y],color='#d5dce3',lw=1.4,zorder=1)
        ax.scatter(stage_x,[y]*6,s=13,color='#d5dce3',zorder=2)
        ax.plot([stage_x[0],stage_x[reached]],[y,y],color=col,lw=2,zorder=2)
        ax.scatter(stage_x[reached],y,s=75,marker='s' if row.get('infrastructure') else 'o',color=col,edgecolor='white',lw=.8,zorder=3)
        for x,flag in zip(flag_x,row['flags']):gate(ax,x,y,flag)
    ax.plot([0,100],[.02,.02],color='#b4bec9',lw=.8)
    ax.plot([72,72],[.02,9.9],color='#d5dce3',lw=.7)

def render(data,png,pdf):
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':'#1e2d3d','pdf.fonttype':42,'ps.fonttype':42})
    fig=plt.figure(figsize=(18,20.5),facecolor='white')
    grid=fig.add_gridspec(3,1,height_ratios=[1.15,1,1],left=.035,right=.965,top=.91,bottom=.12,hspace=.22)
    release_ax=fig.add_subplot(grid[0]);pos=release_ax.get_position();release_ax.set_position([.34,pos.y0,.625,pos.height])
    releases(release_ax,data['release'])
    programme(fig.add_subplot(grid[1]),data['experiment'][:10],'B. Exploratory model changes')
    programme(fig.add_subplot(grid[2]),data['experiment'][10:],'C. Clouds and near-surface air')
    candidate=data['meta']['release_status']=='candidate_preview'
    fig.text(.035,.975,'ReadyESM development scorecard',fontsize=23,fontweight='bold',va='top')
    fig.text(.035,.949,f"Evidence inventory: {data['meta']['as_of']}",fontsize=11,color='#526174')
    version=data['meta']['candidate_version']
    fig.text(.035,.932,f'{version} review preview — release evidence pending' if candidate else f'{version} — release record with disclosed test outcomes',fontsize=11,color='#986b14' if candidate else '#285c8d')
    handles=[Line2D([],[],marker='o',ls='',color='#285c8d',label='tested for this release'),
             Line2D([],[],marker='o',ls='',mfc='white',mec='#687482',label='earlier or unreleased test'),
             Line2D([],[],marker='D',ls='',mfc='white',mec='#b07b13',label='run finished, review pending'),
             Line2D([],[],marker='x',ls='',color='#a83b55',label='source mismatch, excluded')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.03,.087),ncol=4,frameon=False,fontsize=9)
    gates=[Line2D([],[],marker=m,ls='',mfc='white' if f=='.' else FLAG_COLOURS[f],mec=FLAG_COLOURS[f],label=l) for f,m,l in [('P','o','pass'),('F','s','fail'),('~','D','partial / limited'),('.','o','not established'),('-','_','not applicable')]]
    fig.legend(handles=gates,loc='lower left',bbox_to_anchor=(.03,.065),ncol=5,frameon=False,fontsize=9)
    fig.text(.035,.052,'Each row shows the furthest test completed. Items marked in use include optional tools. Climate skill is measured separately.',fontsize=9,color='#526174')
    fig.text(.035,.038,'Off: works when switched off · Resume: continues from a saved state · Balance: energy and water accounted for · Limits: values stay physical · Repeat: matched runs agree · Climate: meets target',fontsize=8.6,color='#526174')
    fig.text(.035,.024,'v0.1.12 tested an optional feature for 31 days but left it off in standard runs. v0.1.13 checked a short restart. v0.1.14/15 allowed stated limit exceptions. Long drift remains open.',fontsize=8.8,color='#526174')
    fig.text(.035,.010,data['meta']['publication_note'],fontsize=8.6,color='#526174')
    png.parent.mkdir(parents=True,exist_ok=True);pdf.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(png,dpi=160,facecolor='white');fig.savefig(pdf,facecolor='white');plt.close(fig)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--png',type=Path,required=True)
    parser.add_argument('--pdf',type=Path,required=True)
    args=parser.parse_args();render(load_data(args.data),args.png,args.pdf)
    print('Rendered',args.png,args.pdf)

if __name__=='__main__':main()
