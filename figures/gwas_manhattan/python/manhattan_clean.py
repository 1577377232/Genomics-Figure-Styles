"""Compact, unannotated SNP/SV/TEV Manhattan renderer. No GWAS refitting."""
from pathlib import Path
import json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
COLORS={'SNP':'#386C8C','SV':'#9689B9','TEV':'#E75A46'}
def draw(snp_chrom,snp_pos,snp_p,sv,layout,title,output,*,regional=None,threshold_p=None,width_mm=183,height_mm=68,outer_padding_bp=12000000,dpi=600):
    """sv must contain chrom, position, P, TE_coverage; unknown coverage stays explicit in audit."""
    plt.rcParams.update({'font.family':'Liberation Sans','font.size':7,'pdf.fonttype':42,'axes.linewidth':.55,'axes.unicode_minus':False})
    layout=layout.copy();layout['chrom']=layout.chrom.astype(str);layout['offset']=np.r_[0,np.cumsum(layout.length.to_numpy()+12000000)[:-1]]
    offsets=layout.set_index('chrom').offset.to_dict();snp_chrom=np.asarray(snp_chrom);snp_pos=np.asarray(snp_pos);snp_p=np.asarray(snp_p)
    valid=np.isfinite(snp_p)&(snp_p>0)&(snp_p<=1);sv=sv.copy();sv['chrom']=sv.chrom.astype(str);sv=sv[sv.chrom.isin(offsets)&np.isfinite(sv.P)&sv.P.gt(0)&sv.P.le(1)]
    if regional:
        chrom,center=regional;chrom=str(chrom);code=40 if chrom=='Z' else int(chrom);lo=max(0,center-2500000);hi=min(int(layout.set_index('chrom').loc[chrom,'length']),center+2500000)
        valid&=(snp_chrom==code)&(snp_pos>=lo)&(snp_pos<=hi);sv=sv[sv.chrom.eq(chrom)&sv.position.between(lo,hi)]
        x=snp_pos[valid]/1e6;sx=sv.position/1e6;xlim=(lo/1e6-.025,hi/1e6+.025);xlabel=f'Chromosome {chrom} position (Mb)'
    else:
        codes=np.array([40 if c=='Z' else int(c) for c in layout.chrom]);valid&=np.isin(snp_chrom,codes); off=np.zeros(42)
        for c in layout.chrom:off[40 if c=='Z' else int(c)]=offsets[c]
        x=(snp_pos[valid].astype(float)+off[snp_chrom[valid]])/1e6;sx=(sv.position+sv.chrom.map(offsets))/1e6
        end=float((layout.offset+layout.length).iloc[-1]/1e6);pad=outer_padding_bp/1e6;xlim=(-pad,end+pad);xlabel='Chromosome'
    pv=snp_p[valid];assert len(pv)>0
    y=-np.log10(pv);sy=-np.log10(sv.P.to_numpy());peak=max(float(y.max()),float(sy.max()) if len(sy) else 0,5);ymax=math.ceil((peak+max(1.,peak*.12))/2)*2
    fig=plt.figure(figsize=(width_mm/25.4,height_mm/25.4));ax=fig.add_axes([.084,.185,.904,.735]);ax.set(xlim=xlim,ylim=(0,ymax));ax.spines[['top','right']].set_visible(False)
    ax.scatter(x,y,s=2.3,c=COLORS['SNP'],alpha=.8,linewidths=0,rasterized=True,zorder=1)
    if threshold_p is not None:ax.axhline(-np.log10(threshold_p),c='#8B9398',lw=.65,ls=(0,(4,4)),zorder=2)
    classes=np.where(sv.TE_coverage.ge(.5),'TEV','SV')
    for cl in ['SV','TEV']:
        keep=classes==cl;ax.scatter(np.asarray(sx)[keep],sy[keep],s=3.1,c=COLORS[cl],alpha=.8,linewidths=0,rasterized=True,zorder=3 if cl=='SV' else 4)
    if not regional:
        mids=(layout.offset+layout.length/2)/1e6;selected=[0];available=(width_mm*(.904))/(xlim[1]-xlim[0])
        for j in range(1,len(layout)-1):
            if (mids.iloc[j]-mids.iloc[selected[-1]])*available>=5.1:selected.append(j)
        if len(selected)>1 and (mids.iloc[-1]-mids.iloc[selected[-1]])*available<5.1:selected.pop()
        selected.append(len(layout)-1);ax.set_xticks(mids.iloc[selected]);ax.set_xticklabels(layout.chrom.iloc[selected])
    ax.tick_params(length=2.5,width=.55,pad=2,labelsize=6.5);ax.set_xlabel(xlabel,fontsize=8,labelpad=3);ax.set_ylabel(r'$-\log_{10}(P)$',fontsize=8,labelpad=3)
    fig.text(.09,.975,title,ha='left',va='top',fontsize=7.2,color='#202830')
    ax.legend(handles=[Line2D([],[],marker='o',ls='',mew=0,ms=3.4,color=c,label=k) for k,c in COLORS.items()],loc='upper right',bbox_to_anchor=(1,1.07),ncol=3,frameon=False,fontsize=6.5,handletextpad=.25,columnspacing=.75,borderaxespad=0)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=dpi,facecolor='white');plt.close(fig)
    return dict(SNP_plotted=len(pv),SV_plotted=int((classes=='SV').sum()),TEV_plotted=int((classes=='TEV').sum()),SV_annotation_unknown=int(sv.TE_coverage.isna().sum()),SNP_minP=float(pv.min()),SV_TEV_minP=float(sv.P.min()),width_mm=width_mm,height_mm=height_mm,dpi=dpi,outer_padding_bp=outer_padding_bp,title=title,no_downsampling=True,no_peak_annotation=True,TEV_definition='TE-only union coverage >= 0.50 of trimmed changed sequence; unknown recorded separately',threshold_P=threshold_p,ymax=ymax)
