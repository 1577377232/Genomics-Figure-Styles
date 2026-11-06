#!/usr/bin/env python3
"""User-approved approximate GGswu(.2) centromeres; keep evidence levels explicit."""
import csv
import re
from collections import defaultdict
import pysam
from prepare_data import OUT, ASSEMBLIES, CHROMS, LENGTHS, union_bp, tsv, write_json
from transfer_centromeres import records


def main():
    target=OUT/'01_reference/v2';target.mkdir(exist_ok=True)
    original={p[0]:(int(p[1]),int(p[2])) for p in (l.split() for l in (OUT/'01_reference/GGswu.centromeres.user_supplied.bed').read_text().splitlines() if not l.startswith('#'))}
    donorlengths={r['donor_chr']:int(r['length_bp']) for r in csv.DictReader((OUT/'01_reference/GGswu.chromosome_accessions.tsv').open(),delimiter='\t')}
    effective={c:(max(0,a),min(b,donorlengths[c])) for c,(a,b) in original.items()}
    tsv(target/'GGswu.centromeres.effective.tsv',['donor_chr','supplied_start','supplied_end','used_start','used_end','clipped_bp'],[[c,*original[c],*effective[c],original[c][1]-effective[c][1]] for c in original])
    rows=list(csv.DictReader((OUT/'01_reference/GGswu_to_GRCg7b.centromere_blocks.tsv').open(),delimiter='\t'));blocks=[]
    for r in rows:blocks.append([r['source_chr'],int(r['source_start']),int(r['source_end']),r['target_chr'],int(r['target_start']),int(r['target_end']),r['strand'],int(r['MAPQ']),r['query'],int(r['flag'])])
    raw=next(r['raw_assembly'] for r in csv.DictReader((ASSEMBLIES/'raw_assemblies_51.tsv').open(),delimiter='\t') if r['sample']=='Ggswu');segments={}
    # Only chr25 was excluded in v1. Reuse all the other measured CIGAR blocks.
    from pathlib import Path
    for acc,desc,seq in records(Path(raw)):
        if re.search(r'chromosome\s+25(?:\s|$)',desc):
            a,b=effective['25']
            for i,m in enumerate(re.finditer('[^N]+',seq)):
                if m.end()>a and m.start()<b:segments[f'{acc}_{i}']=(m.start(),m.end(),max(a,m.start()),min(b,m.end()))
            break
    with pysam.AlignmentFile(str(ASSEMBLIES/'00_bams/Ggswu.sorted.bam'),'rb') as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.mapping_quality<30 or r.query_name not in segments:continue
            offset,segend,a,b=segments[r.query_name];qlen=r.infer_read_length();assert qlen==segend-offset;q=0;t=r.reference_start
            for op,n in r.cigartuples:
                if op in [0,7,8]:
                    sa,sb=(offset+qlen-q-n,offset+qlen-q) if r.is_reverse else (offset+q,offset+q+n);x=max(a,sa);y=min(b,sb)
                    if y>x:
                        ta,tb=(t+sb-y,t+sb-x) if r.is_reverse else (t+x-sa,t+y-sa)
                        blocks.append(['25',x,y,r.reference_name,ta,tb,'-' if r.is_reverse else '+',r.mapping_quality,r.query_name,r.flag])
                    q+=n;t+=n
                elif op in [1,4,5]:q+=n
                elif op in [2,3]:t+=n
            assert q==qlen
    tsv(target/'centromere_alignment_blocks.tsv',['source_chr','source_start','source_end','target_chr','target_start','target_end','strand','MAPQ','query','flag'],blocks)
    candidates=[];audit=[]
    for c,(a,b) in effective.items():
        groups=defaultdict(list)
        for r in blocks:
            if r[0]==c:groups[(r[3],r[6])].append(r)
        total=union_bp([(r[1],r[2]) for rs in groups.values() for r in rs]);best=None
        for (tc,strand),rs in groups.items():
            if tc not in CHROMS:continue
            covered=union_bp([(r[1],r[2]) for r in rs]);ta=min(r[4] for r in rs);tb=max(r[5] for r in rs);multiplicity=sum(r[2]-r[1] for r in rs)/covered
            fraction=covered/(b-a);share=covered/total;ratio=(tb-ta)/(b-a);compression=(tb-ta)/covered
            high=fraction>=.8 and share>=.95 and multiplicity<=1.05 and .25<=ratio<=4
            approximate=fraction>=.15 and share>=.55 and multiplicity<=1.2 and 0<ratio<=2 and .25<=compression<=4 and covered>=15000
            if high or approximate:
                candidate=[tc,ta,tb,c,'high_coverage_transfer' if high else 'approximate_partial_transfer',fraction,share,multiplicity,strand]
                if best is None or (high,fraction)>(best[4]=='high_coverage_transfer',best[5]):best=candidate
        if best:candidates.append(best);audit.append([c,a,b,total,best[0],best[4]])
        else:audit.append([c,a,b,total,'NA','unresolved_not_drawn'])
    # Competing donor intervals must not silently become multiple centromeres on a target chromosome.
    selected={}
    for r in candidates:
        if r[0] not in selected or (r[4]=='high_coverage_transfer',r[5])>(selected[r[0]][4]=='high_coverage_transfer',selected[r[0]][5]):selected[r[0]]=r
    tsv(target/'GRCg7b.centromeres.approximate.tsv',['chrom','start','end','donor_chr','evidence','source_coverage','dominant_target_share','multiplicity','strand'],[selected[c] for c in CHROMS if c in selected])
    tsv(target/'centromere_status.tsv',['donor_chr','start','end','mapped_union_bp','selected_target','status'],audit)
    write_json(OUT/'05_QC/v2/centromere_update.json',{'donor_version':'GCA_024206055.2; confirmed by user','chr25_clip_bp':759,'source_coordinates_retained_as_BED':True,'high_coverage_transfers':sum(r[4]=='high_coverage_transfer' for r in selected.values()),'approximate_partial_transfers':sum(r[4]=='approximate_partial_transfer' for r in selected.values()),'not_active_centromere_validation':True,'input_intervals':len(effective),'drawn_target_chromosomes':[c for c in CHROMS if c in selected],'unresolved_not_fabricated':True})
    print('centromere v2',len(selected),'target chromosomes',flush=True)


if __name__=='__main__':main()
