#!/usr/bin/env python3
"""Supplement centromere markers using measured homology/flanking evidence.

No GGswu genomic coordinate is copied into GRCg7b. Low-coverage and flank-only
positions remain explicitly approximate; unplaced matches stay off the main plot.
"""
from pathlib import Path
from collections import defaultdict
import json
from centromere_anchors_v14 import read_tsv,write_tsv,union_bp

P=Path(__file__).resolve().parents[1]
OUT=P/'01_reference/v14'


def main():
    lengths={x.split()[0]:int(x.split()[1]) for x in (P/'01_reference/main.chrom.sizes').read_text().splitlines()}
    donor_lengths={r['donor_chr']:int(r['length_bp']) for r in read_tsv(P/'01_reference/GGswu.chromosome_accessions.tsv')}
    effective=read_tsv(P/'01_reference/v2/GGswu.centromeres.effective.tsv')
    old={r['donor_chr']:r for r in read_tsv(P/'01_reference/v2/GRCg7b.centromeres.approximate.tsv')}
    prim=read_tsv(OUT/'primary_correspondence.tsv')
    global_rows=read_tsv(OUT/'whole_chromosome_correspondence.tsv')
    blocks=read_tsv(P/'01_reference/v2/centromere_alignment_blocks.tsv')
    flanks=read_tsv(OUT/'centromere_flanking_M_blocks.tsv')
    for rs in [blocks,flanks]:
        for r in rs:
            for k in ['source_start','source_end','target_start','target_end','MAPQ','flag']:r[k]=int(r[k])
    primary={}
    for c in donor_lengths:
        candidates=sorted([r for r in prim if r['source_chr']==c],key=lambda r:-int(r['alignment_bp']))
        if candidates:primary[c]=candidates[0]
    selected={};audits=[];pending=[];unplaced=[]
    # Keep existing results when their target/strand agree with the primary
    # chromosome correspondence. Chr31's previous + repeat-only match fails this.
    for r in effective:
        c=r['donor_chr'];o=old.get(c);pr=primary.get(c)
        if o and pr and (o['chrom'],o['strand'])==(pr['target_chr'],pr['strand']):
            selected[o['chrom']]=dict(chrom=o['chrom'],start=int(o['start']),end=int(o['end']),donor_chr=c,
                evidence=o['evidence'],source_coverage=float(o['source_coverage']),strand=o['strand'],
                method='retained_previous_sequence_transfer',anchor_source_start='NA',anchor_source_end='NA',anchor_target_start='NA',anchor_target_end='NA')
        else:pending.append(r)
    reserved={r['target_chr'] for c,r in primary.items() if c in {e['donor_chr'] for e in effective} and r['target_chr'] in lengths}
    for r in pending:
        c=r['donor_chr'];a=int(r['used_start']);b=int(r['used_end']);dl=donor_lengths[c]
        pr=primary.get(c);target=None;strand=None;selection_note=''
        if pr and pr['target_chr'] in lengths:
            target,strand=pr['target_chr'],pr['strand'];selection_note='primary whole-query chromosome correspondence'
        else:
            # In fragmented microchromosomes prefer the nearest measured unique
            # flanking sequence, not a same-number assumption or a repeat hit.
            choices=[]
            for g in global_rows:
                tc=g['target_chr'];s=g['strand']
                if g['source_chr']!=c or tc not in lengths or tc in selected or tc in reserved or int(g['union_source_bp'])<5000:continue
                xs=[x for x in flanks if x['source_chr']==c and x['target_chr']==tc and x['strand']==s and x['source_end']-x['source_start']>=1000 and (x['source_end']<=a or x['source_start']>=b)]
                if not xs:continue
                distance=min(min(abs(a-x['source_end']),abs(x['source_start']-b)) for x in xs)
                choices.append((distance,-int(g['union_source_bp']),tc,s))
            if choices:
                _,_,target,strand=min(choices);selection_note='nearest unique flank; donor chromosome fragmented across reference sequences'
        if target is None or target in selected:
            projected=[x for x in blocks if pr and x['source_chr']==c and x['target_chr']==pr['target_chr'] and x['strand']==pr['strand']]
            primary_projected=[x for x in projected if not(x['flag'] & 2048)]
            if primary_projected:projected=primary_projected
            ua=min((x['target_start'] for x in projected),default='NA')
            ub=max((x['target_end'] for x in projected),default='NA')
            status='mapped_unplaced_not_in_main_plot' if projected else 'unplaced_or_ambiguous_not_in_main_plot'
            row=dict(donor_chr=c,target_chr=pr['target_chr'] if pr else 'NA',status=status,method=selection_note or 'primary correspondence is unplaced; no unoccupied main-chromosome flank',start=ua,end=ub)
            audits.append(row)
            if projected:unplaced.append(row.copy())
            continue
        xs=[x for x in flanks if x['source_chr']==c and x['target_chr']==target and x['strand']==strand and x['source_end']-x['source_start']>=1000]
        left=sorted([x for x in xs if x['source_end']<=a],key=lambda x:a-x['source_end'])
        right=sorted([x for x in xs if x['source_start']>=b],key=lambda x:x['source_start']-b)
        anchor=None;lo=0;hi=lengths[target]
        if a<=1 and right:
            anchor=right[0]
            if strand=='+':hi=anchor['target_start']
            else:lo=anchor['target_end']
        elif dl-b<=200000 and left:
            anchor=left[0]
            if strand=='+':lo=anchor['target_end']
            else:hi=anchor['target_start']
        if hi<=lo:
            lo=max(0,min(lengths[target]-1,lo));hi=lo+1
        # Existing exact CIGAR projections of the supplied interval. Restrict to
        # a terminal flank-bounded region before clustering repetitive hits.
        core=[x for x in blocks if x['source_chr']==c and x['target_chr']==target and x['strand']==strand and x['target_start']<hi and x['target_end']>lo]
        seeds=sorted(core,key=lambda x:-(x['source_end']-x['source_start']))[:80]
        best=None
        for seed in seeds:
            mid=(seed['target_start']+seed['target_end'])/2
            lower=max(lo,int(mid-100000));upper=min(hi,int(mid+100000))
            rows=[x for x in core if lower<=((x['target_start']+x['target_end'])/2)<=upper]
            if not rows:continue
            source_intervals=[];target_intervals=[]
            for x in rows:
                ta=max(lower,x['target_start']);tb=min(upper,x['target_end'])
                if ta>=tb:continue
                if strand=='+':sa=x['source_start']+ta-x['target_start'];sb=sa+tb-ta
                else:sa=x['source_end']-(tb-x['target_start']);sb=sa+tb-ta
                source_intervals.append((sa,sb));target_intervals.append((ta,tb))
            if not target_intervals:continue
            source_bp=union_bp(source_intervals)
            start=min(x[0] for x in target_intervals);end=max(x[1] for x in target_intervals)
            candidate=(source_bp,-(end-start),start,end,rows)
            if best is None or candidate[:2]>best[:2]:best=candidate
        if best and best[0]>=1000:
            bp,_,start,end,used=best;evidence='approximate_CIGAR_cluster';coverage=bp/(b-a)
            method=selection_note+'; measured centromere-homology cluster, restricted by flank when present'
        elif anchor is not None:
            start,end=lo,hi;evidence='approximate_flank_boundary';coverage=0
            method=selection_note+'; one-sided unique flank bounds terminal candidate; supplied centromere sequence largely absent/unmapped'
        else:
            audits.append(dict(donor_chr=c,target_chr=target,status='no_supported_centromere_location',method=selection_note,start='NA',end='NA'));continue
        assert 0<=start<end<=lengths[target]
        selected[target]=dict(chrom=target,start=start,end=end,donor_chr=c,evidence=evidence,source_coverage=coverage,strand=strand,
            method=method,anchor_source_start=anchor['source_start'] if anchor else 'NA',anchor_source_end=anchor['source_end'] if anchor else 'NA',anchor_target_start=anchor['target_start'] if anchor else 'NA',anchor_target_end=anchor['target_end'] if anchor else 'NA')
    for r in selected.values():
        audits.append(dict(donor_chr=r['donor_chr'],target_chr=r['chrom'],status='drawn_main_chromosome',method=r['method'],start=r['start'],end=r['end']))
    order=list(lengths)
    write_tsv(OUT/'GRCg7b.centromeres.approximate.tsv',[selected[c] for c in order if c in selected])
    changes=[]
    for c in order:
        if c not in selected:continue
        r=selected[c];o=old.get(r['donor_chr'])
        same=o and o['chrom']==c and (int(o['start']),int(o['end']))==(r['start'],r['end'])
        changes.append(dict(donor_chr=r['donor_chr'],target_chr=c,
            V13_start=o['start'] if o else 'NA',V13_end=o['end'] if o else 'NA',
            V14_start=r['start'],V14_end=r['end'],
            change='unchanged' if same else 'revised_using_primary_strand_and_flank' if o else 'added',
            evidence=r['evidence']))
    write_tsv(OUT/'centromere_changes_from_V13.tsv',changes)
    rank={r['donor_chr']:i for i,r in enumerate(effective)};audits.sort(key=lambda x:rank[x['donor_chr']])
    assert len(audits)==len(effective)==35
    write_tsv(OUT/'centromere_status_all35.tsv',audits)
    if unplaced:write_tsv(OUT/'unplaced_centromere_candidates.tsv',unplaced)
    write_tsv(OUT/'main_chromosome_centromere_status.tsv',[dict(chrom=c,drawn=c in selected,donor_chr=selected[c]['donor_chr'] if c in selected else 'NA',reason=selected[c]['evidence'] if c in selected else 'No uniquely supported candidate from supplied donor intervals; not fabricated') for c in order])
    result=dict(input_intervals=35,main_markers=len(selected),unresolved=[x for x in audits if x['status']!='drawn_main_chromosome'],main_chromosomes_without_marker=[c for c in order if c not in selected],coordinates='GRCg7b, BED 0-based half-open; no donor positions copied',approximate_not_CENPA_validated=True,chr25_donor_end_clipped_bp=759)
    (OUT/'centromere_selection_QC.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
