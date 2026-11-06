#!/usr/bin/env python3
"""Inventory GGswu/GRCg7b homology and unique flanking anchors from existing BAM.

Read only original FASTA/BAM. Segment suffixes are gap-split counters, not offsets.
"""
from pathlib import Path
import csv
import re
import json
from collections import defaultdict
import pysam

P = Path(__file__).resolve().parents[1]
OUT = P / '01_reference/v14'
RAW = Path('/media/danding/MyBook/Pan_TE_chicken/03_chicken_assemblies/chicken_assemblies/NCBI/GCA_024206055.2_GGswu_genomic.fna')
BAM = Path('/data/CWD_1/Pan_TE_chicken/01_Pan_TE_Lib/Construct_Pan_TE/Assemblies_51/00_bams/Ggswu.sorted.bam')


def read_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(p, rows):
    with open(p, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def records(path):
    head = None; parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(parts)
                head = line.strip(); parts = []
            else:
                parts.append(line.strip())
        if head is not None:
            yield head, ''.join(parts)


def union_bp(intervals):
    total = 0; end = -1
    for a,b in sorted(intervals):
        total += max(0,b-max(a,end)); end=max(end,b)
    return total


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    supplied = {r['donor_chr']:(int(r['used_start']),int(r['used_end'])) for r in read_tsv(P/'01_reference/v2/GGswu.centromeres.effective.tsv')}
    segments={}; donor_lengths={}
    for head,seq in records(RAW):
        match=re.search(r'chromosome\s+([0-9ZW]+)(?:\s|$)',head)
        if not match:continue
        chrom=match[1]; donor_lengths[chrom]=len(seq);acc=head[1:].split()[0]
        for i,m in enumerate(re.finditer('[^N]+',seq)):
            segments[f'{acc}_{i}']=(chrom,m.start(),m.end())
    print('reconstructed',len(segments),'segments',flush=True)
    groups=defaultdict(list);near=[];alignment_count=0;full=[];primary=[]
    with pysam.AlignmentFile(str(BAM),'rb') as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.mapping_quality<30 or r.query_name not in segments:continue
            c,offset,end=segments[r.query_name];qlen=r.infer_read_length()
            assert qlen==end-offset,(r.query_name,qlen,end-offset)
            aln=sum(n for op,n in r.cigartuples if op in [0,1,2,7,8])
            if not r.is_supplementary:
                primary.append(dict(source_chr=c,query=r.query_name,target_chr=r.reference_name,
                    target_start=r.reference_start,target_end=r.reference_end,
                    strand='-' if r.is_reverse else '+',MAPQ=r.mapping_quality,
                    alignment_bp=aln,identity_NM=1-r.get_tag('NM')/aln))
            if aln<10000 or 1-r.get_tag('NM')/aln<.90:continue
            alignment_count+=1;q=0;t=r.reference_start;strand='-' if r.is_reverse else '+'
            hard_left=r.cigartuples[0][1] if r.cigartuples[0][0]==5 else 0
            qa=r.query_alignment_start+hard_left;qb=r.query_alignment_end+hard_left
            full.append(dict(source_chr=c,query=r.query_name,source_start=offset+(qlen-qb if r.is_reverse else qa),source_end=offset+(qlen-qa if r.is_reverse else qb),target_chr=r.reference_name,target_start=t,target_end=r.reference_end,strand=strand,MAPQ=r.mapping_quality,alignment_bp=aln,flag=r.flag))
            for op,n in r.cigartuples:
                if op in [0,7,8]:
                    a,b=(offset+qlen-q-n,offset+qlen-q) if r.is_reverse else (offset+q,offset+q+n)
                    if n>=200:
                        groups[(c,r.reference_name,strand)].append((a,b))
                        if c in supplied:
                            ca,cb=supplied[c]
                            if b>=ca-2000000 and a<=cb+2000000:
                                near.append(dict(source_chr=c,source_start=a,source_end=b,target_chr=r.reference_name,target_start=t,target_end=t+n,strand=strand,MAPQ=r.mapping_quality,query=r.query_name,alignment_bp=aln,flag=r.flag))
                    q+=n;t+=n
                elif op in [1,4,5]:q+=n
                elif op in [2,3]:t+=n
            assert q==qlen
    summaries=[dict(source_chr=c,target_chr=t,strand=s,union_source_bp=union_bp(v),source_length=donor_lengths[c]) for (c,t,s),v in groups.items()]
    summaries.sort(key=lambda x:(int(x['source_chr']) if x['source_chr'].isdigit() else 100+ord(x['source_chr']),-x['union_source_bp']))
    write_tsv(OUT/'whole_chromosome_correspondence.tsv',summaries)
    write_tsv(OUT/'centromere_flanking_M_blocks.tsv',near)
    write_tsv(OUT/'eligible_alignment_spans.tsv',full)
    write_tsv(OUT/'primary_correspondence.tsv',primary)
    (OUT/'anchor_method.json').write_text(json.dumps(dict(donor=str(RAW),BAM=str(BAM),new_alignment=False,reconstructed_segments=len(segments),eligible_alignments=alignment_count,nearby_M_blocks=len(near),
        flank_filter='nonsecondary MAPQ>=30; alignment>=10kb; identity>=90%; individual M blocks>=200bp; nearest flank selection requires M>=1000bp',
        primary_correspondence_filter='primary nonsecondary MAPQ>=30; NM identity reported, not used to discard primary correspondence; NM includes indels',
        existing_centromere_core_filter='previous exact CIGAR M projections; nonsecondary MAPQ>=30; no additional NM-identity filter',
        query_coordinate_handling='reconstructed gap-split offsets, full hard clips and reverse orientation explicitly handled'),indent=2))
    for c in supplied:
        ranked=[x for x in summaries if x['source_chr']==c][:3]
        print(c,[(x['target_chr'],x['strand'],x['union_source_bp']) for x in ranked],flush=True)


if __name__=='__main__':main()
