#!/usr/bin/env python3
"""Read-only analysis of project inputs; all writes stay below OUTPUT.

Dependencies: Python 3, numpy, pysam. Run individual stages or --stage all.
Coordinates in output tables/BED are 0-based half-open unless POS is explicit.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

OUT = Path(__file__).resolve().parents[1]
PROJECT = OUT.parent
ASSEMBLIES = PROJECT.parent / '01_Pan_TE_Lib/Construct_Pan_TE/Assemblies_51'
SNPBASE = PROJECT / '15_Homologous_SNP_phylogeny'
REF = Path('/data/CWD_1/00_Ref_TERef/7b_genome.fa')
RM = PROJECT / 'GRCg7b/RepeatMasker_out/GRCg7b.fa.out'
DETECT = ASSEMBLIES / '3_TE_detection'
CHROMS = [str(i) for i in range(1,40)] + ['Z','W']
CLASSES = ['LINE','LTR','SINE','DNA']
SIZES = [100_000,500_000,1_000_000]
LENGTHS = {f[0]:int(f[1]) for f in (line.split() for line in Path(str(REF)+'.fai').open())}
OFFSETS = {}
for size in SIZES:
    n=0; OFFSETS[size]={}
    for chrom in CHROMS:
        OFFSETS[size][chrom]=n; n+=math.ceil(LENGTHS[chrom]/size)


def write_json(path, data):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')


def tsv(path, header, rows):
    with path.open('w') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()


def arrays():
    return {s:np.zeros(sum(math.ceil(LENGTHS[c]/s) for c in CHROMS),dtype=np.int64) for s in SIZES}


def add_point(a,chrom,pos):
    if chrom not in CHROMS or not 0<=pos<LENGTHS[chrom]:return
    for s in SIZES:a[s][OFFSETS[s][chrom]+pos//s]+=1


def add_span(a,chrom,start,end):
    for size in SIZES:
        for w in range(start//size,(end-1)//size+1):
            a[size][OFFSETS[size][chrom]+w]+=max(0,min(end,(w+1)*size)-max(start,w*size))


def merged(intervals):
    end=-1;start=-1
    for a,b in sorted(intervals):
        if a>end:
            if end>=0:yield start,end
            start,end=a,b
        else:end=max(end,b)
    if end>=0:yield start,end


def union_bp(intervals):return sum(b-a for a,b in merged(intervals))


def category(raw):
    top=raw.split('/')[0]
    if top in ['RC','Helitron'] or 'helitron' in raw.lower():return 'Helitron_RC'
    if top in CLASSES:return top
    if top=='MITE':return 'DNA'
    return top or 'Unclassified'


def rmrows(path):
    if not path.exists():return
    with path.open() as f:
        for line in f:
            p=line.split()
            if len(p)>=15 and p[5].isdigit() and p[6].isdigit():yield p


def inventory():
    rows=[];vcfs=[];seen=set()
    roots=[PROJECT,ASSEMBLIES,ASSEMBLIES.parent/'Assemblies_51_New',Path('/data/CWD_1/00_Ref_TERef')]
    extensions=('.vcf','.vcf.gz','.bcf','.bed','.bed.gz','.gff','.gff3','.gtf','.gff.gz','.gff3.gz','.gtf.gz','.fai','.fa','.fasta','.fna','.fa.gz','.out','.tbl','.paf','.chain','.chain.gz','.bam','.tsv','.csv','.json','.md')
    for root in roots:
        for directory,dirs,files in os.walk(root,followlinks=False):
            dirs[:]=[x for x in dirs if x not in ['software','softwares','.git','__pycache__','node_modules',OUT.name]]
            for name in files:
                p=Path(directory)/name
                if p in seen or not name.endswith(extensions):continue
                seen.add(p)
                try:
                    st=p.stat();resolved=p.resolve(); status='historical_invalid_exclude' if 'invalid_GRCg6a' in str(p) else 'candidate'
                    rows.append([str(p),str(resolved),st.st_size,st.st_mtime_ns,p.suffix,p.is_symlink(),status])
                    if name.endswith(('.vcf','.vcf.gz','.bcf')):
                        try:
                            with pysam.VariantFile(str(p)) as v:
                                samples=list(v.header.samples); contigs={k:x.length for k,x in v.header.contigs.items()}
                                vcf={'path':str(p),'samples':samples,'sample_count':len(samples),'contigs':contigs,'header_reference':[l for l in str(v.header).splitlines() if l.startswith(('##reference','##source'))]}
                                vcfs.append(vcf)
                        except Exception as e:vcfs.append({'path':str(p),'error':str(e)})
                except OSError as e:rows.append([str(p),'',0,0,'',p.is_symlink(),str(e)])
    tsv(OUT/'00_data_inventory/files.tsv',['path','resolved_path','bytes','mtime_ns','suffix','symlink','status'],rows)
    write_json(OUT/'00_data_inventory/vcf_headers.json',vcfs)
    # File fingerprints are restricted to directly used inputs, not huge unrelated BAMs.
    paths=[REF,Path(str(REF)+'.fai'),RM,Path('/data/CWD_1/00_Ref_TERef/7b_genome.fa.mod.EDTA.TElib.fa'),
           ASSEMBLIES/'2_SV_merge_comparison/merged_jasmine_paper.vcf.gz',ASSEMBLIES/'2_SV_merge_comparison/merged_survivor_previous.vcf.gz',
           SNPBASE/'06_merged/assembly52.union.snps.vcf.gz',SNPBASE/'07_selection/call95.ld02/qc.bim']
    write_json(OUT/'00_data_inventory/input_fingerprints.json',{str(p):{'sha256':sha(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in paths})
    print('inventory',len(rows),'files;',len(vcfs),'VCF headers',flush=True)


def reference():
    a=arrays();nonacgt=arrays();rows=[]
    fasta=pysam.FastaFile(str(REF))
    for chrom in CHROMS:
        length=LENGTHS[chrom];seq=fasta.fetch(chrom).upper().encode()
        for size in SIZES:
            for start in range(0,length,size):
                part=seq[start:start+size];idx=OFFSETS[size][chrom]+start//size
                a[size][idx]=len(part);nonacgt[size][idx]=len(part)-sum(part.count(x) for x in [b'A',b'C',b'G',b'T'])
        rows.append([chrom,'Chr'+chrom,length,'main'])
    rows += [[c,c,n,'MT_excluded' if c=='MT' else 'unplaced_excluded'] for c,n in LENGTHS.items() if c not in CHROMS]
    tsv(OUT/'01_reference/chromosomes.tsv',['chrom','display_label','length_bp','plot_status'],rows)
    with (OUT/'01_reference/main.chrom.sizes').open('w') as f:
        for c in CHROMS:f.write(f'{c}\t{LENGTHS[c]}\n')
    np.savez_compressed(OUT/'01_reference/window_reference.npz',**{f'{k}_{s}':v for k,x in [('effective_bp',a),('nonACGT_bp',nonacgt)] for s,v in x.items()})
    # All dictionaries are checked against every contig length present in each header.
    checks=[]
    for label,path in [('SNP',SNPBASE/'06_merged/assembly52.union.snps.vcf.gz'),('Jasmine',ASSEMBLIES/'2_SV_merge_comparison/merged_jasmine_paper.vcf.gz'),('SURVIVOR',ASSEMBLIES/'2_SV_merge_comparison/merged_survivor_previous.vcf.gz')]:
        with pysam.VariantFile(str(path)) as f:
            bad={c:[x.length,LENGTHS.get(c)] for c,x in f.header.contigs.items() if x.length is not None and x.length!=LENGTHS.get(c)}
            assert not bad,(label,bad)
            checks.append({'dataset':label,'contigs':len(f.header.contigs),'length_mismatches':bad,'missing_main_contigs':[c for c in CHROMS if c not in f.header.contigs]})
    write_json(OUT/'05_QC/reference_dictionary_checks.json',checks)
    print('reference',len(CHROMS),'chromosomes',sum(LENGTHS[c] for c in CHROMS),'bp',flush=True)


def te():
    spans=defaultdict(list);copies=defaultdict(lambda:defaultdict(list));class_raw=Counter();outside=Counter();bad=[]
    with gzip.open(OUT/'02_TE_annotation/RepeatMasker.normalized.bed.gz','wt') as f:
        f.write('#chrom\tstart\tend\tname\tscore\tstrand\tclass\traw_class\tRM_ID\n')
        for p in rmrows(RM):
            chrom=p[4];a=int(p[5])-1;b=int(p[6]);cat=category(p[10]);class_raw[(cat,p[10])]+=1
            if chrom not in LENGTHS or a<0 or b>LENGTHS[chrom] or b<=a:bad.append(p);continue
            if chrom not in CHROMS:outside[cat]+=1;continue
            spans[(chrom,cat)].append((a,b));copies[cat][(chrom,p[14],p[9],p[10])].append((a,b))
            f.write('\t'.join(map(str,[chrom,a,b,p[9],p[0],'-' if p[8]=='C' else '+',cat,p[10],p[14]]))+'\n')
    assert not bad,bad[:5]
    result={};summary=[]
    for cat in sorted(copies):
        covered=arrays();copycount=arrays();fragmentcount=sum(len(v) for (c,k),v in spans.items() if k==cat)
        with gzip.open(OUT/f'02_TE_annotation/{cat}.union.bed.gz','wt') as f:
            for c in CHROMS:
                for a,b in merged(spans[(c,cat)]):
                    add_span(covered,c,a,b);f.write(f'{c}\t{a}\t{b}\n')
        # A reconstructed RepeatMasker ID is an operational copy, not proof of a full-length insertion.
        with gzip.open(OUT/f'02_TE_annotation/{cat}.copies.bed.gz','wt') as f:
            for (c,rid,name,raw),iv in copies[cat].items():
                a=min(x[0] for x in iv);b=max(x[1] for x in iv);mid=(a+b-1)//2;add_point(copycount,c,mid)
                f.write(f'{c}\t{a}\t{b}\t{name}\t{rid}\t{raw}\t{len(iv)}\n')
        for s in SIZES:result[f'{cat}_covered_bp_{s}']=covered[s];result[f'{cat}_copy_count_{s}']=copycount[s]
        summary.append([cat,fragmentcount,len(copies[cat]),sum(covered[100000]),outside[cat]])
    np.savez_compressed(OUT/'02_TE_annotation/window_te.npz',**result)
    tsv(OUT/'02_TE_annotation/class_summary.tsv',['class','main_RM_fragments','main_RM_ID_copies','main_union_bp','excluded_contig_fragments'],summary)
    tsv(OUT/'02_TE_annotation/classification_mapping.tsv',['normalized_class','raw_class','all_contig_fragments'],[[a,b,n] for (a,b),n in sorted(class_raw.items())])
    write_json(OUT/'05_QC/te_annotation_qc.json',{'source':str(RM),'out_of_bounds':len(bad),'copy_definition':'unique (chromosome, RepeatMasker ID, repeat name, raw class), assigned once by span midpoint','coverage_definition':'within-class interval union; inter-class overlaps retained separately','classification':'MITE to DNA; Helitron/RC separate; Unknown and non-TE repeats separate'})
    print('reference TE',summary,flush=True)


def repeats():
    """Reclassify ALL raw RM matches, including records rejected by the old 80% filter.

    Use actually aligned intervals, not the span between OneCode-joined fragments.
    Existing query FASTA includes the VCF anchor; remove query coordinate 0.
    Any DEL with a non-anchor REF mismatch is excluded from a validated TE call.
    """
    report=json.loads((DETECT/'repeat_annotation_summary.json').read_text())
    te_classes={'LINE','LTR','SINE','DNA','MITE','RC','Helitron'}
    nonte={'Simple_repeat','Low_complexity','Satellite'}
    rows=[];counts=Counter();checks=[];fasta=pysam.FastaFile(str(REF));source_fingerprints={}
    for batch in report['results']:
        directory=Path(batch['directory']);rawpath=directory/'repeatmasker_dir/indels.fa.out'
        hits=defaultdict(lambda:defaultdict(list))
        for p in rmrows(rawpath):hits[p[4]][p[10].split('/')[0]].append((int(p[5])-1,int(p[6])))
        source_fingerprints[str(rawpath)]={'sha256':sha(rawpath),'bytes':rawpath.stat().st_size} if rawpath.exists() else {'absent_no_hits':True}
        # Batch query lengths prove no records were omitted by the old filtering.
        qlens={p[0]:int(p[1]) for p in (l.split() for l in (directory/'indels.fa.fai').open())} if (directory/'indels.fa.fai').exists() else {}
        chrom=batch['chrom'];refseq=fasta.fetch(chrom).upper();n=0
        with pysam.VariantFile(str(directory/'genotypes.vcf')) as v:
            for rec in v:
                n+=1;typ=rec.info['SVTYPE'];ref=rec.ref.upper();alt=rec.alts[0].upper();length=abs(len(ref)-len(alt));query=alt if typ=='INS' else ref
                assert typ in ['INS','DEL'] and rec.id in qlens and len(query)==qlens[rec.id]
                assert min(len(ref),len(alt))==1,(rec.id,'not a simple anchored indel')
                assert 0<=rec.start<LENGTHS[chrom]
                actual=refseq[rec.start:rec.start+len(ref)]
                refmatch=(ref==actual);interior_match=(typ=='INS' or ref[1:]==actual[1:])
                counts['REF_exact_match' if refmatch else 'REF_mismatch']+=1
                if not interior_match:counts['DEL_internal_REF_mismatch_excluded']+=1
                cl={c:[(max(1,a)-1,min(len(query),b)-1) for a,b in iv if min(len(query),b)>max(1,a)] for c,iv in hits[rec.id].items()}
                bp={c:union_bp(iv) for c,iv in cl.items()}
                strict_bp=union_bp([x for c,iv in cl.items() if c in te_classes for x in iv])
                other_bp=union_bp([x for c,iv in cl.items() if c in nonte for x in iv])
                expanded_bp=union_bp([x for c,iv in cl.items() if c in te_classes|nonte for x in iv])
                all_bp=union_bp([x for iv in cl.values() for x in iv])
                strict=(strict_bp*2>length and length>50 and interior_match)
                extended=(expanded_bp*2>length and length>50 and interior_match)
                dominant=max((c for c in bp if c in te_classes|nonte),key=lambda c:bp[c],default='None')
                counts['assessed']+=1;counts['strict_gt50']+=strict;counts['extended_gt50']+=extended
                counts['strict_gt80']+=(strict_bp/length>.8 and length>50 and interior_match)
                counts['extended_gt80']+=(expanded_bp/length>.8 and length>50 and interior_match)
                rows.append([chrom,rec.pos,rec.id,typ,length,strict_bp,strict_bp/length,other_bp,expanded_bp,expanded_bp/length,
                             all_bp,dominant,','.join(sorted(bp)),int(strict),int(extended),int(refmatch),int(interior_match),
                             int(expanded_bp>50),batch['batch']])
        assert n==batch['records'],(batch['batch'],n,batch['records'])
        checks.append([batch['batch'],chrom,n,len(qlens),sum(len(v) for v in hits.values())])
    assert counts['assessed']==154228
    with gzip.open(OUT/'03_variants/INS_DEL.repeat_reclassified.tsv.gz','wt') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n')
        w.writerow(['chrom','POS_1based','ID','SVTYPE','event_length','strict_TE_bp','strict_TE_fraction','other_repeat_bp','expanded_repeat_bp','expanded_repeat_fraction','all_RM_bp','dominant_class','matched_classes','strict_TE_gt50','expanded_repeat_gt50','original_REF_exact','DEL_internal_REF_valid','expanded_match_bp_gt50','batch']);w.writerows(rows)
    write_json(OUT/'05_QC/repeat_reclassification.json',{'counts':dict(counts),'query_denominator':'event sequence without VCF anchor','coverage':'union of raw RepeatMasker aligned query intervals; clipped to non-anchor event sequence','TE_classes':sorted(te_classes),'other_repeat_classes':sorted(nonte),'Unknown':'excluded from strict and expanded sets; all_RM_bp retained for sensitivity','threshold':'event_length > 50 bp AND union_coverage > 0.50, exact integer comparison; DEL interior must match GRCg7b','upper_length_limit_bp':100000,'unassessed_types':['DUP:TANDEM','DUP:INT','INV','BND'],'sample_GT_warning':'single ref/1|0 in historical annotation files is synthetic catalog encoding, never population genotypes','source':str(DETECT),'source_raw_fingerprints':source_fingerprints})
    tsv(OUT/'05_QC/repeat_batch_completeness.tsv',['batch','chrom','VCF_records','query_FASTA_records','queries_with_hits'],checks)
    print('repeat reclassification',dict(counts),flush=True)


def sv():
    with gzip.open(OUT/'03_variants/INS_DEL.repeat_reclassified.tsv.gz','rt') as h:
        annotations={r['ID']:r for r in csv.DictReader(h,delimiter='\t')}
    outputs={};summaries={};sample_rows=[];class_counts=Counter();fasta=pysam.FastaFile(str(REF))
    for label,rel in [('Jasmine','2_SV_merge_comparison/merged_jasmine_paper.vcf.gz'),('SURVIVOR','2_SV_merge_comparison/merged_survivor_previous.vcf.gz')]:
        path=ASSEMBLIES/rel;v=pysam.VariantFile(str(path));names=list(v.header.samples);nmiss=np.zeros(len(names),dtype=np.int64);ncall=np.zeros(len(names),dtype=np.int64)
        counts=Counter();seen=set();seenkeys=set();filtered=Counter();tracks={k:arrays() for k in ['all','nonTE_assessed','nonrepeat_assessed','unassessed','TIP_strict','TIP_extended','TIP_strict_INS','TIP_strict_DEL','TIP_extended_INS','TIP_extended_DEL']};bnd=[];refseq=None;lastchrom=None
        with gzip.open(OUT/f'03_variants/{label}.sites.tsv.gz','wt') as h:
            w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['chrom','POS_1based','END_1based_inclusive','midpoint_0based','ID','cluster','SVTYPE','SVLEN','FILTER','SUPP','GT_missing_fraction','annotation_status','strict_TE_fraction','expanded_repeat_fraction','TIP_strict','TIP_extended'])
            for rec in v:
                counts['records']+=1;typ=rec.info.get('SVTYPE','unknown');counts['type_'+typ]+=1
                counts['FILTER_'+','.join(rec.filter.keys())]+=1
                if list(rec.filter.keys())!=['PASS']:filtered['not_PASS']+=1;continue
                if rec.chrom not in CHROMS:filtered['outside_main']+=1;continue
                cluster=rec.info.get('MERGE_CLUSTER',rec.id)
                if cluster in seen:filtered['duplicate_cluster']+=1;continue
                seen.add(cluster)
                # Jasmine/SURVIVOR clusters define loci; exact identical alleles cannot count twice.
                key=(rec.chrom,rec.pos,typ,rec.stop,hashlib.sha256(str(rec.alleles).encode()).hexdigest())
                if key in seenkeys:filtered['duplicate_exact_allele']+=1;continue
                seenkeys.add(key)
                if typ=='BND':bnd.append([rec.chrom,rec.pos,rec.id,str(rec.alts),cluster]);continue
                n=int(rec.info.get('SVLEN',0));n=abs(n)
                if not 50<n<=100000:filtered['length_outside_51_100000']+=1;continue
                assert rec.chrom in LENGTHS and 0<=rec.start<rec.stop<=LENGTHS[rec.chrom],str(rec)[:200]
                # INS uses its breakpoint base. Other event intervals exclude the left VCF anchor.
                midpoint=rec.start if typ=='INS' else (rec.pos+rec.stop-1)//2
                a=annotations.get(rec.id) if label=='Jasmine' else None
                if label=='Jasmine' and typ in ['INS','DEL']:assert a is not None,(rec.id,'missing raw annotation')
                strict=bool(a and a['strict_TE_gt50']=='1');extended=bool(a and a['expanded_repeat_gt50']=='1')
                status='assessed' if a and a['DEL_internal_REF_valid']=='1' else 'REF_mismatch_unresolved' if a else 'unassessed'
                gt=[x.get('GT',()) for x in rec.samples.values()];missing=np.array([not g or any(i is None for i in g) for g in gt]);nmiss+=missing;ncall+=1
                counts['main_unique_nonBND']+=1;counts['main_'+typ]+=1
                add_point(tracks['all'],rec.chrom,midpoint)
                if status=='assessed':
                    if not strict:add_point(tracks['nonTE_assessed'],rec.chrom,midpoint)
                    if not extended:add_point(tracks['nonrepeat_assessed'],rec.chrom,midpoint)
                else:add_point(tracks['unassessed'],rec.chrom,midpoint)
                if strict:
                    add_point(tracks['TIP_strict'],rec.chrom,midpoint);add_point(tracks['TIP_strict_'+typ],rec.chrom,midpoint)
                if extended:
                    add_point(tracks['TIP_extended'],rec.chrom,midpoint);add_point(tracks['TIP_extended_'+typ],rec.chrom,midpoint)
                    class_counts[(typ,a['dominant_class'],'strict_TE' if strict else 'expanded_only')]+=1
                w.writerow([rec.chrom,rec.pos,rec.stop,midpoint,rec.id,cluster,typ,n,','.join(rec.filter.keys()),rec.info.get('SUPP','NA'),float(np.mean(missing)),status,a['strict_TE_fraction'] if a else 'NA',a['expanded_repeat_fraction'] if a else 'NA',int(strict),int(extended)])
        tsv(OUT/f'03_variants/{label}.BND.excluded.tsv',['chrom','POS_1based','ID','ALT','cluster'],bnd)
        for name,a in tracks.items():
            for size in SIZES:outputs[f'{label}_{name}_{size}']=a[size]
        sample_rows += [[label,s,int(n),int(d),float(n/d) if d else 'NA'] for s,n,d in zip(names,nmiss,ncall)]
        summaries[label]={'path':str(path),'assembly_sample_columns':len(names),'counts':dict(counts),'excluded':dict(filtered),'tracks_main_chromosomes':{k:int(sum(a[100000])) for k,a in tracks.items()},'BND_main_excluded':len(bnd),'cluster_definition':'independent merge cluster; exact duplicate alleles removed; multi-allelic/nearby independent clusters retained','genotypes':'discovery genotypes; absent assembly call remains ./., not a proven reference allele'}
        print(label,summaries[label]['tracks_main_chromosomes'],flush=True)
    np.savez_compressed(OUT/'03_variants/window_sv.npz',**outputs)
    write_json(OUT/'05_QC/sv_panels.json',summaries)
    tsv(OUT/'05_QC/SV_discovery_GT_missingness.tsv',['panel','sample','missing_GT','total_main_nonBND_loci','missing_fraction'],sample_rows)
    tsv(OUT/'03_variants/TIP_class_counts.tsv',['SVTYPE','dominant_class','subset','sites'],[[a,b,c,n] for (a,b,c),n in sorted(class_counts.items())])


def snp():
    path=SNPBASE/'06_merged/assembly52.union.snps.vcf.gz';v=pysam.VariantFile(str(path),drop_samples=True)
    raw=arrays();qcarr=arrays();counts=Counter();filters=Counter();chromcounts=Counter();maf_bins=Counter();last=None;posref=None
    qc_hash=hashlib.sha256();bim_hash=hashlib.sha256();fasta=pysam.FastaFile(str(REF));nsamples=52
    with gzip.open(OUT/'03_variants/SNP.PASS.biallelic.call95.maf01.positions.tsv.gz','wt') as h:
        h.write('chrom\tstart\tend\tREF\tALT\tNS\tAN\tAF\tF_MISSING\n')
        for rec in v:
            counts['records']+=1;filters[','.join(rec.filter.keys())]+=1
            key=(rec.chrom,rec.pos)
            if key==last:counts['duplicate_position_records']+=1;continue
            if last and rec.chrom==last[0]:assert rec.pos>last[1]
            last=key
            assert rec.chrom in CHROMS and rec.chrom not in ['Z','W']
            if posref is None or posref[0]!=rec.chrom:posref=(rec.chrom,fasta.fetch(rec.chrom).upper())
            if rec.ref.upper()!=posref[1][rec.start:rec.start+1]:counts['REF_mismatch']+=1;continue
            if list(rec.filter.keys())!=['PASS']:counts['not_PASS']+=1;continue
            is_snp=len(rec.ref)==1 and rec.ref in 'ACGT' and all(len(a)==1 and a in 'ACGT' for a in rec.alts)
            if not is_snp:counts['not_pure_SNP']+=1;continue
            add_point(raw,rec.chrom,rec.start)
            if len(rec.alts)!=1:counts['multiallelic']+=1;continue
            ac=rec.info['AC'];ac=ac[0] if isinstance(ac,tuple) else ac;an=rec.info['AN'];ns=rec.info['NS'];af=ac/an;fm=1-ns/nsamples
            # Integer inequalities avoid rounding at the requested boundaries.
            if ns*100<95*nsamples:counts['callrate_lt95']+=1;continue
            if min(ac,an-ac)*100<an:counts['MAF_lt01']+=1;continue
            counts['QC_sites']+=1;chromcounts[rec.chrom]+=1;add_point(qcarr,rec.chrom,rec.start)
            qc_hash.update(f'{rec.chrom}\t{rec.pos}\n'.encode())
            h.write(f'{rec.chrom}\t{rec.start}\t{rec.pos}\t{rec.ref}\t{rec.alts[0]}\t{ns}\t{an}\t{af:.9g}\t{fm:.9g}\n')
            if counts['records']%2_000_000==0:print('SNP scanned',counts['records'],flush=True)
    n=0
    with (SNPBASE/'07_selection/call95.ld02/qc.bim').open() as h:
        for line in h:
            p=line.split();bim_hash.update(f'{p[0]}\t{p[3]}\n'.encode());n+=1
    assert qc_hash.hexdigest()==bim_hash.hexdigest() and n==counts['QC_sites'],('unpruned SNP QC does not match existing BIM',counts,n)
    assert not counts['REF_mismatch'] and not counts['duplicate_position_records']
    np.savez_compressed(OUT/'03_variants/window_snp.npz',**{f'{label}_{s}':a[s] for label,a in [('SNP_union',raw),('SNP_QC',qcarr)] for s in SIZES})
    # Compute sample missingness directly from the existing UNPRUNED PLINK BED.
    fam=[x.split()[1] for x in (SNPBASE/'07_selection/call95.ld02/qc.fam').read_text().splitlines()]
    assert len(fam)==nsamples
    bedpath=SNPBASE/'07_selection/call95.ld02/qc.bed'
    with bedpath.open('rb') as h:assert h.read(3)==b'\x6c\x1b\x01'
    b=np.memmap(bedpath,dtype='uint8',mode='r',offset=3,shape=(n,math.ceil(nsamples/4)))
    missing=[int(np.count_nonzero(((b[:,j//4]>>(2*(j%4)))&3)==1)) for j in range(nsamples)]
    tsv(OUT/'05_QC/SNP_unpruned_sample_missingness.tsv',['sample','missing_sites','total_sites','missing_fraction'],[[s,m,n,m/n] for s,m in zip(fam,missing)])
    write_json(OUT/'05_QC/snp_qc.json',{'source':str(path),'reference':str(REF),'assemblies_including_reference':nsamples,'assembly_vs_individual':'52 assemblies including reference; some assemblies are haplotypes of the same individual','counts':dict(counts),'FILTER_counts':dict(filters),'main_QC_by_chromosome':dict(chromcounts),'criteria':'PASS; ACGT SNP; exactly 2 alleles; NS/52 >=0.95; min(AC,AN-AC)/AN >=0.01; no LD pruning','match_unpruned_existing_qc_bim':True,'position_sha256':qc_hash.hexdigest(),'missing_chromosomes':['Z','W'],'missing_chromosome_policy':'NA, never zero','quality_note':'assembly alignment filters, not sequencing-depth/GQ filters; QUAL is missing; see original extraction methods'})
    print('SNP',dict(counts),flush=True)


def callability():
    """Population-callable sequence at the same NS>=50/52 threshold as SNP QC."""
    manifest=[l.split() for l in (SNPBASE/'06_merged/callable_masks.tsv').read_text().splitlines()]
    assert len(manifest)==51
    handles=[pysam.TabixFile(p) for _,p in manifest];fasta=pysam.FastaFile(str(REF));call=arrays();summed=arrays();summary=[];checked=0
    selected=gzip.open(OUT/'03_variants/SNP.PASS.biallelic.call95.maf01.positions.tsv.gz','rt');next(selected);line=next(selected,None)
    with gzip.open(OUT/'03_variants/SNP.population_callable_ge95.bed.gz','wt') as bed:
        for chrom in CHROMS[:39]:
            length=LENGTHS[chrom];delta=np.zeros(length+1,dtype=np.int16)
            for handle in handles:
                if chrom not in handle.contigs:continue
                iv=np.array([l.split()[1:3] for l in handle.fetch(chrom)],dtype=np.int64)
                if iv.size:
                    assert np.all(iv[:,0]>=0) and np.all(iv[:,1]<=length)
                    np.add.at(delta,iv[:,0],1);np.add.at(delta,iv[:,1],-1)
            np.cumsum(delta,dtype=np.int16,out=delta);d=delta[:-1]
            assert 0<=np.min(d)<=np.max(d)<=51
            seq=np.frombuffer(fasta.fetch(chrom).upper().encode(),dtype=np.uint8);canonical=(seq==65)|(seq==67)|(seq==71)|(seq==84)
            d+=canonical.astype(np.int16);mask=(d>=50)&canonical
            # Every selected SNP must fall inside the jointly callable sequence.
            while line is not None and line.split('\t',1)[0]==chrom:
                p=line.split('\t',3);assert mask[int(p[1])],(chrom,p[1],'SNP outside >=50/52 callable sequence');checked+=1;line=next(selected,None)
            edges=np.flatnonzero(np.diff(np.r_[np.int8(0),mask.view(np.int8),np.int8(0)]))
            for a,b in edges.reshape(-1,2):bed.write(f'{chrom}\t{a}\t{b}\n')
            for size in SIZES:
                for a in range(0,length,size):
                    idx=OFFSETS[size][chrom]+a//size;call[size][idx]=np.count_nonzero(mask[a:a+size]);summed[size][idx]=np.sum(d[a:a+size],dtype=np.int64)
            summary.append([chrom,length,int(np.count_nonzero(mask)),float(np.mean(d))/52])
            print('callability',chrom,summary[-1][2],flush=True)
            del d,delta,mask,canonical,seq
    assert line is None;selected.close()
    assert checked==json.loads((OUT/'05_QC/snp_qc.json').read_text())['counts']['QC_sites']
    np.savez_compressed(OUT/'03_variants/window_callability.npz',**{f'{k}_{s}':a[s] for k,a in [('SNP_callable95_bp',call),('SNP_called_sample_bp',summed)] for s in SIZES})
    tsv(OUT/'05_QC/SNP_population_callability.tsv',['chrom','length_bp','callable_ge50of52_bp','mean_sample_callable_fraction'],summary)
    write_json(OUT/'00_data_inventory/SNP_callable_mask_fingerprints.json',{sample:{'path':p,'sha256':sha(p),'bytes':Path(p).stat().st_size,'mtime_ns':Path(p).stat().st_mtime_ns} for sample,p in manifest})
    write_json(OUT/'05_QC/SNP_callability_check.json',{'samples':52,'nonreference_masks':51,'criterion':'at least 50/52 assemblies callable, matching NS/52>=.95; reference callable iff A/C/G/T','all_selected_SNPs_within_joint_callable_mask':checked,'Z_W':'unassessed','SNP_display_policy':'no population-callable bases -> NA; passing QC site count=0 retained in SNP_retained_count','mask_sources':dict(manifest)})
    for handle in handles:handle.close()


def windows():
    caches=[np.load(OUT/p) for p in ['01_reference/window_reference.npz','02_TE_annotation/window_te.npz','03_variants/window_snp.npz','03_variants/window_sv.npz','03_variants/window_callability.npz']]
    data={k:z[k] for z in caches for k in z.files};meta=[]
    for size in SIZES:
        rows=[];header=['Chromosome','Start','End','Midpoint','effective_window_bp','nonACGT_bp','ACGT_bp','SNP_callable95_bp','SNP_callable95_fraction','SNP_mean_sample_callable_fraction','SNP_retained_count','SNP_per_callable_Mb']
        for c in CLASSES:header += [f'{c}_covered_bp',f'{c}_coverage_fraction',f'{c}_copy_count',f'{c}_copies_per_Mb']
        variants={'SNP':'SNP_QC','SNP_union':'SNP_union','SV':'Jasmine_all','SV_SURVIVOR':'SURVIVOR_all','SV_nonTE_assessed':'Jasmine_nonTE_assessed','SV_nonrepeat_assessed':'Jasmine_nonrepeat_assessed','SV_TE_unassessed':'Jasmine_unassessed','TIP_strict':'Jasmine_TIP_strict','TIP_extended':'Jasmine_TIP_extended','TIP_strict_INS':'Jasmine_TIP_strict_INS','TIP_strict_DEL':'Jasmine_TIP_strict_DEL','TIP_extended_INS':'Jasmine_TIP_extended_INS','TIP_extended_DEL':'Jasmine_TIP_extended_DEL'}
        for k in variants:header += [f'{k}_count',f'{k}_per_Mb']
        for chrom in CHROMS:
            for start in range(0,LENGTHS[chrom],size):
                i=OFFSETS[size][chrom]+start//size;end=min(start+size,LENGTHS[chrom]);bp=end-start
                cbp=int(data[f'SNP_callable95_bp_{size}'][i]);missing=chrom in ['Z','W'];snpcount=int(data[f'SNP_QC_{size}'][i])
                row=[chrom,start,end,(start+end)/2,bp,int(data[f'nonACGT_bp_{size}'][i]),int(bp-data[f'nonACGT_bp_{size}'][i]),'NA' if missing else cbp,'NA' if missing else cbp/bp,'NA' if missing else int(data[f'SNP_called_sample_bp_{size}'][i])/(52*bp),'NA' if missing else snpcount,'NA' if missing or cbp==0 else snpcount*1e6/cbp]
                for c in CLASSES:
                    cover=int(data[f'{c}_covered_bp_{size}'][i]);copy=int(data[f'{c}_copy_count_{size}'][i]);assert 0<=cover<=bp
                    row += [cover,cover/bp,copy,copy*1e6/bp]
                for k,cachekey in variants.items():
                    n=int(data[f'{cachekey}_{size}'][i]);row += ['NA','NA'] if k.startswith('SNP') and (missing or (k=='SNP' and cbp==0)) else [n,n*1e6/bp]
                rows.append(row)
        tsv(OUT/f'04_window_statistics/windows_{size//1000}kb.tsv',header,rows)
        with (OUT/f'04_window_statistics/windows_{size//1000}kb.bed').open('w') as h:
            for row in rows:h.write('\t'.join(map(str,row[:3]))+'\n')
        for c in CLASSES:
            col=header.index(f'{c}_coverage_fraction');values=np.array([r[col] for r in rows]);meta.append([size,c,len(rows),float(np.mean(values==0)),float(np.median(values)),float(np.quantile(values,.95)),float(max(values))])
        for k in ['SNP','SV','TIP_extended','TIP_strict']:
            col=header.index(k+'_per_Mb');values=np.array([r[col] for r in rows if r[col]!='NA']);meta.append([size,k,len(values),float(np.mean(values==0)),float(np.median(values)),float(np.quantile(values,.95)),float(max(values))])
        print('windows',size,len(rows),flush=True)
    tsv(OUT/'05_QC/window_size_comparison.tsv',['window_bp','track','observed_windows','zero_fraction','median','p95','max'],meta)


def qc():
    files=[list(csv.DictReader((OUT/f'04_window_statistics/windows_{s//1000}kb.tsv').open(),delimiter='\t')) for s in SIZES]
    checks=[]
    for size,rows in zip(SIZES,files):
        for c in CHROMS:
            r=[x for x in rows if x['Chromosome']==c]
            assert int(r[0]['Start'])==0 and int(r[-1]['End'])==LENGTHS[c]
            assert all(int(a['End'])==int(b['Start']) for a,b in zip(r,r[1:]))
        assert sum(int(r['effective_window_bp']) for r in rows)==sum(LENGTHS[c] for c in CHROMS)
        for r in rows:
            if r['Chromosome'] in ['Z','W']:assert r['SNP_count']=='NA'
            if r['SNP_callable95_bp']=='0':assert r['SNP_count']=='NA' and r['SNP_retained_count']=='0'
            assert int(r['TIP_strict_count'])<=int(r['TIP_extended_count'])<=int(r['SV_count'])
            assert int(r['SV_nonTE_assessed_count'])+int(r['TIP_strict_count'])+int(r['SV_TE_unassessed_count'])==int(r['SV_count'])
        for key in ['SNP_count','SV_count','SV_SURVIVOR_count','TIP_strict_count','TIP_extended_count']+[f'{c}_covered_bp' for c in CLASSES]+[f'{c}_copy_count' for c in CLASSES]:
            totals=[sum(int(r[key]) for r in rows if r[key]!='NA') for rows in files];assert len(set(totals))==1,(key,totals)
        checks.append({'window_bp':size,'windows':len(rows),'complete_41_chromosomes':True,'coordinate_partition':True,'cross_resolution_total_conservation':True,'TIP_subset_and_SV_partition':True})
    expected=json.loads((OUT/'05_QC/snp_qc.json').read_text())['counts']['QC_sites']
    assert sum(int(r['SNP_count']) for r in files[-1] if r['SNP_count']!='NA')==expected
    source=json.loads((OUT/'05_QC/sv_panels.json').read_text())
    for column,panel,track in [('SV_count','Jasmine','all'),('SV_SURVIVOR_count','SURVIVOR','all'),('TIP_strict_count','Jasmine','TIP_strict'),('TIP_extended_count','Jasmine','TIP_extended')]:
        assert sum(int(r[column]) for r in files[-1])==source[panel]['tracks_main_chromosomes'][track]
    chromrows=[]
    for c in CHROMS:
        r=[x for x in files[-1] if x['Chromosome']==c];row=[c,LENGTHS[c],sum(int(x['nonACGT_bp']) for x in r)]
        for k in ['SNP_count','SV_count','SV_SURVIVOR_count','TIP_strict_count','TIP_extended_count']:row.append('NA' if k=='SNP_count' and all(x[k]=='NA' for x in r) else sum(int(x[k]) for x in r if x[k]!='NA'))
        chromrows.append(row)
    tsv(OUT/'05_QC/chromosome_track_totals.tsv',['chrom','length_bp','nonACGT_bp','SNP','SV_Jasmine','SV_SURVIVOR','TIP_strict','TIP_extended'],chromrows)
    write_json(OUT/'05_QC/window_integrity.json',checks)
    print('QC passed',flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['inventory','reference','te','repeats','sv','snp','callability','windows','qc','all'],default='all');args=p.parse_args()
    for name in ['inventory','reference','te','repeats','sv','snp','callability','windows','qc']:
        if args.stage in [name,'all']:globals()[name]()


if __name__=='__main__':main()
