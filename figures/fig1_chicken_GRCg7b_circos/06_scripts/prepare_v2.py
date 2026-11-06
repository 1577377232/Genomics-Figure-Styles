#!/usr/bin/env python3
"""Keep population QC, genome-wide discovery, and targeted recovery distinct."""
import csv
import math
from collections import Counter,defaultdict
import numpy as np
import pysam
from prepare_data import OUT,SNPBASE,REF,CHROMS,LENGTHS,SIZES,OFFSETS,arrays,add_point,tsv,write_json
from smallchr_rescue_v2 import TARGETS,DEST


def main():
    dest=OUT/'04_window_statistics/v2';dest.mkdir(exist_ok=True)
    discovery=arrays();ns_sum=arrays();counts=Counter();chromstats=defaultdict(Counter);ns_hist=Counter()
    path=SNPBASE/'06_merged/assembly52.union.snps.vcf.gz'
    with pysam.VariantFile(str(path),drop_samples=True) as v:
        for rec in v:
            counts['source_records']+=1
            if list(rec.filter.keys())!=['PASS'] or len(rec.alleles)!=2 or any(len(a)!=1 or a not in 'ACGT' for a in rec.alleles):continue
            ns=int(rec.info['NS']);assert ns>=2;counts['biallelic_discovery_sites']+=1;chromstats[rec.chrom]['discovery50']+=1;ns_hist[ns]+=1
            add_point(discovery,rec.chrom,rec.start)
            for size in SIZES:ns_sum[size][OFFSETS[size][rec.chrom]+rec.start//size]+=ns
            for n in [2,10,26,47,50,52]:chromstats[rec.chrom][f'NS_ge{n}']+=ns>=n
            if counts['source_records']%5000000==0:print('SNP discovery audit',counts['source_records'],flush=True)
    assert sum(discovery[100000])==counts['biallelic_discovery_sites']
    rescue=arrays();rescue95=arrays()
    with pysam.VariantFile(str(DEST/'smallchr_10kb.biallelic.discovery.vcf.gz'),drop_samples=True) as v:
        for rec in v:
            add_point(rescue,rec.chrom,rec.start)
            if rec.info['NS']>=50:add_point(rescue95,rec.chrom,rec.start)
    call=arrays();call95=arrays();fasta=pysam.FastaFile(str(REF));beds=[pysam.TabixFile(str(p)) for p in sorted(DEST.glob('*.callable.bed.gz'))];assert len(beds)==51
    for c in TARGETS:
        depth=np.zeros(LENGTHS[c],dtype=np.int16)
        for bed in beds:
            if c in bed.contigs:
                for line in bed.fetch(c):
                    _,a,b=line.split('\t');depth[int(a):int(b)]+=1
        s=np.frombuffer(fasta.fetch(c).upper().encode(),dtype=np.uint8);valid=(s==65)|(s==67)|(s==71)|(s==84)
        for size in SIZES:
            for a in range(0,LENGTHS[c],size):
                i=OFFSETS[size][c]+a//size;call[size][i]=np.count_nonzero((depth[a:a+size]>0)&valid[a:a+size]);call95[size][i]=np.count_nonzero((depth[a:a+size]>=49)&valid[a:a+size])
    reference=np.load(OUT/'01_reference/window_reference.npz');callv1=np.load(OUT/'03_variants/window_callability.npz');validstats=[]
    for size in SIZES:
        rows=list(csv.DictReader((OUT/f'04_window_statistics/windows_{size//1000}kb.tsv').open(),delimiter='\t'));header=list(rows[0])+['SNP_discovery50_count','SNP_discovery50_per_Mb','SNP_discovery50_mean_NS','SNP_rescue10_count','SNP_rescue10_per_Mb','SNP_rescue10_call95_count','SNP_rescue10_callable_bp','SNP_rescue10_callable95_bp','SNP_adaptive_discovery_count','SNP_adaptive_discovery_per_Mb','SNP_adaptive_alignment_policy']
        for r in rows:
            c=r['Chromosome'];a=int(r['Start']);i=OFFSETS[size][c]+a//size;bp=int(r['effective_window_bp']);nonref_sum=int(callv1[f'SNP_called_sample_bp_{size}'][i])-int(r['ACGT_bp']);observed=c not in ['Z','W'] and nonref_sum>0;n=int(discovery[size][i]);nr=int(rescue[size][i]);is_rescue=c in TARGETS
            assert n==0 or observed
            r.update(SNP_discovery50_count=n if observed else 'NA',SNP_discovery50_per_Mb=n*1e6/bp if observed else 'NA',SNP_discovery50_mean_NS=int(ns_sum[size][i])/n if n else 'NA',SNP_rescue10_count=nr if is_rescue and call[size][i] else 'NA',SNP_rescue10_per_Mb=nr*1e6/bp if is_rescue and call[size][i] else 'NA',SNP_rescue10_call95_count=int(rescue95[size][i]) if is_rescue and call95[size][i] else 'NA',SNP_rescue10_callable_bp=int(call[size][i]) if is_rescue else 'NA',SNP_rescue10_callable95_bp=int(call95[size][i]) if is_rescue else 'NA')
            adaptive=r['SNP_rescue10_count'] if is_rescue else r['SNP_discovery50_count']
            r.update(SNP_adaptive_discovery_count=adaptive,SNP_adaptive_discovery_per_Mb=int(adaptive)*1e6/bp if adaptive!='NA' else 'NA',SNP_adaptive_alignment_policy='10kb_on_chr_shorter_than_5Mb' if is_rescue else '50kb' if c not in ['Z','W'] else 'unassessed')
        tsv(dest/f'windows_{size//1000}kb.tsv',header,[list(r.values()) for r in rows])
        assert sum(int(r['SNP_discovery50_count']) for r in rows if r['SNP_discovery50_count']!='NA')==counts['biallelic_discovery_sites']
        assert sum(int(r['SNP_rescue10_count']) for r in rows if r['SNP_rescue10_count']!='NA')==int(sum(rescue[size]))
        validstats.append({'window_bp':size,'population95_sites':sum(int(r['SNP_count']) for r in rows if r['SNP_count']!='NA'),'discovery50_sites':counts['biallelic_discovery_sites'],'rescue10_sites':int(sum(rescue[size])),'adaptive_discovery_sites':sum(int(r['SNP_adaptive_discovery_count']) for r in rows if r['SNP_adaptive_discovery_count']!='NA'),'population95_NA_windows':sum(r['SNP_count']=='NA' for r in rows),'discovery50_NA_windows':sum(r['SNP_discovery50_count']=='NA' for r in rows),'population95_nonzero_windows':sum(r['SNP_count']!='NA' and int(r['SNP_count'])>0 for r in rows),'all_windows':len(rows)})
    tsv(OUT/'05_QC/v2/SNP_NS_histogram.tsv',['NS_called_assemblies','biallelic_sites'],sorted(ns_hist.items()))
    tsv(OUT/'05_QC/v2/SNP_callrate_threshold_by_chromosome.tsv',['chrom','discovery50','NS_ge2','NS_ge10','NS_ge26','NS_ge47','NS_ge50','NS_ge52'],[[c,*[chromstats[c][k] for k in ['discovery50','NS_ge2','NS_ge10','NS_ge26','NS_ge47','NS_ge50','NS_ge52']]] for c in CHROMS[:39]])
    assert len({r['adaptive_discovery_sites'] for r in validstats})==1
    write_json(OUT/'05_QC/v2/SNP_density_audit.json',{'counts':dict(counts),'window_statistics':validstats,'physical_reference_main_bp':sum(LENGTHS[c] for c in CHROMS),'population95_percent_of_main_reference':100*9741807/sum(LENGTHS[c] for c in CHROMS),'discovery50_percent_of_main_reference':100*counts['biallelic_discovery_sites']/sum(LENGTHS[c] for c in CHROMS),'panels':'population95 uses original NS>=50/52; discovery50 uses all PASS biallelic assembly SNPs with original 50kb alignment quality; rescue10 is a separate targeted 10kb supplement; adaptive_discovery explicitly uses 10kb for all reference chromosomes shorter than 5Mb and 50kb otherwise','original_uniform_calling_panels_preserved':True,'adaptive_panel_exploratory_and_explicitly_labeled':True,'Z_W':'not called; not fabricated'})
    print('v2 SNP datasets',dict(counts),flush=True)


if __name__=='__main__':main()
