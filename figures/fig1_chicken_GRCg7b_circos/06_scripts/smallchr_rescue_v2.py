#!/usr/bin/env python3
"""Targeted SNP extraction from existing BAMs; original files are read-only.

The sole calling change is min aligned length 50 kb -> 10 kb on the declared
small-chromosome panel. This is an exploratory supplement, never silently
substituted into the original 95%-call-rate or genome-wide density panel.
"""
import concurrent.futures
import csv
import json
import subprocess
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pysam
from prepare_data import OUT, SNPBASE, ASSEMBLIES, REF, LENGTHS, sha, tsv, write_json

TARGETS=['16','22','25']+[str(i) for i in range(29,40)]
DEST=OUT/'03_variants/v2/smallchr_rescue_10kb'
ENV=Path('/home/danding/yes/envs/tegq_longread')


def main():
    DEST.mkdir(parents=True,exist_ok=True);(DEST/'software').mkdir(exist_ok=True)
    source=(SNPBASE/'scripts/extract_assembly_snps.cpp').read_text()
    source=source.replace('for(int c=1;c<=39;++c)', 'for(int c: {16,22,25,29,30,31,32,33,34,35,36,37,38,39})')
    source=source.replace('if(aligned<50000||','if(aligned<10000||').replace('call_mapq20_min_aligned50000','call_mapq20_min_aligned10000').replace('autosomes1-39','small_chromosomes_16_22_25_29-39').replace('assembly_cs_callable_v1','circos_small_chr_10kb_v2')
    cpp=DEST/'software/extract_smallchr_10kb.cpp';cpp.write_text(source);exe=DEST/'software/extract_smallchr_10kb'
    command=['g++','-O3','-std=c++17','-I'+str(ENV/'include'),str(cpp),'-L'+str(ENV/'lib'),'-Wl,-rpath,'+str(ENV/'lib'),'-lhts','-o',str(exe)]
    subprocess.run(command,check=True)
    rows=list(csv.DictReader((ASSEMBLIES/'samples_51.tsv').open(),delimiter='\t'));assert len(rows)==51
    codehash=sha(cpp)
    def run(row):
        name=row['sample'];prefix=DEST/name;done=DEST/f'{name}.done.json';fingerprint={'code':codehash,'BAM':row['bam'],'BAM_mtime_ns':Path(row['bam']).stat().st_mtime_ns}
        if done.exists() and json.loads(done.read_text())['fingerprint']==fingerprint:return name
        cmd=[str(exe),row['bam'],str(REF),name,str(prefix)]
        with (DEST/f'{name}.log').open('w') as h:subprocess.run(cmd,stdout=h,stderr=h,check=True)
        for suffix,preset in [('.snps.vcf.gz','vcf'),('.callable.bed.gz','bed')]:pysam.tabix_index(str(prefix)+suffix,preset=preset,force=True)
        write_json(done,{'fingerprint':fingerprint,'command':cmd,'status':'complete'})
        return name
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for i,name in enumerate(pool.map(run,rows),1):
            if i%10==0 or i==51:print('small-chromosome extraction',i,'/51',flush=True)
    samples=[r['sample'] for r in rows];vfiles=[pysam.VariantFile(str(DEST/f'{s}.snps.vcf.gz')) for s in samples];beds=[pysam.TabixFile(str(DEST/f'{s}.callable.bed.gz')) for s in samples]
    header=pysam.VariantHeader();header.add_line('##source=small_chromosome_10kb_exploratory_recovery')
    header.add_line('##calling=MAPQ20_min_aligned10000_identity0.95; competing_MAPQ5_min10000; unique_reference_coverage; endtrim100; indelflank5; existing_asm5_BAMs')
    header.add_line('##scope=16,22,25,29-39; do not mix with genome-wide 50kb calling without disclosing different ascertainment')
    for c in TARGETS:header.contigs.add(c,length=LENGTHS[c])
    for name,num,typ,desc in [('NS',1,'Integer','Called assemblies including reference'),('AN',1,'Integer','Called haploid alleles'),('AC','A','Integer','Alternate assembly allele count'),('MQMIN',1,'Integer','Minimum mapping quality among ALT-supporting assembly calls')]:header.info.add(name,num,typ,desc)
    header.formats.add('GT',1,'String','Haploid assembly allele');[header.add_sample(s) for s in samples+['GRCg7b']]
    fasta=pysam.FastaFile(str(REF));summary=[];independent_checks=[];total=0
    outfile=DEST/'smallchr_10kb.biallelic.discovery.vcf.gz'
    with pysam.VariantFile(str(outfile),'wz',header=header) as out:
        for chrom in TARGETS:
            length=LENGTHS[chrom];bits=np.zeros(length,dtype=np.uint64);variants=defaultdict(dict)
            for i,(v,bed) in enumerate(zip(vfiles,beds)):
                if chrom in bed.contigs:
                    for line in bed.fetch(chrom):
                        _,a,b=line.split('\t');bits[int(a):int(b)] |= np.uint64(1<<i)
                for rec in v.fetch(chrom):variants[rec.start][i]=(rec.ref,rec.alts[0],int(rec.info['MQ']))
            counts=Counter();positions=np.array(sorted(variants),dtype=np.int64);reference=fasta.fetch(chrom).upper()
            for pos in positions:
                calls=variants[int(pos)];alts={a for r,a,mq in calls.values()};refs={r for r,a,mq in calls.values()};assert refs=={reference[pos]}
                if len(alts)!=1:counts['multiallelic_excluded']+=1;continue
                n=int(bits[pos]).bit_count()+1;ac=len(calls);assert n>ac>=1
                rec=out.new_record(contig=chrom,start=int(pos),alleles=(reference[pos],next(iter(alts))));rec.filter.add('PASS');rec.info['NS']=n;rec.info['AN']=n;rec.info['AC']=(ac,);rec.info['MQMIN']=min(x[2] for x in calls.values())
                for i,s in enumerate(samples):rec.samples[s]['GT']=(1,) if i in calls else (0,) if int(bits[pos])&(1<<i) else (None,)
                rec.samples['GRCg7b']['GT']=(0,);out.write(rec);counts['biallelic_discovery']+=1;counts['callrate_ge95']+=n>=50;total+=1
            # Independent support check: read CIGAR + query bases, rather than the cs parser.
            if chrom in ['35','36','39']:
                for i in range(min(51,len(samples))):
                    chosen=[int(p) for p in positions if i in variants[int(p)]][:4]
                    if not chosen:continue
                    with pysam.AlignmentFile(rows[i]['bam'],'rb') as bam:
                        for pos in chosen:
                            expected=variants[pos][i][1];supported=False
                            for r in bam.fetch(chrom,pos,pos+1):
                                if r.is_secondary or r.is_unmapped or r.mapping_quality<20:continue
                                aligned=sum(n for op,n in r.cigartuples if op in [0,1,2,7,8])
                                if aligned<10000 or 1-r.get_tag('NM')/aligned<.95:continue
                                q=0;t=r.reference_start
                                for op,n in r.cigartuples:
                                    if op in [0,7,8]:
                                        if t<=pos<t+n and r.query_sequence[q+pos-t].upper()==expected:supported=True
                                        q+=n;t+=n
                                    elif op in [1,4]:q+=n
                                    elif op in [2,3]:t+=n
                            assert supported,(samples[i],chrom,pos,expected)
                            independent_checks.append([samples[i],chrom,pos+1,expected,'CIGAR_query_base_support_pass'])
            called=bits!=0;summary.append([chrom,length,int(np.count_nonzero(called)),len(variants),counts['biallelic_discovery'],counts['callrate_ge95'],int(np.max([int(v).bit_count()+1 for v in bits[positions]])) if len(positions) else 0])
            print('rescue SNP',chrom,summary[-1][4],flush=True)
    pysam.tabix_index(str(outfile),preset='vcf',force=True)
    tsv(OUT/'05_QC/v2/smallchr_rescue_summary.tsv',['chrom','length_bp','any_nonreference_callable_bp','union_candidate_positions','biallelic_discovery_sites','callrate_ge95_sites','max_NS_at_candidate_sites'],summary)
    tsv(OUT/'05_QC/v2/rescue_CIGAR_base_validation.tsv',['sample','chrom','POS_1based','ALT','check'],independent_checks)
    write_json(OUT/'05_QC/v2/smallchr_rescue_provenance.json',{'compile_command':command,'source_sha256':codehash,'original_source_sha256':sha(SNPBASE/'scripts/extract_assembly_snps.cpp'),'target_chromosomes':TARGETS,'extraction_threads':4,'whole_genome_alignment_started':False,'existing_BAMs_only':True,'biallelic_discovery_sites':total,'independent_CIGAR_base_checks':len(independent_checks),'calling_change':'aligned block minimum reduced from 50kb to 10kb; all other original quality/uniqueness masks retained','status':'exploratory separate supplement; not a replacement for NS>=50/52 population calls'})


if __name__=='__main__':main()
