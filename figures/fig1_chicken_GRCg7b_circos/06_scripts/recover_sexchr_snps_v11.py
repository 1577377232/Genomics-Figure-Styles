#!/usr/bin/env python3
"""Targeted Z/W discovery from existing assembly BAMs; never realigns assemblies.

Reuses the original autosomal C++ caller without changing any quality filter.
Outputs a sites-only union (not population genotypes / allele frequencies).
Only Z/W adaptive-discovery fields are augmented in versioned window tables.
"""
import concurrent.futures as cf
import csv
import gzip
import hashlib
import heapq
import itertools
import json
from pathlib import Path
import subprocess
import sys
import pysam

PROJECT = Path(__file__).resolve().parents[1]
WORK = PROJECT / '03_variants/v11/sexchr_snps'
BASE = PROJECT.parent
SAMPLES = BASE.parent / '01_Pan_TE_Lib/Construct_Pan_TE/Assemblies_51/samples_51.tsv'
REF = Path('/data/CWD_1/00_Ref_TERef/7b_genome.fa')
ORIGINAL = BASE / '15_Homologous_SNP_phylogeny/scripts/extract_assembly_snps.cpp'
HTS = Path('/home/danding/yes/envs/tegq_longread')


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows):
    with open(path, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run_sample(row):
    prefix = WORK / 'per_sample' / row['sample']
    done = Path(str(prefix) + '.complete.json')
    fingerprint = {k: row[k] for k in ['sample', 'bam']}
    fingerprint.update(bam_size=Path(row['bam']).stat().st_size,
                       bam_mtime_ns=Path(row['bam']).stat().st_mtime_ns,
                       caller_sha256=digest(WORK / 'software/extract_sexchr_snps.cpp'))
    if done.exists() and json.loads(done.read_text()) == fingerprint:
        return row['sample']
    with open(WORK / 'logs' / (row['sample'] + '.log'), 'w') as log:
        subprocess.run([str(WORK / 'software/extract_sexchr_snps'), row['bam'], str(REF),
                        row['sample'], str(prefix)], stdout=log, stderr=log, check=True)
    pysam.tabix_index(str(prefix) + '.snps.vcf.gz', preset='vcf', force=True)
    pysam.tabix_index(str(prefix) + '.callable.bed.gz', preset='bed', force=True)
    done.write_text(json.dumps(fingerprint, indent=2))
    return row['sample']


def merge_intervals(iterable):
    start = end = None
    for a, b in iterable:
        if start is None:
            start, end = a, b
        elif a <= end:
            end = max(end, b)
        else:
            yield start, end
            start, end = a, b
    if start is not None:
        yield start, end


def main():
    for d in ['software', 'per_sample', 'logs']:
        (WORK / d).mkdir(parents=True, exist_ok=True)
    original = ORIGINAL.read_text()
    old = 'for(int c=1;c<=39;++c) {\n        string chr=to_string(c);'
    assert original.count(old) == 2
    source = original.replace(old, 'for (const string& chr: vector<string>{"Z", "W"}) {\n       ')
    source = source.replace('autosomes1-39;', 'sex_chromosomes_Z_W;')
    source = source.replace('assembly_cs_callable_v1', 'assembly_cs_callable_sexchr_v11')
    source = source.replace('Missing autosome', 'Missing sex chromosome')
    cpp = WORK / 'software/extract_sexchr_snps.cpp'
    cpp.write_text(source)
    command = ['g++', '-O3', '-std=c++17', '-I' + str(HTS / 'include'), str(cpp),
               '-L' + str(HTS / 'lib'), '-Wl,-rpath,' + str(HTS / 'lib'), '-lhts',
               '-o', str(WORK / 'software/extract_sexchr_snps')]
    subprocess.run(command, check=True)
    rows = read_tsv(SAMPLES)
    fasta = pysam.FastaFile(str(REF))
    lengths = {c: fasta.get_reference_length(c) for c in ['Z', 'W']}
    inventory = []
    for row in rows:
        with pysam.AlignmentFile(row['bam'], 'rb') as bam:
            idx = {s.contig: s.mapped for s in bam.get_index_statistics()}
            for c in lengths:
                assert bam.get_reference_length(c) == lengths[c]
                inventory.append(dict(sample=row['sample'], chromosome=c, bam=row['bam'],
                                      reference_bp=lengths[c], mapped_records=idx[c]))
    write_tsv(WORK / 'BAM_ZW_inventory.tsv', inventory)
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for sample in pool.map(run_sample, rows):
            print('completed', sample, flush=True)

    sample_qc = [q for row in rows for q in read_tsv(WORK / 'per_sample' / (row['sample'] + '.chromosome_qc.tsv'))]
    write_tsv(WORK / 'per_sample_chromosome_qc.tsv', sample_qc)
    vcfs = [pysam.TabixFile(str(WORK / 'per_sample' / (r['sample'] + '.snps.vcf.gz'))) for r in rows]
    beds = [pysam.TabixFile(str(WORK / 'per_sample' / (r['sample'] + '.callable.bed.gz'))) for r in rows]
    union = {}
    sites = {}
    summary = []
    vcf_path = WORK / 'assembly51.ZW.discovery.biallelic.sites.vcf.gz'
    bed_path = WORK / 'assembly51.ZW.callable_union.bed.gz'
    with pysam.BGZFile(str(vcf_path), 'w') as out, pysam.BGZFile(str(bed_path), 'w') as bout:
        header = ('##fileformat=VCFv4.2\n##source=V11_existing_assembly_cs_ZW_discovery\n'
                  '##reference=GRCg7b\n##FILTER=<ID=PASS,Description="Original callable assembly filters; unique biallelic union">\n'
                  '##INFO=<ID=ALT_ASM,Number=1,Type=Integer,Description="Number of assemblies supporting the alternate allele; not genotype call rate">\n')
        header += ''.join(f'##contig=<ID={c},length={n}>\n' for c, n in lengths.items())
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        out.write(header.encode())
        for chrom, length in lengths.items():
            def intervals(bed):
                if chrom not in bed.contigs:
                    return
                for line in bed.fetch(chrom):
                    _, a, b = line.split('\t')
                    yield int(a), int(b)
            union[chrom] = list(merge_intervals(heapq.merge(*(intervals(b) for b in beds))))
            for a, b in union[chrom]:
                bout.write(f'{chrom}\t{a}\t{b}\n'.encode())
            def variants(vcf, sample_index):
                if chrom not in vcf.contigs:
                    return
                for line in vcf.fetch(chrom):
                    z = line.split('\t')
                    assert z[6] == 'PASS' and z[9] == '1'
                    yield int(z[1]), z[3], z[4], sample_index
            merged = heapq.merge(*(variants(v, i) for i, v in enumerate(vcfs)))
            kept = []; multi = 0
            reference = fasta.fetch(chrom).upper()
            callable_index = 0
            for pos, group in itertools.groupby(merged, key=lambda x: x[0]):
                records = list(group)
                refs = {x[1] for x in records}; alts = {x[2] for x in records}
                assert len(refs) == 1 and next(iter(refs)) == reference[pos-1]
                if len(alts) != 1:
                    multi += 1
                    continue
                ref, alt = next(iter(refs)), next(iter(alts))
                assert ref in 'ACGT' and alt in 'ACGT' and ref != alt
                while union[chrom][callable_index][1] <= pos-1:
                    callable_index += 1
                assert union[chrom][callable_index][0] <= pos-1 < union[chrom][callable_index][1]
                support = len({x[3] for x in records})
                assert support == len(records)
                out.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tALT_ASM={support}\n'.encode())
                kept.append(pos-1)
            sites[chrom] = kept
            qc = [q for q in sample_qc if q['chrom'] == chrom]
            summary.append(dict(chromosome=chrom, length_bp=length,
                assemblies_screened=len(rows), assemblies_with_eligible_alignment=sum(int(q['eligible_alignments']) > 0 for q in qc),
                assemblies_with_callable_bases=sum(int(q['callable_bp']) > 0 for q in qc),
                callable_union_bp=sum(b-a for a,b in union[chrom]),
                unique_PASS_biallelic_SNPs=len(kept), excluded_multiallelic_positions=multi))
    for f in vcfs + beds:
        f.close()
    pysam.tabix_index(str(vcf_path), preset='vcf', force=True)
    pysam.tabix_index(str(bed_path), preset='bed', force=True)
    window_qc = []
    for kb in [100, 500, 1000]:
        data = read_tsv(PROJECT / f'04_window_statistics/v2/windows_{kb}kb.tsv')
        for d in data:
            d['SNP_v11_ZW_callable_bp'] = 'NA'
            if d['Chromosome'] not in lengths:
                continue
            chrom = d['Chromosome']; a = int(d['Start']); b = int(d['End'])
            import bisect
            n = bisect.bisect_left(sites[chrom], b) - bisect.bisect_left(sites[chrom], a)
            callable_bp = sum(max(0, min(b,y)-max(a,x)) for x,y in union[chrom] if x < b and y > a)
            d['SNP_adaptive_discovery_count'] = str(n) if callable_bp else 'NA'
            d['SNP_adaptive_discovery_per_Mb'] = str(n*1e6/(b-a)) if callable_bp else 'NA'
            d['SNP_adaptive_alignment_policy'] = '50kb_ZW_v11_existing_BAM'
            d['SNP_v11_ZW_callable_bp'] = str(callable_bp)
            window_qc.append(dict(window_kb=kb, chromosome=chrom, start=a, end=b,
                                 SNP_count=n if callable_bp else 'NA', callable_bp=callable_bp,
                                 state='observed' if callable_bp else 'unassessed'))
        outdir = PROJECT / '04_window_statistics/v11'
        outdir.mkdir(parents=True, exist_ok=True)
        write_tsv(outdir / f'windows_{kb}kb.tsv', data)
    write_tsv(WORK / 'chromosome_summary.tsv', summary)
    write_tsv(WORK / 'ZW_window_QC.tsv', window_qc)
    metadata = dict(reference=str(REF), reference_fai_sha256=digest(str(REF)+'.fai'),
        original_caller=str(ORIGINAL), original_caller_sha256=digest(ORIGINAL), modified_caller_sha256=digest(cpp),
        sample_manifest=str(SAMPLES), sample_manifest_sha256=digest(SAMPLES),
        filters='unchanged: competing MAPQ>=5/alignment>=10kb; eligible MAPQ>=20/alignment>=50kb/identity>=95%; unique reference coverage; 100bp end trim; indel flank5bp; ACGT; PASS; union biallelic',
        site_unit='unique reference position; reference allele compared to GRCg7b; sites-only union; no joint genotypes or allele frequencies',
        missingness='any callable bp required per window; otherwise NA; density denominator is physical window bp, not callable bp',
        sex_caveat='Assembly discovery, not sex-balanced population genotyping. W alignments may include homologous/misassembled sequence; no sex metadata or sample AF is inferred.',
        outputs_sha256={str(f.relative_to(PROJECT)):digest(f) for f in [vcf_path,bed_path]})
    (WORK / 'method.json').write_text(json.dumps(metadata, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
