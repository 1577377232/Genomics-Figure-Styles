#!/usr/bin/env python3
"""Independent CIGAR/query-sequence check of real SNPs; does not parse cs tags."""
import bisect
import csv
import gzip
import json
from pathlib import Path
import pysam

PROJECT = Path(__file__).resolve().parents[1]
WORK = PROJECT / '03_variants/v11/sexchr_snps'
MANIFEST = PROJECT.parent.parent / '01_Pan_TE_Lib/Construct_Pan_TE/Assemblies_51/samples_51.tsv'


def main():
    with open(MANIFEST) as f:
        paths = {r['sample']: r['bam'] for r in csv.DictReader(f, delimiter='\t')}
    reports = []
    for sample in ['Asil', 'GG1F', 'GRCg7w', 'Silkie1', 'SLWX2', 'Ggswu']:
        variants = {'Z': [], 'W': []}
        with gzip.open(WORK / 'per_sample' / (sample + '.snps.vcf.gz'), 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                z = line.split('\t')
                variants[z[0]].append((int(z[1])-1, z[3], z[4]))
        with pysam.AlignmentFile(paths[sample], 'rb') as bam:
            for c, all_sites in variants.items():
                step = max(1, len(all_sites)//250)
                selected = all_sites[::step]
                positions = [v[0] for v in selected]
                expected = {v[0]: v[2] for v in selected}
                observed = {}
                for a in bam.fetch(c):
                    if a.flag & (4|256|512|1024) or a.mapping_quality < 20:
                        continue
                    aligned = sum(n for op, n in a.cigartuples if op in [0, 1, 2, 7, 8])
                    if aligned < 50000 or 1-a.get_tag('NM')/aligned < .95:
                        continue
                    left = bisect.bisect_left(positions, a.reference_start)
                    right = bisect.bisect_left(positions, a.reference_end)
                    if left == right:
                        continue
                    sequence = a.query_sequence.upper()
                    ref = a.reference_start; query = 0
                    for op, n in a.cigartuples:
                        if op in [0, 7, 8]:
                            i = bisect.bisect_left(positions, ref)
                            j = bisect.bisect_left(positions, ref+n)
                            for pos in positions[i:j]:
                                assert pos not in observed, (sample, c, pos, 'multiple eligible alignments')
                                observed[pos] = sequence[query+pos-ref]
                            ref += n; query += n
                        elif op in [1, 4]:
                            query += n
                        elif op in [2, 3]:
                            ref += n
                assert observed == expected, (sample, c, len(observed), len(expected))
                reports.append(dict(sample=sample, chromosome=c, available_SNPs=len(all_sites),
                                    deterministic_spread_subset=len(selected), CIGAR_ALT_matches=len(observed)))
    result = dict(status='PASS', method='Independent BAM CIGAR and query base; no cs parsing; evenly spaced real SNP subset', checks=reports)
    (WORK / 'independent_CIGAR_validation.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
