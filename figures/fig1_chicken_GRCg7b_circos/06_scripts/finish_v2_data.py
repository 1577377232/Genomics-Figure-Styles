#!/usr/bin/env python3
"""Make the explicitly labeled exploratory view without altering source panels."""
import csv
from prepare_data import OUT,CHROMS,LENGTHS,SIZES,tsv,write_json
from smallchr_rescue_v2 import TARGETS
import json


def main():
    assert set(TARGETS)=={c for c in CHROMS[:39] if LENGTHS[c]<5000000}
    totals={}
    for size in SIZES:
        p=OUT/f'04_window_statistics/v2/windows_{size//1000}kb.tsv';rows=list(csv.DictReader(p.open(),delimiter='\t'));header=list(rows[0])
        cols=['SNP_adaptive_discovery_count','SNP_adaptive_discovery_per_Mb','SNP_adaptive_alignment_policy']
        for col in cols:
            if col not in header:header.append(col)
        for r in rows:
            c=r['Chromosome'];n=r['SNP_rescue10_count'] if c in TARGETS else r['SNP_discovery50_count'];bp=int(r['effective_window_bp'])
            r.update(SNP_adaptive_discovery_count=n,SNP_adaptive_discovery_per_Mb=int(n)*1e6/bp if n!='NA' else 'NA',SNP_adaptive_alignment_policy='10kb_on_chr_shorter_than_5Mb' if c in TARGETS else '50kb' if c not in ['Z','W'] else 'unassessed')
        tsv(p,header,[[r[k] for k in header] for r in rows]);totals[size]=sum(int(r[cols[0]]) for r in rows if r[cols[0]]!='NA')
    assert len(set(totals.values()))==1
    p=OUT/'05_QC/v2/SNP_density_audit.json';audit=json.loads(p.read_text());audit.pop('no_cohort_or_calling_filter_mixture',None)
    for r in audit['window_statistics']:r['adaptive_discovery_sites']=totals[r['window_bp']]
    audit.update(original_uniform_calling_panels_preserved=True,adaptive_panel_exploratory_and_explicitly_labeled=True,adaptive_policy='All reference autosomes <5Mb use the targeted 10kb calling panel; other autosomes retain the 50kb panel. Site call rate is not required for this exploratory discovery view.',adaptive_source_assignment={c:'targeted_10kb' if c in TARGETS else 'original_50kb' if c not in ['Z','W'] else 'unassessed' for c in CHROMS})
    write_json(p,audit);print('adaptive exploratory view',totals[1000000],'sites',flush=True)


if __name__=='__main__':main()
