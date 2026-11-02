#!/usr/bin/env python3
"""Validate joint Manhattan delivery and build a local workbook/gallery."""
import argparse,csv,hashlib,html,json,math,re,struct,subprocess
from datetime import datetime,timezone
from pathlib import Path

def rows(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def run(*a):return subprocess.check_output(a,text=True)
def png_info(p):
    with Path(p).open('rb') as f:
        assert f.read(8)==b'\x89PNG\r\n\x1a\n';size=None;dpi=None
        while True:
            n,tag=struct.unpack('>I4s',f.read(8))
            if tag in [b'IHDR',b'pHYs']:
                b=f.read(n)
                if tag==b'IHDR':size=struct.unpack('>II',b[:8])
                else:x,y,u=struct.unpack('>IIB',b);assert u==1;dpi=(x*.0254,y*.0254)
            else:f.seek(n,1)
            f.read(4)
            if tag==b'IEND':break
    assert size and dpi
    return dict(width=size[0],height=size[1],dpi_x=dpi[0],dpi_y=dpi[1])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,required=True);ap.add_argument('--script',type=Path,default=Path(__file__).with_name('plot_joint_manhattan.R'));a=ap.parse_args();out=a.outdir.resolve()
    cfg=json.loads((out/'config/analysis.json').read_text());style=json.loads((out/cfg['display_style_file']).read_text())
    traits=rows(out/'config/traits.tsv');sources=rows(out/'config/sources.tsv')
    assert len(traits)==cfg['expected_count'] and len({x['trait_id'] for x in traits})==len(traits)
    expected={'single':cfg['expected_singles'],'module':cfg['expected_modules']}
    qcs=[];source_qc=[];chrom_qc=[]
    if cfg.get('Z_analysis')=='joint':
        assert len(sources)==4*len(traits) and {s['panel'] for s in sources}=={'A','Z_joint'}
    for t in traits:
        tid=t['trait_id'];q=json.loads((out/'audit'/f'{tid}__PASS.json').read_text());assert q['status']=='PASS'
        assert q.get('Z_analysis')==cfg.get('Z_analysis','unspecified')
        assert q['analysis_kind']==t['kind'] and q['p_column']==('p_score' if t['kind']=='module' else 'p_wald')
        assert not q['chr38_plotted'] and not q['threshold_label'] and not q['TE_type_distinction'] and not q['Z_sex_distinction']
        assert q['SNP_point_size_mm']==q['TE_point_size_mm']==style['snp_size_mm']
        assert q['SNP_display_shape']==q['TE_display_shape']==16
        assert q['threshold_line_count']==1 and q['display_threshold_p']==style['threshold_p']
        assert q['SNP_colour']==style['colours']['SNP'] and q['TE_colour']==style['colours']['TE']
        assert q['clipped_peaks']==0 and q['P_zero']==0 and not q['downsampling']
        ss=rows(out/'audit'/f'{tid}__sources.tsv');cc=rows(out/'audit'/f'{tid}__chromosome_counts.tsv')
        assert all(r['chr'] in {str(i) for i in range(1,38)}|{'Z'} for r in cc)
        for marker in ['SNP','TE']:
            s=[x for x in ss if x['marker']==marker]
            if cfg.get('Z_analysis')=='joint':assert len(s)==2 and {x['panel'] for x in s}=={'A','Z_joint'} and len({x['N'] for x in s})==1
            native=sum(int(x['valid']) for x in s);drawn=sum(int(x['plotted']) for x in s)
            dropped=sum(int(x['excluded_chr38']) for x in s)
            assert native==q['M_'+marker] and drawn==q[marker+'_plotted'] and native==drawn+dropped
            assert dropped==q[marker+'_excluded_chr38']
            assert drawn==sum(int(x['valid_tests']) for x in cc if x['marker']==marker)
            assert math.isclose(min(float(x['minP']) for x in s),q['minP_'+marker],rel_tol=1e-13)
            assert math.isclose(min(float(x['minP_plotted']) for x in s),q['minP_plotted_'+marker],rel_tol=1e-13)
        snapshot={(r['marker'],r['panel']):r for r in sources if r['trait_id']==tid}
        for r in ss:
            z=snapshot[(r['marker'],r['panel'])]
            assert r['source_sha256']==z['expected_sha256'] and int(r['valid'])==int(z['expected_valid'])
        for ext in ['png','pdf']:assert sha(q[ext])==q[ext+'_sha256']
        image=png_info(q['png']);assert abs(image['dpi_x']-600)<.05 and abs(image['dpi_y']-600)<.05
        assert abs(image['width']/600*25.4-style['width_mm'])<.1
        assert abs(image['height']/600*25.4-style['height_mm'])<.1
        text=run('pdftotext','-layout',q['pdf'],'-')
        assert tid in text and 'Chromosome' in text and re.search(r'\bZ\b',text)
        assert not re.search(r'P\s*=',text) and 'female' not in text and 'male' not in text
        assert 'TE insertion' not in text and 'TE deletion' not in text
        fonts=run('pdffonts',q['pdf']);assert 'LiberationSans' in fonts and 'yes' in fonts
        images=run('pdfimages','-list',q['pdf'])
        assert sum(bool(re.match(r'^\s*\d+\s+\d+\s+image\s',line)) for line in images.splitlines())==1
        qcs.append(q);source_qc.extend(ss);chrom_qc.extend(cc)
    for kind,n in expected.items():
        assert sum(q['analysis_kind']==kind for q in qcs)==n
        folder=out/('modules' if kind=='module' else 'traits')
        assert len(list(folder.glob('*.pdf')))==len(list(folder.glob('*.png')))==n
        label='modules' if kind=='module' else 'quantitative_traits'
        book=out/'overview'/f'SNP_TE_{n}_{label}_booklet.pdf'
        pages=int(re.search(r'^Pages:\s+(\d+)',run('pdfinfo',str(book)),re.M)[1]);assert pages==n
    fail=list((out/'audit').glob('*__FAIL.json'));assert not fail
    # Optional workbook; all result values and input hashes remain local.
    import pandas as pd
    frame=pd.DataFrame(qcs)
    with pd.ExcelWriter(out/'plot_catalog.xlsx',engine='openpyxl') as w:
        frame.to_excel(w,sheet_name='joint_plots',index=False)
        pd.DataFrame(source_qc).to_excel(w,sheet_name='source_QC',index=False)
        pd.DataFrame(chrom_qc).to_excel(w,sheet_name='chromosome_counts',index=False)
        if (out/'audit/Z_coverage.tsv').exists():pd.read_csv(out/'audit/Z_coverage.tsv',sep='\t').to_excel(w,sheet_name='Z_coverage',index=False)
        for s in w.book:s.freeze_panes='A2';s.auto_filter.ref=s.dimensions
    cards=[]
    for q in qcs:
        png=Path(q['png']).relative_to(out);pdf=Path(q['pdf']).relative_to(out)
        name=html.escape(q['trait_id']+' '+q['display_name_zh'])
        cards.append(f'<article data-kind="{q["analysis_kind"]}"><h3>{name}</h3><a href="{png}"><img loading="lazy" src="{png}" alt="{name}"></a><a href="{pdf}">PDF</a></article>')
    page='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>SNP + TE GEMMA</title><style>body{font:16px sans-serif;margin:24px;color:#202830}article{margin:28px 0;border-top:1px solid #ddd}img{width:100%;height:auto}input{padding:8px;min-width:300px}button{padding:8px;margin:8px}</style><h1>SNP + TE：89数量性状与7模块</h1><p>蓝色SNP、珊瑚色TE；Z采用雌雄联合GWAS，每个位点一个P值。1–37号和Z，单条参考线无数值标签。</p><p><a href="overview/SNP_TE_89_quantitative_traits_booklet.pdf">89性状图册</a> · <a href="overview/SNP_TE_7_modules_booklet.pdf">7模块图册</a> · <a href="plot_catalog.xlsx">完整目录</a></p><input id="query" placeholder="按ID或性状名称筛选"><button onclick="kind='' ;filter()">全部</button><button onclick="kind='single';filter()">数量性状</button><button onclick="kind='module';filter()">模块</button>'''+''.join(cards)+'''<script>let kind='';function filter(){let q=document.getElementById('query').value.toLowerCase();document.querySelectorAll('article').forEach(e=>e.style.display=(!kind||e.dataset.kind===kind)&&e.textContent.toLowerCase().includes(q)?'':'none')}document.getElementById('query').oninput=filter</script></html>'''
    (out/'index.html').write_text(page)
    records=[]
    for folder in ['traits','modules','overview','scripts','config']:
        records.extend(p for p in (out/folder).rglob('*') if p.is_file())
    with (out/'audit/delivery_SHA256.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['file','bytes','sha256'])
        for p in sorted(records):w.writerow([str(p.relative_to(out)),p.stat().st_size,sha(p)])
    result=dict(status='PASS',verified_utc=datetime.now(timezone.utc).isoformat(),joint_plots=len(qcs),
                single_traits=expected['single'],multivariate_modules=expected['module'],source_files=len(source_qc),
                PNG_dpi=600,point_size_mm=style['snp_size_mm'],colours=style['colours'],
                Z_sex_display=cfg.get('Z_display','identical aesthetics'),Z_analysis=cfg.get('Z_analysis','unspecified'),
                chr38_plotted=False,threshold_label=False,threshold_line_count=1,no_downsampling=True,
                individual_PDFs=len(qcs),individual_PNGs=len(qcs),script_sha256=sha(a.script))
    (out/'audit/FINAL.PASS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
