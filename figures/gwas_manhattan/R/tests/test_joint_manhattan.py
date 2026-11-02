#!/usr/bin/env python3
"""End-to-end TEST-ONLY fixtures. Never previews or biological results."""
import csv,hashlib,json,os,subprocess,tempfile,unittest
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
SCRIPT=HERE/'plot_joint_manhattan.R'
STYLE=HERE/'user_style.two_colour.json'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(p,rows,fields):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def save(p,x):Path(p).write_text(json.dumps(x)+'\n')

def fixture(out,zero=False,bad_order=False):
    (out/'config').mkdir();(out/'TEST_ONLY_inputs').mkdir()
    # Tiny, artificial contigs and association values; no actual project data.
    fai=out/'TEST_ONLY.fai';fai.write_text('1\t1000\t0\t50\t51\n2\t800\t0\t50\t51\n38\t500\t0\t50\t51\nZ\t400\t0\t50\t51\n')
    style=json.loads(STYLE.read_text());save(out/'config/user_style.json',style)
    cfg=dict(expected_count=2,expected_singles=1,expected_modules=1,display_style_file='config/user_style.json',fai=str(fai),fai_sha256=sha(fai),
        chromosomes=['2','1','Z'] if bad_order else ['1','2','Z'],chromosome_aliases={},excluded_plot_chromosomes=['38'],
        alpha=.05,reference_build='TEST_ONLY_NOT_A_REAL_GENOME',style_commit='TEST_ONLY',Z_analysis='joint')
    save(out/'config/analysis.json',cfg)
    traits=[];sources=[]
    for tid,kind in [('TEST_ONLY_WALD','single'),('TEST_ONLY_MODULE','module')]:
        traits.append(dict(trait_id=tid,display_name_zh='TEST ONLY',display_name_en='TEST ONLY '+kind,category='TEST_ONLY',wave='TEST_ONLY',
            preview_order=0,BW42_source_flag=False,existing_QC_HOLD_TE='',kind=kind,module_dimensions=2 if kind=='module' else 1,module_traits=''))
        for marker in ['SNP','TE']:
            for panel in ['A','Z_joint']:
                pc='p_wald' if kind=='single' else 'p_score'
                values=[('1',10,.1),('1',20,1e-7),('2',30,.01),('38',5,1e-8),('1',50,'NA')] if panel=='A' else [('Z',200,5e-4)]
                if zero and tid=='TEST_ONLY_WALD' and marker=='SNP' and panel=='A':values[0]=('1',10,0.)
                rows=[dict(chr=c,rs='TEST_ONLY.INS.'+str(i),ps=bp,**{pc:p}) for i,(c,bp,p) in enumerate(values)]
                if marker=='TE' and kind=='single':
                    for r in rows:r.update(beta=.1,se=.05)
                file=out/'TEST_ONLY_inputs'/f'{tid}.{marker}.{panel}.tsv';table(file,rows,list(rows[0]))
                good=[float(r[pc]) for r in rows if r[pc]!='NA']
                sources.append(dict(trait_id=tid,marker=marker,panel=panel,N=500,p_column=pc,file=str(file),receipt='',
                    expected_reported=len(rows),expected_valid=len(good),expected_minP=min(good),expected_P_zero=int(0. in good),expected_sha256=sha(file),method='TEST_ONLY',model='TEST_ONLY'))
    table(out/'config/traits.tsv',traits,list(traits[0]));table(out/'config/sources.tsv',sources,list(sources[0]))

def render(out,selection='all',threshold=None):
    env=os.environ.copy()
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
    cmd=['Rscript',str(SCRIPT),'--outdir',str(out),'--jobs','1','--trait-id',selection,'--no-overview']
    if threshold:cmd+=['--threshold-p',str(threshold)]
    return subprocess.run(cmd,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

class PlotTests(unittest.TestCase):
    def test_real_render_counts_coordinates_p_values_and_resume(self):
        with tempfile.TemporaryDirectory(prefix='GWAS_TEST_ONLY_') as t:
            out=Path(t);fixture(out);r=render(out);self.assertEqual(r.returncode,0,r.stdout)
            for tid,pc in [('TEST_ONLY_WALD','p_wald'),('TEST_ONLY_MODULE','p_score')]:
                q=json.loads((out/'audit'/f'{tid}__PASS.json').read_text())
                self.assertEqual(q['p_column'],pc);self.assertEqual(q['Z_analysis'],'joint')
                for marker in ['SNP','TE']:
                    self.assertEqual(q['M_'+marker],5);self.assertEqual(q[marker+'_plotted'],4)
                    self.assertEqual(q[marker+'_excluded_chr38'],1);self.assertEqual(q[marker+'_Z_tests'],1)
                    self.assertAlmostEqual(q['minP_'+marker],1e-8);self.assertAlmostEqual(q['minP_plotted_'+marker],1e-7)
                    self.assertEqual(q[marker+'_point_size_mm'],1.1);self.assertEqual(q[marker+'_display_shape'],16)
                self.assertEqual(q['threshold_line_count'],1);self.assertEqual(q['display_threshold_p'],1e-5)
                self.assertEqual(q['Bonferroni_SNP'],.01);self.assertEqual(q['Bonferroni_TE'],.01)
                self.assertFalse(q['threshold_label']);self.assertFalse(q['downsampling']);self.assertFalse(q['chr38_plotted'])
                self.assertFalse(q['TE_type_distinction']);self.assertFalse(q['Z_sex_distinction']);self.assertGreater(q['ymax'],7)
            with (out/'audit/chromosome_layout.tsv').open() as f:layout=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual([r['chr'] for r in layout],['1','2','Z']);self.assertEqual(float(layout[1]['offset_bp']),12001000)
            again=render(out);self.assertEqual(again.returncode,0,again.stdout);self.assertIn('REUSE PASS',again.stdout)
            changed=render(out,'TEST_ONLY_WALD',1e-4);self.assertEqual(changed.returncode,0,changed.stdout);self.assertNotIn('REUSE PASS',changed.stdout)
            self.assertEqual(json.loads((out/'audit/TEST_ONLY_WALD__PASS.json').read_text())['display_threshold_p'],1e-4)

    def test_p_zero_rejected_without_arbitrary_floor(self):
        with tempfile.TemporaryDirectory(prefix='GWAS_TEST_ONLY_') as t:
            out=Path(t);fixture(out,zero=True);r=render(out,'TEST_ONLY_WALD')
            self.assertNotEqual(r.returncode,0);self.assertIn('P=0 underflow',r.stdout)
            self.assertFalse((out/'audit/TEST_ONLY_WALD__PASS.json').exists())

    def test_nonbiological_chromosome_order_rejected(self):
        with tempfile.TemporaryDirectory(prefix='GWAS_TEST_ONLY_') as t:
            out=Path(t);fixture(out,bad_order=True);r=render(out)
            self.assertNotEqual(r.returncode,0);self.assertIn('biological_order',r.stdout)

if __name__=='__main__':unittest.main(verbosity=2)
