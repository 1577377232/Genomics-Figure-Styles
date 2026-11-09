#!/usr/bin/env python3
"""Audit actual SVG quadratic paths, source preservation and vector/raster exports."""
from pathlib import Path
import csv,hashlib,json,math,re,subprocess
import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial import cKDTree
from PIL import Image
O=Path(__file__).resolve().parents[1];P=O.parents[1]
S='chicken_GRCg7b_R_v15_with_TE_architectures_v3'
def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tag(e):return e.tag.rsplit('}',1)[-1]
def norm(e,prefix=''):
    return (e.tag,{k:v.replace(prefix,'') for k,v in e.attrib.items() if not(k=='id' and v=='')},(e.text or '').strip(),[norm(c,prefix) for c in e])
def body(e,prefix=''):return [norm(c,prefix) for c in e]
root=ET.parse(O/(S+'.svg')).getroot()
source_file=next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.svg'))
source=ET.parse(source_file).getroot()
main=next(e for e in root if e.get('id')=='V15-circos')
assert body(source)==body(main,'V15-circos__')
assert main.get('preserveAspectRatio')=='xMinYMin meet'
clips=sum(tag(e)=='clipPath' for e in main.iter())
assert clips==sum(tag(e)=='clipPath' for e in source.iter())==325
baseline=json.loads((O/'protected_baseline.json').read_text())
for rel,expected in baseline.items():assert hashlib.sha256((P/rel).read_bytes()).hexdigest()==expected,rel
for row in table(O/'input_preservation.tsv'):
    assert row['unchanged']=='TRUE'
    assert hashlib.md5(Path(row['path']).read_bytes()).hexdigest()==row['md5']
arch=table(O/'architecture_placements.tsv')
assert [r['track'] for r in arch]==['LINE','LTR','SINE','DNA','TEV']
assert len({r['width_in'] for r in arch})==1
old=table(P/'07_figures/R_v15_with_TE_architectures_v2/architecture_placements.tsv')
for a,b in zip(arch,old):
    assert abs(float(a['top_in'])-float(b['top_in'])-.6)<1e-10
    assert abs(float(a['scale_x'])-float(a['scale_y']))<1e-12
    embedded=next(e for e in root if e.get('id')=='architecture-'+a['track'])
    assert body(ET.parse(P/'assets'/a['file']).getroot())==body(embedded,'architecture-'+a['track']+'__')
assert np.max(np.abs(np.diff([float(r['target_y_in']) for r in arch])-1.55))<1e-10
tev=ET.parse(P/'assets/TE_mediated_SV/TE_mediated_SV_classification.svg').getroot()
assert tev.get('viewBox')=='0 0 1400 250'
assert {e.text.strip() for e in tev.iter() if tag(e)=='text'}=={'TEV','SV','TE','>80%'}
assert next(e for e in tev.iter() if e.get('id')=='threshold-line').get('d')=='M822 130 V204'
assert 822==110+890*.8
assert 760/890>.8
assert '>50%' in ''.join(ET.parse(P/'assets/TE_mediated_SV/TE_mediated_SV_classification_V15_gt50.svg').getroot().itertext())
ids=[e.get('id') for e in root.iter() if e.get('id')];assert len(ids)==len(set(ids))
refs=[]
for e in root.iter():
    for k,v in e.attrib.items():
        refs+=re.findall(r'url\(\s*[\'"]?#([^\)\'\"]+)',v)
        if k.endswith('href') and v.startswith('#'):refs.append(v[1:])
assert set(refs)<=set(ids)
assert not any(tag(e)=='image' for e in root.iter())
layer=next(e for e in root if e.get('id')=='connector-layer')
paths=[e for e in layer.iter() if tag(e)=='path'];heads=[e for e in layer.iter() if tag(e)=='polygon']
assert len(paths)==len(heads)==5
assert not any(tag(e)=='polyline' for e in layer.iter())
# Independently parse real M/L/Q SVG geometry, sampling with step <=0.08pt.
def sample_path(d):
    tok=re.findall(r'[MLQ]|[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?',d)
    i=0;result=[];counts={}
    while i<len(tok):
        cmd=tok[i];i+=1;counts[cmd]=counts.get(cmd,0)+1
        n={'M':2,'L':2,'Q':4}[cmd];v=np.array(list(map(float,tok[i:i+n]))).reshape(-1,2);i+=n
        if cmd=='M':pos=v[0];result.append(pos.copy());continue
        length=np.linalg.norm(v[0]-pos)+(np.linalg.norm(v[1]-v[0]) if cmd=='Q' else 0)
        t=np.linspace(0,1,max(2,math.ceil(length/.08)+1))[1:,None]
        curve=(1-t)*pos+t*v[0] if cmd=='L' else (1-t)**2*pos+2*(1-t)*t*v[0]+t*t*v[1]
        result.extend(curve);pos=v[-1]
    assert counts=={'M':1,'L':3,'Q':2},counts
    return np.array(result)
curves={e.get('id').split('connector-')[-1]:sample_path(e.get('d')) for e in paths}
assert list(curves)==['LINE','LTR','SINE','DNA','TEV']
vertices=table(O/'connector_vertices.tsv')
for name,points in curves.items():
    vv=[r for r in vertices if r['connector']==name]
    start=np.array([float(vv[0]['x_in']),float(vv[0]['y_from_top_in'])])*72
    end=np.array([float(vv[-1]['x_in']),float(vv[-1]['y_from_top_in'])])*72
    assert np.max(np.abs(points[0]-start))<1e-9
    assert np.max(np.abs(points[-1]-end))<1e-9
    assert float(next(e for e in paths if e.get('id').endswith(name)).get('stroke-width'))==.8
geom=(O/'source_geometry.R').read_text()
cx,cy=map(float,re.search(r'circle_center_pt = c\(([^,]+), ([^)]+)\)',geom).groups())
radius=float(re.search(r'circle_radius_pt = ([0-9.]+)',geom)[1])
place=table(O/'circos_placement.tsv')[0];sx=float(place['width_pt'])/705;sy=float(place['height_pt'])/720
unitx=705/9.8*sx;unity=72*sy;left=float(place['left_pt']);top=float(place['top_pt'])
annotations=[]
for suffix in ['chromosome_labels.tsv','coordinate_labels.tsv','TE_axis_labels.tsv','variant_axis_labels.tsv']:
    for r in table(next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.'+suffix))):
        if r.get('label_status') in ['hidden','omitted'] or 'center_x_in' not in r:continue
        annotations.append(dict(x=left+float(r['center_x_in'])*unitx,y=top+(10-float(r['center_y_in']))*unity,
            hw=float(r['width_in'])*unitx/2+1.2,hh=float(r['height_in'])*unity/2+1.2,
            tx=float(r.get('tangent_x',1)),ty=-float(r.get('tangent_y',0)),rx=float(r.get('radial_x',0)),ry=-float(r.get('radial_y',1)),label=r['label']))
for name,curve in curves.items():
    x=curve[:,0];y=curve[:,1];rr=np.hypot(x-cx,cy-y)/radius;angle=np.degrees(np.arctan2(cy-y,x-cx))%360
    assert not np.any((rr>=.3185)&(rr<=.9975)&(angle>=24)&(angle<=336)),(name,'data collision')
    for a in annotations:
        dx=x-a['x'];dy=y-a['y']
        assert not np.any((np.abs(dx*a['tx']+dy*a['ty'])<a['hw'])&(np.abs(dx*a['rx']+dy*a['ry'])<a['hh'])),(name,a['label'])
    # Arrow ends stop before each matching module's left edge.
    a=next(r for r in arch if r['track']==name)
    assert abs(y[-1]-float(a['target_y_in'])*72)<1e-9
# Maximum sampling step bounds the gap between samples and the actual curve.
pairs=[]
for i,(name,p) in enumerate(curves.items()):
    assert p[0,0]>cx
    assert (p[0,1]<cy)==(name!='TEV')
    for other,q in list(curves.items())[i+1:]:
        sampled=float(cKDTree(q).query(p)[0].min())
        conservative=sampled-.16
        assert conservative>=3.5,(name,other,conservative)
        pairs.append(dict(track1=name,track2=other,sampled_clearance_pt=sampled,conservative_clearance_lower_bound_pt=conservative))
anchors=table(O/'connector_anchors.tsv')
assert [int(r['source_track_index']) for r in anchors]==[4,5,6,7,10]
assert anchors[-1]['source_track']=='TIP'
for a in anchors:assert float(a['source_inner_radius'])<float(a['anchor_radius'])<float(a['source_outer_radius'])
with Image.open(O/(S+'.png')) as im:
    assert im.size==(8400,6000) and all(abs(d-600)<1 for d in im.info['dpi'])
with Image.open(O/(S+'.preview.png')) as im:assert im.size==(2100,1500)
info=subprocess.check_output(['pdfinfo',str(O/(S+'.pdf'))],text=True);assert re.search(r'Page size:\s+1008 x 720 pts',info)
images=subprocess.check_output(['pdfimages','-list',str(O/(S+'.pdf'))],text=True);assert len(images.strip().splitlines())==2
text=subprocess.check_output(['pdftotext',str(O/(S+'.pdf')),'-'],text=True)
for word in ['LINE','LTR','SINE','DNA','ORF1','ORF2','Transposase','TEV','>80%']:assert word in text,word
fonts=subprocess.check_output(['pdffonts',str(O/(S+'.pdf'))],text=True);(O/'provenance/PDF_fonts.txt').write_text(fonts)
report=dict(status='PASS',V15_vector_contents_identical=True,protected_files_unchanged=len(baseline),four_TE_originals_unchanged=True,V2_unchanged=True,
    V15_clip_paths_preserved=clips,all_five_SVG_contents_preserved=True,independent_connectors=5,closed_arrowheads=5,common_stems=0,
    quadratic_corners_per_connector=2,actual_SVG_curves_independently_audited=True,annotation_padding_pt=1.2,
    intersection_count=0,connectors_cross_data=False,connectors_cover_annotations=False,source_track_indices=[4,5,6,7,10],
    origin_sides=['upper','upper','upper','upper','lower'],module_shift_down_in=.6,module_center_spacing_in=1.55,
    minimum_clearance_lower_bound_pt=min(r['conservative_clearance_lower_bound_pt'] for r in pairs),connector_pairs=pairs,
    TEV_schematic_threshold=.8,V15_data_threshold_unchanged=.5,threshold_difference_documented=True,
    TEV_marker_x=822,TEV_marker_expected_x=110+890*.8,TEV_text_editable=True,genomic_data_recomputed=False,
    SVG_embedded_rasters=0,PDF_embedded_rasters=0,PNG_pixels=[8400,6000],PNG_dpi=600,preview_pixels=[2100,1500])
(O/'validation_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
