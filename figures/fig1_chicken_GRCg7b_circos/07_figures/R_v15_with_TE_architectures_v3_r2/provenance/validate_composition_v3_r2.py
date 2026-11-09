#!/usr/bin/env python3
"""Independent audit of V3 R2 actual cubic paths, true track origins and vectors."""
from pathlib import Path
import csv,hashlib,json,math,re,subprocess
import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial import cKDTree
from PIL import Image
O=Path(__file__).resolve().parents[1];P=O.parents[1]
S='chicken_GRCg7b_R_v15_with_TE_architectures_v3_r2'
def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tag(e):return e.tag.rsplit('}',1)[-1]
def normalized(e,prefix=''):
    return (e.tag,{k:v.replace(prefix,'') for k,v in e.attrib.items() if not(k=='id' and v=='')},(e.text or '').strip(),[normalized(c,prefix) for c in e])
def body(e,prefix=''):return [normalized(c,prefix) for c in e]
root=ET.parse(O/(S+'.svg')).getroot()
main=next(e for e in root if e.get('id')=='V15-circos')
original=next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.svg'))
source=ET.parse(original).getroot()
assert body(source)==body(main,'V15-circos__')
clips=sum(tag(e)=='clipPath' for e in main.iter());assert clips==325
for rel,sha in json.loads((O/'protected_baseline.json').read_text()).items():assert hashlib.sha256((P/rel).read_bytes()).hexdigest()==sha,rel
for r in table(O/'input_preservation.tsv'):
    assert r['unchanged']=='TRUE' and hashlib.md5(Path(r['path']).read_bytes()).hexdigest()==r['md5']
prior=P/'07_figures/R_v15_with_TE_architectures_v3'
placement=table(O/'circos_placement.tsv')[0]
assert placement==table(prior/'circos_placement.tsv')[0],'Left Circos scale or position changed'
arch=table(O/'architecture_placements.tsv');old=table(prior/'architecture_placements.tsv')
assert [r['track'] for r in arch]==['LINE','LTR','SINE','DNA','TEV']
for a,b in zip(arch,old):
    assert abs(float(a['width_in'])/float(b['width_in'])-1.25)<1e-12
    assert abs(float(a['scale_x'])-float(a['scale_y']))<1e-12
    asset=ET.parse(P/'assets'/a['file']).getroot()
    embedded=next(e for e in root if e.get('id')=='architecture-'+a['track'])
    assert body(asset)==body(embedded,'architecture-'+a['track']+'__')
    texts=[e for e in asset.iter() if tag(e)=='text']
    heading=next(e for e in texts if e.text==a['track'])
    assert heading.get('font-size')=='58'
assert np.max(np.abs(np.diff([float(r['target_y_in']) for r in arch])-1.45))<1e-12
tev=ET.parse(P/'assets/TE_mediated_SV/TE_mediated_SV_classification_v3_r2.svg').getroot()
assert tev.get('viewBox')=='0 0 1400 250'
assert {e.text.strip() for e in tev.iter() if tag(e)=='text'}=={'SV','TE','TEV','>80%'}
assert sum(tag(e)=='text' and e.text=='TEV' for e in tev.iter())==1
assert not any(e.get('id') in ['classification-arrow','tev-result','tev-result-label'] for e in tev.iter())
assert not any(tag(e)=='polygon' for e in tev.iter())
assert next(e for e in tev.iter() if e.get('id')=='threshold-line').get('d')=='M1054 130 V204'
assert 1054==110+1180*.8 and 1008/1180>.8
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
# Independently parse exported cubic geometry. Bound step using max derivative.
def sample_path(d):
    tok=re.findall(r'[MLC]|[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?',d);i=0;result=[];counts={}
    while i<len(tok):
        cmd=tok[i];i+=1;counts[cmd]=counts.get(cmd,0)+1;n={'M':2,'L':2,'C':6}[cmd]
        v=np.array(list(map(float,tok[i:i+n]))).reshape(-1,2);i+=n
        if cmd=='M':pos=v[0];result.append(pos.copy());continue
        bound=np.linalg.norm(v[0]-pos) if cmd=='L' else 3*max(np.linalg.norm(v[0]-pos),np.linalg.norm(v[1]-v[0]),np.linalg.norm(v[2]-v[1]))
        t=np.linspace(0,1,max(2,math.ceil(bound/.08)+1))[1:,None]
        curve=(1-t)*pos+t*v[0] if cmd=='L' else (1-t)**3*pos+3*(1-t)**2*t*v[0]+3*(1-t)*t*t*v[1]+t**3*v[2]
        result.extend(curve);pos=v[-1]
    assert counts=={'M':1,'C':2,'L':2},counts
    return np.array(result)
curves={e.get('id').split('connector-')[-1]:sample_path(e.get('d')) for e in paths}
assert list(curves)==['LINE','LTR','SINE','DNA','TEV']
geom=(O/'source_geometry.R').read_text()
cx,cy=map(float,re.search(r'circle_center_pt = c\(([^,]+), ([^)]+)\)',geom).groups())
radius=float(re.search(r'circle_radius_pt = ([0-9.]+)',geom)[1])
sx=float(placement['scale_x']);sy=float(placement['scale_y']);left=float(placement['left_pt']);top=float(placement['top_pt'])
unitx=705/9.8*sx;unity=72*sy
annotations=[]
for suffix in ['chromosome_labels.tsv','coordinate_labels.tsv','TE_axis_labels.tsv','variant_axis_labels.tsv']:
    for r in table(next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.'+suffix))):
        if r.get('label_status') in ['hidden','omitted'] or 'center_x_in' not in r:continue
        annotations.append(dict(x=left+float(r['center_x_in'])*unitx,y=top+(10-float(r['center_y_in']))*unity,hw=float(r['width_in'])*unitx/2+1.2,hh=float(r['height_in'])*unity/2+1.2,tx=float(r.get('tangent_x',1)),ty=-float(r.get('tangent_y',0)),rx=float(r.get('radial_x',0)),ry=-float(r.get('radial_y',1)),label=r['label']))
anchors=table(O/'connector_anchors.tsv');panels=table(next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.sector_panels.tsv')))
assert [int(r['source_track_index']) for r in anchors]==[4,5,6,7,10]
for a in anchors:
    assert a['source_chromosome']=='1' and a['source_position_bp']=='0'
    assert a['anchor_side']=='upper' and a['exact_sector_endpoint']=='TRUE'
    panel=next(r for r in panels if r['chromosome']=='1' and r['track']==a['source_track'])
    theta=float(panel['start_degrees']);rad=float(a['anchor_radius'])
    assert float(a['anchor_degrees'])==theta==24
    assert float(a['source_inner_radius'])<rad<float(a['source_outer_radius'])
    expected=np.array([cx+radius*rad*math.cos(math.radians(theta)),cy-radius*rad*math.sin(math.radians(theta))])
    observed=curves[a['track']][0]
    assert np.max(np.abs(expected-observed))<1e-9
    assert observed[0]>cx and observed[1]<cy
    curve=curves[a['track']];x=curve[:,0];y=curve[:,1]
    rr=np.hypot(x-cx,cy-y)/radius;angle=np.degrees(np.arctan2(cy-y,x-cx))%360
    # Boundary contact at the true source is intentional; interior crossings fail.
    assert not np.any((rr>=.3185)&(rr<=.9975)&(angle>24+1e-6)&(angle<336-1e-6)),a['track']
    for b in annotations:
        dx=x-b['x'];dy=y-b['y']
        assert not np.any((np.abs(dx*b['tx']+dy*b['ty'])<b['hw'])&(np.abs(dx*b['rx']+dy*b['ry'])<b['hh'])),(a['track'],b['label'])
    dest=next(r for r in arch if r['track']==a['track'])
    assert np.max(np.abs(curve[-1]-np.array([float(dest['target_x_in']),float(dest['target_y_in'])])*72))<1e-9
pairs=[]
for i,(name,p) in enumerate(curves.items()):
    for other,q in list(curves.items())[i+1:]:
        distance=float(cKDTree(q).query(p)[0].min());lower=distance-.16
        assert lower>=3.5,(name,other,lower)
        pairs.append(dict(track1=name,track2=other,sampled_clearance_pt=distance,conservative_clearance_lower_bound_pt=lower))
# Verify the source edge against an actual coral polygon in the original V15 SVG.
# A witness 0.05 degrees INSIDE Chr1 must lie in that polygon; the anchor must lie
# on its actual boundary. This guards against relying on a visually guessed radius.
a=anchors[-1];rad=float(a['anchor_radius']);ang=24.05
witness=np.array([(cx+radius*rad*math.cos(math.radians(ang))-left)/sx,(cy-radius*rad*math.sin(math.radians(ang))-top)/sy])
origin=(curves['TEV'][0]-np.array([left,top]))/np.array([sx,sy])
def inside(poly,point):
    x,y=point;state=False
    for (ax,ay),(bx,by) in zip(poly,poly[1:]+poly[:1]):
        if (ay>y)!=(by>y) and x<(bx-ax)*(y-ay)/(by-ay)+ax:state=not state
    return state
def point_segment(p,a,b):
    p=np.array(p);a=np.array(a);b=np.array(b);v=b-a;vv=float(v@v)
    t=max(0,min(1,float((p-a)@v)/vv)) if vv else 0
    return float(np.linalg.norm(p-(a+t*v)))
matches=[]
for i,e in enumerate(source.iter()):
    style=e.get('style','');d=e.get('d','')
    if 'fill:rgb(' not in style or 'stroke:none' not in style or not d or re.search(r'[CQAHV]',d):continue
    rgb=[float(q.rstrip('%')) for q in style.split('fill:rgb(')[1].split(')')[0].split(',')]
    if rgb[0]<=rgb[1]+5 or rgb[0]<=rgb[2]+5:continue
    for part in re.findall(r'M([^M]+)',d):
        nums=list(map(float,re.findall(r'-?\d+(?:\.\d+)?',part)));poly=list(zip(nums[::2],nums[1::2]))
        if len(poly)>2 and inside(poly,witness):
            gap=min(point_segment(origin,a,b) for a,b in zip(poly,poly[1:]+poly[:1]))
            matches.append(dict(source_element_index=i,coral_rgb_percent=rgb,anchor_boundary_distance_source_pt=gap,witness_source_point=witness.tolist()))
assert len(matches)==1,matches
assert matches[0]['anchor_boundary_distance_source_pt']<.005,matches
# DEBUG markers must have exactly the same centers as actual connector origins.
debug=ET.parse(O/(S+'.DEBUG.svg')).getroot();markers=next(e for e in debug if e.get('id')=='DEBUG-anchor-markers')
circles=[e for e in markers if tag(e)=='circle'];assert len(circles)==5
for a,node in zip(anchors,circles):
    assert abs(float(node.get('cx'))-float(a['x_in'])*72)<1e-9
    assert abs(float(node.get('cy'))-float(a['y_from_top_in'])*72)<1e-9
assert not any(e.get('id')=='DEBUG-anchor-markers' for e in root.iter())
for suffix,size,dpi in [('.png',(10200,6000),600),('.preview.png',(2550,1500),150),('.DEBUG.preview.png',(2550,1500),150),('.DEBUG.anchor_zoom.png',(2130,1380),300)]:
    with Image.open(O/(S+suffix)) as im:
        assert im.size==size,(suffix,im.size)
        assert all(abs(d-dpi)<1 for d in im.info['dpi'])
info=subprocess.check_output(['pdfinfo',str(O/(S+'.pdf'))],text=True);assert re.search(r'Page size:\s+1224 x 720 pts',info)
images=subprocess.check_output(['pdfimages','-list',str(O/(S+'.pdf'))],text=True);assert len(images.strip().splitlines())==2
text=subprocess.check_output(['pdftotext',str(O/(S+'.pdf')),'-'],text=True)
for word in ['LINE','LTR','SINE','DNA','TEV','ORF1','ORF2','Transposase','>80%']:assert word in text,word
fonts=subprocess.check_output(['pdffonts',str(O/(S+'.pdf'))],text=True);(O/'provenance/PDF_fonts.txt').write_text(fonts)
report=dict(status='PASS',revision='V3 revision 2',V15_vector_contents_identical=True,V15_placement_and_scale_identical=True,V15_clip_paths_preserved=325,
    protected_old_files_unchanged=len(json.loads((O/'protected_baseline.json').read_text())),V2_and_previous_V3_unchanged=True,four_original_TE_SVGs_unchanged=True,
    module_uniform_scale_factor=1.25,module_width_in=4.9375,module_title_font_source_px=58,module_title_final_pt=58*4.9375*72/1400,module_center_spacing_in=1.45,
    TEV_internal_arrow_removed=True,TEV_result_box_removed=True,TEV_visible_texts=['TEV','SV','TE','>80%'],TEV_threshold_marker_x=1054,
    independent_connectors=5,closed_arrowheads=5,shared_stems=0,connector_intersections=0,actual_cubic_paths_audited=True,connectors_cross_data=False,connectors_cover_annotations=False,
    all_origins_exact_upper_Chr1_endpoints=True,source_track_indices=[4,5,6,7,10],TEV_source_track='TIP / TEV density',TEV_anchor_radius=rad,TEV_anchor_degrees=24,TEV_source_chromosome='1',TEV_source_position_bp=0,
    TEV_origin_verified_against_original_coral_polygon=matches[0],debug_markers_verified=True,
    minimum_curve_clearance_lower_bound_pt=min(r['conservative_clearance_lower_bound_pt'] for r in pairs),connector_pairs=pairs,
    schematic_TE_threshold=.8,V15_data_TE_threshold_unchanged=.5,genomic_data_recomputed=False,
    SVG_embedded_rasters=0,PDF_embedded_rasters=0,PDF_inches=[17,10],PNG_pixels=[10200,6000],PNG_dpi=600,preview_pixels=[2550,1500],DEBUG_zoom_pixels=[2130,1380])
(O/'validation_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
