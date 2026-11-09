#!/usr/bin/env python3
"""Audit V3_R3 actual straight geometry, normal gaps, source assets and exports."""
from pathlib import Path
import csv,hashlib,json,math,re,subprocess
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
O=Path(__file__).resolve().parents[1];P=O.parents[1];S='chicken_GRCg7b_R_v15_with_TE_architectures_v3_r3'
def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tag(e):return e.tag.rsplit('}',1)[-1]
def norm(e,prefix=''):
    return (e.tag,{k:v.replace(prefix,'') for k,v in e.attrib.items() if not(k=='id' and v=='')},(e.text or '').strip(),[norm(c,prefix) for c in e])
def body(e,prefix=''):return [norm(c,prefix) for c in e]
root=ET.parse(O/(S+'.svg')).getroot();main=next(e for e in root if e.get('id')=='V15-circos')
source_file=next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.svg'));source=ET.parse(source_file).getroot()
assert body(source)==body(main,'V15-circos__')
assert sum(tag(e)=='clipPath' for e in main.iter())==325
baseline=json.loads((O/'protected_baseline.json').read_text())
for rel,sha in baseline.items():assert hashlib.sha256((P/rel).read_bytes()).hexdigest()==sha,rel
for r in table(O/'input_preservation.tsv'):assert r['unchanged']=='TRUE' and hashlib.md5(Path(r['path']).read_bytes()).hexdigest()==r['md5']
prior=P/'07_figures/R_v15_with_TE_architectures_v3_r2';placement=table(O/'circos_placement.tsv')[0]
assert placement==table(prior/'circos_placement.tsv')[0]
geo={k:float(v) for k,v in table(O/'layout_geometry.tsv')[0].items()}
cx=geo['circle_center_x_in']*72;cy=geo['circle_center_y_from_top_in']*72;R=geo['circle_radius_in']*72
sx=float(placement['scale_x']);sy=float(placement['scale_y']);left=float(placement['left_pt']);top=float(placement['top_pt'])
unitx=705/9.8*sx;unity=72*sy
annotations=[]
for suffix in ['chromosome_labels.tsv','coordinate_labels.tsv','TE_axis_labels.tsv','variant_axis_labels.tsv']:
    for r in table(next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.'+suffix))):
        if r.get('label_status') in ['hidden','omitted'] or 'center_x_in' not in r:continue
        annotations.append(dict(x=left+float(r['center_x_in'])*unitx,y=top+(10-float(r['center_y_in']))*unity,
            hw=float(r['width_in'])*unitx/2,hh=float(r['height_in'])*unity/2,
            tx=float(r.get('tangent_x',1)),ty=-float(r.get('tangent_y',0)),rx=float(r.get('radial_x',0)),ry=-float(r.get('radial_y',1)),label=r['label']))
map_file=next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.chromosome_display_mapping.tsv'))
track_file=next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.track_geometry.tsv'))
mapping=table(map_file);tracks=table(track_file);outer=max(float(r['outer_radius']) for r in tracks)
cosines=[];sines=[]
for row in mapping:
    a=float(row['start_degrees']);b=float(row['end_degrees']);angles=[a,b]+[i for i in range(-720,811,90) if a<=i<=b]
    cosines+=list(map(lambda q:math.cos(math.radians(q)),angles));sines+=list(map(lambda q:math.sin(math.radians(q)),angles))
sector_right=cx+R*outer*max(cosines)
label_right=max(a['x']+abs(a['tx'])*a['hw']+abs(a['rx'])*a['hh'] for a in annotations)
bound=(max(sector_right,label_right)+.5)/72
assert abs(bound-geo['circos_right_bound_in'])<1e-10
assert abs(geo['module_left_in']-bound-.70)<1e-10
assert abs(geo['page_width_in']-geo['module_left_in']-7.40625-.20)<1e-10
circos_top=min(cy-R*outer*max(sines),min(a['y']-abs(a['ty'])*a['hw']-abs(a['ry'])*a['hh'] for a in annotations))/72
circos_bottom=max(cy-R*outer*min(sines),max(a['y']+abs(a['ty'])*a['hw']+abs(a['ry'])*a['hh'] for a in annotations))/72
assert abs(circos_bottom-circos_top-geo['circos_height_in'])<1e-10
arch=table(O/'architecture_placements.tsv');old_arch=table(prior/'architecture_placements.tsv')
expected=['TE_architectures/LINE_CR1_composition_1_5x.svg','TE_architectures/LTR_ERV_composition_1_5x.svg','TE_architectures/SINE_generic_composition_1_5x.svg','TE_architectures/DNA_TIR_composition_1_5x.svg','TE_mediated_SV/TEV_classification_composition_1_5x.svg']
assert [r['file'] for r in arch]==expected
heading_ratios=[]
for a,b in zip(arch,old_arch):
    assert float(a['width_in'])==7.40625
    assert abs(float(a['width_in'])/float(b['width_in'])-1.5)<1e-12
    assert abs(float(a['scale_x'])-float(a['scale_y']))<1e-12
    new=ET.parse(P/'assets'/a['file']).getroot();old=ET.parse(P/'assets'/b['file']).getroot()
    heading=next(e for e in new.iter() if tag(e)=='text' and e.text==a['track']);assert heading.get('font-size')=='39'
    old_heading=next(e for e in old.iter() if tag(e)=='text' and e.text==a['track']);assert old_heading.get('font-size')=='58'
    heading_ratios.append(39/58*1.5)
    old_heading.set('font-size','39')
    assert body(old)==body(new),'Supplied asset changes exceed heading compensation: '+a['track']
    embedded=next(e for e in root if e.get('id')=='architecture-'+a['track']);assert body(new)==body(embedded,'architecture-'+a['track']+'__')
heights=[float(a['height_pt'])/72 for a in arch];tops=[float(a['top_in']) for a in arch]
for i in range(4):assert abs(tops[i+1]-tops[i]-heights[i]-.08)<1e-10
column_height=tops[-1]+heights[-1]-tops[0]
assert abs(column_height-geo['column_height_in'])<1e-10
assert column_height<(circos_bottom-circos_top)*.92
assert abs((tops[0]+tops[-1]+heights[-1])/2-(circos_top+circos_bottom)/2)<1e-10
tev=ET.parse(P/'assets'/expected[-1]).getroot()
assert {e.text.strip() for e in tev.iter() if tag(e)=='text'}=={'SV','TE','TEV','>80%'}
assert not any(e.get('id') in ['classification-arrow','tev-result','tev-result-label'] for e in tev.iter())
assert next(e for e in tev.iter() if e.get('id')=='threshold-line').get('d')=='M1054 130 V204'
ids=[e.get('id') for e in root.iter() if e.get('id')];assert len(ids)==len(set(ids))
refs=[]
for e in root.iter():
    for k,v in e.attrib.items():
        refs+=re.findall(r'url\(\s*[\'"]?#([^\)\'\"]+)',v)
        if k.endswith('href') and v.startswith('#'):refs.append(v[1:])
assert set(refs)<=set(ids)
assert not any(tag(e)=='image' for e in root.iter())
layer=next(e for e in root if e.get('id')=='connector-layer')
paths=[e for e in layer.iter() if tag(e)=='path'];heads=[e for e in layer.iter() if tag(e)=='polygon'];assert len(paths)==len(heads)==5
lines={}
for node in paths:
    d=node.get('d');assert re.findall(r'[A-Za-z]',d)==['M','L','L','L']
    n=list(map(float,re.findall(r'-?\d+(?:\.\d+)?',d)));p=np.array(n).reshape(-1,2);assert p.shape==(4,2)
    assert float(node.get('stroke-width'))==.85 and node.get('stroke-linejoin')=='round'
    assert abs(p[-1,1]-p[-2,1])<1e-9
    lines[node.get('id').split('connector-')[-1]]=p
assert list(lines)==['LINE','LTR','SINE','DNA','TEV']
vertices=table(O/'connector_vertices.tsv');anchors=table(O/'connector_anchors.tsv')
assert [int(a['source_track_index']) for a in anchors]==[4,5,6,7,10]
all_points=[]
# Sample segments for label/panel/canvas tests; pair clearance uses exact segments.
def samples(p):
    return np.vstack([(1-t)*a+t*b for a,b in zip(p,p[1:]) for t in np.linspace(0,1,max(2,math.ceil(np.linalg.norm(b-a)/.06)+1))])
def no_annotations(arr):
    x=arr[:,0];y=arr[:,1]
    for a in annotations:
        dx=x-a['x'];dy=y-a['y']
        assert not np.any((np.abs(dx*a['tx']+dy*a['ty'])<a['hw']+1.2)&(np.abs(dx*a['rx']+dy*a['ry'])<a['hh']+1.2)),a['label']
normal=np.array([math.sin(math.radians(24)),math.cos(math.radians(24))])
for a in anchors:
    name=a['track'];p=lines[name];all_points.extend(map(tuple,p))
    vv=np.array([[float(r['x_in'])*72,float(r['y_from_top_in'])*72] for r in vertices if r['connector']==name]);assert np.max(np.abs(vv-p))<1e-9
    assert a['source_chromosome']=='1' and a['source_position_bp']=='0' and a['anchor_side']=='upper'
    r=float(a['anchor_radius']);true=np.array([cx+R*r*math.cos(math.radians(24)),cy-R*r*math.sin(math.radians(24))])
    assert np.max(np.abs(true-np.array([float(a['true_x_in']),float(a['true_y_from_top_in'])])*72))<1e-9
    assert np.max(np.abs(p[0]-true-normal*.14*72))<1e-9
    assert float(a['source_inner_radius'])<r<float(a['source_outer_radius'])
    assert a['clearance_direction']=='end_face_normal_into_opening' and p[0,1]<cy
    arr=samples(p);no_annotations(arr)
    gap=samples(np.array([true,p[0]]));no_annotations(gap)
    rr=np.hypot(arr[:,0]-cx,cy-arr[:,1])/R;ang=np.degrees(np.arctan2(cy-arr[:,1],arr[:,0]-cx))%360
    assert not np.any((rr>=.3185)&(rr<=.9975)&(ang>24+1e-6)&(ang<336-1e-6)),name
    for other in arch:
        if other['track']==name:continue
        x=arr[:,0]/72;y=arr[:,1]/72
        assert not np.any((x>=float(other['left_in']))&(x<=float(other['left_in'])+float(other['width_in']))&(y>=float(other['top_in']))&(y<=float(other['top_in'])+float(other['height_pt'])/72)),(name,other['track'])
    arow=next(q for q in arch if q['track']==name)
    asset=ET.parse(P/'assets'/arow['file']).getroot()
    rects=[e for e in asset.iter() if tag(e)=='rect' and e.get('stroke')==arow['color']]
    visible_x=min(float(e.get('x'))-float(e.get('stroke-width'))/2 for e in rects)
    visible_edge=float(arow['left_in'])+visible_x*7.40625/1400
    assert abs(visible_edge-p[-1,0]/72-.09)<1e-10
assert len(set(all_points))==20
# The second elbows are staggered increasingly left-to-right in module order.
assert all(lines[a][2,0]<lines[b][2,0] for a,b in zip(list(lines),list(lines)[1:]))
def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def intersects(a,b,c,d):
    return max(min(a[0],b[0]),min(c[0],d[0]))<=min(max(a[0],b[0]),max(c[0],d[0])) and max(min(a[1],b[1]),min(c[1],d[1]))<=min(max(a[1],b[1]),max(c[1],d[1])) and orient(a,b,c)*orient(a,b,d)<=0 and orient(c,d,a)*orient(c,d,b)<=0
def pdistance(p,a,b):
    v=b-a;t=max(0,min(1,float((p-a)@v)/float(v@v)));return float(np.linalg.norm(p-a-t*v))
def sdist(a,b,c,d):
    assert not intersects(a,b,c,d),'Connector intersection'
    return min(pdistance(a,c,d),pdistance(b,c,d),pdistance(c,a,b),pdistance(d,a,b))
pairs=[]
for i,(name,p) in enumerate(lines.items()):
    for other,q in list(lines.items())[i+1:]:
        dist=min(sdist(a,b,c,d) for a,b in zip(p,p[1:]) for c,d in zip(q,q[1:]));assert dist>=3.5
        pairs.append(dict(track1=name,track2=other,minimum_clearance_pt=dist))
recorded=table(O/'connector_pair_QC.tsv')
for a,b in zip(recorded,pairs):assert abs(float(a['minimum_clearance_pt'])-b['minimum_clearance_pt'])<1e-9
# Preserve the R2 proof that true TEV endpoint meets an actual source coral polygon.
proof=json.loads((prior/'validation_report.json').read_text())['TEV_origin_verified_against_original_coral_polygon']
assert float(anchors[-1]['anchor_radius'])==float(table(prior/'connector_anchors.tsv')[-1]['anchor_radius'])
assert proof['anchor_boundary_distance_source_pt']<.005
# DEBUG distinguishes true anchor and offset start, never part of formal output.
debug=ET.parse(O/(S+'.DEBUG.svg')).getroot();markers=next(e for e in debug if e.get('id')=='DEBUG-anchor-markers')
for cls,xx,yy in [('true-anchor','true_x_in','true_y_from_top_in'),('connector-start','connector_start_x_in','connector_start_y_from_top_in')]:
    nodes=[e for e in markers if e.get('class')==cls];assert len(nodes)==5
    for a,n in zip(anchors,nodes):assert abs(float(n.get('cx'))-float(a[xx])*72)<1e-9 and abs(float(n.get('cy'))-float(a[yy])*72)<1e-9
assert not any(e.get('id')=='DEBUG-anchor-markers' for e in root.iter())
W=geo['page_width_in'];H=geo['page_height_in']
for suffix,dpi in [('.png',600),('.preview.png',150),('.DEBUG.preview.png',150)]:
    with Image.open(O/(S+suffix)) as im:assert im.size==(round(W*dpi),round(H*dpi)) and all(abs(d-dpi)<1 for d in im.info['dpi'])
with Image.open(O/(S+'.DEBUG.anchor_zoom.png')) as im:assert im.size==(2130,1380)
info=subprocess.check_output(['pdfinfo',str(O/(S+'.pdf'))],text=True);wh=re.search(r'Page size:\s+([0-9.]+) x ([0-9.]+) pts',info)
assert abs(float(wh[1])-math.ceil(W*72))<.02 and abs(float(wh[2])-math.ceil(H*72))<.02
pdf_geo={k:float(v) for k,v in table(O/'PDF_geometry.tsv')[0].items()}
assert pdf_geo['content_scale']==1 and 0<=pdf_geo['right_padding_pt']<1
images=subprocess.check_output(['pdfimages','-list',str(O/(S+'.pdf'))],text=True);assert len(images.strip().splitlines())==2
text=subprocess.check_output(['pdftotext',str(O/(S+'.pdf')),'-'],text=True)
for word in ['LINE','LTR','SINE','DNA','TEV','ORF1','ORF2','Transposase','>80%']:assert word in text,word
fonts=subprocess.check_output(['pdffonts',str(O/(S+'.pdf'))],text=True);(O/'provenance/PDF_fonts.txt').write_text(fonts)
report=dict(status='PASS',version='V3_R3 compact column',source_asset_revision='06e8112',new_five_SVG_assets_used=expected,all_inputs_unchanged=True,protected_old_files_unchanged=len(baseline),V15_vector_contents_identical=True,V15_scale_and_placement_identical=True,V15_clip_paths_preserved=325,
    module_scale_vs_R2=1.5,module_width_in=7.40625,internal_text_scale_vs_R2=1.5,heading_font_source_px=39,heading_final_pt=39*7.40625*72/1400,heading_visual_ratio_vs_R2=heading_ratios[0],
    module_gap_from_Circos_bound_in=.70,module_canvas_vertical_gap_in=.08,column_height_in=column_height,circos_height_in=circos_bottom-circos_top,column_height_ratio=column_height/(circos_bottom-circos_top),column_vertically_centered=True,
    start_clearance_in=.14,start_clearance_direction='endpoint outward normal into gap (user-confirmed)',arrow_tip_to_visible_module_edge_in=.09,independent_connectors=5,straight_segments_each=3,commands_each=['M','L','L','L'],Bezier_commands=0,shared_vertices=0,shared_segments=0,intersections=0,
    minimum_connector_clearance_pt=min(p['minimum_clearance_pt'] for p in pairs),connector_pairs=pairs,connectors_cover_labels=False,connectors_cross_main_data=False,connectors_cross_other_modules=False,
    true_sources_upper_Chr1_endpoints=True,source_track_indices=[4,5,6,7,10],TEV_source='TIP / innermost TEV density',TEV_true_anchor_radius=float(anchors[-1]['anchor_radius']),TEV_true_anchor_angle_degrees=24,TEV_original_coral_boundary_proof=proof,DEBUG_true_anchor_and_start_verified=True,
    TEV_internal_arrow_absent=True,TEV_result_box_absent=True,schematic_threshold=.8,V15_data_threshold_unchanged=.5,genomic_statistics_recalculated=False,SVG_embedded_rasters=0,PDF_embedded_rasters=0,PDF_inches=[float(wh[1])/72,float(wh[2])/72],PDF_right_padding_pt=pdf_geo["right_padding_pt"],PDF_content_scale=1,PNG_pixels=[round(W*600),round(H*600)],PNG_dpi=600,preview_pixels=[round(W*150),round(H*150)])
(O/'validation_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
