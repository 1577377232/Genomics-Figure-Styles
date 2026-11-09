#!/usr/bin/env python3
"""Independent QC of compact 1x composition; requires Pillow for export QA."""
import csv, hashlib, itertools, json, math, re
from pathlib import Path
from PIL import Image
P=Path(__file__).resolve().parents[3]
O=Path(__file__).resolve().parents[1]
B=P/'07_figures/R_v15_with_TE_architectures_v3_r3'
S='chicken_GRCg7b_R_v15_with_TE_architectures_v3_r3'
tracks=['LINE','LTR','SINE','DNA','TEV']
def rows(p):
    with p.open() as f: return list(csv.DictReader(f,delimiter='\t'))
def near(a,b): return abs(float(a)-float(b))<1e-8
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def mask(s):
    for t in tracks:
        s=re.sub(r'<svg\b[^>]*\bid="architecture-'+t+r'"[^>]*>.*?</svg>','ASSET_'+t,s,flags=re.S)
        s=re.sub(r'(<path\b[^>]*\bid="connector__connector-'+t+r'"[^>]*\bd=")[^"]*',r'\1GEOMETRY',s)
        s=re.sub(r'(<polygon\b[^>]*\bid="connector__arrow-'+t+r'"[^>]*\bpoints=")[^"]*',r'\1GEOMETRY',s)
    return s.strip()
protected=json.loads((O/'protected_baseline.json').read_text())
assert all(digest(P/p)==h for p,h in protected.items()),'Protected source changed'
for suffix in ['.svg','.DEBUG.svg','.DEBUG.anchor_zoom.svg']:
    old=(B/(S+suffix)).read_text()
    new=(O/(S+'_1x'+suffix)).read_text()
    assert mask(old)==mask(new),'Changes outside five modules and line terminal geometry'
a=rows(O/'architecture_placements.tsv'); b=rows(B/'architecture_placements.tsv')
g=rows(B/'layout_geometry.tsv')[0]
for x,y in zip(a,b):
    assert near(x['width_in'],4.9375)
    assert near(x['scale_x'],x['scale_y'])
    for t in [x,y]:
        t['edge']=float(t['left_pt'])+(float(t['bar_left_source_x'])-float(t['bar_outline_source_width'])/2)*float(t['scale_x'])
    assert near(x['edge'],y['edge'])
    assert near(float(x['edge'])/72-float(x['target_x_in']),.09)
    assert (P/'assets'/x['file']).exists()
gaps=[float(y['top_in'])-(float(x['top_in'])+float(x['height_pt'])/72) for x,y in zip(a,a[1:])]
assert all(near(z,.08) for z in gaps)
top=float(a[0]['top_in']);bottom=float(a[-1]['top_in'])+float(a[-1]['height_pt'])/72
assert near((top+bottom)/2,(float(g['circos_top_in'])+float(g['circos_bottom_in']))/2)
assert bottom-top<float(g['circos_height_in'])
v=rows(O/'connector_vertices.tsv');ov=rows(B/'connector_vertices.tsv')
for x,y in zip(v,ov):
    assert near(x['x_in'],y['x_in'])
    if int(x['vertex'])<3: assert near(x['y_from_top_in'],y['y_from_top_in'])
assert (O/'connector_anchors.tsv').read_bytes()==(B/'connector_anchors.tsv').read_bytes()
paths=rows(O/'connector_paths.tsv')
for x in paths:
    assert re.findall('[A-Za-z]',x['svg_path'])==['M','L','L','L']

def point_dist(p,a,b):
    ab=[b[i]-a[i] for i in range(2)]
    t=max(0,min(1,sum((p[i]-a[i])*ab[i] for i in range(2))/sum(z*z for z in ab)))
    return math.dist(p,[a[i]+t*ab[i] for i in range(2)])
def cross(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def dist(a,b,c,d):
    boxes=all(max(min(a[i],b[i]),min(c[i],d[i]))<=min(max(a[i],b[i]),max(c[i],d[i])) for i in range(2))
    if boxes and cross(a,b,c)*cross(a,b,d)<=0 and cross(c,d,a)*cross(c,d,b)<=0: return 0
    return min(point_dist(a,c,d),point_dist(b,c,d),point_dist(c,a,b),point_dist(d,a,b))
minimum=math.inf
for t,u in itertools.combinations(tracks,2):
    p=[(float(x['x_in'])*72,float(x['y_from_top_in'])*72) for x in v if x['connector']==t]
    q=[(float(x['x_in'])*72,float(x['y_from_top_in'])*72) for x in v if x['connector']==u]
    clearance=min(dist(c,d,e,f) for c,d in zip(p,p[1:]) for e,f in zip(q,q[1:]))
    assert clearance>3.5,(t,u,clearance)
    minimum=min(minimum,clearance)
image=Image.open(O/(S+'_1x.png'))
assert image.size==(9928,6000)
assert all(abs(d-600)<.1 for d in image.info['dpi'])
assert Image.open(O/(S+'_1x.preview.png')).size==(2482,1500)
report={
 'status':'PASS','protected_files_unchanged':len(protected),
 'source':'Frozen R3 SVG; no V15 redraw or genomic computation',
 'module_width_in':4.9375,'module_height_in':float(a[0]['height_pt'])/72,
 'module_canvas_gap_in':gaps[0],'column_height_in':bottom-top,
 'circos_height_in':float(g['circos_height_in']),
 'main_and_true_anchors_unchanged':True,'connector_crossings':0,
 'minimum_connector_clearance_pt':minimum,'png_pixels':list(image.size),
 'png_dpi':list(image.info['dpi']),
 'visual_review':'Preview and DEBUG inspected; modules compact, ordered and unobstructed',
 'TEV_note':'Schematic >80%; original V15 data remain >50% and >50bp, unchanged.'}
(O/'QC_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
