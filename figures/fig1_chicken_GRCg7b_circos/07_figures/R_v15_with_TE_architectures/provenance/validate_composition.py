#!/usr/bin/env python3
"""Independent vector-content, source-preservation and export audit."""
from pathlib import Path
import collections,csv,hashlib,json,math,re,subprocess
import xml.etree.ElementTree as ET
from PIL import Image

OUT=Path(__file__).resolve().parents[1]
P=OUT.parents[1]
S='chicken_GRCg7b_R_v15_with_TE_architectures'
SVG='{http://www.w3.org/2000/svg}'
def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tag(e):return e.tag.rsplit('}',1)[-1]
def normalized(e,prefix=''):
    attrs={k:v.replace(prefix,'') if prefix else v for k,v in e.attrib.items() if not(k=='id' and v=='')}
    return (e.tag,attrs,(e.text or '').strip(),[normalized(c,prefix) for c in e])
def body(e,prefix=''):return [normalized(c,prefix) for c in e]
root=ET.parse(OUT/(S+'.svg')).getroot()
main=next(e for e in root if e.get('id')=='V15-circos')
original=next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.svg'))
source=ET.parse(original).getroot()
assert body(source)==body(main,'V15-circos__'),'Original V15 vector contents changed'
assert main.get('preserveAspectRatio')=='xMinYMin meet'
placements=table(OUT/'architecture_placements.tsv')
for r in placements:
    element=next(e for e in root if e.get('id')=='architecture-'+r['track'])
    asset=ET.parse(P/'assets/TE_architectures'/r['file']).getroot()
    assert body(asset)==body(element,'architecture-'+r['track']+'__'),r['track']
    assert abs(float(r['scale_x'])-float(r['scale_y']))<1e-12
ids=[e.get('id') for e in root.iter() if e.get('id')]
assert len(set(ids))==len(ids),'Duplicate SVG IDs'
refs=[]
for e in root.iter():
    for k,v in e.attrib.items():
        refs.extend(re.findall(r'url\(\s*[\'"]?#([^\)\'\"]+)',v))
        if k.endswith('href') and v.startswith('#'):refs.append(v[1:])
assert set(refs)<=set(ids),set(refs)-set(ids)
assert not any(tag(e)=='image' for e in root.iter()),'Raster image embedded into SVG'
source_clips=sum(tag(e)=='clipPath' for e in source.iter())
main_clips=sum(tag(e)=='clipPath' for e in main.iter())
assert source_clips==main_clips and source_clips>0
assert not any(tag(e)=='text' and e.text and 'SVA' in e.text for e in root.iter())
connectors=next(e for e in root if e.get('id')=='connector-layer')
assert sum(tag(e)=='polygon' for e in connectors.iter())==4,'Expected four arrowheads'
assert sum(tag(e)=='polyline' for e in connectors.iter())==5,'Expected one stem and four branches'
for r in table(OUT/'input_preservation.tsv'):
    f=Path(r['path']);assert f.is_file(),f
    assert hashlib.md5(f.read_bytes()).hexdigest()==r['md5'],f
    assert r['unchanged']=='TRUE'
rows=table(OUT/'connector_vertices.tsv')
groups=collections.defaultdict(list)
for r in rows:groups[r['connector']].append((float(r['x_in'])*72,float(r['y_from_top_in'])*72))
# Independent geometric check against all sector/track panels and annotation boxes.
geom=(OUT/'source_geometry.R').read_text()
cx,cy=map(float,re.search(r'circle_center_pt = c\(([^,]+), ([^)]+)\)',geom).groups())
radius=float(re.search(r'circle_radius_pt = ([0-9.]+)',geom)[1])
positions=table(OUT/'circos_placement.tsv')[0]
sx=float(positions['width_pt'])/705;sy=float(positions['height_pt'])/720
left=float(positions['left_pt']);top=float(positions['top_pt'])
unitx=705/9.8*sx;unity=72*sy
annotations=[]
for suffix in ['chromosome_labels.tsv','coordinate_labels.tsv','TE_axis_labels.tsv','variant_axis_labels.tsv']:
    file=next((P/'07_figures/R_v15').glob('chicken_GRCg7b_R_v15_1000kb*.'+suffix))
    for r in table(file):
        if r.get('label_status') in ['hidden','omitted']:continue
        if 'center_x_in' not in r:continue
        a=dict(cx=left+float(r['center_x_in'])*unitx,cy=top+(10-float(r['center_y_in']))*unity,
            hw=float(r['width_in'])*unitx/2+1.2,hh=float(r['height_in'])*unity/2+1.2,
            tx=float(r.get('tangent_x',1)),ty=-float(r.get('tangent_y',0)),
            rx=float(r.get('radial_x',0)),ry=-float(r.get('radial_y',1)),label=r['label'])
        annotations.append(a)
for name,vertices in groups.items():
    for a,b in zip(vertices,vertices[1:]):
        for i in range(501):
            t=i/500;x=a[0]+t*(b[0]-a[0]);y=a[1]+t*(b[1]-a[1])
            rr=math.hypot(x-cx,cy-y)/radius;angle=math.degrees(math.atan2(cy-y,x-cx))%360
            assert not(.3185<=rr<=.8895 and 24<=angle<=336),(name,'crosses data')
            for ann in annotations:
                dx=x-ann['cx'];dy=y-ann['cy']
                assert not(abs(dx*ann['tx']+dy*ann['ty'])<ann['hw'] and
                           abs(dx*ann['rx']+dy*ann['ry'])<ann['hh']),(name,ann['label'],'annotation collision')
assert groups['stem'][0][1]>cy,'Origin must lie below circle center'
with Image.open(OUT/(S+'.png')) as im:
    assert im.size==(8400,6000)
    assert all(abs(d-600)<1 for d in im.info['dpi'])
with Image.open(OUT/(S+'.preview.png')) as im:assert im.size==(2100,1500)
info=subprocess.check_output(['pdfinfo',str(OUT/(S+'.pdf'))],text=True)
assert re.search(r'Page size:\s+1008 x 720 pts',info),info
images=subprocess.check_output(['pdfimages','-list',str(OUT/(S+'.pdf'))],text=True)
assert len(images.strip().splitlines())==2,'Embedded raster image found in PDF'
text=subprocess.check_output(['pdftotext',str(OUT/(S+'.pdf')),'-'],text=True)
for word in ['LINE','LTR','SINE','DNA','ORF1','ORF2','Transposase','Env']:
    assert word in text,word
fonts=subprocess.check_output(['pdffonts',str(OUT/(S+'.pdf'))],text=True)
assert 'LiberationSans' in fonts
(OUT/'provenance/PDF_fonts.txt').write_text(fonts)
report=dict(status='PASS',V15_vector_contents_identical=True,all_four_architecture_vector_contents_identical=True,
    original_inputs_unchanged=len(table(OUT/'input_preservation.tsv')),V15_clip_paths_preserved=source_clips,
    duplicate_ids=0,unresolved_SVG_references=0,arrowheads=4,stems=1,branches=4,
    connectors_cross_data=False,connectors_cross_annotation_boxes=False,origin_below_circle_center=True,
    original_asset_text_nodes_retained=sum(tag(e)=='text' for e in root.iter()),
    images_embedded_in_SVG=0,images_embedded_in_PDF=0,PDF_inches=[14,10],PNG_pixels=[8400,6000],PNG_dpi=600,
    preview_pixels=[2100,1500],genomic_data_recalculated=False)
(OUT/'validation_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))

# Refresh the final artifact index after this independent QC writes its report.
paths=sorted(f for f in OUT.iterdir() if f.is_file() and f.name!='output_manifest.tsv'
    and f.suffix not in ['.log','.zip'])
with (OUT/'output_manifest.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['path','bytes','md5'])
    for p in paths:w.writerow([p.name,p.stat().st_size,hashlib.md5(p.read_bytes()).hexdigest()])
