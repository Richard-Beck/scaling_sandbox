#!/usr/bin/env python3
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np,tifffile
from PIL import Image,ImageDraw
from skimage.segmentation import find_boundaries
import json
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/lfct_mcf10a_fld102'
def main():
    m=json.loads((DATA/'source_manifest.json').read_text());doc=ET.parse(ROOT/m['xml']['path'])
    canvas=Image.new('RGB',(2448,1808),'white')
    for row,f in enumerate([192,384]):
        spots=doc.find(f'Model/AllSpots/SpotsInFrame[@frame="{f}"]')
        for col,ch in enumerate(['phase','rfp']):
            a=tifffile.imread(ROOT/m['frames'][f][ch]['path']);lo,hi=np.percentile(a,[1,99.5]);a=(np.clip((a.astype(float)-lo)/(hi-lo),0,1)*255).astype('uint8');rgb=np.stack([a]*3,-1)
            masks=tifffile.imread(DATA/'cpsam'/ch/f'{f:04d}.tif');rgb[find_boundaries(masks)]=[0,255,100]
            im=Image.fromarray(rgb);draw=ImageDraw.Draw(im)
            labs=[]
            for s in spots:
                x,y=float(s.get('POSITION_X')),float(s.get('POSITION_Y'));labs.append(int(masks[round(y),round(x)]));draw.ellipse((x-3,y-3,x+3,y+3),fill='red')
            print(f,ch,'masks',int(masks.max()),'tracked',len(spots),'inside',sum(l>0 for l in labs),'unique inside',len(set(labs)-{0}))
            draw.text((12,12),f'{ch} XML frame {f} green=mask red=tracked centroid',fill='yellow');canvas.paste(im,(col*1224,row*904))
    canvas.resize((1468,1084)).save(DATA/'pilot_overlay.jpg',quality=94)
if __name__=='__main__':main()
