import sys,tempfile,unittest
from pathlib import Path
import numpy as np,tifffile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from geometry import mask_objects,assign_masks,pixel_square_gap,clearances
from analyze import frame_analysis,eligible_runs
class Core(unittest.TestCase):
    def test_source_mapping(self):
        from prepare import source_indices
        self.assertEqual(source_indices('test',3,{0,1,2}),([0,1,2],[]))
        self.assertEqual(source_indices('test',3,{1,2,3}),([1,2,3],[]))
        with self.assertRaises(ValueError):source_indices('test',3,{0,1,2,3})
        indices,skip=source_indices('u87/FLD_1',288,set(range(289)))
        self.assertEqual([indices[i] for i in [0,1,144,287]],[0,2,145,288]);self.assertEqual(skip,[1])
    def test_touching_squares_zero(self):
        a=np.zeros((30,30),np.uint16);a[5:10,5:10]=1;a[5:10,10:15]=2
        objs=mask_objects(a);self.assertEqual(pixel_square_gap(*objs),0)
    def test_diagonal_touch_zero(self):
        a=np.zeros((30,30),np.uint16);a[5:10,5:10]=1;a[10:15,10:15]=2
        self.assertEqual(pixel_square_gap(*mask_objects(a)),0)
    def test_known_empty_gap(self):
        a=np.zeros((60,60),np.uint16);a[5:10,5:10]=1;a[5:10,20:25]=2
        self.assertEqual(pixel_square_gap(*mask_objects(a)),10)
    def test_gap_matches_brute_force(self):
        # Irregular, disjoint masks: no reliance on bounding boxes for final distance.
        a=np.zeros((30,40),np.uint16);a[4:12,4:6]=1;a[10:12,4:14]=1;a[15:22,17:23]=2
        objs=mask_objects(a);p=np.argwhere(a==1);q=np.argwhere(a==2);delta=np.maximum(np.abs(p[:,None,:]-q[None,:,:])-1,0)
        expected=float(np.sqrt((delta*delta).sum(2)).min());self.assertEqual(pixel_square_gap(*objs),expected)
    def test_shared_mask_never_duplicated(self):
        a=np.zeros((30,30),np.uint16);a[5:20,5:20]=1
        m=assign_masks(a,mask_objects(a),[{'x':8,'y':8},{'x':15,'y':15}])
        self.assertTrue(all(r['status']=='shared_mask' and r['label'] is None for r in m))
    def test_all_masks_are_neighbors_even_untracked(self):
        a=np.zeros((100,100),np.uint16);a[30:40,30:40]=1;a[30:40,42:50]=2
        s=[dict(id=5,frame=0,x=35,y=35,track=1,name='Track',branch=5,original=[])]
        rows,_=frame_analysis(a,np.zeros_like(a),s)
        self.assertEqual(rows[0]['clearance_px'],2);self.assertFalse(rows[0]['eligible']);self.assertIn('neighbor_too_close',rows[0]['exclusions'])
    def test_same_track_id_still_neighbor(self):
        a=np.zeros((100,100),np.uint16);a[30:40,30:40]=1;a[30:40,42:52]=2
        s=[dict(id=i,frame=0,x=x,y=35,track=7,name='Shared track',branch=i,original=[]) for i,x in [(1,35),(2,47)]]
        rows,_=frame_analysis(a,np.zeros_like(a),s);self.assertTrue(all(not r['eligible'] for r in rows))
    def test_missing_nucleus_alone_does_not_reject(self):
        a=np.zeros((150,150),np.uint16);a[60:80,60:80]=1
        s=[dict(id=1,frame=0,x=70,y=70,track=7,name='Track',branch=1,original=[])]
        rows,_=frame_analysis(a,np.zeros_like(a),s);self.assertTrue(rows[0]['eligible'])
    def test_unmatched_nearby_centroid_is_uncertainty(self):
        a=np.zeros((150,150),np.uint16);a[60:80,60:80]=1
        s=[dict(id=i,frame=0,x=x,y=70,track=i,name='Track',branch=i,original=[]) for i,x in [(1,70),(2,105)]]
        rows,_=frame_analysis(a,np.zeros_like(a),s);self.assertIn('uncertain_neighbor',rows[0]['exclusions'])
    def test_orphan_nucleus_is_uncertainty(self):
        a=np.zeros((150,150),np.uint16);a[60:80,60:80]=1;n=np.zeros_like(a);n[60:70,105:115]=1
        s=[dict(id=1,frame=0,x=70,y=70,track=1,name='Track',branch=1,original=[])]
        rows,_=frame_analysis(a,n,s);self.assertIn('uncertain_neighbor',rows[0]['exclusions'])
    def test_edge_rejected(self):
        a=np.zeros((100,100),np.uint16);a[30:40,0:10]=1
        s=[dict(id=1,frame=0,x=5,y=35,track=1,name='Track',branch=1,original=[])]
        rows,_=frame_analysis(a,np.zeros_like(a),s);self.assertIn('field_edge',rows[0]['exclusions'])
    def make_rows(self,frames):
        return [dict(frame=f,branch_start_spot_id=1,track_id=8,track_name='Track_123',spot_id=100+f,clearance_px=20.,eligible=True) for f in frames]
    def test_ten_consecutive_required(self):
        runs,_=eligible_runs(self.make_rows(range(10)),10);self.assertEqual(runs[0]['n_frames'],10)
        runs,_=eligible_runs(self.make_rows([0,1,2,3,4,6,7,8,9,10]),10);self.assertEqual(runs,[])
    def test_invalid_frame_breaks_run(self):
        rows=self.make_rows(range(12));rows[6]['eligible']=False;runs,longest=eligible_runs(rows,10)
        self.assertEqual(runs,[]);self.assertEqual(longest,6)
    def test_branches_never_combined(self):
        rows=self.make_rows(range(12))
        for r in rows[6:]:r['branch_start_spot_id']=2
        self.assertEqual(eligible_runs(rows,10)[0],[])
    def test_compressed_stack_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'masks.tif';a=np.arange(30*40,dtype=np.uint16).reshape(30,40)
            with tifffile.TiffWriter(p) as w:
                for f in [1,13]:w.write(a+f,compression='deflate',description=str(f),metadata=None)
            with tifffile.TiffFile(p) as t:
                self.assertEqual([int(x.description) for x in t.pages],[1,13]);np.testing.assert_array_equal(t.pages[1].asarray(),a+13)
if __name__=='__main__':unittest.main()
