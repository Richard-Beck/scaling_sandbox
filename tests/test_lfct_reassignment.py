"""Behavioral checks for avoiding duplicated or ambiguous segmentation assignments."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from build_lfct_mcf10a_cpsam import assign_masks,mask_objects,add_nuclear_overlap_matches
class Assignments(unittest.TestCase):
    def assign(self,m,xy):
        return assign_masks(m,mask_objects(m),[{'x':x,'y':y} for x,y in xy])
    def test_contained_preserves_mask_label(self):
        m=np.zeros((50,50),np.uint16);m[10:20,10:20]=7
        a=self.assign(m,[(15.25,14.75)])[0]
        self.assertEqual((a['label'],a['status']),(7,'contained'))
    def test_shared_mask_not_duplicated(self):
        m=np.zeros((50,50),np.uint16);m[10:30,10:30]=3
        a=self.assign(m,[(15,15),(25,25)])
        self.assertTrue(all(x['label'] is None and x['status']=='shared_mask' for x in a))
    def test_nearby_accepts_only_close_mask(self):
        m=np.zeros((70,70),np.uint16);m[10:20,10:20]=1
        a=self.assign(m,[(23,15),(65,65)])
        self.assertEqual(a[0]['status'],'nearby');self.assertEqual(a[0]['gap'],4)
        self.assertIsNone(a[1]['label'])
    def test_equal_distance_ambiguous(self):
        m=np.zeros((50,50),np.uint16);m[10:20,10:20]=1;m[10:20,30:40]=2
        a=self.assign(m,[(24.5,15)])[0]
        self.assertEqual(a['status'],'ambiguous_nearby');self.assertIsNone(a['label'])
    def test_nearby_does_not_steal_mask(self):
        m=np.zeros((50,50),np.uint16);m[10:20,10:20]=1
        a=self.assign(m,[(15,15),(23,15)])
        self.assertEqual(a[0]['label'],1);self.assertEqual(a[1]['status'],'nearest_mask_occupied')
    def test_competing_nearby_left_unassigned(self):
        m=np.zeros((50,50),np.uint16);m[10:20,10:20]=1
        a=self.assign(m,[(22,15),(15,22)])
        self.assertTrue(all(x['label'] is None and x['status']=='competing_centroids' for x in a))
    def test_overlap_fallback_and_ambiguous_nuclei(self):
        cell=np.zeros((100,100),np.uint16);cell[5:80,5:80]=1
        nuc=np.zeros_like(cell);nuc[60:70,60:70]=4
        cm=self.assign(cell,[(15,15)]);nm=self.assign(nuc,[(15,15)])
        add_nuclear_overlap_matches(cell,nuc,mask_objects(nuc),cm,nm)
        self.assertEqual(nm[0]['status'],'cell_overlap');self.assertEqual(nm[0]['label'],4)
        nuc[40:50,40:50]=5;nm=self.assign(nuc,[(15,15)])
        add_nuclear_overlap_matches(cell,nuc,mask_objects(nuc),cm,nm)
        self.assertIsNone(nm[0]['label'])
    def test_out_of_bounds_never_clamped_into_cell(self):
        m=np.ones((20,20),np.uint16)
        a=self.assign(m,[(-1,10)])[0]
        self.assertEqual(a['status'],'out_of_bounds');self.assertIsNone(a['label'])
if __name__=='__main__':unittest.main()
