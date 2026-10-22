from pathlib import Path
import sys
import unittest
import tempfile
import gzip
import numpy as np
from PIL import Image
from pypdf import PdfReader,PdfWriter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
from signature_figure4_atlas_20261020 import physical_faces,selection_features,choose,FEATURES
from signature_figure4_atlas_pdf_20261020 import add_exact_png
from signature_figure4_atlas_csv_20261020 import WideShards
from signature_figure4_atlas_final_20261020 import fullstate_identity,forcing_identity,run_identities,observation_equivalence,validate_contexts,production_resource_choice,LAW_KEYS


class AtlasTests(unittest.TestCase):
    def test_physical_faces_preserve_native_spring_circumference_extremum(self):
        x,t,rows=physical_faces([0.,2.,4.],[2.,4.,3.],4)
        np.testing.assert_allclose(t[:,0],[0,1,2,3,4])
        np.testing.assert_allclose(x[:,-1],[2,3,4,3.5,3])
        np.testing.assert_array_equal(rows,[0,1,1,2])
        native=np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]])
        np.testing.assert_array_equal(native[rows][1],native[1])
        self.assertEqual(set(rows),{0,1,2})

    def test_feature_medoids_use_unit_invariant_scaled_distances(self):
        matrix=np.array([[0]*7,[1]*7,[2]*7],float)
        first=choose(matrix,['a','b','c'],'condition',123)
        changed=choose(matrix*np.arange(1,8)*100,['a','b','c'],'condition',123)
        self.assertEqual(first['medoid'],1)
        self.assertEqual(first['medoid'],changed['medoid'])
        self.assertNotEqual(first['random'],first['medoid'])
        self.assertEqual(first,choose(matrix,['a','b','c'],'condition',123))

    def test_zero_iqr_nonconstant_feature_uses_range_and_singleton_has_no_random(self):
        matrix=np.zeros((5,7));matrix[-1,0]=5.
        result=choose(matrix,['e','d','c','b','a'],'condition',123)
        self.assertEqual(result['scales'][0],5.)
        self.assertEqual(result['medoid'],3)
        self.assertIsNone(choose(np.zeros((1,7)),['one'],'condition',123)['random'])

    def test_last_half_native_feature_weights_and_qualification(self):
        t=np.arange(5,dtype=float)
        out=selection_features(t,np.array([2,2,2,.05,.05]),np.ones(5,dtype=complex),np.array([1,1,2,1,1]),1.,4.)
        self.assertAlmostEqual(out['multiple_fraction'],.5)
        self.assertAlmostEqual(out['low_contrast_fraction'],.5)
        self.assertAlmostEqual(out['mean_qualified_count'],1.)
        self.assertAlmostEqual(out['P_coverage_01'],1.)
        self.assertEqual(set(FEATURES)-set(out),set())

    def test_pdf_native_png_is_lossless_flate_with_no_quantization(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);png=root/'source.png';pdf=root/'native.pdf'
            rgba=np.array([[[1,2,3,255],[253,254,255,255]],[[8,7,6,255],[200,100,0,255]]],dtype=np.uint8)
            Image.fromarray(rgba).save(png)
            writer=PdfWriter();record=add_exact_png(writer,png);writer.write(pdf)
            image=PdfReader(pdf).pages[0]['/Resources']['/XObject']['/NativePNG'].get_object()
            self.assertEqual(image['/Filter'],'/FlateDecode')
            self.assertEqual(image.get_data(),rgba[:,:,:3].tobytes())
            self.assertTrue(record['opaque_alpha_255'])
            rgba[0,0,3]=254;Image.fromarray(rgba).save(png)
            with self.assertRaisesRegex(ValueError,'alpha may not be discarded'):
                add_exact_png(PdfWriter(),png)

    def test_wide_native_csv_shards_preserve_double_bits_and_headers(self):
        with tempfile.TemporaryDirectory() as directory:
            values=np.array([[i,np.nextafter(1.+i,2.+i),-0.] for i in range(12)],dtype=float)
            writer=WideShards(directory,'native','frame,active0,active1',limit=140)
            np.savetxt(writer,values,fmt=['%d','%.17g','%.17g'],delimiter=',');writer.close()
            self.assertGreater(len(writer.files),1)
            arrays=[]
            for part in writer.files:
                with gzip.open(part['path'],'rt') as stream:
                    self.assertEqual(stream.readline().strip(),'frame,active0,active1')
                    arrays.append(np.loadtxt(stream,delimiter=',',ndmin=2))
            self.assertEqual(np.vstack(arrays).tobytes(),values.tobytes())

    def test_complete_initial_state_includes_inactive_pools_but_not_zip_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);x=np.arange(4)+.5;state=np.array([[1.,2,3,4],[5.,6,7,8]])
            np.savez(root/'a.npz',state=state,x=x,label='first container')
            np.savez_compressed(root/'b.npz',state=state,x=x,label='other metadata')
            a=fullstate_identity(root/'a.npz','wave_pinning',4.)
            self.assertEqual(a,fullstate_identity(root/'b.npz','wave_pinning',4.))
            state[1,0]+=1.;np.savez(root/'c.npz',state=state,x=x)
            self.assertNotEqual(a['fullstate_and_geometry_sha256'],fullstate_identity(root/'c.npz','wave_pinning',4.)['fullstate_and_geometry_sha256'])

    def test_spring_restart_identity_retains_geometry_and_both_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);state=np.r_[np.ones(4),np.ones(4)*2,3.]
            np.savez(root/'a.npz',state=state,material_positions=(np.arange(4)+.5)/4)
            a=fullstate_identity(root/'a.npz','spring_ring',3.)
            state[-1]=4.;np.savez(root/'b.npz',state=state,material_positions=(np.arange(4)+.5)/4)
            self.assertNotEqual(a['fullstate_and_geometry_sha256'],fullstate_identity(root/'b.npz','spring_ring',3.)['fullstate_and_geometry_sha256'])
            np.savez(root/'missing.npz',state=np.ones(4),material_positions=(np.arange(4)+.5)/4)
            with self.assertRaisesRegex(ValueError,'both fields'):fullstate_identity(root/'missing.npz','spring_ring',3.)

    def test_sigma_pairing_and_seed_independence_are_distinct_identities(self):
        law={k:'declared' for k in LAW_KEYS};law.update(gaussian_knots_sha256='actual_knots',seed=11,sigma=.5)
        case=dict(kind='noise',model='goryachev',C_label=20.,f=.1,sigma=.5,seed=11,native_horizon=32000.,profile_cadence=5.,origin='E1')
        initial=dict(fullstate_and_geometry_sha256='all_native_species_and_geometry')
        forcing=forcing_identity(case,dict(job=dict(mean=.2)),law)
        one=run_identities(case,initial,forcing,law)
        other=dict(case,sigma=2.,profile_cadence=40.,origin='other')
        otherlaw=dict(law,sigma=2.)
        two=run_identities(other,initial,forcing_identity(other,dict(job=dict(mean=.2)),otherlaw),otherlaw)
        self.assertNotEqual(one[0],two[0]);self.assertEqual(one[1],two[1])
        companion=dict(case,profile_cadence=40.,origin='E0')
        self.assertEqual(one,run_identities(companion,initial,forcing,law))
        fresh=dict(case,seed=12);freshlaw=dict(law,seed=12)
        three=run_identities(fresh,initial,forcing_identity(fresh,dict(job=dict(mean=.2)),freshlaw),freshlaw)
        self.assertNotEqual(one[0],three[0]);self.assertNotEqual(one[1],three[1])
        extended=dict(case,native_horizon=64000.)
        self.assertNotEqual(one[0],run_identities(extended,initial,forcing,law)[0])
        changedpath=dict(law,gaussian_knots_sha256='different_actual_knots')
        self.assertNotEqual(one[0],run_identities(case,initial,forcing,changedpath)[0])
        self.assertNotEqual(one[1],run_identities(case,initial,forcing,changedpath)[1])

    def test_deterministic_controls_have_no_effective_stochastic_seed(self):
        case=dict(kind='mean',sigma=1.,seed=19)
        forcing=forcing_identity(case,dict(job=dict(mean=.2)),{})
        self.assertIsNone(forcing['effective_stochastic_seed']);self.assertEqual(forcing['sigma'],0.)
        self.assertEqual(forcing,forcing_identity(dict(case,seed=200),dict(job=dict(mean=.2)),{}))

    def test_native_companion_checks_shared_times_and_inactive_arrays(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);field=np.arange(12.).reshape(3,4);inactive=field+100
            np.savez(root/'old.npz',time=np.array([0.,2.,4.]),primary=field,inactive=inactive,x=np.arange(4)+.5)
            fine=np.zeros((5,4));fine[::2]=field;v=np.zeros_like(fine);v[::2]=inactive
            np.savez(root/'new.npz',time=np.arange(5.),primary=fine,inactive=v,x=np.arange(4)+.5)
            checks=observation_equivalence(root/'old.npz',root/'new.npz');self.assertTrue(all(c['bitwise_equal'] for c in checks))
            v[2,0]+=1.;np.savez(root/'wrong.npz',time=np.arange(5.),primary=fine,inactive=v,x=np.arange(4)+.5)
            with self.assertRaises(AssertionError):observation_equivalence(root/'old.npz',root/'wrong.npz')

    def test_matching_unavailable_scalar_values_do_not_create_nan_audit_receipts(self):
        import json
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('old','new'):np.savez(root/(name+'.npz'),time=[0.,1.],heading=[0.,np.nan])
            checks=observation_equivalence(root/'old.npz',root/'new.npz')
            json.dumps(checks,allow_nan=False)
            self.assertEqual(next(c for c in checks if c['field']=='heading')['matching_nonfinite_values'],1)

    def test_final_coverage_cannot_be_inferred_from_completed_prototype(self):
        with self.assertRaisesRegex(ValueError,'three frozen final contexts'):validate_contexts({'E0v2':[]})

    def test_post_pilot_resource_override_retains_measured_margin(self):
        records=[dict(slot=0,worker_sha256='verified_native_worker',resource=dict(peak_resident_memory_KiB=2*1024**2))]
        result=production_resource_choice(16,records,4)
        self.assertEqual(result['production_allocation_GiB'],4)
        self.assertTrue(result['scientific_inputs_and_native_values_unchanged'])
        with self.assertRaisesRegex(ValueError,'measured pilot margin'):production_resource_choice(16,records,3)
        with self.assertRaisesRegex(ValueError,'positive native pilot'):production_resource_choice(16,[])


if __name__=='__main__':unittest.main()
