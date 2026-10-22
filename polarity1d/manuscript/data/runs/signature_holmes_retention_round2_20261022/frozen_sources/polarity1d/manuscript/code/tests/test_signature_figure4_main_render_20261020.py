"""Native main examples retain measured float64 values, identity and geometry."""
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import matplotlib.pyplot as plt
from PIL import Image

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure4_main_render_20261020 as M


class MainNativeTests(unittest.TestCase):
    def fixture(self,root,size=20.):
        root.mkdir(parents=True,exist_ok=True)
        time=np.asarray([0.,5.,10.]);C=np.asarray([size,size*1.1,size])
        field=np.asarray([[np.nextafter(.1,1.),.2],[.3,.4],[.5,np.nextafter(.6,0.)]])
        native=root/"field.csv.gz"
        with gzip.open(native,"wt",newline="") as stream:
            writer=csv.writer(stream);writer.writerow(["case_index","native_frame","time","actualC","active_cell_0000","active_cell_0001"])
            for i in range(3):writer.writerow([1,i,format(time[i],".17g"),format(C[i],".17g"),*[format(v,".17g") for v in field[i]]])
        selection=root/"selection.csv"
        record=dict(case_index=1,case_id=f"seed_C{size:g}",active_field="native active",model="wave_pinning",C_label=size,f=.03,sigma=1.,kind="noise",native_horizon=10.,profile_cadence=5.)
        with selection.open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(record));writer.writeheader();writer.writerow(record)
        q=root/"q.csv.gz"
        with gzip.open(q,"wt",newline="") as stream:
            writer=csv.writer(stream);writer.writerow(["case_index","native_cell_index","q"]);writer.writerows([[1,0,.25],[1,1,.75]])
        file_record=lambda path:dict(path=str(path),sha256=M.P.sha(path))
        page=dict(slot=0,passed=True,selected_cases=[record["case_id"]],condition=record,decimal_significant_digits=17,
            float64_roundtrip_equal=True,sourcefield_resampling=False,sourcefield_quantization=False,
            selection_CSV=file_record(selection),companions=dict(q_grid=file_record(q)),
            expected_native_binary_hashes={"1":dict(field_sha256=M.binary_sha(field),time_sha256=M.binary_sha(time),actualC_sha256=M.binary_sha(C),rows=3,cells=2)},
            native_field_CSV_shards=[dict(file_record(native),rows=3,roundtrip_float64_equal=True)])
        index=root/"index.json";index.write_text(json.dumps(dict(summary=dict(passed=True),pages=[page])))
        return file_record(index),record["case_id"],field,C

    def test_roundtrip_exact_values_and_native_circumference_extrema(self):
        with tempfile.TemporaryDirectory() as temp:
            record,case,field,C=self.fixture(Path(temp))
            movie=M.selected_native(record,0,case,Path(temp))
            np.testing.assert_array_equal(movie["field"],field)
            x,t,strips=M.Atlas.physical_faces(movie["time"],movie["C"],2)
            self.assertEqual(float(x.max()),float(C.max()))
            self.assertEqual(set(strips),{0,1,2})
            self.assertEqual(movie["evidence"]["native_values"],6)

    def test_rehashed_csv_cannot_override_independent_native_hashes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);record,case,_,_=self.fixture(root)
            with gzip.open(root/"field.csv.gz","rt") as stream:text=stream.read()
            with gzip.open(root/"field.csv.gz","wt") as stream:stream.write(text.replace("0.40000000000000002","0.41"))
            index=json.loads((root/"index.json").read_text());index["pages"][0]["native_field_CSV_shards"][0]["sha256"]=M.P.sha(root/"field.csv.gz")
            (root/"index.json").write_text(json.dumps(index));record["sha256"]=M.P.sha(root/"index.json")
            with self.assertRaisesRegex(ValueError,"independent native binary hashes"):
                M.selected_native(record,0,case,root)

    def test_shared_limits_require_matched_model_field_input_and_actual_sizes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);a,aid,_,_=self.fixture(root/"small",20);b,bid,_,_=self.fixture(root/"large",60)
            small=M.selected_native(a,0,aid,root);large=M.selected_native(b,0,bid,root)
            limits=M.paired_limits(small,large);self.assertEqual(limits["vmin"],0.)
            self.assertEqual(limits["vmax"],float(small["field"].max()))
            self.assertEqual(len(limits["native_binary_source_hashes"]),2)
            self.assertEqual(limits['physical_arc_extent'],[0.,66.])
            self.assertEqual(limits['native_time_extent'],[0.,10.])
            large['time'][-1]=20.
            with self.assertRaisesRegex(ValueError,'no hidden clock rescaling'):
                M.paired_limits(small,large)
            large['time'][-1]=10.
            large["evidence"]["selection"]["f"]=".1"
            with self.assertRaisesRegex(ValueError,"same model/field/noise/input"):
                M.paired_limits(small,large)

    def test_incomplete_final_statistics_refuses_before_rendering(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);producer=root/"producer.json";producer.write_text(json.dumps(dict(passed=True,complete=False)))
            spec=dict(schema="figure4_main_candidate_v1",final_statistics=dict(producer=dict(path=str(producer),sha256=M.P.sha(producer))))
            with self.assertRaisesRegex(ValueError,"predicate failed"):
                M.validated_inputs(spec,root)

    def test_native_panel_preserves_shared_limits_and_fits_actual_text(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);a,aid,_,_=self.fixture(root/"small",20);b,bid,_,_=self.fixture(root/"large",60)
            small=M.selected_native(a,0,aid,root);large=M.selected_native(b,0,bid,root)
            pdf,records,bounds=M.native_panel([(small,large)],root/"native.pdf")
            self.assertTrue(bounds["passed"])
            self.assertEqual(records[0]["shared_color_limits"],records[1]["shared_color_limits"])
            self.assertFalse(records[0]["concentration_interpolated"])
            self.assertTrue(Path(pdf["path"]).exists())
            clipped,cliprecords,_=M.native_panel([(small,large)],root/'clipped.pdf',{'wave_pinning':[5.,10.]})
            self.assertEqual(cliprecords[0]['full_native_time_extent'],[0.,10.])
            self.assertEqual(cliprecords[0]['display_native_time_extent'],[5.,10.])
            np.testing.assert_array_equal(small['field'],M.selected_native(a,0,aid,root)['field'])

    def test_actual_composition_keeps_native_artwork_under_transparent_label_layer(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);a,aid,_,_=self.fixture(root/'small',20);b,bid,_,_=self.fixture(root/'large',60)
            measured=root/'measured.pdf';fig,ax=plt.subplots(figsize=(12.,2.7));ax.set_facecolor('#205090');fig.savefig(measured);plt.close(fig)
            panels=[dict(id=letter,pdf=dict(path=str(measured),sha256=M.P.sha(measured))) for letter in 'BCD']
            spec=dict(native_pairs=[dict(small=dict(atlas_index=a,slot=0,case_id=aid),large=dict(atlas_index=b,slot=0,case_id=bid))],statistical_panels=panels,framing='SYNTHETIC layout regression',no_clean_mapping=True)
            path=root/'spec.json';path.write_text(json.dumps(spec))
            with patch.object(M,'validated_inputs',return_value={}):result=M.render(path,root/'output','composition')
            pixels=np.asarray(Image.open(result['main']['png']['path']).convert('RGB'))
            self.assertGreater(np.mean(np.max(pixels,axis=2)-np.min(pixels,axis=2)>40),.05)
            slots=result['composition']['slots']
            self.assertAlmostEqual(result['composition']['height_points'],
                slots[0]['rectangle'][3]+slots[1]['rectangle'][3]+max(row['rectangle'][3] for row in slots[2:])+120.)
            self.assertTrue(result['composition']['aspect_ratio_preserved'])


if __name__=="__main__":unittest.main()
