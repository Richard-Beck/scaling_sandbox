"""Protocol-grid and phase-application controls without native integration."""
from pathlib import Path
from types import SimpleNamespace
import collections
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_campaign_20261020 as campaign
import signature_noise_campaign_20261020 as noise_campaign
from signature_inputs_20261020 import AreaCue, make_train


class NoiseCampaignControls(unittest.TestCase):
    def test_full_grid_has_all_noise_and_matched_controls(self):
        rows=noise_campaign.condition_grid()
        self.assertEqual(collections.Counter(row["kind"] for row in rows),
                         {"noise":432,"mean":36,"zero":18})
        self.assertEqual(len({row["seed"] for row in rows if row["kind"]=="noise"}),432)
        self.assertEqual({row["requested_reference_length"] for row in rows},{19,29,39})

    def test_pilot_is_exact_full_condition_subset(self):
        full=noise_campaign.condition_grid();pilot=noise_campaign.condition_grid("pilot")
        def key(row):return tuple((name,value) for name,value in row.items() if name!="task")
        self.assertEqual(len(pilot),12)
        self.assertTrue({key(row) for row in pilot}<={key(row) for row in full})
        self.assertEqual(collections.Counter(row["kind"] for row in pilot),
                         {"noise":4,"mean":4,"zero":4})

    def test_bounded_pilot_covers_longest_two_model_conditions(self):
        rows=noise_campaign.condition_grid("bounded_pilot")
        self.assertEqual(len(rows),6)
        self.assertEqual({row["requested_reference_length"] for row in rows},{39})
        self.assertEqual({row["model"] for row in rows},{"wave_pinning","debelly"})
        self.assertEqual(collections.Counter(row["kind"] for row in rows),
                         {"noise":2,"mean":2,"zero":2})
        self.assertEqual({row["gap_multiple_of_pulse"] for row in rows},{16.})
        config=noise_campaign.protocol("bounded_tanh")
        self.assertEqual(config["hard_maximum_fraction_of_native_peak"],.04)
        self.assertEqual(config["input_transform"],"bounded_tanh")

    def test_absolute_clock_noise_is_only_applied_during_gaps(self):
        schedule=make_train(4.,1.5,4,4.)
        x=(np.arange(201)+.5)*19./201
        class Captured(Exception):pass
        def spy(model,length,state,phases,output_dir=None,cells=None,max_step=None):
            self.assertEqual(max_step,1.25)
            for phase in phases:
                value=phase["cue"](.25,x,19.)
                if phase["name"].endswith("_gap"):
                    index=int(phase["name"].split("_")[0]);start=schedule.segments[index][1]
                    np.testing.assert_allclose(value,start+.25)
                elif phase["name"]=="recovery":np.testing.assert_array_equal(value,0.)
                else:self.assertAlmostEqual(float(value.sum()*19./201),1.)
            raise Captured
        native=SimpleNamespace(run_phases=spy)
        forcing=lambda time,x,length:np.full_like(x,time)
        case=dict(model="wave_pinning",native_length=19.,cells=201,cue_area=1.,cue_width=3.8,recovery_duration=10.)
        with tempfile.TemporaryDirectory() as name:
            with self.assertRaises(Captured):
                campaign._run_trial(native,SimpleNamespace(AreaCue=AreaCue),None,case,None,None,
                    Path(name)/"trial",schedule=schedule,campaign_root=Path(name),
                    gap_input=forcing,max_step=1.25)

    def test_native_guards_preserve_conservation_and_signed_mechanics(self):
        with tempfile.TemporaryDirectory() as name:
            path=Path(name)/"trace.npz"
            active=np.array([[1.]*5,[2.]*5,[1.]*5])
            np.savez(path,field_active=active,field_inactive=3.-active)
            case=dict(model="wave_pinning",native_length=19.,cells=5)
            result=noise_campaign.physical_diagnostics(path,case)
            self.assertTrue(result["valid"])
            self.assertEqual(result["conservation"]["maximum_relative_drift"],0.)
            np.savez(path,**{"field_"+field:np.ones((3,5)) for field in
                ("rac","rho","mca_bound","mca_unbound")},field_displacement=-active)
            self.assertTrue(noise_campaign.physical_diagnostics(path,dict(case,model="debelly"))["valid"])

    def test_applied_input_distribution_reports_strong_draws_without_clipping(self):
        class RampedInput:
            cells=5;length=19.;knot=.5;kind="noise"
            def __call__(self,time,x,length):return np.full_like(x,time)
        forcing=RampedInput();schedule=make_train(4.,1.5,4,4.)
        result=noise_campaign.sampled_gap_distribution(forcing,schedule,1.)
        gaps=[(a,b) for kind,a,b in schedule.segments if kind=="gap"]
        expected=sum((b-a)*(a+b)/2 for a,b in gaps)/sum(b-a for a,b in gaps)
        self.assertAlmostEqual(result["time_averaged_spatial_mean"],expected)
        self.assertEqual(result["maximum"],gaps[-1][1])
        self.assertGreater(result["maximum_over_native_reference_peak"],1.)
        audited=noise_campaign.AppliedInputAudit(forcing,1.)
        for time in (1.,4.):np.testing.assert_array_equal(audited(time,np.ones(5),19.),time)
        self.assertEqual(audited.audit()["maximum"],4.)
        self.assertFalse(audited.audit()["clipping_or_realized_normalization"])


if __name__=="__main__":unittest.main()
