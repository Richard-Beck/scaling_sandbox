"""Prospective design, observation clocks and retained-input prefix controls."""
from pathlib import Path
import collections
import copy
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure4_campaign_20261020 as campaign
import signature_figure4_basal_20261020 as intermediate_basal


class Figure4CampaignControls(unittest.TestCase):
    def test_e1_requires_fixed_strengths_with_complete_existing_companions(self):
        with self.assertRaises(ValueError):campaign.grid("E1")
        rows=campaign.grid("E1",strengths={model:[.01,.1] for model in campaign.MODELS})
        self.assertEqual(len(rows),864)
        self.assertEqual({row["sigma"] for row in rows},{.25,.5,2.})
        self.assertEqual(len({row["seed"] for row in rows}),288)
        grouped=collections.defaultdict(set)
        for row in rows:grouped[(row["model"],row["size"],row["fraction"],row["replicate"])].add(row["seed"])
        self.assertTrue(all(len(seeds)==1 for seeds in grouped.values()))
        for model in campaign.MODELS:
            self.assertEqual({row["size"] for row in rows if row["model"]==model},set(campaign.sizes(model)))
        pilot=campaign.grid("E1_pilot",strengths={model:.1 for model in campaign.MODELS})
        self.assertEqual(len(pilot),8)
        self.assertEqual({row["sigma"] for row in pilot},{1.,2.})

    def test_e2_and_e3_include_the_required_distinct_controls(self):
        e2=campaign.grid("E2");e3=campaign.grid("E3")
        self.assertEqual(collections.Counter(r["kind"] for r in e2),{"noise":432,"mean":36})
        self.assertEqual(collections.Counter(r["kind"] for r in e3),{"noise":480,"mean":40,"basal":8})
        self.assertEqual({r["fraction"] for r in e3 if r["kind"]=="noise"},{.001,.003,.01,.03,.1})
        self.assertEqual(len({r["seed"] for r in e3 if r["kind"]=="noise"}),480)

    def test_e1_companion_batches_are_disjoint_and_keep_all_sigma_seed_cells(self):
        choices={model:[.03,.1] for model in campaign.MODELS};choices["goryachev"]=[.02,.03]
        full=campaign.grid("E1",strengths=choices)
        keys=sorted({(j["model"],j["size"],j["fraction"]) for j in full})
        late=[dict(model=m,size=s,fraction=f) for m,s,f in keys if m=="goryachev" and f==.02 and s<60]
        ready=[dict(model=m,size=s,fraction=f) for m,s,f in keys if dict(model=m,size=s,fraction=f) not in late]
        first,deferred=campaign.select_conditions(full,ready);second,other=campaign.select_conditions(full,late)
        self.assertEqual((len(first),len(second),len(deferred)),(792,72,72))
        self.assertEqual({j["full_stage_slot"] for j in first}|{j["full_stage_slot"] for j in second},set(range(864)))
        self.assertFalse({j["full_stage_slot"] for j in first}&{j["full_stage_slot"] for j in second})
        self.assertEqual(set(collections.Counter((j["model"],j["size"],j["fraction"]) for j in first).values()),{36})
        with self.assertRaises(ValueError):campaign.select_conditions(full,late+late)

    def test_e4_uses_every_trigger_with_matching_extended_controls(self):
        trigger=dict(model="wave_pinning",size=60.,fraction=.001,heading_cosine=.9,longest_usable_lag=8000.)
        rows=campaign.grid("E4",triggers=[trigger])
        self.assertEqual(collections.Counter(r["kind"] for r in rows),{"noise":12,"mean":1,"basal":1})
        for value in (.8,float("nan")):
            with self.assertRaises(ValueError):campaign.grid("E4",triggers=[dict(trigger,heading_cosine=value)])
        with self.assertRaises(ValueError):campaign.grid("E4",triggers=[trigger,trigger])

    def test_dense_profiles_and_extended_windows_preserve_native_solver_and_noise(self):
        baseline=campaign.read(campaign.SOURCE/"manifest.json")["native_protocols"]
        original=copy.deepcopy(baseline)
        for platform in ("chemical","spring"):
            e2=campaign.protocol(platform,"E2",baseline);e4=campaign.protocol(platform,"E4",baseline)
            self.assertEqual(e2["profile_interval"],e2["scalar_interval"])
            self.assertEqual(e4["profile_interval"],e4["scalar_interval"])
            self.assertEqual(e4["horizon"],4*baseline[platform]["horizon"])
            for key in ("solver","noise"):
                self.assertEqual(e2[key],baseline[platform][key]);self.assertEqual(e4[key],baseline[platform][key])
        self.assertEqual(baseline,original)

    def test_spring_checkpoint_restart_replays_the_native_original_endpoint(self):
        from scipy.integrate import solve_ivp
        solver=dict(rtol=1e-8,atol=1e-10,max_step=.125)
        times=np.arange(41)*.25
        rhs=lambda t,y:np.array([-y[0]+.2*np.sin(t)])
        jac=lambda t,y:np.array([[-1.]])
        initial=np.array([1.])
        old=solve_ivp(rhs,(0.,4.),initial,t_eval=times[times<=4.],method="BDF",jac=jac,**solver)
        replay=campaign.spring_native_solve(rhs,jac,initial,times,solver,4.)
        self.assertTrue(replay.success)
        np.testing.assert_array_equal(replay.t,times)
        np.testing.assert_array_equal(replay.y[:,:len(old.t)],old.y)
        # Global time is retained in the nonautonomous RHS after the restart.
        reference=solve_ivp(rhs,(0.,10.),initial,t_eval=times,method="BDF",jac=jac,**solver)
        np.testing.assert_allclose(replay.y,reference.y,rtol=2e-7,atol=2e-9)
        direct=campaign.spring_native_solve(rhs,jac,initial,times,solver)
        np.testing.assert_array_equal(direct.y,reference.y)
        with self.assertRaisesRegex(ValueError,"exact internal"):
            campaign.spring_native_solve(rhs,jac,initial,times,solver,4.1)

    def test_chemical_bdf_restart_retains_the_native_observation_closure(self):
        from scipy.integrate import BDF
        kwargs=dict(rtol=1e-8,atol=1e-10,max_step=.125)
        rhs=lambda t,y:-y+.2*np.sin(t)
        old=BDF(rhs,0.,np.array([1.]),4.,**kwargs)
        original_steps=[]
        while old.status=="running":old.step();original_steps.append((old.t,old.y.copy()))
        repaired=campaign.OriginalCapBDF(BDF,4.,rhs,0.,np.array([1.]),10.,**kwargs)
        replay=[];times=[]
        while repaired.status=="running":
            repaired.step();times.append(repaired.t)
            if repaired.t<=4.:replay.append((repaired.t,repaired.y.copy()))
        np.testing.assert_array_equal([v[0] for v in replay],[v[0] for v in original_steps])
        np.testing.assert_array_equal([v[1] for v in replay],[v[1] for v in original_steps])
        self.assertEqual(times.count(4.),1)
        self.assertEqual(repaired.t,10.)
        self.assertGreater(repaired.nfev,old.nfev)

    def test_e0_interim_gate_precedes_any_new_payload(self):
        with tempfile.TemporaryDirectory() as temporary:
            gate=Path(temporary)/"gate.json";campaign.write(gate,dict(e0_interim_delivered=False))
            out=campaign.MS/"data/runs/figure4_test_must_not_be_created"
            with self.assertRaisesRegex(ValueError,"Deliver E0"):
                campaign.prepare(out,"E2",gate)
            self.assertFalse(out.exists())

    def test_intermediate_spring_basal_waits_for_full_e1_e2_native_gates(self):
        with tempfile.TemporaryDirectory() as temporary:
            gate=Path(temporary)/"gate.json";campaign.write(gate,dict(e1_e2_complete=False))
            out=campaign.MS/"data/runs/figure4_intermediate_basal_must_not_be_created"
            with self.assertRaisesRegex(ValueError,"All E1/E2"):
                intermediate_basal.prepare(out,gate)
            self.assertFalse(out.exists())

    def test_forcing_preflight_requires_every_disjoint_worker_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            out=Path(temporary);jobs=[dict(slot=i,model="wave_pinning",size=20.,fraction=.01,sigma=1.,seed=i,input_sha256="basal") for i in range(2)]
            campaign.write(out/"manifest.json",dict(jobs=jobs,source_sha256={},input_sha256={}))
            digest=campaign.sha(out/"manifest.json");campaign.write(out/"freeze.json",dict(manifest_sha256=digest))
            for shard in range(2):
                proof=out/"results/forcing_preflight"/f"{shard:04d}.json";campaign.write(proof,jobs[shard])
                campaign.write(out/"results/forcing_preflight_shards"/f"shard_{shard:03d}.json",dict(passed=True,manifest_sha256=digest,
                    selected_tasks=[shard],shards=2,shard=shard,workers=[dict(slot=shard,receipt=str(proof),receipt_sha256=campaign.sha(proof))]))
                result=campaign.collect_forcing_preflight(out,2)
                self.assertEqual(result["passed"],shard==1)
            proof.write_text("changed forcing identity")
            with self.assertRaisesRegex(ValueError,"receipt changed"):campaign.collect_forcing_preflight(out,2)

    def test_extended_native_reference_is_hashed_before_the_long_solve(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            for name in ("trace.npz","primary.npz","final_state.npz"):(folder/name).write_bytes(b"native fields")
            status=folder/"status.json";campaign.write(status,dict(output_sha256={"final_state.npz":campaign.sha(folder/"final_state.npz")}))
            old={key:str(folder/name) for key,name in (("trace","trace.npz"),("profile","primary.npz"))}
            old.update({key+"_sha256":campaign.sha(value) for key,value in old.copy().items()})
            job=dict(retained_prefix=dict(original_case=old,original_status=str(status),original_status_sha256=campaign.sha(status)))
            self.assertEqual(len(campaign.check_retained_reference_artifacts(job)),3)
            (folder/"final_state.npz").write_bytes(b"changed native coupled checkpoint")
            with self.assertRaisesRegex(ValueError,"reference artifact changed"):campaign.check_retained_reference_artifacts(job)

    def test_held_pilot_needs_a_matching_concrete_delivery_release(self):
        with tempfile.TemporaryDirectory() as temporary:
            out=Path(temporary);campaign.write(out/"manifest.json",dict(stage="E1_pilot"))
            manifest=dict(native_production_held_for_E0_interim=True)
            with self.assertRaisesRegex(RuntimeError,"held until concrete E0"):
                campaign.check_launch_gate(out,manifest)
            interim=out/"delivered_interim.json";campaign.write(interim,dict(delivered=True))
            release=dict(e0_interim_delivered=True,manifest_sha256=campaign.sha(out/"manifest.json"),
                interim_artifacts={str(interim):campaign.sha(interim)})
            campaign.write(out/"results/e0_interim_launch_release.json",release)
            checked=campaign.check_launch_gate(out,manifest)
            self.assertEqual(checked["e0_delivery_release_sha256"],campaign.sha(out/"results/e0_interim_launch_release.json"))
            campaign.write(interim,dict(delivered=False))
            with self.assertRaisesRegex(ValueError,"artifact changed"):
                campaign.check_launch_gate(out,manifest)

    def test_execution_collection_retains_and_verifies_native_pattern_histories(self):
        with tempfile.TemporaryDirectory() as temporary:
            out=Path(temporary);folder=out/"workers/0000";folder.mkdir(parents=True)
            job=dict(model="wave_pinning",size=20.,fraction=.02,kind="noise",seed=7,sigma=1.,
                independent_noise_seed=True,new_gaussian_path=True,input_sha256="basal",slot=0,
                protocol=dict(horizon=32000.,scalar_interval=5.,profile_interval=5.))
            manifest=dict(stage="E2",jobs=[job],source_sha256={},input_sha256={},counts=dict(noise=1),parent_manifest_sha256="parent")
            campaign.write(out/"manifest.json",manifest);digest=campaign.sha(out/"manifest.json")
            campaign.write(out/"freeze.json",dict(manifest_sha256=digest))
            for name in ("trace.npz","primary.npz","local_patterns.jsonl"):(folder/name).write_text("scientific artifact\n")
            artifacts={p.name:campaign.sha(p) for p in folder.iterdir()}
            campaign.write(folder/"result.json",dict(job=job,manifest_sha256=digest,artifact_sha256=artifacts,
                wall_seconds=1.,peak_resident_memory_KiB=100))
            campaign.write(folder/"status.json",dict(complete=True,output_sha256={"result.json":campaign.sha(folder/"result.json")}))
            self.assertTrue(campaign.collect(out)["complete"])
            (folder/"local_patterns.jsonl").write_text("changed native pattern history\n")
            with self.assertRaisesRegex(ValueError,"artifact changed"):campaign.collect(out)

    def test_horizon_extension_retains_native_gaussian_prefixes(self):
        native=campaign.native_module(campaign.SOURCE,"chemical")
        grid=native.ring.RingGrid(20.,40)
        short=native.NoiseInput(grid,20.,883,.01,sigma=1.)
        extended=native.NoiseInput(grid,80.,883,.01,sigma=1.)
        np.testing.assert_array_equal(short.knots,extended.knots[:len(short.knots)])
        high_variance=native.NoiseInput(grid,20.,883,.01,sigma=2.)
        np.testing.assert_array_equal(short.knots,high_variance.knots)
        self.assertAlmostEqual(high_variance.audit()["theoretical_input_variance"],.0001*(np.exp(4)-1))
        spring=campaign.native_module(campaign.SOURCE,"spring")
        first=spring.MaterialNoise(20,1.0414,10.,884,.01,knot=.125,coupling=.125)
        longer=spring.MaterialNoise(20,1.0414,40.,884,.01,knot=.125,coupling=.125)
        np.testing.assert_array_equal(first.knots,longer.knots[:len(first.knots)])

    def test_realized_input_statistics_and_query_audit_do_not_normalize_the_path(self):
        native=campaign.native_module(campaign.SOURCE,"chemical")
        field=native.NoiseInput(native.ring.RingGrid(20.,40),20.,883,.01,sigma=2.)
        job=dict(sigma=2.,kind="noise",mean=.01,fraction=.1)
        result=campaign.input_statistics(field,job,dict(horizon=20.,noise=dict(knot=5.)))
        self.assertAlmostEqual(result["theoretical_input_mean"],.01)
        self.assertNotEqual(result["realized_input_windows"]["postburn"]["mean"],.01)
        self.assertIn("0.999",result["realized_input_windows"]["postburn"]["quantiles"])
        queried=campaign.QueriedInput(field)
        expected=field(2.5).copy()
        np.testing.assert_array_equal(queried(2.5),expected)
        self.assertEqual(queried.report(.1)["maximum"],float(expected.max()))
        self.assertFalse(queried.report(.1)["clipping_or_pathwise_normalization"])

    def test_sigma1_equivalence_checks_the_full_history_not_only_aggregates(self):
        with tempfile.TemporaryDirectory() as temporary:
            source=Path(temporary)/"source";fresh=Path(temporary)/"fresh"
            source.mkdir();fresh.mkdir()
            for folder in (source,fresh):
                np.savez(folder/"trace.npz",time=np.arange(4.),R=[1.,2.,1.,1.])
                np.savez(folder/"primary.npz",time=[0.,3.],primary=np.ones((2,3)))
                np.savez(folder/"final_state.npz",time=3.,state=np.ones((2,3)))
            job=dict(sigma=1.,paired_sigma1_source=dict(original_case=dict(
                trace=str(source/"trace.npz"),profile=str(source/"primary.npz"))))
            campaign.compare_sigma1_evolution(fresh,job)
            self.assertTrue(campaign.read(fresh/"sigma1_evolution_equivalence.json")["passed"])
            np.savez(fresh/"primary.npz",time=np.arange(4.),primary=np.ones((4,3)))
            campaign.compare_sigma1_evolution(fresh,job)
            self.assertTrue(campaign.read(fresh/"sigma1_evolution_equivalence.json")["passed"])
            np.savez(fresh/"trace.npz",time=np.arange(4.),R=[1.,1.,2.,1.])
            with self.assertRaises(AssertionError):campaign.compare_sigma1_evolution(fresh,job)

    def test_long_native_prefix_checks_interior_profiles_and_coupled_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            source=Path(temporary)/"source";fresh=Path(temporary)/"fresh";source.mkdir();fresh.mkdir()
            np.savez(source/"trace.npz",time=np.arange(4.),P_real=np.ones(4))
            np.savez(source/"primary.npz",time=[0.,1.,3.],primary=np.ones((3,3)))
            np.savez(source/"final_state.npz",state=np.ones((2,3)),time=3.)
            campaign.write(source/"status.json",dict(output_sha256={"final_state.npz":campaign.sha(source/"final_state.npz")}))
            np.savez(fresh/"trace.npz",time=np.arange(7.),P_real=np.ones(7))
            fields=np.ones((7,3));fields[4:]=2.
            np.savez(fresh/"primary.npz",time=np.arange(7.),primary=fields)
            np.savez(fresh/"full_states.npz",time=[0.,3.,6.],states=np.ones((3,2,3)))
            old=dict(trace=str(source/"trace.npz"),trace_sha256=campaign.sha(source/"trace.npz"),
                profile=str(source/"primary.npz"),profile_sha256=campaign.sha(source/"primary.npz"))
            job=dict(retained_prefix=dict(original_case=old,horizon=3.,original_status=str(source/"status.json"),
                original_status_sha256=campaign.sha(source/"status.json")))
            campaign.check_retained_evolution_prefix(fresh,job)
            self.assertTrue(campaign.read(fresh/"retained_native_prefix_equivalence.json")["passed"])
            fields[1]=2.;np.savez(fresh/"primary.npz",time=np.arange(7.),primary=fields)
            with self.assertRaises(AssertionError):campaign.check_retained_evolution_prefix(fresh,job)
            fields[1]=1.;np.savez(fresh/"primary.npz",time=np.arange(7.),primary=fields)
            states=np.ones((3,2,3));states[1,1,0]=2.
            np.savez(fresh/"full_states.npz",time=[0.,3.,6.],states=states)
            with self.assertRaises(AssertionError):campaign.check_retained_evolution_prefix(fresh,job)


if __name__=="__main__":unittest.main()
