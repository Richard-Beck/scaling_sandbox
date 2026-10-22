"""Independent physical moment and periodic-region observer checks."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure4_native_audit_20261020 as audit


class NativeObservationAudit(unittest.TestCase):
    def test_activity_reference_magnitude_is_distinct_from_coherence(self):
        q=(np.arange(16)+.5)/16
        first=2+np.cos(2*np.pi*q)
        rows=audit.profile_observations(np.stack((first,2*first)),q,20.)
        np.testing.assert_allclose(rows["P_real"],[1.,2.],atol=1e-14)
        np.testing.assert_allclose(rows["q1"],[.25,.25],atol=1e-14)

    def test_actual_circumference_converts_diluted_spring_density_to_amount(self):
        q=(np.arange(16)+.5)/16;first=2+np.cos(2*np.pi*q)
        rows=audit.profile_observations(np.stack((first,first/2)),q,[1.,2.])
        np.testing.assert_allclose(rows["M0"],[2.,2.],atol=1e-14)
        np.testing.assert_allclose(rows["P_real"],[1.,1.],atol=1e-14)

    def test_periodic_component_crossing_zero_is_one_region(self):
        q=(np.arange(8)+.5)/8;profile=np.array([2.,0,0,0,0,0,0,2.])
        rows=audit.profile_observations(profile[None,:],q,8.)
        np.testing.assert_array_equal(rows["current_half_contrast_regions"],[1])
        with self.assertRaises(AssertionError):audit._compare([2.],[1.],"current_half_contrast_regions")

    def test_post_hash_scientific_file_change_invalidates_receipt(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"native.npz";path.write_bytes(b"retained native fields")
            snapshots={str(path):audit.file_snapshot(path)}
            audit.check_snapshots(snapshots)
            path.write_bytes(b"replacement native fields")
            with self.assertRaises(ValueError):audit.check_snapshots(snapshots)

    def test_exact_artifact_reuse_checks_actual_execution_identity_and_seed(self):
        import copy
        C=audit.C
        with tempfile.TemporaryDirectory() as folder:
            base=Path(folder);original=base/"original";selected=base/"selected"
            old=original/"workers/0000";fresh=selected/"workers/0361"
            old.mkdir(parents=True);fresh.mkdir(parents=True)
            job=dict(slot=0,full_stage_slot=361,model="spring_ring",seed=123,
                native_solver_restart_at_original_horizon=3200.,protocol=dict(horizon=12800.))
            manifest=dict(jobs=[job],source_sha256={},input_sha256={})
            C.write(original/"manifest.json",manifest);digest=C.sha(original/"manifest.json")
            C.write(original/"freeze.json",dict(manifest_sha256=digest))
            (old/"trace.npz").write_bytes(b"exact actual native arrays")
            artifacts={"trace.npz":C.sha(old/"trace.npz")}
            C.write(old/"result.json",dict(job=job,manifest_sha256=digest,artifact_sha256=artifacts))
            C.write(old/"status.json",dict(complete=True,output_sha256={"result.json":C.sha(old/"result.json")}))
            audit_root=original/"results/independent_native_audit_20261020";prior=audit_root/"workers/0000.json"
            C.write(prior,dict(passed=True,native_manifest_sha256=digest,native_result_sha256=C.sha(old/"result.json"),
                native_status_sha256=C.sha(old/"status.json"),hashed_native_file_snapshots={str(old/"trace.npz"):audit.file_snapshot(old/"trace.npz")}))
            C.write(audit_root/"manifest.json",dict(run=str(original),native_manifest_sha256=digest,source_sha256={},tasks=[0]))
            audit_digest=C.sha(audit_root/"manifest.json");C.write(audit_root/"freeze.json",dict(manifest_sha256=audit_digest))
            C.write(audit_root/"status.json",dict(full_native_run_covered=True,native_manifest_sha256=digest,audit_manifest_sha256=audit_digest,
                receipts=[dict(slot=0,path=str(prior),sha256=C.sha(prior),passed=True)]))
            (fresh/"trace.npz").symlink_to(old/"trace.npz")
            new=copy.deepcopy(job);new.pop("full_stage_slot");new["slot"]=361
            new["native_artifact_reuse"]=dict(run=str(original),slot=0,native_manifest_sha256=digest,
                result_sha256=C.sha(old/"result.json"),status_sha256=C.sha(old/"status.json"),artifact_sha256=artifacts,
                worker_audit=str(prior),worker_audit_sha256=C.sha(prior),source_audit_manifest=str(audit_root/"manifest.json"),
                source_audit_manifest_sha256=audit_digest,source_audit_status=str(audit_root/"status.json"),source_audit_status_sha256=C.sha(audit_root/"status.json"))
            result=audit.check_native_artifact_reuse(selected,new,manifest)
            self.assertEqual(result["actual_native_slot"],0)
            self.assertFalse(result["native_reintegration_performed"])
            new["seed"]=124
            with self.assertRaisesRegex(ValueError,"scientific condition"):audit.check_native_artifact_reuse(selected,new,manifest)
            new["seed"]=123;(old/"trace.npz").write_bytes(b"changed retained array")
            with self.assertRaisesRegex(ValueError,"changed"):audit.check_native_artifact_reuse(selected,new,manifest)
            (audit_root/"status.json").write_text((audit_root/"status.json").read_text()+" ")
            with self.assertRaisesRegex(ValueError,"full source-audit"):audit.check_native_artifact_reuse(selected,new,manifest)
            with self.assertRaises(ValueError):audit.check_snapshots(result["source_provenance_snapshots"])

    def test_selected_source_gate_binds_roles_paths_recipe_and_frozen_selection(self):
        import copy
        C=audit.C
        with tempfile.TemporaryDirectory() as temporary:
            base=Path(temporary);original=base/"original";selected=base/"selected"
            job=dict(slot=0,model="spring_ring",seed=123,retained_prefix=dict(horizon=3200.))
            C.write(original/"manifest.json",dict(jobs=[job],source_sha256={},input_sha256={}))
            digest=C.sha(original/"manifest.json");C.write(original/"freeze.json",dict(manifest_sha256=digest))
            status=original/"workers/0000/status.json";C.write(status,dict(complete=True,slurm_job_id="actual-original"))
            driver=selected/"frozen_sources/figure4_campaign_20261020.py";driver.parent.mkdir(parents=True);driver.write_text("frozen native adapter")
            recipe=base/"recipe.json";C.write(recipe,dict(original_manifest_sha256=digest,prefix_rule=C.ODE_PREFIX_RULE,native_adapter_sha256=C.sha(driver)))
            reuse=dict(run=str(original),slot=0,worker_audit="actual-audit.json",worker_audit_sha256="audit-hash")
            row=dict(original_slot=0,role="passing_original",actual=reuse,original_status=dict(path=str(status),sha256=C.sha(status)))
            native=dict(original_attempt=dict(run=str(original),manifest_sha256=digest),
                recovery_recipe=dict(path=str(recipe),sha256=C.sha(recipe)),jobs=[dict(job,native_artifact_reuse=reuse)],actual_native_selection=[row])
            def save(value,rows=None):
                C.write(selected/"manifest.json",value)
                C.write(selected/"results/complete_original_slot_selection_authority.json",dict(
                    native_manifest_sha256=C.sha(selected/"manifest.json"),all_original_slots=1,
                    original_manifest_sha256=digest,recipe_sha256=C.sha(recipe),response_outcomes_used_for_selection=False,
                    selection=value["actual_native_selection"] if rows is None else rows))
            save(native);self.assertEqual(audit.selected_native_source_gates(selected,native)[0]["role"],"passing_original")
            changed=copy.deepcopy(native);changed["actual_native_selection"][0]["role"]="qualified_old_cap_restart";save(changed)
            with self.assertRaisesRegex(ValueError,"passing original"):audit.selected_native_source_gates(selected,changed)
            changed=copy.deepcopy(native);changed["actual_native_selection"][0]["original_status"]["path"]=str(base/"other_status.json");save(changed)
            with self.assertRaisesRegex(ValueError,"path/hash"):audit.selected_native_source_gates(selected,changed)
            save(native,[])
            with self.assertRaisesRegex(ValueError,"differs from the frozen"):audit.selected_native_source_gates(selected,native)

    def test_e4_original_population_requires_a_retained_prefix(self):
        self.assertEqual(audit.e4_population(dict(new_gaussian_path=False,retained_prefix={})),"original_seed_continuation")
        self.assertEqual(audit.e4_population(dict(new_gaussian_path=True)),"additional_independent_seed")
        for invalid in ({},{'new_gaussian_path':False},{'new_gaussian_path':True,'retained_prefix':{}}):
            with self.assertRaises(ValueError):audit.e4_population(invalid)


if __name__=="__main__":unittest.main()
