"""Native recovery cannot select outcomes or rerun passing/live conditions."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure4_recovery_20261020 as recovery


class RecoveryAuthorityControls(unittest.TestCase):
    def frozen(self,root,jobs):
        C=recovery.C
        C.write(root/"manifest.json",dict(stage="E4",jobs=jobs,source_sha256={},input_sha256={}))
        C.write(root/"freeze.json",dict(manifest_sha256=C.sha(root/"manifest.json")))

    def test_actual_failed_native_prefix_is_required_before_a_retry(self):
        C=recovery.C
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/"original";out=Path(temporary)/"not_created"
            self.frozen(root,[dict(slot=0,retained_prefix=dict(horizon=3200.))])
            C.write(root/"workers/0000/status.json",dict(status="complete",native_integration_complete=True))
            with self.assertRaisesRegex(ValueError,"preserved native-prefix"):
                recovery.prepare_pilot(root,out,[0])
            self.assertFalse(out.exists())

    def test_full_selection_waits_for_every_original_slot_to_be_terminal(self):
        C=recovery.C
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/"original";out=Path(temporary)/"not_created"
            self.frozen(root,[dict(slot=i) for i in range(688)])
            C.write(root/"workers/0000/status.json",dict(status="running",complete=False))
            recipe=Path(temporary)/"recipe.json"
            C.write(recipe,dict(original_manifest_sha256=C.sha(root/"manifest.json"),
                native_adapter_sha256=C.sha(C.__file__),prefix_rule=C.ODE_PREFIX_RULE))
            with self.assertRaisesRegex(ValueError,"remains live"):
                recovery.prepare_selection(root,out,[],recipe)
            self.assertFalse(out.exists())

    def test_platform_qualification_ignores_only_the_unrelated_adapter(self):
        C=recovery.C
        with tempfile.TemporaryDirectory() as temporary:
            alternate=Path(temporary)/"spring_only.py"
            alternate.write_text(Path(C.__file__).read_text().replace("class OriginalCapBDF:","class UnrelatedChemicalAdapter:"))
            self.assertEqual(recovery.platform_adapter_sha(alternate,"spring"),recovery.platform_adapter_sha(C.__file__,"spring"))
            with self.assertRaises(StopIteration):recovery.platform_adapter_sha(alternate,"chemical")

    def test_frozen_recovery_cli_has_its_canonical_import_dependencies(self):
        from unittest.mock import patch
        C=recovery.C
        with tempfile.TemporaryDirectory() as temporary:
            base=Path(temporary);original=base/"original";driver=original/"frozen_sources/figure4_campaign_20261020.py"
            driver.parent.mkdir(parents=True);driver.write_text("original native driver")
            C.write(original/"manifest.json",dict(jobs=[],input_sha256={},source_sha256={str(driver.relative_to(original)):C.sha(driver)}))
            C.write(original/"freeze.json",dict(manifest_sha256=C.sha(original/"manifest.json")))
            ms=base/"manuscript";out=ms/"data/runs/recovery"
            with patch.object(C,"MS",ms):manifest=recovery.freeze_base(original,out)
            for name in ("signature_figure4_campaign_20261020.py","signature_figure4_native_audit_20261020.py"):
                relative="frozen_sources/"+name
                self.assertTrue((out/relative).is_file())
                self.assertEqual(C.sha(out/relative),manifest["source_sha256"][relative])


if __name__=="__main__":unittest.main()
