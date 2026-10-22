"""Controls for explicitly bounded extensions of a native calibration assay."""
from pathlib import Path
import sys
import unittest
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_campaign_20261020 as campaign


class CalibrationHorizonControls(unittest.TestCase):
    def test_extension_retains_native_protocol_and_finite_scan(self):
        baseline = campaign.config()
        extended = campaign.config(calibration_horizon_overrides={"holmes_model3": 100000.})
        self.assertEqual(extended["calibration_horizon"]["holmes_model3"], 100000.)
        self.assertEqual(extended["calibration"]["scans"]["holmes_model3"],
                         [1., 10., 100., 1000., 3000., 6000., 10000., 30000., 60000., 100000.])
        for key in ("recovery", "native_peak", "cue", "post_withdrawal_replacement_threshold", "field_tail"):
            self.assertEqual(extended[key], baseline[key])
        self.assertEqual(campaign.config(), baseline)
        self.assertEqual(extended["calibration_horizon"]["wave_pinning"],
                         baseline["calibration_horizon"]["wave_pinning"])

    def test_invalid_extension_is_rejected(self):
        for model, horizon in (("holmes_model3", 9999.), ("spring", float("inf")),
                               ("unknown", 100000.), ("debelly", float("nan"))):
            with self.subTest(model=model, horizon=horizon), self.assertRaises(ValueError):
                campaign.config(calibration_horizon_overrides={model: horizon})

    def test_cross_node_nfs_mount_number_preserves_hashed_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            campaign.write(root/"manifest.json", {"cases": [{"task": 0}]})
            scientific=root/"workers/0000/trace.npz"
            scientific.parent.mkdir(parents=True)
            scientific.write_bytes(b"frozen native fields")
            campaign.write(scientific.parent/"result.json", {"complete": True})
            paths=[scientific,scientific.parent/"result.json"]
            snapshots={}
            for path in paths:
                fingerprint=campaign.stat_fingerprint(path)
                fingerprint["device"] += 1
                snapshots[str(path.relative_to(root))] = dict(
                    fingerprint=fingerprint, sha256=campaign.sha(path))
            shards=root/"results/artifact_shards"
            campaign.write(shards/"shard_000.json", dict(
                index=0,shards=1,complete=True,
                manifest_sha256=campaign.sha(root/"manifest.json"),
                workers=[dict(task=0,artifact_snapshots=snapshots)]))
            self.assertTrue(campaign.verified_shard_receipts(root,shards)[0]["hash_proof_reused_from_shard"])
            scientific.write_bytes(b"changed native fields")
            with self.assertRaisesRegex(ValueError,"stat identity changed"):
                campaign.verified_shard_receipts(root,shards)


if __name__ == "__main__":
    unittest.main()
