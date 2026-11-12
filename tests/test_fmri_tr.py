"""Verify TR selection, logging and workflow wiring without running FSL."""

import argparse
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from common.fmri_metadata import positive_tr, read_fmri_metadata, resolve_tr
from common import script_logging


class TRSelectionTests(unittest.TestCase):
    def select(self, metadata, cli=None):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = resolve_tr(metadata, cli)
        return result, output.getvalue()

    def test_cli_wins_without_inspecting_json_tr(self):
        for metadata in ({"RepetitionTime": 3.0}, {"RepetitionTime": "invalid"}, {}):
            value, log = self.select(metadata, 2.0)
            self.assertEqual(value, 2.0)
            self.assertEqual(log, "Using TR = 2.0 s (source: cli-option)\n")

    def test_json_when_cli_is_absent(self):
        value, log = self.select({"RepetitionTime": 2.5})
        self.assertEqual(value, 2.5)
        self.assertIn("Using TR = 2.5 s (source: JSON)", log)

    def test_missing_or_invalid_json_tr_uses_fallback(self):
        for metadata in ({}, *({"RepetitionTime": value} for value in (
            None, True, 0, -1, "2.0", [], {}, float("nan"), float("inf"),
        ))):
            with self.subTest(metadata=metadata):
                value, log = self.select(metadata)
                self.assertEqual(value, 1.42)
                self.assertIn("Using TR = 1.42 s (source: fallback)", log)

    def test_legacy_json_converts_milliseconds_only_without_cli(self):
        metadata = {"RepetitionTime": 1420, "ObjOrderList": [1],
                    "n_slices": 1, "costum_timings": [0]}
        value, log = self.select(metadata)
        self.assertEqual(value, 1.42)
        self.assertIn("milliseconds to seconds", log)
        self.assertEqual(self.select(metadata, 2.0)[0], 2.0)
        self.assertEqual(self.select({"RepetitionTime": 1420})[0], 1420)

    def test_cli_validation(self):
        for value in ("nan", "inf", "-inf", "0", "-1", "bad", True):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                positive_tr(value)
        self.assertEqual(positive_tr("2.0"), 2.0)

    def test_sidecar_reading_and_failures(self):
        with tempfile.TemporaryDirectory() as folder:
            for extension in (".nii", ".nii.gz"):
                image = Path(folder) / ("scan" + extension)
                sidecar = Path(folder) / "scan.json"
                sidecar.write_text('{"RepetitionTime": 2.0}')
                self.assertEqual(read_fmri_metadata(image), {"RepetitionTime": 2.0})
                for contents in ("{broken", "[]", "null"):
                    sidecar.write_text(contents)
                    with contextlib.redirect_stdout(io.StringIO()):
                        self.assertEqual(read_fmri_metadata(image), {})
                sidecar.unlink()
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(read_fmri_metadata(image), {})

    def test_all_sources_are_written_to_processing_log(self):
        log_files = []
        with tempfile.TemporaryDirectory() as folder:
            try:
                with patch.object(script_logging, "_LOG_FILES", log_files), patch.dict(os.environ, {
                    "AIDAMRI_DISABLE_SCRIPT_LOG": "0", "AIDAMRI_DISABLE_PROCESS_LOG": "0",
                }), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    script_logging.setup_script_logging(folder, "process.log")
                    resolve_tr({}, 2.0)
                    resolve_tr({"RepetitionTime": 2.5})
                    resolve_tr({})
                log = (Path(folder) / "process.log").read_text()
                for value, source in ((2.0, "cli-option"), (2.5, "JSON"), (1.42, "fallback")):
                    self.assertIn(f"Using TR = {value} s (source: {source})", log)
            finally:
                for log_file in log_files:
                    log_file.close()


class ProcessingTRTests(unittest.TestCase):
    def setUp(self):
        # Execute real CLI/TR code, replacing only imaging dependencies.
        dependencies = {name: MagicMock() for name in (
            "nipype", "nipype.interfaces", "nipype.interfaces.fsl", "nibabel",
            "numpy", "regress", "getSingleRegTable", "cv2", "create_seed_rois",
            "fsl_mean_ts", "common.bet",
        )}
        spec = importlib.util.spec_from_file_location(
            "process_fmri_tr_test", BIN_DIR / "3.3_fMRIActivity" / "process_fMRI.py"
        )
        self.module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, dependencies):
            spec.loader.exec_module(self.module)

    def run_cli(self, metadata, options):
        with tempfile.TemporaryDirectory() as folder, contextlib.ExitStack() as stack:
            image = Path(folder) / "scan.nii.gz"
            image.touch()
            if metadata is not None:
                contents = metadata if isinstance(metadata, str) else json.dumps(metadata)
                (Path(folder) / "scan.json").write_text(contents)
            for name in ("start_output_tracking", "setup_script_logging", "startProcess",
                         "copyAtlasOfData", "create_txt_file", "delete_txt_file"):
                stack.enter_context(patch.object(self.module, name))
            self.module.regress.startRegression.return_value = ("rgr", "srgr", "sfrgr")
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            try:
                self.module.main(["-i", str(image), *options])
            except SystemExit:
                self.module.startProcess.assert_not_called()
                raise
            self.module.startProcess.assert_called_once()
            return self.module.regress.startRegression.call_args.args

    def test_cli_aliases_override_json_in_regression(self):
        for alias in ("--tr", "--TR", "-t"):
            with self.subTest(alias=alias):
                args = self.run_cli({"RepetitionTime": 3.0}, [alias, "2.0"])
                self.assertEqual(args[3], 2.0)
                self.assertFalse(args[4])

    def test_json_and_fallback_without_stc(self):
        for metadata, expected in (({"RepetitionTime": 2.5}, 2.5), (None, 1.42), ("{broken", 1.42)):
            with self.subTest(metadata=metadata):
                self.assertEqual(self.run_cli(metadata, [])[3], expected)

    def test_stc_preserves_cli_json_and_fallback(self):
        metadata = {"ObjOrderList": [1], "costum_timings": [0], "n_slices": 1}
        for value, options, expected in (
            (1420, ["--tr", "2.0"], 2.0), (2500, [], 2.5), (None, [], 1.42),
        ):
            with self.subTest(value=value, options=options):
                args = self.run_cli(dict(metadata, RepetitionTime=value), ["-stc", *options])
                self.assertEqual(args[3], expected)
                self.assertTrue(args[4])

    def test_invalid_cli_and_missing_stc_metadata_fail_before_processing(self):
        for options in (["--tr", "nan"], ["--tr", "0"], ["-stc"]):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as exc:
                    self.run_cli(None, options)
                self.assertEqual(exc.exception.code, 2)

    def test_unusable_metadata_skips_physio_lookup(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(self.module.copyRawPhysioData("/tmp/scan.nii.gz", "/tmp", metadata={}), [])


class BatchTRTests(unittest.TestCase):
    def test_batch_forwards_explicit_tr_only(self):
        with patch.dict(sys.modules, {"tqdm": MagicMock()}):
            import batchProc
        original_cwd = os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as folder:
                (Path(folder) / "test_EPI.nii.gz").touch()
                for cfg in ({}, {"func_tr": 2.0}):
                    with patch.object(batchProc, "run_subprocess", return_value=0) as run:
                        batchProc.executeScripts(folder, "func", "process", cfg)
                    command = run.call_args.args[0]
                    self.assertEqual("--tr 2.0" in command, "func_tr" in cfg)
        finally:
            os.chdir(original_cwd)


if __name__ == "__main__":
    unittest.main()
