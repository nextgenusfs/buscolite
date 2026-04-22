"""
Tests for the busco module.
"""

import os
import subprocess
import tempfile
from unittest.mock import patch

import pytest  # noqa: F401 - Needed for pytest fixtures

from buscolite.busco import (
    check_lineage,
    load_config,
    load_cutoffs,
    predict_and_validate,
    process_busco_jobs,
)


def test_load_config(mock_busco_lineage):
    """Test loading configuration from dataset.cfg."""
    config = load_config(mock_busco_lineage)

    assert isinstance(config, dict)
    assert config["name"] == "test_lineage"
    assert config["species"] == "test_species"
    assert config["domain"] == "test_domain"


def test_load_cutoffs(mock_busco_lineage):
    """Test loading cutoffs from scores_cutoff and lengths_cutoff."""
    cutoffs = load_cutoffs(mock_busco_lineage)

    assert isinstance(cutoffs, dict)
    assert "busco1" in cutoffs
    assert "busco2" in cutoffs

    # Check busco1 values
    assert cutoffs["busco1"]["score"] == 100.0
    assert cutoffs["busco1"]["sigma"] == 1.5
    assert cutoffs["busco1"]["length"] == 300

    # Check busco2 values
    assert cutoffs["busco2"]["score"] == 200.0
    assert cutoffs["busco2"]["sigma"] == 1.0  # Should be 1 because sigma was 0.0
    assert cutoffs["busco2"]["length"] == 400


def test_check_lineage_valid(mock_busco_lineage):
    """Test checking a valid lineage directory."""
    valid, message = check_lineage(mock_busco_lineage)

    assert valid is True
    assert message == ""


def test_check_lineage_invalid_dir(temp_dir):
    """Test checking an invalid directory."""
    invalid_dir = os.path.join(temp_dir, "nonexistent")
    valid, message = check_lineage(invalid_dir)

    assert valid is False
    assert "is not a directory" in message


def test_check_lineage_missing_dir(mock_busco_lineage):
    """Test checking a lineage with a missing directory."""
    # Create a new lineage directory without the hmms directory
    import shutil
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy everything except hmms directory
        new_lineage = os.path.join(tmpdir, "lineage_no_hmms")
        os.makedirs(new_lineage)

        # Copy files
        for f in [
            "dataset.cfg",
            "scores_cutoff",
            "lengths_cutoff",
            "ancestral",
            "ancestral_variants",
        ]:
            shutil.copy(os.path.join(mock_busco_lineage, f), os.path.join(new_lineage, f))

        # Create prfl directory
        os.makedirs(os.path.join(new_lineage, "prfl"))

        # Test the lineage without hmms directory
        valid, message = check_lineage(new_lineage)

        assert valid is False
        assert "hmms directory was not found" in message


def test_check_lineage_missing_file(mock_busco_lineage):
    """Test checking a lineage with a missing file."""
    # Remove the scores_cutoff file
    os.remove(os.path.join(mock_busco_lineage, "scores_cutoff"))

    valid, message = check_lineage(mock_busco_lineage)

    assert valid is False
    assert "scores_cutoff file is missing" in message


def _minimal_fadict():
    return {"contig1": "ACGT" * 25}


def test_predict_and_validate_saves_debug_artifacts_on_failure(tmp_path):
    """On augustus CalledProcessError, debug artifacts are persisted and
    the exception carries a .buscolite_debug_dir path."""
    prfl = str(tmp_path / "busco42.prfl")
    open(prfl, "w").close()

    cmd = ["augustus", "--fake", "arg"]
    err = subprocess.CalledProcessError(-4, cmd, output=None, stderr="illegal instruction\n")

    with patch("buscolite.busco.proteinprofile", side_effect=err):
        with pytest.raises(subprocess.CalledProcessError) as excinfo:
            predict_and_validate(
                _minimal_fadict(),
                "contig1",
                prfl,
                cutoffs={"busco42": {"score": 100.0}},
                species="human",
                start=0,
                end=20,
                strand="+",
                configpath="/nope",
                blast_score=1.0,
            )

    debug_dir = getattr(excinfo.value, "buscolite_debug_dir", None)
    assert debug_dir is not None and os.path.isdir(debug_dir)
    assert os.path.isfile(os.path.join(debug_dir, "input.fasta"))
    with open(os.path.join(debug_dir, "input.fasta")) as f:
        body = f.read()
    assert ">contig1" in body
    with open(os.path.join(debug_dir, "cmd.txt")) as f:
        assert f.read().splitlines() == cmd
    with open(os.path.join(debug_dir, "returncode.txt")) as f:
        assert f.read().strip() == "-4"
    with open(os.path.join(debug_dir, "stderr.txt")) as f:
        assert "illegal instruction" in f.read()


def test_predict_and_validate_sigsegv_saves_artifacts(tmp_path):
    """Any negative returncode (e.g. -11 SIGSEGV) is treated as a signal
    death: artifacts are persisted and the exception is re-raised."""
    prfl = str(tmp_path / "busco42.prfl")
    open(prfl, "w").close()

    cmd = ["augustus", "--fake", "arg"]
    err = subprocess.CalledProcessError(-11, cmd, output=None, stderr="segfault\n")

    with patch("buscolite.busco.proteinprofile", side_effect=err):
        with pytest.raises(subprocess.CalledProcessError) as excinfo:
            predict_and_validate(
                _minimal_fadict(),
                "contig1",
                prfl,
                cutoffs={"busco42": {"score": 100.0}},
                species="human",
                start=0,
                end=20,
                strand="+",
                configpath="/nope",
                blast_score=1.0,
            )

    debug_dir = getattr(excinfo.value, "buscolite_debug_dir", None)
    assert debug_dir is not None and os.path.isdir(debug_dir)
    with open(os.path.join(debug_dir, "returncode.txt")) as f:
        assert f.read().strip() == "-11"


def test_predict_and_validate_benign_nonzero_exit_is_nohit(tmp_path):
    """rc>0 CalledProcessError (augustus' normal "No feasible path found in
    HMM" exit) returns False silently — no debug dir, no re-raise."""
    prfl = str(tmp_path / "busco42.prfl")
    open(prfl, "w").close()

    debug_root = os.path.join(tempfile.gettempdir(), "buscolite_debug")
    before = set(os.listdir(debug_root)) if os.path.isdir(debug_root) else set()

    cmd = ["augustus", "--fake", "arg"]
    err = subprocess.CalledProcessError(
        1, cmd, output=None, stderr="augustus: ERROR\n\tNo feasible path found in HMM\n"
    )

    with patch("buscolite.busco.proteinprofile", side_effect=err):
        result = predict_and_validate(
            _minimal_fadict(),
            "contig1",
            prfl,
            cutoffs={"busco42": {"score": 100.0}},
            species="human",
            start=0,
            end=20,
            strand="+",
            configpath="/nope",
            blast_score=1.0,
        )

    assert result is False
    after = set(os.listdir(debug_root)) if os.path.isdir(debug_root) else set()
    assert after == before, "benign rc>0 must not leave debug dirs behind"


def test_predict_and_validate_cleans_up_on_success(tmp_path):
    """When augustus returns no predictions, no debug directory is left behind."""
    prfl = str(tmp_path / "busco42.prfl")
    open(prfl, "w").close()

    with patch("buscolite.busco.proteinprofile", return_value={}):
        result = predict_and_validate(
            _minimal_fadict(),
            "contig1",
            prfl,
            cutoffs={"busco42": {"score": 100.0}},
            species="human",
            start=0,
            end=20,
            strand="+",
            configpath="/nope",
            blast_score=1.0,
        )

    assert result is False


class _FakeLogger:
    def __init__(self):
        self.warnings = []

    def warning(self, msg):
        self.warnings.append(msg)


def _make_job(i):
    return {
        "contig": "contig1",
        "prfl": "/tmp/fake_{}.prfl".format(i),
        "start": i * 100,
        "end": (i + 1) * 100,
        "strand": "both",
        "score": 42.0,
        "priority": 1,
    }


def test_process_busco_jobs_resilient_one_failure():
    """One augustus failure should not abort processing of remaining jobs."""
    logger = _FakeLogger()
    err = subprocess.CalledProcessError(-4, ["augustus"], stderr="boom")
    err.buscolite_debug_dir = "/tmp/fake_debug"
    success = ("bx", {"status": "fragmented", "contig": "c"})

    with patch("buscolite.busco.predict_and_validate", side_effect=[err, success, False]):
        busco_id, results, submitted, skipped = process_busco_jobs(
            "bx", [_make_job(0), _make_job(1), _make_job(2)], {}, {}, "human", logger
        )

    assert busco_id == "bx"
    assert submitted == 3
    assert len(results) == 1
    assert len(logger.warnings) == 1
    assert "/tmp/fake_debug" in logger.warnings[0]
    assert "bx" in logger.warnings[0]


def test_process_busco_jobs_all_fail():
    """If every job raises, function returns empty results with all warnings logged."""
    logger = _FakeLogger()
    err = subprocess.CalledProcessError(-4, ["augustus"], stderr="boom")

    with patch("buscolite.busco.predict_and_validate", side_effect=[err, err, err]):
        busco_id, results, submitted, skipped = process_busco_jobs(
            "bx", [_make_job(0), _make_job(1), _make_job(2)], {}, {}, "human", logger
        )

    assert results == []
    assert submitted == 3
    assert skipped == 0
    assert len(logger.warnings) == 3


def test_process_busco_jobs_early_exit_on_complete():
    """Finding a complete result should stop submission of remaining jobs."""
    logger = _FakeLogger()
    complete = ("bx", {"status": "complete", "contig": "c"})

    with patch("buscolite.busco.predict_and_validate", side_effect=[complete]) as m:
        busco_id, results, submitted, skipped = process_busco_jobs(
            "bx", [_make_job(0), _make_job(1), _make_job(2)], {}, {}, "human", logger
        )

    assert submitted == 1
    assert skipped == 2
    assert len(results) == 1
    assert m.call_count == 1
