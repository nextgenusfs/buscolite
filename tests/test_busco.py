"""
Tests for the busco module.
"""

import os
import pytest
from buscolite.busco import load_config, load_cutoffs, check_lineage


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
    import tempfile
    import shutil

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
            shutil.copy(
                os.path.join(mock_busco_lineage, f), os.path.join(new_lineage, f)
            )

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
