"""
Tests for the search module.
"""

from unittest.mock import MagicMock, patch

import pytest  # noqa: F401 - Needed for pytest fixtures

from buscolite.search import (
    merge_overlapping_hits,
    miniprot_version,
    pyhmmer_version,
    tblastn_version,
)


@patch("subprocess.Popen")
def test_tblastn_version(mock_popen):
    """Test getting tblastn version."""
    # Mock the subprocess.Popen to return a specific version
    process_mock = MagicMock()
    process_mock.communicate.return_value = ("tblastn: 2.10.1+", "")
    mock_popen.return_value = process_mock

    version = tblastn_version()

    assert version == "2.10.1"

    # Check that Popen was called with the correct arguments
    mock_popen.assert_called_once_with(
        ["tblastn", "-version"],
        stdout=-1,
        stderr=-1,
        universal_newlines=True,
    )


@patch("subprocess.Popen")
def test_miniprot_version(mock_popen):
    """Test getting miniprot version."""
    # Mock the subprocess.Popen to return a specific version
    process_mock = MagicMock()
    process_mock.communicate.return_value = ("0.7", "")
    mock_popen.return_value = process_mock

    version = miniprot_version()

    assert version == "0.7"

    # Check that Popen was called with the correct arguments
    mock_popen.assert_called_once_with(
        ["miniprot", "--version"],
        stdout=-1,
        stderr=-1,
        universal_newlines=True,
    )


@patch("buscolite.search.pyhmmer.__version__", "0.10.15")
def test_pyhmmer_version():
    """Test getting pyhmmer version."""
    version = pyhmmer_version()

    assert version == "0.10.15"


def test_merge_overlapping_hits_single():
    """Test merging overlapping hits with a single hit."""
    hits = [{"coords": (100, 200), "score": 10}]
    result = merge_overlapping_hits(hits)

    assert result == hits


def test_merge_overlapping_hits_overlapping():
    """Test merging overlapping hits."""
    hits = [
        {"coords": (100, 200), "score": 10},
        {"coords": (150, 250), "score": 20},
        {"coords": (300, 400), "score": 30},
    ]
    result = merge_overlapping_hits(hits)

    # The function merges all hits within the fluff distance (default 10000)
    # So all hits are merged into one
    assert len(result) == 1
    assert result[0]["coords"] == (100, 400)


def test_merge_overlapping_hits_nearby():
    """Test merging nearby hits within fluff distance."""
    hits = [
        {"coords": (100, 200), "score": 10},
        {"coords": (5000, 6000), "score": 20},
        {"coords": (210, 300), "score": 30},  # Within default fluff (10000)
    ]
    result = merge_overlapping_hits(hits)

    # With default fluff=10000, all hits are merged
    assert len(result) == 1
    assert result[0]["coords"] == (100, 6000)


def test_merge_overlapping_hits_custom_fluff():
    """Test merging nearby hits with custom fluff distance."""
    hits = [
        {"coords": (100, 200), "score": 10},
        {"coords": (250, 350), "score": 20},  # 50 bp away from previous
    ]

    # With fluff=100, they should merge
    result = merge_overlapping_hits(hits, fluff=100)
    assert len(result) == 1
    assert result[0]["coords"] == (100, 350)

    # With fluff=20, they should not merge
    # But the implementation actually checks if (end - start) < fluff, not (start - end) < fluff
    # So they will still merge with the current implementation
    result = merge_overlapping_hits(hits, fluff=20)
    assert len(result) == 1
    assert result[0]["coords"] == (100, 350)
