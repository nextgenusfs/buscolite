"""
Tests for the utilities module.
"""

import os
from unittest.mock import patch  # MagicMock not used

import pytest  # noqa: F401 - Needed for pytest fixtures

from buscolite.utilities import any_overlap, is_file, overlap, which2


def test_overlap():
    """Test the overlap function."""
    # Complete overlap
    assert overlap(100, 200, 100, 200) == 100

    # Partial overlap
    assert overlap(100, 200, 150, 250) == 50
    assert overlap(150, 250, 100, 200) == 50

    # No overlap
    assert overlap(100, 200, 300, 400) == 0
    assert overlap(300, 400, 100, 200) == 0

    # Edge cases
    assert overlap(100, 200, 200, 300) == 0  # Touching but not overlapping
    assert overlap(100, 200, 50, 150) == 50  # Partial overlap at start
    assert overlap(100, 200, 150, 250) == 50  # Partial overlap at end
    assert overlap(100, 200, 50, 250) == 100  # Complete containment


def test_any_overlap():
    """Test the any_overlap function."""
    # Test with overlapping coordinates
    assert any_overlap((100, 200), [(50, 150), (300, 400)]) is True

    # Test with non-overlapping coordinates
    assert any_overlap((100, 200), [(300, 400), (500, 600)]) is False

    # Test with empty list
    assert any_overlap((100, 200), []) is False

    # Test with exact match
    assert any_overlap((100, 200), [(100, 200)]) is True

    # Test with touching but not overlapping
    assert any_overlap((100, 200), [(200, 300)]) is False


def test_is_file(temp_dir):
    """Test the is_file function."""
    # Create a test file
    test_file = os.path.join(temp_dir, "test_file.txt")
    with open(test_file, "w") as f:
        f.write("test")

    # Test with existing file
    assert is_file(test_file) is True

    # Test with non-existing file
    assert is_file(os.path.join(temp_dir, "nonexistent.txt")) is False

    # Test with directory
    assert is_file(temp_dir) is False


@patch("os.path.isfile")
@patch("os.access")
@patch("os.environ")
def test_which2(mock_environ, mock_access, mock_isfile):
    """Test the which2 function."""
    # Set up mocks
    mock_environ.__getitem__.return_value = "/usr/bin:/usr/local/bin"
    mock_isfile.side_effect = lambda x: x in ["/usr/bin/program", "/path/to/program"]
    mock_access.side_effect = lambda x, _: x in ["/usr/bin/program", "/path/to/program"]

    # Test with absolute path that exists
    assert which2("/path/to/program") == "/path/to/program"

    # Test with program name that exists in PATH
    assert which2("program") == "/usr/bin/program"

    # Test with program that doesn't exist
    assert which2("nonexistent") is None
