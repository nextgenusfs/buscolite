"""
Tests for the main module.
"""

# import os  # Not used
# import sys  # Not used directly, but needed for patching
from unittest.mock import patch  # MagicMock not used

import pytest

from buscolite.__main__ import parse_args


def test_parse_args_required():
    """Test parse_args with required arguments."""
    # Test with all required arguments
    test_args = [
        "-i",
        "input.fasta",
        "-o",
        "output_name",
        "-m",
        "genome",
        "-l",
        "/path/to/lineage",
    ]

    args = parse_args(test_args)

    assert args.input == "input.fasta"
    assert args.out == "output_name"
    assert args.mode == "genome"
    assert args.lineage == "/path/to/lineage"
    assert args.cpus == 1  # Default value
    assert args.species == "anidulans"  # Default value
    assert args.flanks == 2000  # Default value


def test_parse_args_all_options():
    """Test parse_args with all arguments."""
    # Test with all arguments
    test_args = [
        "-i",
        "input.fasta",
        "-o",
        "output_name",
        "-m",
        "proteins",
        "-l",
        "/path/to/lineage",
        "-c",
        "8",
        "-s",
        "human",
        "-f",
        "5000",
    ]

    args = parse_args(test_args)

    assert args.input == "input.fasta"
    assert args.out == "output_name"
    assert args.mode == "proteins"
    assert args.lineage == "/path/to/lineage"
    assert args.cpus == 8
    assert args.species == "human"
    assert args.flanks == 5000


def test_parse_args_invalid_mode():
    """Test parse_args with invalid mode."""
    # Test with invalid mode
    test_args = [
        "-i",
        "input.fasta",
        "-o",
        "output_name",
        "-m",
        "invalid",
        "-l",
        "/path/to/lineage",
    ]

    with pytest.raises(SystemExit):
        parse_args(test_args)


@patch("sys.stderr")
def test_parse_args_no_args(_):
    """Test parse_args with no arguments."""
    # Test with no arguments
    test_args = []

    with pytest.raises(SystemExit):
        parse_args(test_args)
