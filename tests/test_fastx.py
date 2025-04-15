"""
Tests for the fastx module.
"""

# import os  # Not used in this file

# import pytest  # Not used directly, but needed for pytest fixtures

from buscolite.fastx import (
    RevComp,
    dict2stats,
    fasta2dict,
    fasta2headers,
    fasta2lengths,
    getSeqRegions,
    softwrap,
    translate,
)


def test_revcomp():
    """Test the RevComp function."""
    # Test with simple sequence
    assert RevComp("ATGC") == "GCAT"

    # Test with ambiguous nucleotides
    assert RevComp("ATGCRYSWKMBDHVN") == "NBDHVKMWSRYGCAT"

    # Test with lowercase (should be converted to uppercase)
    assert RevComp("atgc") == "GCAT"

    # Test with empty string
    assert RevComp("") == ""


def test_translate():
    """Test the translate function."""
    # Test forward strand translation
    assert translate("ATGAAATAG", "+", 0) == "MK*"

    # Test reverse strand translation
    assert translate("CTATTTCAT", "-", 0) == "MK*"

    # Test with phase 1
    assert translate("CATGAAATAG", "+", 1) == "MK*"

    # Test with phase 2
    assert translate("ACATGAAATAG", "+", 2) == "MK*"

    # Test with incomplete codon at end
    assert translate("ATGAAA", "+", 0) == "MK"

    # Test with alternative start codon (table 11)
    assert translate("TTGAAATAG", "+", 0, table=11) == "MK*"


def test_fasta2dict(mock_fasta_file):
    """Test the fasta2dict function."""
    # Test with default parameters
    result = fasta2dict(mock_fasta_file)

    assert isinstance(result, dict)
    assert len(result) == 2
    assert "seq1" in result
    assert "seq2" in result
    assert result["seq1"] == "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    assert result["seq2"] == "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"


def test_fasta2headers(mock_fasta_file):
    """Test the fasta2headers function."""
    # Test with default parameters
    result = fasta2headers(mock_fasta_file)

    assert isinstance(result, set)
    assert len(result) == 2
    assert "seq1" in result
    assert "seq2" in result


def test_fasta2lengths(mock_fasta_file):
    """Test the fasta2lengths function."""
    # Test with default parameters
    result = fasta2lengths(mock_fasta_file)

    assert isinstance(result, dict)
    assert len(result) == 2
    assert "seq1" in result
    assert "seq2" in result
    assert result["seq1"] == 52
    assert result["seq2"] == 52


def test_getseqregions():
    """Test the getSeqRegions function."""
    # Create a test dictionary
    seqs = {"test": "ATGCATGCATGCATGCATGC"}

    # Test with single region
    result = getSeqRegions(seqs, "test", [(1, 10)])
    assert result == "ATGCATGCAT"

    # Test with multiple regions
    result = getSeqRegions(seqs, "test", [(1, 5), (11, 15)])
    assert result == "ATGCAGCATG"

    # Test with unsorted regions (should be sorted)
    result = getSeqRegions(seqs, "test", [(11, 15), (1, 5)])
    assert result == "ATGCAGCATG"


def test_dict2stats():
    """Test the dict2stats function."""
    # Create a test dictionary
    fadict = {"seq1": "ATGCATGCAT", "seq2": "GCATGCATGC", "seq3": "ATGCATGCATGCATGCAT"}

    # Test the function
    result = dict2stats(fadict)

    assert isinstance(result, dict)
    assert result["n"] == 3
    assert result["size"] == 38


def test_softwrap():
    """Test the softwrap function."""
    # Create a test string
    test_str = "A" * 100

    # Test with default wrap (80)
    result = softwrap(test_str)
    assert result == "A" * 80 + "\n" + "A" * 20

    # Test with custom wrap
    result = softwrap(test_str, every=50)
    assert result == "A" * 50 + "\n" + "A" * 50

    # Test with string shorter than wrap
    result = softwrap("ATGC", every=10)
    assert result == "ATGC"
