"""
Tests for the new filtering functions.
"""

from buscolite.utilities import filter_low_scoring_matches, remove_duplicate_gene_matches


def test_filter_low_scoring_matches():
    """Test the 85% bitscore threshold filter."""
    # Test data: BUSCO with multiple matches
    busco_results = {
        "BUSCO1": [
            {"hit": "gene1", "bitscore": 100.0},
            {"hit": "gene2", "bitscore": 90.0},  # 90% of top - should pass
            {"hit": "gene3", "bitscore": 80.0},  # 80% of top - should fail
        ],
        "BUSCO2": [
            {"hit": "gene4", "bitscore": 200.0},
            {"hit": "gene5", "bitscore": 170.0},  # 85% of top - should pass
            {"hit": "gene6", "bitscore": 160.0},  # 80% of top - should fail
        ],
    }

    filtered = filter_low_scoring_matches(busco_results, threshold=0.85)

    # Check BUSCO1: should keep gene1 and gene2, remove gene3
    assert len(filtered["BUSCO1"]) == 2
    assert filtered["BUSCO1"][0]["hit"] == "gene1"
    assert filtered["BUSCO1"][1]["hit"] == "gene2"

    # Check BUSCO2: should keep gene4 and gene5, remove gene6
    assert len(filtered["BUSCO2"]) == 2
    assert filtered["BUSCO2"][0]["hit"] == "gene4"
    assert filtered["BUSCO2"][1]["hit"] == "gene5"


def test_filter_low_scoring_matches_single():
    """Test filtering with single match per BUSCO."""
    busco_results = {
        "BUSCO1": [{"hit": "gene1", "bitscore": 100.0}],
    }

    filtered = filter_low_scoring_matches(busco_results, threshold=0.85)

    # Single match should always be kept
    assert len(filtered["BUSCO1"]) == 1
    assert filtered["BUSCO1"][0]["hit"] == "gene1"


def test_filter_low_scoring_matches_empty():
    """Test filtering with empty results."""
    busco_results = {}

    filtered = filter_low_scoring_matches(busco_results, threshold=0.85)

    assert len(filtered) == 0


def test_remove_duplicate_gene_matches():
    """Test duplicate gene removal across BUSCOs."""
    # Test data: gene1 matches both BUSCO1 and BUSCO2
    busco_results = {
        "BUSCO1": [
            {"hit": "gene1", "bitscore": 100.0},
            {"hit": "gene2", "bitscore": 90.0},
        ],
        "BUSCO2": [
            {"hit": "gene1", "bitscore": 80.0},  # Lower score - should be removed
            {"hit": "gene3", "bitscore": 85.0},
        ],
    }

    filtered = remove_duplicate_gene_matches(busco_results, score_key="bitscore")

    # gene1 should only appear in BUSCO1 (higher score)
    assert len(filtered["BUSCO1"]) == 2
    assert filtered["BUSCO1"][0]["hit"] == "gene1"
    assert filtered["BUSCO1"][1]["hit"] == "gene2"

    # gene1 should be removed from BUSCO2
    assert len(filtered["BUSCO2"]) == 1
    assert filtered["BUSCO2"][0]["hit"] == "gene3"


def test_remove_duplicate_gene_matches_no_duplicates():
    """Test duplicate removal when there are no duplicates."""
    busco_results = {
        "BUSCO1": [{"hit": "gene1", "bitscore": 100.0}],
        "BUSCO2": [{"hit": "gene2", "bitscore": 90.0}],
    }

    filtered = remove_duplicate_gene_matches(busco_results, score_key="bitscore")

    # Nothing should be removed
    assert len(filtered["BUSCO1"]) == 1
    assert len(filtered["BUSCO2"]) == 1
    assert filtered["BUSCO1"][0]["hit"] == "gene1"
    assert filtered["BUSCO2"][0]["hit"] == "gene2"


def test_remove_duplicate_gene_matches_genome_mode():
    """Test duplicate removal for genome mode (with contig and location)."""
    busco_results = {
        "BUSCO1": [
            {"contig": "chr1", "location": [100, 200], "bitscore": 100.0},
        ],
        "BUSCO2": [
            {
                "contig": "chr1",
                "location": [100, 200],
                "bitscore": 80.0,
            },  # Same location, lower score
        ],
    }

    filtered = remove_duplicate_gene_matches(busco_results, score_key="bitscore")

    # The same genomic location should only appear once (in BUSCO1)
    assert len(filtered["BUSCO1"]) == 1
    assert "BUSCO2" not in filtered or len(filtered["BUSCO2"]) == 0
