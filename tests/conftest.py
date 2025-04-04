"""
Configuration file for pytest.
"""
import os
import pytest
import tempfile


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def mock_busco_lineage(temp_dir):
    """Create a mock BUSCO lineage directory with required files and structure."""
    lineage_dir = os.path.join(temp_dir, "mock_lineage")
    
    # Create required directories
    os.makedirs(os.path.join(lineage_dir, "hmms"), exist_ok=True)
    os.makedirs(os.path.join(lineage_dir, "prfl"), exist_ok=True)
    
    # Create required files
    with open(os.path.join(lineage_dir, "dataset.cfg"), "w") as f:
        f.write("name=test_lineage\n")
        f.write("species=test_species\n")
        f.write("domain=test_domain\n")
    
    with open(os.path.join(lineage_dir, "scores_cutoff"), "w") as f:
        f.write("busco1\t100.0\n")
        f.write("busco2\t200.0\n")
    
    with open(os.path.join(lineage_dir, "lengths_cutoff"), "w") as f:
        f.write("busco1\t0\t1.5\t300\n")
        f.write("busco2\t0\t0.0\t400\n")
    
    # Create empty files for the remaining required files
    open(os.path.join(lineage_dir, "ancestral"), "w").close()
    open(os.path.join(lineage_dir, "ancestral_variants"), "w").close()
    
    # Create a mock HMM file
    with open(os.path.join(lineage_dir, "hmms", "busco1.hmm"), "w") as f:
        f.write("HMMER3/f [3.1b2 | February 2015]\n")
        f.write("NAME  busco1\n")
        f.write("ACC   PF00001.1\n")
        f.write("DESC  Test BUSCO HMM\n")
        f.write("LENG  300\n")
        f.write("ALPH  amino\n")
        f.write("//\n")
    
    return lineage_dir


@pytest.fixture
def mock_fasta_file(temp_dir):
    """Create a mock FASTA file for testing."""
    fasta_path = os.path.join(temp_dir, "test.fasta")
    
    with open(fasta_path, "w") as f:
        f.write(">seq1\n")
        f.write("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")
        f.write(">seq2\n")
        f.write("GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n")
    
    return fasta_path


@pytest.fixture
def mock_protein_fasta(temp_dir):
    """Create a mock protein FASTA file for testing."""
    fasta_path = os.path.join(temp_dir, "test_proteins.fasta")
    
    with open(fasta_path, "w") as f:
        f.write(">protein1\n")
        f.write("MVNLKPTSAGRTWLKTIIIGVISAIILVVVIVIILIITSRRLNR\n")
        f.write(">protein2\n")
        f.write("MEEAKQKVVDFLNSKGYKGSTFHRVIPSFYVKETTDLPAKGLVD\n")
    
    return fasta_path
