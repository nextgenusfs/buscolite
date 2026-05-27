import os
import signal
import subprocess
import sys

from .utilities import execute_timeout


def augustus_version():
    p1 = subprocess.Popen(
        ["augustus", "--version"],
        stderr=subprocess.STDOUT,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    ).communicate()
    stdout, stderr = p1
    if isinstance(stdout, str):
        try:
            stdout = stdout.decode("ascii", "ignore").encode("ascii")
        except AttributeError:
            pass
    version = stdout.split(" is ")[0]
    if "(" in version:
        version = version.split("(")[-1]
    if ")" in version:
        version = version.split(")")[0]
    return version


def augustus_functional():
    """
    Run augustus with --proteinprofile on a small bundled example and check
    whether the binary is usable on this host. Returns True if augustus
    produced the expected gene prediction marker, False otherwise.

    On failure, a diagnostic line is written to stderr distinguishing:
      * augustus binary missing from PATH
      * augustus killed by a signal (typically SIGILL — the binary contains
        CPU instructions this host cannot run, e.g. a linux/amd64 bioconda
        build under Apple Silicon / Rosetta 2, or pre-AVX hardware)
      * augustus exited non-zero with its own "augustus: ERROR" message
      * augustus exited cleanly but produced no recognizable output
    """
    datadir = os.path.join(os.path.dirname(__file__), "data")
    cmd = [
        "augustus",
        "--species=human",
        "--softmasking=1",
        "--gff3=on",
        "--UTR=off",
        "--proteinprofile={}".format(os.path.join(datadir, "HsDHC.prfl")),
        "--stopCodonExcludedFromCDS=False",
        os.path.join(datadir, "example.fa"),
    ]
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=120,
        )
    except FileNotFoundError:
        sys.stderr.write("augustus binary not found on PATH\n")
        return False
    except subprocess.TimeoutExpired:
        sys.stderr.write("augustus --proteinprofile self-test timed out after 120s\n")
        return False

    stdout = (proc.stdout or "").strip()
    stderr = (proc.stderr or "").strip()
    rc = proc.returncode

    if rc is not None and rc < 0:
        sig = -rc
        try:
            sig_name = signal.Signals(sig).name
        except (ValueError, AttributeError):
            sig_name = ""
        sys.stderr.write(
            "augustus --proteinprofile self-test was killed by signal {} ({}). "
            "This typically means the augustus binary contains CPU "
            "instructions your host does not support — most commonly seen "
            "when running a linux/amd64 docker image on Apple Silicon via "
            "Rosetta 2, or on older x86_64 hardware lacking AVX/AVX2/BMI2. "
            "Run on real x86_64 hardware or use a native install.\n".format(sig, sig_name)
        )
        return False

    if stderr.startswith("augustus: ERROR"):
        sys.stderr.write(stderr + "\n")
        return False

    for line in stdout.split("\n"):
        if line.strip().startswith("# start gene g1"):
            return True

    sys.stderr.write(
        "augustus --proteinprofile self-test exited rc={} but produced no "
        "gene prediction. stderr was: {!r}\n".format(rc, stderr[:500])
    )
    return False


def proteinprofile(
    fasta,
    prfl,
    species="anidulans",
    start=False,
    end=False,
    strand="both",
    configpath=False,
):
    if not configpath:
        try:
            configpath = os.environ["AUGUSTUS_CONFIG_PATH"]
        except KeyError:
            pass
    cmd = [
        "augustus",
        "--species={}".format(species),
        "--softmasking=1",
        "--gff3=on",
        "--UTR=off",
        "--genemodel=exactlyone",
        "--proteinprofile={}".format(prfl),
        "--stopCodonExcludedFromCDS=False",
        "--strand={}".format(strand),
    ]
    if configpath:
        cmd.append("--AUGUSTUS_CONFIG_PATH={}".format(configpath))
    cmd.append(fasta)
    preds = {}
    ID = None
    for line in execute_timeout(cmd):
        if "\tAUGUSTUS\tgene\t" in line:
            cols = line.rstrip().split("\t")
            ID = cols[8].replace("ID=", "")
            preds[ID] = {
                "contig": cols[0],
                "strand": cols[6],
                "location": (int(cols[3]) + start, int(cols[4]) + start),
                "coords": [],
                "phase": [],
            }
        elif "\tAUGUSTUS\tCDS\t" in line:
            cols = line.rstrip().split("\t")
            ID = cols[8].split("Parent=")[-1].split(".")[0]
            preds[ID]["coords"].append((int(cols[3]) + start, int(cols[4]) + start))
            preds[ID]["phase"].append(int(cols[7]))
    return preds
