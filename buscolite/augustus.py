import sys
import os
import subprocess
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
    functional = False
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    stdout, stderr = proc.communicate()
    stderr = stderr.strip()
    if isinstance(stdout, str):
        try:
            stdout = stdout.decode("ascii", "ignore").encode("ascii")
        except AttributeError:
            pass
    stdout = stdout.strip().split("\n")
    if stderr.startswith("augustus: ERROR"):
        sys.stderr.write(stderr + "\n")
        return functional
    else:
        for line in stdout:
            line = line.strip()
            if line.startswith("# start gene g1"):
                functional = True
    return functional


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
