# import importlib.metadata  # Unused import
import os
import shutil
import subprocess
import uuid
from urllib.parse import unquote

import pyhmmer
from packaging.version import parse as parse_version

from .fastx import fasta2dict, fasta2headers
from .gff import validate_models
from .utilities import execute, runprocess


def tblastn_version():
    """
    Get the version of tblastn installed on the system.

    Returns
    -------
    str
        The version string of tblastn
    """
    p1 = subprocess.Popen(
        ["tblastn", "-version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    vers = p1.communicate()[0].split("+")[0]
    vers = vers.split(" ")[-1]
    return vers


def miniprot_version():
    """
    Get the version of miniprot installed on the system.

    Returns
    -------
    str
        The version string of miniprot
    """
    p1 = subprocess.Popen(
        ["miniprot", "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    vers = p1.communicate()[0].rstrip()
    return vers


def pyhmmer_version():
    """
    Get the version of pyhmmer installed on the system.

    Returns
    -------
    str
        The version string of pyhmmer
    """
    return pyhmmer.__version__


def miniprot_prefilter(
    input, query, cutoffs, tmpdir=False, cpus=1, maxhits=3, maxintron=10000, buscodb="."
):
    """
    Prefilter alignments using miniprot for BUSCO analysis.

    Parameters
    ----------
    input : str
        Path to the input genome FASTA file
    query : str
        Path to the query protein FASTA file
    cutoffs : dict
        Dictionary containing score cutoffs for each BUSCO model
    tmpdir : str or bool, optional
        Path to temporary directory, or False to use system default
    cpus : int, optional
        Number of CPU threads to use (default: 1)
    maxhits : int, optional
        Maximum number of hits to return per query (default: 3)
    maxintron : int, optional
        Maximum intron size (default: 10000)
    buscodb : str, optional
        Path to the BUSCO database directory (default: ".")

    Returns
    -------
    dict
        Dictionary containing filtered alignment results
    """
    # Continue with function implementation
    # setup tmpdir
    if not tmpdir:
        tmpdir = os.path.join("/tmp", str(uuid.uuid4()))
        os.makedirs(tmpdir)
    # check fasta headers
    headers = fasta2headers(input)
    # setup miniprot command
    cmd = [
        "miniprot",
        "-G",
        str(maxintron),
        "-N",
        str(maxhits),
        "-t",
        str(cpus),
        "--gff",
        os.path.abspath(input),
        os.path.abspath(query),
    ]
    results = {}
    Genes = {}
    for line in execute(cmd):
        if line.startswith("##PAF"):
            (
                pname,
                plen,
                pstart,
                pend,
                strand,
                cname,
                _,  # clen (not used)
                cstart,
                cend,
                _,  # match (not used)
                _,  # naln (not used)
                _,  # mapq (not used)
            ) = line.split("\t")[1:13]
            score = 0
            escore = 0.0
            for x in line.split("\t")[13:]:
                if x.startswith("AS:"):
                    score = int(x.split(":")[-1])
            qcov = (int(pend) - int(pstart)) / float(int(plen))
            if cname not in headers:
                if cname.startswith("gb|"):
                    cname = cname.replace("gb|", "").rstrip("|")
            if pname not in results:
                results[pname] = {cname: []}
            if cname not in results[pname]:
                results[pname][cname] = [
                    {
                        "evalue": escore,
                        "score": score,
                        "qcov": qcov,
                        "qlen": int(plen),
                        "cov": (int(pstart), int(pend)),
                        "coords": (int(cstart), int(cend)),
                        "strand": strand,
                    }
                ]
            else:
                results[pname][cname].append(
                    {
                        "evalue": escore,
                        "score": score,
                        "qcov": qcov,
                        "qlen": int(plen),
                        "cov": (int(pstart), int(pend)),
                        "coords": (int(cstart), int(cend)),
                        "strand": strand,
                    }
                )
        elif line.startswith(("##gff-version", "\n")):
            continue
        else:
            line = line.rstrip()
            # skip lines that aren't 9 columns
            if not line.count("\t") == 8:
                continue
            (
                contig,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = line.split("\t")
            if feature not in [
                "mRNA",
                "CDS",
                "stop_codon",
            ]:
                continue
            attributes = unquote(attributes)
            source = unquote(source)
            feature = unquote(feature)
            start = int(start)
            end = int(end)
            ID = None
            Parent = None
            Name = None
            # These variables are not used but kept for future compatibility
            # Product = None
            # GeneFeature = None
            # gbkey = None
            info = {}
            for field in attributes.split(";"):
                try:
                    k, v = field.split("=", 1)
                    info[k] = v.strip()
                except (IndexError, ValueError):
                    pass
            # now can lookup in info dict for values
            ID = info.get("ID", None)
            Parent = info.get("Parent", None)
            Target = info.get("Target", None)
            Identity = info.get("Identity", None)
            Positive = info.get("Positive", None)
            Rank = info.get("Rank", None)
            if Target:
                Name = Target.split()[0]
            # now we can do add to dictionary these parsed values
            # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
            if feature in ["mRNA"]:
                if ID not in Genes:
                    Genes[ID] = {
                        "name": Name,
                        "type": ["mRNA"],
                        "transcript": [],
                        "cds_transcript": [],
                        "protein": [],
                        "5UTR": [[]],
                        "3UTR": [[]],
                        "gene_synonym": [],
                        "codon_start": [],
                        "ids": [ID],
                        "CDS": [[]],
                        "mRNA": [[]],
                        "strand": strand,
                        "EC_number": [[]],
                        "location": (start, end),
                        "contig": contig,
                        "product": ["miniprot alignment"],
                        "source": source,
                        "phase": [[]],
                        "db_xref": [[]],
                        "go_terms": [[]],
                        "note": [
                            [
                                f"IDENTITY:{Identity}",
                                f"RANK:{Rank}",
                                f"Positive:{Positive}",
                            ]
                        ],
                        "partialStart": [],
                        "partialStop": [],
                        "pseudo": False,
                        "score": [],
                        "target": [],
                    }
                else:
                    if start < Genes[ID]["location"][0]:
                        Genes[ID]["location"] = (start, Genes[ID]["location"][1])
                    if end > Genes[ID]["location"][1]:
                        Genes[ID]["location"] = (Genes[ID]["location"][0], end)
                    Genes[ID]["Note"].append(
                        [
                            f"IDENTITY:{Identity}",
                            f"RANK:{Rank}",
                            f"Positive:{Positive}",
                        ]
                    )
            elif feature in ["CDS"]:
                if Parent:
                    # determine which transcript this is get index from id
                    i = Genes[Parent]["ids"].index(Parent)
                    Genes[Parent]["CDS"][i].append((start, end))
                    Genes[Parent]["mRNA"][i].append((start, end))
                    Genes[Parent]["score"].append(round(float(Identity) * 100, 2))
                    Genes[Parent]["target"].append(Target)
                    # add phase
                    try:
                        Genes[Parent]["phase"][i].append(int(phase))
                    except ValueError:
                        Genes[Parent]["phase"][i].append("?")
            elif feature in ["stop_codon"]:
                if Parent:
                    i = Genes[Parent]["ids"].index(Parent)
                    # here we need to extend the last CDS if + and first if - strand
                    if strand == "+":
                        lt = Genes[Parent]["CDS"][i][-1]
                        nLT = (lt[0], end)
                        Genes[Parent]["CDS"][i][-1] = nLT
                        Genes[Parent]["mRNA"][i][-1] = nLT
                    elif strand == "-":
                        ft = Genes[Parent]["CDS"][i][0]
                        nFT = (start, ft[1])
                        Genes[Parent]["CDS"][i][0] = nFT
                        Genes[Parent]["mRNA"][i][0] = nFT

    # we can then process the Genes dictionary for full length genes
    SeqRecords = fasta2dict(input)
    annotation = validate_models(Genes, SeqRecords, table=1)
    g = {}
    for k, v in annotation.items():
        if "mRNA" in v["type"]:
            if False in v["partialStart"] and False in v["partialStop"] and not v["pseudo"]:
                # check complete models with hmmmer to validate and reformat
                hmmfile = os.path.join(buscodb, "hmms", "{}.hmm".format(v["name"]))
                if not os.path.isfile(hmmfile):
                    continue
                # now we can check via hmmer
                hmm_result = hmmer_search_single(hmmfile, v["protein"][0].rstrip("*"))
                if len(hmm_result) > 0:
                    if hmm_result[0]["bitscore"] > cutoffs[v["name"]]["score"]:
                        final = {
                            "contig": v["contig"],
                            "strand": v["strand"],
                            "location": v["location"],
                            "coords": v["CDS"][0],
                            "phase": v["phase"][0],
                            "transcript": v["cds_transcript"][0],
                            "translation": v["protein"][0],
                            "status": "complete",
                            "hmmer": hmm_result[0],
                            "miniprot_score": results[v["name"]][v["contig"]][0]["score"],
                        }
                        if v["name"] not in g:
                            g[v["name"]] = [final]
                        else:
                            g[v["name"]].append(final)
    # want to return a dictionary of busco hits
    final = {}
    njobs = 0
    for k, v in sorted(results.items()):
        if k not in g:
            for x, z in sorted(v.items()):
                merged = merge_overlapping_hits(z, fluff=maxintron)
                for w in merged:
                    njobs += 1
                    if k not in final:
                        final[k] = [
                            {
                                "contig": x,
                                "coords": w["coords"],
                                "evalue": w["evalue"],
                                "score": w["score"],
                                "strand": w["strand"],
                            }
                        ]
                    else:
                        final[k].append(
                            {
                                "contig": x,
                                "coords": w["coords"],
                                "evalue": w["evalue"],
                                "score": w["score"],
                                "strand": w["strand"],
                            }
                        )
    return g, final, njobs


def blast_prefilter(
    input,
    query,
    logger=None,
    tmpdir=False,
    evalue=1e-3,
    cpus=1,
    maxhits=3,
    maxintron=10000,
    # logger parameter is kept for API compatibility but not used
):
    """
    Prefilters alignments for augustus using tblastn.

    Parameters
    ----------
    input : str
        Path to the input genome FASTA file
    query : str
        Path to the query protein FASTA file
    logger : object
        Logger object for logging messages
    tmpdir : str or bool, optional
        Path to temporary directory, or False to use system default
    evalue : float, optional
        E-value threshold for BLAST hits (default: 1e-3)
    cpus : int, optional
        Number of CPU threads to use (default: 1)
    maxhits : int, optional
        Maximum number of hits to return per query (default: 3)
    maxintron : int, optional
        Maximum intron size (default: 10000)

    Returns
    -------
    tuple
        (dict, int) - Dictionary containing filtered alignment results and number of jobs
    """
    # figure out cpus
    if cpus > 1:  # need to see if tblastn is safe multithreading
        if parse_version("2.2.32") < parse_version(tblastn_version()) < parse_version("2.10.2"):
            cpus = 1
    # setup tmpdir
    if not tmpdir:
        tmpdir = os.path.join("/tmp", str(uuid.uuid4()))
        os.makedirs(tmpdir)
    # check fasta headers
    headers = fasta2headers(input)
    # start by formatting blast db/dustmasker filtered format
    cmd = [
        "dustmasker",
        "-in",
        os.path.abspath(input),
        "-infmt",
        "fasta",
        "-parse_seqids",
        "-outfmt",
        "maskinfo_asn1_bin",
        "-out",
        "genome_dust.asnb",
    ]
    runprocess(cmd, cwd=tmpdir)
    cmd = [
        "makeblastdb",
        "-in",
        os.path.abspath(input),
        "-dbtype",
        "nucl",
        "-parse_seqids",
        "-mask_data",
        "genome_dust.asnb",
        "-out",
        "genome",
    ]
    runprocess(cmd, cwd=tmpdir)
    cmd = [
        "tblastn",
        "-num_threads",
        str(cpus),
        "-db",
        os.path.join(tmpdir, "genome"),
        "-query",
        query,
        "-max_target_seqs",
        str(maxhits),
        "-db_soft_mask",
        "11",
        "-threshold",
        "999",
        "-max_intron_length",
        str(maxintron),
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 sseqid slen sstart send qseqid qlen qstart qend pident length evalue score",
    ]
    results = {}
    for line in execute(cmd):
        cols = line.rstrip().split("\t")
        qcov = (int(cols[7]) - int(cols[6])) / float(int(cols[5]))
        score = float(cols[11])
        escore = float(cols[10])
        contig = cols[0]
        if contig not in headers:
            if contig.startswith("gb|"):
                contig = contig.replace("gb|", "").rstrip("|")
        if int(cols[3]) > int(cols[2]):
            start = int(cols[2])
            end = int(cols[3])
        else:
            start = int(cols[3])
            end = int(cols[2])
        if cols[4] not in results:
            results[cols[4]] = {contig: []}
        if contig not in results[cols[4]]:
            results[cols[4]][contig] = [
                {
                    "evalue": escore,
                    "score": score,
                    "qcov": qcov,
                    "qlen": int(cols[5]),
                    "cov": (int(cols[6]), int(cols[7])),
                    "coords": (start, end),
                }
            ]
        else:
            results[cols[4]][contig].append(
                {
                    "evalue": escore,
                    "score": score,
                    "qcov": qcov,
                    "qlen": int(cols[5]),
                    "cov": (int(cols[6]), int(cols[7])),
                    "coords": (start, end),
                }
            )
    # clean up intermediate files
    shutil.rmtree(tmpdir)

    # now we can parse and filter results
    # want to return a dictionary of busco hits
    final = {}
    njobs = 0
    for k, v in sorted(results.items()):
        for x, z in sorted(v.items()):
            merged = merge_overlapping_hits(z, fluff=maxintron)
            for w in merged:
                njobs += 1
                if k not in final:
                    final[k] = [
                        {
                            "contig": x,
                            "coords": w["coords"],
                            "evalue": w["evalue"],
                            "score": w["score"],
                        }
                    ]
                else:
                    final[k].append(
                        {
                            "contig": x,
                            "coords": w["coords"],
                            "evalue": w["evalue"],
                            "score": w["score"],
                        }
                    )
    return final, njobs


def merge_overlapping_hits(queryList, fluff=10000):
    """
    Merge overlapping or nearby hits from BLAST or miniprot searches.

    Parameters
    ----------
    queryList : list
        List of dictionaries containing hit information with 'coords' key
    fluff : int, optional
        Maximum distance between hits to consider them for merging (default: 10000)

    Returns
    -------
    list
        List of merged hits
    """
    if len(queryList) < 2:
        return queryList
    # sort by coords and try to merge if within fluff
    queryList.sort(key=lambda x: x["coords"][0])
    result = [queryList[0]]
    for x in queryList[1:]:
        if (
            x["coords"][0] <= result[-1]["coords"][1]
            or (result[-1]["coords"][1] - x["coords"][0]) < fluff
        ):
            result[-1]["coords"] = (
                result[-1]["coords"][0],
                max(result[-1]["coords"][1], x["coords"][1]),
            )
        else:
            result.append(x)
    return result


def hmmer_search_single(hmmfile, seq):
    """
    Search a single protein sequence against an HMM profile using pyhmmer.

    Parameters
    ----------
    hmmfile : str
        Path to the HMM profile file
    seq : str
        Protein sequence to search

    Returns
    -------
    list
        List of dictionaries containing search results with the following keys:
        - name: Name of the HMM profile
        - bitscore: Bit score of the hit
        - evalue: E-value of the hit
        - domains: List of domain hits with coordinates and scores
    """
    hmm = next(pyhmmer.plan7.HMMFile(hmmfile))
    prot = pyhmmer.easel.TextSequence(sequence=seq)
    alphabet = pyhmmer.easel.Alphabet.amino()
    results = []
    for top_hits in pyhmmer.hmmsearch([hmm], [prot.digitize(alphabet)]):
        for hit in top_hits:
            cog = hit.best_domain.alignment.hmm_name.decode()
            domains = []
            for h in hit.domains:
                domains.append(
                    {
                        "hmm_from": h.alignment.hmm_from,
                        "hmm_to": h.alignment.hmm_to,
                        "env_from": h.env_from,
                        "env_to": h.env_to,
                        "score": h.score,
                    }
                )
            results.append(
                {
                    "name": cog,
                    "bitscore": hit.score,
                    "evalue": hit.evalue,
                    "domains": domains,
                }
            )
    return results


def hmmer_search(hmmfile, sequences):
    """
    Search multiple protein sequences against an HMM profile using pyhmmer.

    Parameters
    ----------
    hmmfile : str
        Path to the HMM profile file
    sequences : list
        List of digitized protein sequences to search

    Returns
    -------
    list
        List of dictionaries containing search results with the following keys:
        - name: Name of the HMM profile
        - hit: Name of the sequence that matched
        - bitscore: Bit score of the hit
        - evalue: E-value of the hit
        - domains: List of domain hits with coordinates and scores
    """
    hmm = next(pyhmmer.plan7.HMMFile(hmmfile))
    results = []
    for top_hits in pyhmmer.hmmsearch([hmm], sequences):
        for hit in top_hits:
            cog = hit.best_domain.alignment.hmm_name.decode()
            domains = []
            for h in hit.domains:
                domains.append(
                    {
                        "hmm_from": h.alignment.hmm_from,
                        "hmm_to": h.alignment.hmm_to,
                        "env_from": h.env_from,
                        "env_to": h.env_to,
                        "score": h.score,
                    }
                )
            results.append(
                {
                    "name": cog,
                    "hit": hit.name.decode(),
                    "bitscore": hit.score,
                    "evalue": hit.evalue,
                    "domains": domains,
                }
            )
    return results
