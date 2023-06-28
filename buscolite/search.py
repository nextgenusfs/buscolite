import subprocess
from pkg_resources import parse_version
import uuid
import os
import shutil
import pyhmmer
from .utilities import runprocess, execute
from .fastx import fasta2headers, fasta2dict
from .gff import validate_models
from urllib.parse import unquote


def tblastn_version():
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
    p1 = subprocess.Popen(
        ["miniprot", "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    vers = p1.communicate()[0].rstrip()
    return vers


def pyhmmer_version():
    return pyhmmer.__version__


def miniprot_prefilter(
    input, query, cutoffs, tmpdir=False, cpus=1, maxhits=3, maxintron=10000, buscodb="."
):
    """
    Prefilters alignments for augustus using miniprot.

    Parameters
    ----------
    input : filename
        DNA (genome) sequence as nucleotides
    query : filename
        Protein (amino acids) sequences
    logger : int
        phase to start translation [0,1,2]

    Returns
    -------
    results : dict
        query seq id keyed dictionary

    """
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
                clen,
                cstart,
                cend,
                match,
                naln,
                mapq,
            ) = line.split("\t")[1:13]
            score = 0
            escore = 0.0
            for x in line.split("\t")[13:]:
                if x.startswith("AS:"):
                    score = int(x.split(":")[-1])
            qcov = (int(pend) - int(pstart)) / float(int(plen))
            if not cname in headers:
                if cname.startswith("gb|"):
                    cname = cname.replace("gb|", "").rstrip("|")
            if not pname in results:
                results[pname] = {cname: []}
            if not cname in results[pname]:
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
            Product = None
            GeneFeature = None
            gbkey = None
            info = {}
            for field in attributes.split(";"):
                try:
                    k, v = field.split("=", 1)
                    info[k] = v.strip()
                except (IndexError, ValueError) as E:
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
                if not ID in Genes:
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
            if (
                False in v["partialStart"]
                and False in v["partialStop"]
                and v["pseudo"] == False
            ):
                # check complete models with hmmmer to validate and reformat
                hmmfile = os.path.join(buscodb, "hmms", "{}.hmm".format(v["name"]))
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
                            "miniprot_score": results[v["name"]][v["contig"]][0][
                                "score"
                            ],
                        }
                        if not v["name"] in g:
                            g[v["name"]] = [final]
                        else:
                            g[v["name"]].append(final)
    # want to return a dictionary of busco hits
    final = {}
    njobs = 0
    for k, v in sorted(results.items()):
        if not k in g:
            for x, z in sorted(v.items()):
                merged = merge_overlapping_hits(z, fluff=maxintron)
                for w in merged:
                    njobs += 1
                    if not k in final:
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
    input, query, logger, tmpdir=False, evalue=1e-3, cpus=1, maxhits=3, maxintron=10000
):
    """
    Prefilters alignments for augustus using tblastn.

    Parameters
    ----------
    input : filename
        DNA (genome) sequence as nucleotides
    query : filename
        Protein (amino acids) sequences
    logger : int
        phase to start translation [0,1,2]

    Returns
    -------
    results : dict
        query seq id keyed dictionary

    """
    # figure out cpus
    if cpus > 1:  # need to see if tblastn is safe multithreading
        if (
            parse_version("2.2.32")
            < parse_version(tblastn_version())
            < parse_version("2.10.2")
        ):
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
        if not contig in headers:
            if contig.startswith("gb|"):
                contig = contig.replace("gb|", "").rstrip("|")
        if int(cols[3]) > int(cols[2]):
            start = int(cols[2])
            end = int(cols[3])
        else:
            start = int(cols[3])
            end = int(cols[2])
        if not cols[4] in results:
            results[cols[4]] = {contig: []}
        if not contig in results[cols[4]]:
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
                if not k in final:
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
