import os
import tempfile
import concurrent.futures
from natsort import natsorted
import pyhmmer
import pyfastx
from .__init__ import __version__
from .log import startLogging
from .search import (
    hmmer_search_single,
    hmmer_search,
    miniprot_prefilter,
    miniprot_version,
    pyhmmer_version,
)
from .augustus import proteinprofile, augustus_version, augustus_functional
from .fastx import (
    softwrap,
    getSeqRegions,
    translate,
    fasta2dict,
    fasta2lengths,
    dict2stats,
)
from .utilities import any_overlap


def load_config(lineage):
    config = {}
    with open(os.path.join(lineage, "dataset.cfg"), "r") as infile:
        for line in infile:
            key, value = line.rstrip().split("=")
            config[key] = value
    return config


def load_cutoffs(lineage):
    cutoffs = {}
    with open(os.path.join(lineage, "scores_cutoff"), "r") as infile:
        for line in infile:
            busco, score = line.rstrip().split("\t")
            if busco not in cutoffs:
                cutoffs[busco] = {"score": float(score)}
            else:
                cutoffs[busco]["score"] = float(score)
    with open(os.path.join(lineage, "lengths_cutoff"), "r") as infile:
        for line in infile:
            busco, _, sigma, length = line.rstrip().split("\t")
            if float(sigma) == 0.0:
                sigma = 1
            if busco not in cutoffs:
                cutoffs[busco] = {"sigma": float(sigma), "length": int(float(length))}
            else:
                cutoffs[busco]["sigma"] = float(sigma)
                cutoffs[busco]["length"] = int(float(length))
    return cutoffs


def check_lineage(lineage):
    lineage = os.path.abspath(lineage)
    if not os.path.isdir(lineage):
        return False, "{} is not a directory".format(lineage)
    dirs = ["prfl", "hmms"]
    files = [
        "ancestral",
        "ancestral_variants",
        "dataset.cfg",
        "lengths_cutoff",
        "scores_cutoff",
    ]
    for d in dirs:
        if not os.path.isdir(os.path.join(lineage, d)):
            return False, "{} directory was not found in {}".format(d, lineage)
    for f in files:
        if not os.path.isfile(os.path.join(lineage, f)):
            return False, "{} file is missing from {}".format(f, lineage)
    return True, ""


def predict_and_validate(
    fadict,
    contig,
    prfl,
    cutoffs,
    species,
    start,
    end,
    strand,
    configpath,
    blast_score,
):
    """
    Predict and validate protein coding regions using Augustus and HMMER.

    Args:
        fadict (dict): Dictionary containing contig sequences.
        contig (str): Name of the contig.
        prfl (str): Path to the protein profile file.
        cutoffs (dict): Dictionary of cutoff scores for each protein.
        species (str): Species name for Augustus prediction.
        start (int): Start position of the region.
        end (int): End position of the region.
        strand (str): Strand of the region ('+' or '-') for prediction.
        configpath (str): Path to the configuration file for Augustus.
        blast_score (float): BLAST score threshold for validation.

    Returns:
        tuple: A tuple containing the BUSCO name and the final prediction dictionary if a valid prediction is found, otherwise False.
    """
    # augustus no longer accepts fasta from stdin, so write tempfile
    t_fasta = tempfile.NamedTemporaryFile(suffix=".fasta")
    with open(t_fasta.name, "w") as f:
        f.write(">{}\n{}\n".format(contig, fadict[contig][start:end]))

    # run augustus on the regions
    aug_preds = proteinprofile(
        t_fasta.name,
        prfl,
        species=species,
        start=start,
        end=end,
        strand=strand,
        configpath=configpath,
    )
    t_fasta.close()
    busco_name = os.path.basename(prfl).split(".")[0]
    if len(aug_preds) > 0:
        hmmfile = os.path.join(
            os.path.dirname(os.path.dirname(prfl)), "hmms", "{}.hmm".format(busco_name)
        )
        for k, v in aug_preds.items():
            transcript = getSeqRegions(fadict, v["contig"], v["coords"])
            if v["strand"] == "+":
                protein = translate(transcript, v["strand"], v["phase"][0])
            else:
                protein = translate(transcript, v["strand"], v["phase"][-1])
            if protein.startswith("M") and protein.endswith("*"):
                status = "complete"
            else:
                status = "fragmented"
            # now we can check via hmmer
            hmm_result = hmmer_search_single(hmmfile, protein.rstrip("*"))
            if len(hmm_result) > 0:
                if hmm_result[0]["bitscore"] > cutoffs[busco_name]["score"]:
                    final = {
                        "contig": v["contig"],
                        "strand": v["strand"],
                        "location": v["location"],
                        "coords": v["coords"],
                        "phase": v["phase"],
                        "transcript": transcript,
                        "translation": protein,
                        "status": status,
                        "hmmer": hmm_result[0],
                        "miniprot_score": blast_score,
                    }
                    return (busco_name, final)
    return False


def runbusco(
    input,
    lineage,
    mode="genome",
    species="anidulans",
    cpus=1,
    offset=2000,
    verbosity=3,
    logger=False,
    check_augustus=True,
):
    """
    Run BUSCO analysis on genome or protein sequences.

    Parameters:
    - input (str): Path to the input genome or protein sequences.
    - lineage (str): Path to the BUSCO lineage directory.
    - mode (str, optional): Analysis mode, either 'genome' or 'proteins'. Defaults to 'genome'.
    - species (str, optional): Species name for Augustus prediction. Defaults to 'anidulans'.
    - cpus (int, optional): Number of CPUs to use. Defaults to 1.
    - offset (int, optional): Offset value for sequence extraction. Defaults to 2000.
    - verbosity (int, optional): Level of verbosity for logging. Defaults to 3.
    - logger (bool or logger object, optional): Custom logger object. Defaults to False.
    - check_augustus (bool, optional): Flag to check Augustus functionality. Defaults to True.

    Returns:
    - b_final (dict): Final BUSCO results.
    - missing (list): List of BUSCOs not found.
    - stats (dict): Statistics of BUSCO analysis.
    - Config (dict): Parsed configuration details.
    """
    if not logger:
        logger = startLogging()
    # check that the lineage path/files all present
    check, msg = check_lineage(lineage)
    if not check:
        logger.error("Error: {}".format(msg))
        raise SystemExit(1)

    # parse the configs and cutoffs
    Config = load_config(lineage)
    CutOffs = load_cutoffs(lineage)

    if mode == "genome":
        # check augustus functionality
        aug_version = augustus_version()
        if verbosity >= 2:
            logger.info(
                "BUSCOlite v{}; Augustus v{}; miniprot v{}; pyhmmer v{}; pyfastx v{}".format(
                    __version__,
                    aug_version,
                    miniprot_version(),
                    pyhmmer_version(),
                    pyfastx.__version__,
                )
            )
        if check_augustus:
            if not augustus_functional():
                logger.error(
                    "Augustus PPX (--proteinprofile) is non-functional. Usually caused by compilation errors."
                )
                raise SystemExit(1)
        if verbosity >= 1:
            logger.info(
                "{} lineage contains {} BUSCO models".format(
                    Config["name"], len(CutOffs)
                )
            )

        # load genome into dictionary
        seq_records = fasta2dict(input)
        fa_stats = dict2stats(seq_records)
        if verbosity >= 2:
            logger.info(
                "Input genome is {} MB and {} contigs".format(
                    round(fa_stats["size"] / 1e6, 2), fa_stats["n"]
                )
            )

        # run miniprot filter from ancesteral proteins
        query = os.path.join(lineage, "ancestral")
        if verbosity >= 1:
            logger.info(
                "Prefiltering predictions using miniprot of ancestral sequences"
            )
        complete, coords, njobs = miniprot_prefilter(
            input, query, CutOffs, cpus=cpus, buscodb=lineage
        )
        if verbosity >= 1:
            logger.info(
                "Found {} complete models from miniprot, now launching {} augustus/pyhmmer [species={}] jobs for {} BUSCO models".format(
                    len(complete), njobs, species, len(coords)
                )
            )
        # run busco analysis using threadpool, limit io as much as possible
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus + 2) as executor:
            for k, v in coords.items():
                busco_prlf = os.path.join(lineage, "prfl", "{}.prfl".format(k))
                for i in v:
                    contig = i["contig"]
                    start = i["coords"][0] - offset
                    if start < 0:
                        start = 0
                    end = i["coords"][1] + offset
                    if end > len(seq_records[contig]):
                        end = len(seq_records[contig])
                    if i["strand"] == "+":
                        aug_strand = "forward"
                    elif i["strand"] == "-":
                        aug_strand = "reverse"
                    else:
                        aug_strand = "both"
                    # send to executor
                    results.append(
                        executor.submit(
                            predict_and_validate,
                            seq_records,
                            contig,
                            busco_prlf,
                            CutOffs,
                            species,
                            start,
                            end,
                            aug_strand,
                            False,
                            i["score"],
                        )
                    )

        # process mt results
        if len(complete) > 0:
            b_results = complete
        else:
            b_results = {}
        for r in results:
            if isinstance(r.result(), tuple):
                b, res = r.result()
                if b not in b_results:
                    b_results[b] = [res]
                else:
                    b_results[b].append(res)

        # first run is finished, now see which buscos are missing
        missing = []
        for b in CutOffs.keys():
            if b not in b_results:
                missing.append(b)
        if verbosity >= 1:
            logger.info(
                "Found {} BUSCOs in first pass, trying harder to find remaining {}".format(
                    len(b_results), len(missing)
                )
            )
        filt_variants = tempfile.NamedTemporaryFile(suffix=".fasta")
        avariants = fasta2dict(
            os.path.join(lineage, "ancestral_variants"), full_header=True
        )
        seen = set()
        with open(filt_variants.name, "w") as outfile:
            for title, seq in avariants.items():
                try:
                    z, num = title.split(" ")
                except ValueError:
                    if "_" in title:
                        z, num = title.rsplit("_", 1)
                    else:
                        z = title
                if z in missing:
                    seen.add(z)
                    outfile.write(">{} {}\n{}\n".format(z, num, softwrap(seq)))
        if verbosity >= 2:
            logger.info(
                "Trying to use ancestral variants to recover {} BUSCOs".format(
                    len(seen)
                )
            )
        complete2, coords2, njobs2 = miniprot_prefilter(
            input, filt_variants.name, CutOffs, cpus=cpus, buscodb=lineage
        )
        filt_variants.close()
        if verbosity >= 1:
            logger.info(
                "Found {} from miniprot, now launching {} augustus/pyhmmer jobs for {} BUSCO models".format(
                    len(complete2), njobs2, len(coords2)
                )
            )
        if len(coords2) > len(seen):
            logger.error("Something is wrong with parsing the acenstral variants.....")
            print(coords2.keys())
            raise SystemExit(1)

        if len(complete2) > 0:
            for k, v in complete2.items():
                if k not in b_results:
                    b_results[k] = v
                else:
                    b_results[k] += v

        # try thread pool
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus + 2) as executor:
            for k, v in coords2.items():
                busco_prlf = os.path.join(lineage, "prfl", "{}.prfl".format(k))
                for i in v:
                    contig = i["contig"]
                    start = i["coords"][0] - offset
                    if start < 0:
                        start = 0
                    end = i["coords"][1] + offset
                    if end > len(seq_records[contig]):
                        end = len(seq_records[contig])
                    if i["strand"] == "+":
                        aug_strand = "forward"
                    elif i["strand"] == "-":
                        aug_strand = "reverse"
                    else:
                        aug_strand = "both"
                    # send to executor
                    results.append(
                        executor.submit(
                            predict_and_validate,
                            seq_records,
                            contig,
                            busco_prlf,
                            CutOffs,
                            species,
                            start,
                            end,
                            aug_strand,
                            False,
                            i["score"],
                        )
                    )
        # process results; screen for false duplicates here
        for r in results:
            if isinstance(r.result(), tuple):
                b, res = r.result()
                if b not in b_results:
                    b_results[b] = [res]
                else:
                    b_results[b].append(res)

        # finally loop through results, classify and reorganize
        b_final = {}
        missing = []
        for b in CutOffs.keys():
            if b not in b_results:
                missing.append(b)
        stats = {
            "total": 0,
            "single-copy": 0,
            "fragmented": 0,
            "duplicated": 0,
            "missing": len(missing),
        }
        for k, v in natsorted(b_results.items()):
            stats["total"] += 1
            if (
                len(v) > 1
            ):  # duplicate processing here, overlapping regions are not duplicates
                dp = []
                for i, x in enumerate(
                    sorted(v, key=lambda y: y["hmmer"]["bitscore"], reverse=True)
                ):
                    # if first item in list then highest score
                    if i == 0:
                        dp.append(x)
                    else:
                        if x["contig"] in [
                            a["contig"] for a in dp
                        ]:  # means contig already there
                            if not any_overlap(
                                x["location"], [a["location"] for a in dp]
                            ):
                                dp.append(x)
                # these are actually duplicated
                if len(dp) > 1:
                    for y, z in enumerate(dp):
                        z["status"] = "duplicated"
                        if y > 0:
                            name = "{}_{}".format(k, y)
                        else:
                            name = k
                        stats["duplicated"] += 1
                        b_final[name] = z
                else:
                    if dp[0]["status"] == "fragmented":
                        stats["fragmented"] += 1
                    else:
                        stats["single-copy"] += 1
                    b_final[k] = dp[0]
            else:
                if v[0]["status"] == "fragmented":
                    stats["fragmented"] += 1
                else:
                    stats["single-copy"] += 1
                b_final[k] = v[0]
        return b_final, missing, stats, Config

    elif mode == "proteins":
        if verbosity >= 2:
            logger.info(
                "BUSCOlite v{}; pyhmmer v{}".format(__version__, pyhmmer.__version__)
            )
        if verbosity >= 1:
            logger.info(
                "{} lineage contains {} BUSCO models".format(
                    Config["name"], len(CutOffs)
                )
            )
        # load proteome into easel digitized sequence to pass to pyhmmer
        falengths = fasta2lengths(input)
        alphabet = pyhmmer.easel.Alphabet.amino()
        sequences = []
        with pyhmmer.easel.SequenceFile(
            input, digital=True, alphabet=alphabet
        ) as seq_file:
            sequences = list(seq_file)
        if verbosity >= 1:
            logger.info(
                "Loaded {} protein sequences from {}".format(len(sequences), input)
            )
        # now we can loop over the hmms in the lineage and run hmmer on each
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus + 2) as executor:
            for f in os.listdir(os.path.join(lineage, "hmms")):
                if f.endswith(".hmm"):
                    hmmfile = os.path.join(lineage, "hmms", f)
                    results.append(executor.submit(hmmer_search, hmmfile, sequences))
        # process mt results
        b_results = {}
        for r in results:
            if isinstance(r.result(), list):
                for x in r.result():
                    if x["bitscore"] > CutOffs[x["name"]]["score"]:
                        if x["name"] not in b_results:
                            b_results[x["name"]] = [x]
                        else:
                            b_results[x["name"]].append(x)
        b_final = {}
        missing = []
        for b in CutOffs.keys():
            if b not in b_results:
                missing.append(b)
        stats = {
            "total": len(CutOffs),
            "single-copy": 0,
            "fragmented": 0,
            "duplicated": 0,
            "missing": len(missing),
        }
        for k, v in natsorted(b_results.items()):
            if len(v) > 1:  # duplicates
                for i, x in enumerate(
                    sorted(v, key=lambda y: y["bitscore"], reverse=True)
                ):
                    x["status"] = "duplicated"
                    x["length"] = falengths[x["hit"]]
                    if i > 0:
                        name = "{}_{}".format(k, i)
                    else:
                        name = k
                        stats["duplicated"] += 1
                    b_final[name] = x
            else:
                x = v[0]
                x["length"] = falengths[x["hit"]]
                if "status" in x and x["status"] == "fragmented":
                    stats["fragmented"] += 1
                    x["status"] = "fragmented"
                else:
                    x["status"] = "complete"
                    stats["single-copy"] += 1
                b_final[k] = x

        return b_final, missing, stats, Config
