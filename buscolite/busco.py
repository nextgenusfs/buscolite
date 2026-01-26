import concurrent.futures
import os
import tempfile

import pyfastx
import pyhmmer
from natsort import natsorted

from .__init__ import __version__
from .augustus import augustus_functional, augustus_version, proteinprofile
from .fastx import dict2stats, fasta2dict, fasta2lengths, getSeqRegions, softwrap, translate
from .log import startLogging
from .search import (
    hmmer_search,
    hmmer_search_single,
    miniprot_prefilter,
    miniprot_version,
    pyhmmer_version,
)
from .utilities import any_overlap, filter_low_scoring_matches, remove_duplicate_gene_matches


def load_config(lineage):
    """
    Load the BUSCO dataset configuration file.

    Parameters
    ----------
    lineage : str
        Path to the BUSCO lineage directory containing the dataset.cfg file

    Returns
    -------
    dict
        Dictionary containing the configuration parameters from dataset.cfg
    """
    config = {}
    with open(os.path.join(lineage, "dataset.cfg"), "r") as infile:
        for line in infile:
            key, value = line.rstrip().split("=")
            config[key] = value
    return config


def load_cutoffs(lineage):
    """
    Load the BUSCO score and length cutoffs from the lineage directory.

    Parameters
    ----------
    lineage : str
        Path to the BUSCO lineage directory containing the cutoff files

    Returns
    -------
    dict
        Dictionary containing the score and length cutoffs for each BUSCO model
        Format: {busco_id: {"score": float, "sigma": float, "length": int}}
    """
    cutoffs = {}
    with open(os.path.join(lineage, "scores_cutoff"), "r") as infile:
        for line in infile:
            busco, score = line.rstrip().split("\t")
            if busco not in cutoffs:
                cutoffs[busco] = {"score": float(score)}
            else:
                cutoffs[busco]["score"] = float(score)
    # this file is not there anymore in odb12
    length_cutoffs = os.path.join(lineage, "lengths_cutoff")
    if os.path.isfile(length_cutoffs):
        with open(length_cutoffs, "r") as infile:
            for line in infile:
                busco, _, sigma, length = line.rstrip().split("\t")
                if float(sigma) == 0.0:
                    sigma = 1
                if busco not in cutoffs:
                    cutoffs[busco] = {
                        "sigma": float(sigma),
                        "length": int(float(length)),
                    }
                else:
                    cutoffs[busco]["sigma"] = float(sigma)
                    cutoffs[busco]["length"] = int(float(length))
    return cutoffs


def check_lineage(lineage):
    """
    Verify that the BUSCO lineage directory contains all required files and directories.

    Parameters
    ----------
    lineage : str
        Path to the BUSCO lineage directory

    Returns
    -------
    tuple
        (bool, str) - Boolean indicating if the lineage is valid, and an error message if not
    """
    lineage = os.path.abspath(lineage)
    if not os.path.isdir(lineage):
        return False, "{} is not a directory".format(lineage)
    dirs = ["prfl", "hmms"]
    files = [
        "ancestral",
        "ancestral_variants",
        "dataset.cfg",
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
            logger.info("{} lineage contains {} BUSCO models".format(Config["name"], len(CutOffs)))

        # load genome into dictionary
        seq_records = fasta2dict(input)
        fa_stats = dict2stats(seq_records)
        if verbosity >= 2:
            logger.info(
                "Input genome is {} MB and {} contigs".format(
                    round(fa_stats["size"] / 1e6, 2), fa_stats["n"]
                )
            )

        # run miniprot filter from ancestral proteins + variants (combined for efficiency)
        # We need to track which hits come from ancestral vs variants for prioritization
        # Priority order: ancestral complete > ancestral incomplete > variant incomplete

        # First, run miniprot on ancestral sequences
        ancestral = os.path.join(lineage, "ancestral")
        if verbosity >= 1:
            logger.info("Prefiltering predictions using miniprot of ancestral sequences")
        complete_ancestral, coords_ancestral, njobs_ancestral = miniprot_prefilter(
            input, ancestral, CutOffs, cpus=cpus, buscodb=lineage
        )

        # Tag ancestral coords with priority=1 (highest)
        for busco_id in coords_ancestral:
            for region in coords_ancestral[busco_id]:
                region["priority"] = 1
                region["source"] = "ancestral"

        # Second, run miniprot on ancestral_variants for BUSCOs not yet complete
        ancestral_variants = os.path.join(lineage, "ancestral_variants")
        complete_variants = {}
        coords_variants = {}
        njobs_variants = 0

        if os.path.isfile(ancestral_variants):
            # Only search for BUSCOs not already complete from ancestral
            missing = [b for b in CutOffs.keys() if b not in complete_ancestral]

            if len(missing) > 0:
                # Create filtered variants file with only missing BUSCOs
                filt_variants = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
                avariants = fasta2dict(ancestral_variants, full_header=True)
                seen = set()
                for title, seq in avariants.items():
                    # Parse the BUSCO ID from the variant header
                    try:
                        busco_id, variant_num = title.split(" ")
                    except ValueError:
                        if "_" in title:
                            busco_id, variant_num = title.rsplit("_", 1)
                        else:
                            busco_id = title
                            variant_num = "0"
                    if busco_id in missing:
                        seen.add(busco_id)
                        # Write with normalized header
                        filt_variants.write(
                            ">{} {}\n{}\n".format(busco_id, variant_num, softwrap(seq))
                        )
                filt_variants.close()

                if verbosity >= 2:
                    logger.info(
                        "Trying ancestral variants for {} BUSCOs not complete from ancestral".format(
                            len(seen)
                        )
                    )

                complete_variants, coords_variants, njobs_variants = miniprot_prefilter(
                    input, filt_variants.name, CutOffs, cpus=cpus, buscodb=lineage
                )

                # Clean up temp file
                os.unlink(filt_variants.name)

                # Tag variant coords with priority=2 (lower than ancestral)
                for busco_id in coords_variants:
                    for region in coords_variants[busco_id]:
                        region["priority"] = 2
                        region["source"] = "ancestral_variants"

        # Combine complete models from both sources
        complete = complete_ancestral.copy()
        for busco_id, models in complete_variants.items():
            if busco_id not in complete:
                complete[busco_id] = models
            else:
                complete[busco_id].extend(models)

        # Combine coords, merging and sorting by priority
        coords = {}
        for busco_id in set(list(coords_ancestral.keys()) + list(coords_variants.keys())):
            coords[busco_id] = []
            if busco_id in coords_ancestral:
                coords[busco_id].extend(coords_ancestral[busco_id])
            if busco_id in coords_variants:
                coords[busco_id].extend(coords_variants[busco_id])
            # Sort by priority (1=ancestral first, 2=variants second)
            coords[busco_id].sort(key=lambda x: x["priority"])

        njobs = njobs_ancestral + njobs_variants

        if verbosity >= 1:
            logger.info(
                "Found {} complete models from miniprot, now launching {} augustus/pyhmmer [species={}] jobs for {} BUSCO models".format(
                    len(complete), njobs, species, len(coords)
                )
            )

        # Helper function to process all jobs for a single BUSCO with early exit
        def process_busco_jobs(busco_id, job_list, seq_records, cutoffs, species):
            """
            Process all jobs for a single BUSCO in priority order.
            Returns early if a complete model is found.

            Returns:
                tuple: (busco_id, results_list, jobs_submitted, jobs_skipped)
            """
            results = []
            jobs_submitted = 0
            jobs_skipped = 0

            for job in job_list:
                # Run Augustus + HMMER validation
                result = predict_and_validate(
                    seq_records,
                    job["contig"],
                    job["prfl"],
                    cutoffs,
                    species,
                    job["start"],
                    job["end"],
                    job["strand"],
                    False,
                    job["score"],
                )
                jobs_submitted += 1

                if isinstance(result, tuple):
                    b, res = result
                    results.append(res)

                    # Early exit: if we found a COMPLETE model, skip remaining jobs
                    if res.get("status") == "complete":
                        remaining_jobs = len(job_list) - (job_list.index(job) + 1)
                        jobs_skipped = remaining_jobs
                        break  # Exit early!

            return (busco_id, results, jobs_submitted, jobs_skipped)

        # run busco analysis using threadpool with early exit optimization
        # Process multiple BUSCOs in parallel, each with its own early exit logic
        if len(complete) > 0:
            b_results = complete
        else:
            b_results = {}

        # Track which BUSCOs are already complete from miniprot
        found_complete = set(complete.keys())

        # Group jobs by BUSCO ID and priority for early exit
        jobs_by_busco = {}
        for k, v in coords.items():
            busco_prlf = os.path.join(lineage, "prfl", "{}.prfl".format(k))
            if not os.path.isfile(busco_prlf):
                continue
            jobs_by_busco[k] = []
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
                jobs_by_busco[k].append(
                    {
                        "contig": contig,
                        "prfl": busco_prlf,
                        "start": start,
                        "end": end,
                        "strand": aug_strand,
                        "score": i["score"],
                        "priority": i["priority"],
                    }
                )

        # Submit one task per BUSCO (each task processes multiple jobs with early exit)
        # This allows parallel processing across BUSCOs while maintaining early exit within each
        total_jobs_submitted = 0
        total_jobs_skipped = 0
        futures = []

        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
            for busco_id, job_list in jobs_by_busco.items():
                # Skip if already complete from miniprot
                if busco_id in found_complete:
                    total_jobs_skipped += len(job_list)
                    continue

                # Submit this BUSCO's job processor to the executor
                future = executor.submit(
                    process_busco_jobs,
                    busco_id,
                    job_list,
                    seq_records,
                    CutOffs,
                    species,
                )
                futures.append(future)

            # Collect results as they complete
            for future in concurrent.futures.as_completed(futures):
                busco_id, results, jobs_submitted, jobs_skipped = future.result()
                total_jobs_submitted += jobs_submitted
                total_jobs_skipped += jobs_skipped

                # Add results to b_results
                if busco_id not in b_results:
                    b_results[busco_id] = results
                else:
                    b_results[busco_id].extend(results)

        if verbosity >= 2:
            logger.info(
                "Augustus completed: {} jobs submitted, {} jobs skipped via early exit".format(
                    total_jobs_submitted, total_jobs_skipped
                )
            )

        # Apply BUSCO v6-style filtering
        # NOTE: In genome mode, we skip the duplicate gene filter because:
        # 1. We already handle overlapping genes later (line 584)
        # 2. Using gene span for duplicate detection is too strict - different BUSCOs
        #    can legitimately predict overlapping genes
        # 3. BUSCO v6 uses exon-level overlap checking, which is more sophisticated
        if verbosity >= 2:
            logger.info(
                "Found {} matches from Augustus/HMMER".format(
                    sum(len(v) for v in b_results.values())
                )
            )

        # Step 2: Remove matches scoring less than 85% of top bitscore for each BUSCO
        # NOTE: Temporarily disabled - this filter may be too aggressive for genome mode
        # TODO: Investigate if BUSCO v6 actually applies this in genome mode with Augustus
        # For genome mode, the bitscore is in the nested 'hmmer' dict
        # Keep reference to original matches before flattening
        # b_results_original = {}
        # b_results_for_filter = {}
        # for busco_id, matches in b_results.items():
        #     b_results_original[busco_id] = matches  # Keep original
        #     b_results_for_filter[busco_id] = []
        #     for match in matches:
        #         # Create a flattened version for filtering
        #         flat_match = match.copy()
        #         if "hmmer" in match:
        #             flat_match["bitscore"] = match["hmmer"]["bitscore"]
        #         b_results_for_filter[busco_id].append(flat_match)

        # filtered = filter_low_scoring_matches(b_results_for_filter, threshold=0.85)

        # # Map back to original structure
        # b_results = {}
        # for busco_id, matches in filtered.items():
        #     b_results[busco_id] = []
        #     for match in matches:
        #         # Find the corresponding original match using object identity
        #         # The filtered matches are the same objects from b_results_for_filter
        #         for i, flat_match in enumerate(b_results_for_filter[busco_id]):
        #             if flat_match is match:
        #                 # Use the original match from b_results_original
        #                 b_results[busco_id].append(b_results_original[busco_id][i])
        #                 break

        # matches_after_85_filter = sum(len(v) for v in b_results.values())
        # if verbosity >= 2:
        #     logger.info(
        #         "After 85% bitscore filter: {} matches remain ({} removed)".format(
        #             matches_after_85_filter, matches_after_dup_filter - matches_after_85_filter
        #         )
        #     )

        # finally loop through results, classify and reorganize
        b_final = {}
        missing = []
        for b in CutOffs.keys():
            # Count as missing if not in results OR if results list is empty (filtered out)
            if b not in b_results or len(b_results[b]) == 0:
                missing.append(b)
        stats = {
            "total": len(CutOffs),
            "single-copy": 0,
            "fragmented": 0,
            "duplicated": 0,
            "missing": len(missing),
        }
        for k, v in natsorted(b_results.items()):
            # Skip empty lists (can happen after filtering)
            if len(v) == 0:
                continue
            elif len(v) > 1:  # duplicate processing here, overlapping regions are not duplicates
                dp = []
                for i, x in enumerate(
                    sorted(v, key=lambda y: y["hmmer"]["bitscore"], reverse=True)
                ):
                    # if first item in list then highest score
                    if i == 0:
                        dp.append(x)
                    else:
                        if x["contig"] in [a["contig"] for a in dp]:  # means contig already there
                            if not any_overlap(x["location"], [a["location"] for a in dp]):
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
                elif len(dp) == 1:
                    if dp[0]["status"] == "fragmented":
                        stats["fragmented"] += 1
                    else:
                        stats["single-copy"] += 1
                    b_final[k] = dp[0]
                # else: dp is empty, skip this BUSCO
            else:  # len(v) == 1
                if v[0]["status"] == "fragmented":
                    stats["fragmented"] += 1
                else:
                    stats["single-copy"] += 1
                b_final[k] = v[0]
        return b_final, missing, stats, Config

    elif mode == "proteins":
        if verbosity >= 2:
            logger.info("BUSCOlite v{}; pyhmmer v{}".format(__version__, pyhmmer.__version__))
        if verbosity >= 1:
            logger.info("{} lineage contains {} BUSCO models".format(Config["name"], len(CutOffs)))
        # load proteome into easel digitized sequence to pass to pyhmmer
        falengths = fasta2lengths(input)
        alphabet = pyhmmer.easel.Alphabet.amino()
        sequences = []
        with pyhmmer.easel.SequenceFile(input, digital=True, alphabet=alphabet) as seq_file:
            sequences = list(seq_file)
        if verbosity >= 1:
            logger.info("Loaded {} protein sequences from {}".format(len(sequences), input))
        # now we can loop over the hmms in the lineage and run hmmer on each
        results = []
        hmmfilelookup = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus + 2) as executor:
            for f in os.listdir(os.path.join(lineage, "hmms")):
                if f.endswith(".hmm"):
                    if "at" in f:
                        hmmfilelookup[f.split("at")[0]] = f.split(".hmm")[0]
                    hmmfile = os.path.join(lineage, "hmms", f)
                    results.append(executor.submit(hmmer_search, hmmfile, sequences))
        # process mt results
        # danger with odb12, the hmm name does not match the file name in scores_cutoffs.... #dumb
        b_results = {}
        for r in results:
            if isinstance(r.result(), list):
                for x in r.result():
                    if x["name"] not in CutOffs:
                        if x["name"] in hmmfilelookup:
                            modelName = hmmfilelookup[x["name"]]
                        else:
                            modelName = x["name"]
                    else:
                        modelName = x["name"]
                    if x["bitscore"] > CutOffs[modelName]["score"]:
                        if modelName not in b_results:
                            b_results[modelName] = [x]
                        else:
                            b_results[modelName].append(x)

        # Apply BUSCO v6-style filtering
        # Step 1: Remove matches scoring less than 85% of top bitscore for each BUSCO
        if verbosity >= 2:
            logger.info(
                "Filtering {} initial matches using 85% bitscore threshold".format(
                    sum(len(v) for v in b_results.values())
                )
            )
        b_results = filter_low_scoring_matches(b_results, threshold=0.85)

        # Step 2: Remove duplicate gene matches (same gene matching multiple BUSCOs)
        b_results = remove_duplicate_gene_matches(b_results, score_key="bitscore")
        if verbosity >= 2:
            logger.info(
                "After filtering: {} matches remain".format(sum(len(v) for v in b_results.values()))
            )

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
                for i, x in enumerate(sorted(v, key=lambda y: y["bitscore"], reverse=True)):
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
