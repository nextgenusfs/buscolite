import errno
import os

# import signal  # Used in open_pipe function
import subprocess
import sys

from natsort import natsorted

from .__init__ import __version__


def summary_writer(result, missing, cmd, cfg, handle, mode="genome"):
    for x in missing:
        result[x] = {"status": "missing"}
    handle.write("# BUSCOlite v{}\n".format(__version__))
    handle.write(
        "# The lineage dataset is: {} (Creation date: {}, number of species: {}, number of BUSCOs: {})\n".format(
            cfg["name"],
            cfg["creation_date"],
            cfg["number_of_species"],
            cfg["number_of_BUSCOs"],
        )
    )
    handle.write("# To reproduce this run: {}\n".format(" ".join(cmd)))
    handle.write("#\n")
    if mode == "genome":
        handle.write("# Busco id\tStatus\tContig\tStart\tEnd\tHMM_score\tLength\tminiprot_score\n")
        for k, v in natsorted(result.items()):
            if "_" in k:
                busco = k.rsplit("_", 1)[0]
            else:
                busco = k
            if v["status"] == "missing":
                handle.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        busco,
                        v["status"],
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                    )
                )
            else:
                handle.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        busco,
                        v["status"],
                        v["contig"],
                        min(v["location"]),
                        max(v["location"]),
                        round(v["hmmer"]["bitscore"], 2),
                        len(v["translation"].rstrip("*")),
                        v["miniprot_score"],
                    )
                )

    elif mode == "proteins":
        handle.write("# Busco id\tStatus\tSequence\tScore\tLength\n")
        for k, v in result.items():
            if "_" in k:
                busco = k.rsplit("_", 1)[0]
            else:
                busco = k
            if v["status"] == "missing":
                handle.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        busco, v["status"], "", "", "", "", "", "", ""
                    )
                )
            else:
                handle.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(
                        busco,
                        v["status"],
                        v["hit"],
                        round(v["bitscore"], 2),
                        v["length"],
                    )
                )


def overlap(start1, end1, start2, end2):
    """how much does the range (start1, end1) overlap with (start2, end2)"""
    return max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)


def any_overlap(query, lcoords):
    # if any list of coords overlaps with query
    v = False
    for x in lcoords:
        if overlap(query[0], query[1], x[0], x[1]) > 0:
            v = True
    return v


def filter_low_scoring_matches(busco_results, threshold=0.85):
    """
    Filter out matches that score less than a threshold percentage of the top bitscore
    for each BUSCO. This implements the BUSCO v6 filtering logic.

    Parameters
    ----------
    busco_results : dict
        Dictionary of BUSCO results where keys are BUSCO IDs and values are lists
        of match dictionaries containing 'bitscore' keys
    threshold : float, optional
        Minimum score threshold as fraction of top score (default: 0.85)

    Returns
    -------
    dict
        Filtered dictionary with low-scoring matches removed
    """
    filtered = {}
    for busco_id, matches in busco_results.items():
        if len(matches) == 0:
            continue
        elif len(matches) == 1:
            filtered[busco_id] = matches
        else:
            # Find the maximum bitscore for this BUSCO
            max_bitscore = max(m["bitscore"] for m in matches)
            cutoff = threshold * max_bitscore
            # Keep only matches above the cutoff
            kept_matches = [m for m in matches if m["bitscore"] >= cutoff]
            if len(kept_matches) > 0:
                filtered[busco_id] = kept_matches
    return filtered


def remove_duplicate_gene_matches(busco_results, score_key="bitscore"):
    """
    When the same gene/sequence matches multiple BUSCOs, keep only the highest
    scoring match. This prevents a single gene from being counted multiple times.

    Parameters
    ----------
    busco_results : dict
        Dictionary of BUSCO results where keys are BUSCO IDs and values are lists
        of match dictionaries
    score_key : str, optional
        Key to use for scoring. Can be a simple key like 'bitscore' or a nested
        key like 'hmmer.bitscore' (default: 'bitscore')

    Returns
    -------
    dict
        Filtered dictionary with duplicate gene matches removed
    """
    # Build a reverse mapping: gene_id -> list of (busco_id, match_dict, score)
    gene_to_buscos = {}
    for busco_id, matches in busco_results.items():
        for match in matches:
            # For protein mode, use 'hit' as gene identifier
            # For genome mode, use combination of contig and location
            if "hit" in match:
                gene_id = match["hit"]
            elif "contig" in match and "location" in match:
                gene_id = "{}_{}:{}".format(
                    match["contig"], match["location"][0], match["location"][1]
                )
            else:
                continue

            # Handle nested score keys (e.g., 'hmmer.bitscore')
            if "." in score_key:
                keys = score_key.split(".")
                score = match
                for key in keys:
                    score = score.get(key, {})
                if not isinstance(score, (int, float)):
                    score = 0
            else:
                score = match.get(score_key, 0)

            if gene_id not in gene_to_buscos:
                gene_to_buscos[gene_id] = []
            gene_to_buscos[gene_id].append((busco_id, match, score))

    # Find genes that match multiple BUSCOs
    genes_to_remove = {}
    for gene_id, busco_matches in gene_to_buscos.items():
        if len(busco_matches) > 1:
            # Sort by score descending
            busco_matches.sort(key=lambda x: x[2], reverse=True)
            # Keep the best match, mark others for removal
            for busco_id, match, score in busco_matches[1:]:
                if busco_id not in genes_to_remove:
                    genes_to_remove[busco_id] = []
                genes_to_remove[busco_id].append(gene_id)

    # Remove the duplicate matches
    filtered = {}
    for busco_id, matches in busco_results.items():
        if busco_id in genes_to_remove:
            genes_to_skip = set(genes_to_remove[busco_id])
            kept_matches = []
            for match in matches:
                if "hit" in match:
                    gene_id = match["hit"]
                elif "contig" in match and "location" in match:
                    gene_id = "{}_{}:{}".format(
                        match["contig"], match["location"][0], match["location"][1]
                    )
                else:
                    gene_id = None

                if gene_id not in genes_to_skip:
                    kept_matches.append(match)

            if len(kept_matches) > 0:
                filtered[busco_id] = kept_matches
        else:
            filtered[busco_id] = matches

    return filtered


def runprocess(cmd, stdout=False, stderr=False, cwd=".", debug=False):
    if not stdout and not stderr:
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif stdout and not stderr:
        with open(stdout, "w") as outfile:
            proc = subprocess.Popen(cmd, cwd=cwd, stdout=outfile, stderr=subprocess.PIPE)
    elif not stdout and stderr:
        with open(stderr, "w") as outfile:
            proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=outfile)
    elif stdout and stderr:
        if stdout == stderr:
            with open(stdout, "w") as outfile:
                proc = subprocess.Popen(cmd, cwd=cwd, stdout=outfile, stderr=outfile)
        else:
            with open(stdout, "w") as outfile1:
                with open(stderr, "w") as outfile2:
                    proc = subprocess.Popen(cmd, cwd=cwd, stdout=outfile1, stderr=outfile2)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        sys.stderr.write("CMD ERROR: {}".format(" ".join(cmd)))
        if stdout:
            sys.stderr.write(stdout.decode("utf-8"))
        if stderr:
            sys.stderr.write(stderr.decode("utf-8"))
        sys.exit(1)
    if debug:
        if stdout:
            sys.stderr.write(stdout.decode("utf-8"))
        if stderr:
            sys.stderr.write(stderr.decode("utf-8"))


def execute(cmd):
    DEVNULL = open(os.devnull, "w")
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, stderr=DEVNULL)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def execute_timeout(cmd, timeout=120):
    # stream execute but add a timeout; capture stderr so callers can see
    # why a subprocess failed instead of discarding it to /dev/null
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    try:
        stdout_data, stderr_data = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.terminate()
        try:
            p.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
        return
    if p.returncode:
        raise subprocess.CalledProcessError(
            p.returncode, cmd, output=stdout_data, stderr=stderr_data
        )
    for stdout_line in stdout_data.splitlines(keepends=True):
        yield stdout_line


def check_inputs(inputs):
    for filename in inputs:
        if not is_file(filename):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)


def is_file(f):
    if os.path.isfile(f):
        return True
    else:
        return False


def which2(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)  # fname is not used
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def open_pipe(command, mode="r", buff=1024 * 1024):
    import signal  # Import signal locally where it's used

    # import subprocess  # Already imported at the module level

    if "r" in mode:
        return subprocess.Popen(
            command,
            shell=True,
            bufsize=buff,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL),
        ).stdout
    elif "w" in mode:
        return subprocess.Popen(
            command,
            shell=True,
            bufsize=buff,
            universal_newlines=True,
            stdin=subprocess.PIPE,
        ).stdin
    return None


NORMAL = 0
PROCESS = 1
PARALLEL = 2

WHICH_BZIP2 = which2("bzip2")
WHICH_PBZIP2 = which2("pbzip2")


def open_bz2(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import bz2

        return bz2.BZ2File(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_BZIP2:
            return open_bz2(filename, mode, buff, NORMAL)
        if "r" in mode:
            return open_pipe("bzip2 -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("bzip2 >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PBZIP2:
            return open_bz2(filename, mode, buff, PROCESS)
        if "r" in mode:
            return open_pipe("pbzip2 -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("pbzip2 >" + filename, mode, buff)
    return None


WHICH_GZIP = which2("gzip")
WHICH_PIGZ = which2("pigz")


def open_gz(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import gzip

        return gzip.GzipFile(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_GZIP:
            return open_gz(filename, mode, buff, NORMAL)
        if "r" in mode:
            return open_pipe("gzip -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("gzip >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PIGZ:
            return open_gz(filename, mode, buff, PROCESS)
        if "r" in mode:
            return open_pipe("pigz -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("pigz >" + filename, mode, buff)
    return None


WHICH_XZ = which2("xz")


def open_xz(
    filename, mode="r", buff=1024 * 1024, external=PARALLEL
):  # external is kept for API compatibility
    if WHICH_XZ:
        if "r" in mode:
            return open_pipe("xz -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("xz >" + filename, mode, buff)
    return None


def zopen(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    """
    Open pipe, zipped, or unzipped file automagically

    # external == 0: normal zip libraries
    # external == 1: (zcat, gzip) or (bzcat, bzip2)
    # external == 2: (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)
    """
    if "r" in mode and "w" in mode:
        return None
    if filename.startswith("!"):
        return open_pipe(filename[1:], mode, buff)
    elif filename.endswith(".bz2"):
        return open_bz2(filename, mode, buff, external)
    elif filename.endswith(".gz"):
        return open_gz(filename, mode, buff, external)
    elif filename.endswith(".xz"):
        return open_xz(filename, mode, buff, external)
    else:
        return open(filename, mode, buff)
    # This line is unreachable, but kept for consistency with other functions
    # return None
