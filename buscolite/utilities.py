import subprocess
import signal
import sys
import os
import errno
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
        handle.write(
            "# Busco id\tStatus\tContig\tStart\tEnd\tHMM_score\tLength\tminiprot_score\n"
        )
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
    return max(
        max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0
    )


def any_overlap(query, lcoords):
    # if any list of coords overlaps with query
    v = False
    for x in lcoords:
        if overlap(query[0], query[1], x[0], x[1]) > 0:
            v = True
    return v


def runprocess(cmd, stdout=False, stderr=False, cwd=".", debug=False):
    if not stdout and not stderr:
        proc = subprocess.Popen(
            cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    elif stdout and not stderr:
        with open(stdout, "w") as outfile:
            proc = subprocess.Popen(
                cmd, cwd=cwd, stdout=outfile, stderr=subprocess.PIPE
            )
    elif not stdout and stderr:
        with open(stderr, "w") as outfile:
            proc = subprocess.Popen(
                cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=outfile
            )
    elif stdout and stderr:
        if stdout == stderr:
            with open(stdout, "w") as outfile:
                proc = subprocess.Popen(cmd, cwd=cwd, stdout=outfile, stderr=outfile)
        else:
            with open(stdout, "w") as outfile1:
                with open(stderr, "w") as outfile2:
                    proc = subprocess.Popen(
                        cmd, cwd=cwd, stdout=outfile1, stderr=outfile2
                    )
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
    popen = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, universal_newlines=True, stderr=DEVNULL
    )
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def execute_timeout(cmd, timeout=120):
    # stream execute but add a timeout
    DEVNULL = open(os.devnull, "w")
    try:
        p = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, universal_newlines=True, stderr=DEVNULL
        )
        return_code = p.wait(timeout=timeout)
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)
        else:
            for stdout_line in iter(p.stdout.readline, ""):
                yield stdout_line
            p.stdout.close()
    except subprocess.TimeoutExpired:
        # print(f"Timeout for {cmd} ({timeout}s) expired", file=sys.stderr)
        p.terminate()
        return ""


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

    fpath, fname = os.path.split(program)
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
    import subprocess
    import signal

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


def open_xz(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
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
    return None
