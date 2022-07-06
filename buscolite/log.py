import logging

GREEN = "\033[92m"
END = "\033[0m"


def startLogging(logfile=False):
    stdoutformat = logging.Formatter(
        GREEN + "%(asctime)s" + END + ": %(message)s", datefmt="[%b %d %I:%M %p]"
    )
    fileformat = logging.Formatter("%(asctime)s: %(message)s", datefmt="[%x %H:%M:%S]")
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    logger.addHandler(sth)
    if logfile:
        fhnd = logging.FileHandler(logfile)
        fhnd.setLevel(logging.DEBUG)
        fhnd.setFormatter(fileformat)
        logger.addHandler(fhnd)
    return logger
