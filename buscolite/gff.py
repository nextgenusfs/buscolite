from urllib.parse import unquote
from collections import OrderedDict
import re
import sys
import concurrent.futures
from natsort import natsorted
from .fastx import translate, codon_table, getSeqRegions, RevComp


def gffwriter(result, handle):
    cleaned = {}
    for k, v in result.items():
        if v["status"] != "missing":
            cleaned[k] = v
    # sort results by contig and location
    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    sGenes = natsorted(iter(cleaned.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # write GFF3 output from runbusco dictionary result
    handle.write("##gff-version 3\n")
    for k, v in sortedGenes.items():
        extraAnnots = (
            "Note=status={},hmmer_score={},hmmer_evalue={},miniprot_score={},".format(
                v["status"],
                round(v["hmmer"]["bitscore"], 2),
                round(v["hmmer"]["evalue"], 2),
                v["miniprot_score"],
            )
        )

        handle.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tID={};{}\n".format(
                v["contig"],
                "buscolite",
                "gene",
                min(v["location"]),
                max(v["location"]),
                round(v["hmmer"]["bitscore"], 2),
                v["strand"],
                ".",
                k,
                extraAnnots,
            )
        )
        handle.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tID={}-T1;Parent={}\n".format(
                v["contig"],
                "buscolite",
                "mRNA",
                min(v["location"]),
                max(v["location"]),
                round(v["hmmer"]["bitscore"], 2),
                v["strand"],
                ".",
                k,
                k,
            )
        )
        if v["strand"] == "+":
            sortedCDS = sorted(v["coords"], key=lambda tup: tup[0])
            current_phase = v["phase"][0]
        else:
            sortedCDS = sorted(v["coords"], key=lambda tup: tup[0], reverse=True)
            current_phase = v["phase"][-1]
        num_cds = len(sortedCDS)
        for y in range(0, num_cds):
            handle.write(
                "{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.exon{:};Parent={:}-T1;\n".format(
                    v["contig"],
                    "buscolite",
                    sortedCDS[y][0],
                    sortedCDS[y][1],
                    v["strand"],
                    current_phase,
                    k,
                    y + 1,
                    k,
                )
            )
            handle.write(
                "{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:}-T1;\n".format(
                    v["contig"],
                    "buscolite",
                    sortedCDS[y][0],
                    sortedCDS[y][1],
                    v["strand"],
                    current_phase,
                    k,
                    k,
                )
            )
            current_phase = (
                current_phase - (int(sortedCDS[y][1]) - int(sortedCDS[y][0]) + 1)
            ) % 3
            if current_phase == 3:
                current_phase = 0


def miniprot_gff_parser(line, Genes):
    if line.startswith("\n") or line.startswith("#"):
        return Genes
    line = line.rstrip()
    # skip lines that aren't 9 columns
    if not line.count("\t") == 8:
        return Genes
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
        return Genes
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
    return Genes


def validate_models(
    annotation, fadict, logger=sys.stderr.write, table=1, gap_filter=False
):
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    # try to use multithreading here, not sure its necessary but new syntax for me with concurrent
    results = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for k, v in list(annotation.items()):
            results.append(
                executor.submit(
                    validate_and_translate_models,
                    k,
                    v,
                    fadict,
                    {"gap_filter": gap_filter, "table": table, "logger": logger},
                )
            )
            # check, update = validate_and_translate_models(v, SeqRecords, gap_filter=gap_filter, table=table, logger=logger)
    # pass the updates back to the main dictionary
    for r in results:
        gene, update = r.result()
        for key, value in update.items():
            annotation[gene][key] = value
        # some assertion statements here to ensure parsing is correct
        assert_lengths_fail = []
        for z in [
            "type",
            "mRNA",
            "CDS",
            "codon_start",
            "phase",
            "5UTR",
            "3UTR",
            "protein",
            "transcript",
            "cds_transcript",
            "partialStart",
            "partialStop",
            "product",
        ]:
            if len(annotation[gene]["ids"]) != len(annotation[gene][z]):
                assert_lengths_fail.append(
                    (z, annotation[gene][z], len(annotation[gene][z]))
                )
        if len(assert_lengths_fail) > 0:
            logger(
                "ERROR in parsing gene {}\n{}\n{}\n".format(
                    gene, assert_lengths_fail, annotation[gene]
                )
            )
            raise SystemExit(1)
    return annotation


def validate_and_translate_models(
    k, v, SeqRecords, gap_filter=False, table=1, logger=sys.stderr.write
):
    # take the Genes dictionary and validate gene models
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    # return sorted mRNA, sorted CDS, transcripts, translations, proper phase
    results = {
        "gene_synonym": [],
        "location": v["location"],
        "mRNA": [],
        "CDS": [],
        "protein": [],
        "transcript": [],
        "cds_transcript": [],
        "codon_start": [],
        "partialStart": [],
        "partialStop": [],
        "type": [],
    }
    assert len(v["ids"]) == len(v["type"])
    for i in range(0, len(v["ids"])):
        if v["type"][i] in ["mRNA", "tRNA", "ncRNA", "rRNA", "transcript"]:
            if v["strand"] == "+":
                sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0], reverse=True)
            mrnaSeq = getSeqRegions(SeqRecords, v["contig"], sortedExons)
            results["mRNA"].append(sortedExons)
            results["transcript"].append(mrnaSeq)
        if v["type"][i] in ["mRNA", "transcript"]:
            if not v["CDS"][i]:
                results["CDS"].append([])
                results["type"].append("ncRNA")
                results["protein"].append(None)
                results["cds_transcript"].append(None)
                results["codon_start"].append(None)
                results["partialStart"].append(None)
                results["partialStop"].append(None)
            else:
                results["type"].append("mRNA")
                if v["strand"] == "+":
                    sortedCDS = sorted(v["CDS"][i], key=lambda tup: tup[0])
                else:
                    sortedCDS = sorted(
                        v["CDS"][i], key=lambda tup: tup[0], reverse=True
                    )
                # get the codon_start by getting first CDS phase + 1
                indexStart = [
                    x for x, y in enumerate(v["CDS"][i]) if y[0] == sortedCDS[0][0]
                ]
                cdsSeq = getSeqRegions(SeqRecords, v["contig"], sortedCDS)
                protSeq, codon_start = (None,) * 2
                if (
                    "?" in v["phase"][i]
                ):  # dont know the phase -- malformed GFF3, try to find best CDS
                    translateResults = []
                    for y in [1, 2, 3]:
                        protSeq = translate(cdsSeq, v["strand"], y - 1, table=table)
                        numStops = protSeq.count("*")
                        if protSeq[-1] == "*":
                            numStops -= 1
                        translateResults.append((y, numStops, protSeq))
                    sortedResults = sorted(translateResults, key=lambda tup: tup[1])
                    codon_start = sortedResults[0][0]
                    protSeq = sortedResults[0][2]
                    v["phase"][i] = codon_start - 1
                else:
                    try:
                        codon_start = int(v["phase"][i][indexStart[0]]) + 1
                    except IndexError:
                        pass
                    # translate and get protein sequence
                    protSeq = translate(
                        cdsSeq, v["strand"], codon_start - 1, table=table
                    )
                results["codon_start"].append(codon_start)
                if codon_start > 1:
                    if v["strand"] == "+":
                        cdsSeq = cdsSeq[codon_start - 1 :]
                    elif v["strand"] == "-":
                        endTrunc = len(cdsSeq) - codon_start - 1
                        cdsSeq = cdsSeq[0:endTrunc]
                    else:
                        logger(
                            "ERROR nonsensical strand ({}) for gene {}\n".format(
                                v["strand"], v["ids"][i]
                            )
                        )
                results["cds_transcript"].append(cdsSeq)
                results["CDS"].append(sortedCDS)
                results["protein"].append(protSeq)
                if protSeq:
                    if protSeq.endswith("*"):
                        results["partialStop"].append(False)
                    else:
                        results["partialStop"].append(True)
                    if codon_start == 1 and protSeq.startswith("M"):
                        results["partialStart"].append(False)
                    else:
                        results["partialStart"].append(True)
                    if protSeq.rstrip("*").count("*") > 0:
                        results["pseudo"] = True
        else:
            results["CDS"].append([])
            results["type"].append(v["type"][i])
            results["codon_start"].append(None)
            results["partialStart"].append(None)
            results["partialStop"].append(None)
            results["protein"].append(None)
            results["cds_transcript"].append(None)
        # since its possible updated the mRNA/CDS fields, double check that gene coordinates are ok
        all_mRNA_coords = [item for sublist in results["mRNA"] for item in sublist]
        try:
            results["location"] = (
                min(all_mRNA_coords, key=lambda item: item[0])[0],
                max(all_mRNA_coords, key=lambda item: item[1])[1],
            )
        except ValueError:
            continue
        # clean up any repeated synonym
        if len(v["gene_synonym"]) > 1:
            uniqueSynonyms = set(v["gene_synonym"])
            results["gene_synonym"] = list(uniqueSynonyms)
    return (k, results)


def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i : i + every])
    return "\n".join(lines)


def longest_orf(annot, fadict, minlen=50, table=1):
    # this is for alignment phasing, where we go through each gene
    # and see if can find full length coding region
    Clean = {}
    for k, v in annot.items():
        if v["strand"] == "+":
            sortedExons = sorted(v["mRNA"][0], key=lambda tup: tup[0])
        else:
            sortedExons = sorted(v["mRNA"][0], key=lambda tup: tup[0], reverse=True)
        mrnaSeq = getSeqRegions(fadict, v["contig"], sortedExons)
        if v["strand"] == "-":
            searchSeq = RevComp(mrnaSeq).upper()
        else:
            searchSeq = mrnaSeq.upper()
        # get all possible ORFs from the mRNA sequence
        valid_starts = "|".join(codon_table[table]["start"])
        allORFS = re.findall(
            rf"((?:{valid_starts})(?:\S{{3}})*?T(?:AG|AA|GA))", searchSeq
        )
        longestORF = None
        # if you found ORFs, then we'll get the longest one and see if meets criterea
        if len(allORFS) > 0:
            longestORF = max(allORFS, key=len)
            if len(longestORF) > minlen * 3:  # than at least min prot length
                # prot = translate(longestORF, "+", 0, table=table)
                # now find out where this ORF is in the coords, sort all back to left --> right
                if v["strand"] == "-":
                    longestORF = RevComp(longestORF)
                    sortedExons = sorted(v["mRNA"][0], key=lambda tup: tup[0])
                # find left most position
                left_pos = mrnaSeq.upper().index(longestORF)
                # map the coords
                CDS = []
                cov = 0
                lenOrf = len(longestORF)
                c = 0
                for i, (s, e) in enumerate(sortedExons):
                    gap_len = e - s + 1
                    last = c
                    c += gap_len
                    if c > left_pos:
                        start_pos = s + left_pos - last
                        break

                for i, (s, f) in enumerate(sortedExons):
                    if s <= start_pos <= f:  # means should start in this exon
                        if (start_pos + lenOrf - 1) <= f:  # then single coord
                            CDS.append((start_pos, (start_pos + lenOrf - 1)))
                            break
                        else:  # spans multiple coords
                            CDS.append((start_pos, f))
                            cov += f - start_pos + 1
                    elif len(CDS) > 0:
                        if cov <= lenOrf:
                            remainder = lenOrf - cov
                            if (f - s) < remainder:
                                CDS.append((s, f))
                                cov += f - s + 1
                            else:
                                CDS.append((s, s + remainder - 1))
                                break
                # see if this is correct
                cdsSeq = getSeqRegions(fadict, v["contig"], CDS)
                try:
                    assert cdsSeq.upper() == longestORF.upper()
                    # okay, then we can add CDS and change type
                    v["type"] = ["mRNA"]
                    if v["strand"] == "+":
                        v["CDS"] = [sorted(CDS, key=lambda tup: tup[0])]
                    else:
                        v["CDS"] = [sorted(CDS, key=lambda tup: tup[0], reverse=True)]
                    v["phase"] = ["?"]
                    Clean[k] = v
                except AssertionError:
                    Clean[k] = v
            else:
                Clean[k] = v
        else:  # did not find orfs so add as is
            Clean[k] = v
    return Clean
