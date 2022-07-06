
import subprocess
from pkg_resources import parse_version
import uuid
import os
import shutil
import pyhmmer
from .utilities import runprocess, execute
from .fasta import fasta2headers

def tblastn_version():
    p1 = subprocess.Popen(['tblastn', '-version'],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)
    vers = p1.communicate()[0].split('+')[0]
    vers = vers.split(' ')[-1]
    return vers


def blast_prefilter(input, query, logger, tmpdir=False,
                    evalue=1e-3, cpus=1, maxhits=3,
                    maxintron=10000):
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
        if parse_version('2.2.32') < parse_version(tblastn_version()) < parse_version('2.10.2'):
            cpus = 1
    # setup tmpdir
    if not tmpdir:
        tmpdir = os.path.join('/tmp', str(uuid.uuid4()))
        os.makedirs(tmpdir)
    # check fasta headers
    headers = fasta2headers(input)
    # start by formatting blast db/dustmasker filtered format
    cmd = ['dustmasker', '-in', os.path.abspath(input), '-infmt', 'fasta', '-parse_seqids',
           '-outfmt', 'maskinfo_asn1_bin', '-out', 'genome_dust.asnb']
    runprocess(cmd, cwd=tmpdir)
    cmd = ['makeblastdb', '-in', os.path.abspath(input), '-dbtype', 'nucl',
           '-parse_seqids', '-mask_data', 'genome_dust.asnb', '-out', 'genome']
    runprocess(cmd, cwd=tmpdir)
    cmd = ['tblastn', '-num_threads', str(cpus), '-db', os.path.join(tmpdir, 'genome'),
           '-query', query, '-max_target_seqs', str(maxhits),
           '-db_soft_mask', '11',
           '-threshold', '999', '-max_intron_length', str(maxintron),
           '-evalue', str(evalue),
           '-outfmt', '6 sseqid slen sstart send qseqid qlen qstart qend pident length evalue score']
    results = {}
    for line in execute(cmd):
        cols = line.rstrip().split('\t')
        qcov = (int(cols[7]) - int(cols[6])) / float(int(cols[5]))
        score = float(cols[11])
        escore = float(cols[10])
        contig = cols[0]
        if not contig in headers:
            if contig.startswith('gb|'):
                contig = contig.replace('gb|', '').rstrip('|')
        if int(cols[3]) >int(cols[2]):
            start = int(cols[2])
            end = int(cols[3])
        else:
            start = int(cols[3])
            end = int(cols[2])
        if not cols[4] in results:
            results[cols[4]] = {contig: []}
        if not contig in results[cols[4]]:
            results[cols[4]][contig] = [{'evalue': escore, 'score': score,
                                         'qcov': qcov, 'qlen': int(cols[5]),
                                         'cov': (int(cols[6]), int(cols[7])),
                                         'coords': (start, end)}]
        else:
            results[cols[4]][contig].append({'evalue': escore, 'score': score,
                                            'qcov': qcov, 'qlen': int(cols[5]),
                                            'cov': (int(cols[6]), int(cols[7])),
                                            'coords': (start, end)})
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
                    final[k] = [{'contig': x, 'coords': w['coords'], 'evalue': w['evalue'], 'score': w['score']}]
                else:
                    final[k].append({'contig': x, 'coords': w['coords'],'evalue': w['evalue'], 'score': w['score']})
    return final, njobs


def merge_overlapping_hits(queryList, fluff=10000):
    if len(queryList) < 2:
        return queryList
    # sort by coords and try to merge if within fluff
    queryList.sort(key=lambda x: x['coords'][0])
    result = [queryList[0]]
    for x in queryList[1:]:
        if x['coords'][0] <= result[-1]['coords'][1] or (result[-1]['coords'][1] - x['coords'][0]) < fluff:
            result[-1]['coords'] = (result[-1]['coords'][0], max(result[-1]['coords'][1], x['coords'][1]))
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
                domains.append({'hmm_from': h.alignment.hmm_from,
                                'hmm_to': h.alignment.hmm_to,
                                'env_from': h.env_from,
                                'env_to': h.env_to,
                                'score': h.score})
            results.append({'name': cog,
                            'bitscore': hit.score, 'evalue': hit.evalue,
                            'domains': domains})
    return results


def hmmer_search(hmmfile, sequences):
    hmm = next(pyhmmer.plan7.HMMFile(hmmfile))
    results = []
    for top_hits in pyhmmer.hmmsearch([hmm], sequences):
        for hit in top_hits:
            cog = hit.best_domain.alignment.hmm_name.decode()
            domains = []
            for h in hit.domains:
                domains.append({'hmm_from': h.alignment.hmm_from,
                                'hmm_to': h.alignment.hmm_to,
                                'env_from': h.env_from,
                                'env_to': h.env_to,
                                'score': h.score})
            results.append({'name': cog, 'hit': hit.name.decode(),
                            'bitscore': hit.score, 'evalue': hit.evalue,
                            'domains': domains})
    return results