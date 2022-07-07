import sys
import os
import tempfile
import concurrent.futures
from natsort import natsorted
import pyhmmer
from .__version__ import __version__
from .log import startLogging
from .search import blast_prefilter, hmmer_search_single, hmmer_search, tblastn_version
from .augustus import proteinprofile, augustus_version, augustus_functional
from .fasta import softwrap, getSeqRegions, translate, fasta2dict, fasta2lengths, dict2stats


def load_config(lineage):
    config = {}
    with open(os.path.join(lineage, 'dataset.cfg'), 'r') as infile:
        for line in infile:
            key, value = line.rstrip().split('=')
            config[key] = value
    return config


def load_cutoffs(lineage):
    cutoffs = {}
    with open(os.path.join(lineage, 'scores_cutoff'), 'r') as infile:
        for line in infile:
            busco, score = line.rstrip().split('\t')
            if not busco in cutoffs:
                cutoffs[busco] = {'score': float(score)}
            else:
                cutoffs[busco]['score'] = float(score)
    with open(os.path.join(lineage, 'lengths_cutoff'), 'r') as infile:
        for line in infile:
            busco, _, sigma, length = line.rstrip().split('\t')
            if float(sigma) == 0.0:
                sigma = 1
            if not busco in cutoffs:
                cutoffs[busco] = {'sigma': float(sigma), 'length': int(float(length))}
            else:
                cutoffs[busco]['sigma'] = float(sigma)
                cutoffs[busco]['length'] = int(float(length))
    return cutoffs


def check_lineage(lineage):
    lineage = os.path.abspath(lineage)
    if not os.path.isdir(lineage):
        return False, '{} is not a directory'.format(lineage)
    dirs = ['prfl', 'hmms']
    files = ['ancestral', 'ancestral_variants', 'dataset.cfg',
             'lengths_cutoff', 'scores_cutoff']
    for d in dirs:
        if not os.path.isdir(os.path.join(lineage, d)):
            return False, '{} directory was not found in {}'.format(d, lineage)
    for f in files:
        if not os.path.isfile(os.path.join(lineage, f)):
            return False, '{} file is missing from {}'.format(f, lineage)
    return True, ''


def predict_and_validate(fadict, contig, prfl, cutoffs,
                         species, start,
                         end, configpath,
                         blast_evalue, blast_score):
    # augustus no longer accepts fasta from stdin, so write tempfile
    t_fasta = tempfile.NamedTemporaryFile(suffix=".fasta")
    with open(t_fasta.name, 'w') as f:
        f.write('>{}\n{}\n'.format(contig, fadict[contig][start:end]))
    # run augustus on the regions
    aug_preds = proteinprofile(t_fasta.name, prfl,
                               species=species, start=start,
                               end=end, configpath=configpath)
    t_fasta.close()
    busco_name = os.path.basename(prfl).split('.')[0]
    if len(aug_preds) > 0:
        hmmfile = os.path.join(os.path.dirname(os.path.dirname(prfl)), 'hmms', '{}.hmm'.format(busco_name))
        for k, v in aug_preds.items():
            transcript = getSeqRegions(fadict, v['contig'], v['coords'])
            if v['strand'] == '+':
                protein = translate(transcript, v['strand'], v['phase'][0])
            else:
                protein = translate(transcript, v['strand'], v['phase'][-1])
            if protein.startswith('M') and protein.endswith('*'):
                status = 'complete'
            else:
                status = 'fragmented'
            # now we can check via hmmer
            hmm_result = hmmer_search_single(hmmfile, protein.rstrip('*'))
            if len(hmm_result) > 0:
                if hmm_result[0]['bitscore'] > cutoffs[busco_name]['score']:
                    final = {'contig': v['contig'], 'strand': v['strand'],
                            'location': v['location'], 'coords': v['coords'],
                            'phase': v['phase'], 'transcript': transcript,
                            'translation': protein, 'status': status, 'hmmer': hmm_result[0],
                            'blast_evalue': blast_evalue, 'blast_score': blast_score
                            }
                    return (busco_name, final)
    return False


def runbusco(input, lineage, mode='genome', species='anidulans',
             cpus=1, evalue=1e-50, offset=2000, silent=False, logger=False):
    if not logger:
        logger = startLogging()
    # check that the lineage path/files all present
    check, msg = check_lineage(lineage)
    if not check:
        logger.error('Error: {}'.format(msg))
        sys.exit(1)

    # parse the configs and cutoffs
    Config = load_config(lineage)
    CutOffs = load_cutoffs(lineage)

    if mode == 'genome':
        # check augustus functionality
        aug_version = augustus_version()
        blast_version = tblastn_version()
        logger.info('BUSCOlite v{}; Augustus v{}; tblastn v{}'.format(__version__, aug_version, blast_version))
        if not augustus_functional():
            logger.error('Augustus PPX (--proteinprofile) is non-functional. Usually caused by compilation errors.')
            sys.exit(1)
        if not silent:
            logger.info('{} lineage contains {} BUSCO models'.format(Config['name'], len(CutOffs)))

        # load genome into dictionary
        seq_records = fasta2dict(input)
        fa_stats = dict2stats(seq_records)
        if not silent:
            logger.info('Input genome is {} MB and {} contigs'.format(
                round(fa_stats['size'] / 1e6, 2), fa_stats['n']
            ))

        # run blast filter from ancesteral proteins
        query = os.path.join(lineage, 'ancestral')
        if not silent:
            logger.info('Prefiltering predictions using tblastn of ancestral sequences')
        coords, njobs = blast_prefilter(input, query, logger, evalue=evalue, cpus=cpus)
        if not silent:
            logger.info('Now launching {} augustus/phymmer jobs for {} BUSCO models'.format(njobs, len(coords)))

        # run busco analysis using threadpool, limit io as much as possible
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus+4) as executor:
            for k, v in coords.items():
                busco_prlf = os.path.join(lineage, 'prfl', '{}.prfl'.format(k))
                for i in v:
                    contig = i['contig']
                    start = i['coords'][0] - offset
                    if start < 0:
                        start = 0
                    end = i['coords'][1] + offset
                    if end > len(seq_records[contig]):
                        end = len(seq_records[contig])
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
                            False,
                            i['evalue'],
                            i['score']
                        )
                    )
        # process mt results
        b_results = {}
        for r in results:
            if isinstance(r.result(), tuple):
                b, res = r.result()
                if not b in b_results:
                    b_results[b] = [res]
                else:
                    b_results[b].append(res)
       # first run is finished, now see which buscos are missing
        missing = []
        for b in CutOffs.keys():
            if b not in b_results:
                missing.append(b)
        if not silent:
            logger.info('Found {} BUSCOs in first pass, trying harder to find remaining {}'.format(
                len(b_results), len(missing)))
        filt_variants = tempfile.NamedTemporaryFile(suffix=".fasta")
        avariants = fasta2dict(os.path.join(lineage, 'ancestral_variants'), full_header=True)
        seen = set()
        with open(filt_variants.name, 'w') as outfile:
            for title, seq in avariants.items():
                try:
                    z, num = title.split(' ')
                except ValueError:
                    if '_' in title:
                        z, num = title.rsplit('_', 1)
                    else:
                        z = title
                if z in missing:
                    seen.add(z)
                    outfile.write('>{} {}\n{}\n'.format(z, num, softwrap(seq)))
        logger.info('Trying to use ancestral variants to recover {} BUSCOs'.format(len(seen)))
        coords2, njobs2 = blast_prefilter(input, filt_variants.name, logger, cpus=cpus, evalue=1e-5)
        filt_variants.close()
        if not silent:
            logger.info('Now launching {} augustus/pyhmmer jobs for {} BUSCO models'.format(njobs2, len(coords2)))
        if len(coords2) > len(seen):
            logger.error('Something is wrong with parsing the acenstral variants.....')
            print(coords2.keys())
            sys.exit(1)
        # try thread pool
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus+4) as executor:
            for k, v in coords2.items():
                busco_prlf = os.path.join(lineage, 'prfl', '{}.prfl'.format(k))
                for i in v:
                    contig = i['contig']
                    start = i['coords'][0] - offset
                    if start < 0:
                        start = 0
                    end = i['coords'][1] + offset
                    if end > len(seq_records[contig]):
                        end = len(seq_records[contig])
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
                            False,
                            i['evalue'],
                            i['score']
                        )
                    )
        # process mt results
        for r in results:
            if isinstance(r.result(), tuple):
                b, res = r.result()
                if not b in b_results:
                    b_results[b] = [res]
                else:
                    b_results[b].append(res)

        # finally loop through results, classify and reorganize
        b_final = {}
        missing = []
        for b in CutOffs.keys():
            if b not in b_results:
                missing.append(b)
        stats = {'total': 0, 'single-copy': 0, 'fragmented': 0, 'duplicated': 0}
        for k, v in natsorted(b_results.items()):
            stats['total'] += 1
            if len(v) > 1:  # duplicates
                for i, x in enumerate(sorted(v, key=lambda y: y['hmmer']['bitscore'], reverse=True)):
                    x['status'] = 'duplicated'
                    if i > 0:
                        name = '{}_{}'.format(k, i)
                    else:
                        name = k
                        stats['duplicated'] += 1
                    b_final[name] = x
            else:
                if v[0]['status'] == 'fragmented':
                    stats['fragmented'] += 1
                else:
                    stats['single-copy'] += 1
                b_final[k] = v[0]
        return b_final, missing, stats, Config

    elif mode == 'proteins':
        if not silent:
            logger.info('BUSCOlite v{}; pyhmmer v{}'.format(__version__, pyhmmer.__version__))
            logger.info('{} lineage contains {} BUSCO models'.format(Config['name'], len(CutOffs)))
        # load proteome into easel digitized sequence to pass to pyhmmer
        falengths = fasta2lengths(input)
        alphabet = pyhmmer.easel.Alphabet.amino()
        sequences = []
        with pyhmmer.easel.SequenceFile(input, digital=True, alphabet=alphabet) as seq_file:
            sequences = list(seq_file)
        if not silent:
            logger.info('Loaded {} protein sequences from {}'.format(len(sequences), input))
        # now we can loop over the hmms in the lineage and run hmmer on each
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus+4) as executor:
            for f in os.listdir(os.path.join(lineage, 'hmms')):
                if f.endswith('.hmm'):
                    hmmfile = os.path.join(lineage, 'hmms', f)
                    results.append(
                        executor.submit(
                            hmmer_search,
                            hmmfile,
                            sequences
                        )
                    )
        # process mt results
        b_results = {}
        for r in results:
            if isinstance(r.result(), list):
                for x in r.result():
                    if x['bitscore'] > CutOffs[x['name']]['score']:
                        if not x['name'] in b_results:
                            b_results[x['name']] = [x]
                        else:
                            b_results[x['name']].append(x)
        b_final = {}
        missing = []
        for b in CutOffs.keys():
            if b not in b_results:
                missing.append(b)
        stats = {'total': 0, 'single-copy': 0, 'fragmented': 0, 'duplicated': 0}
        for k, v in natsorted(b_results.items()):
            stats['total'] += 1
            if len(v) > 1:  # duplicates
                for i, x in enumerate(sorted(v, key=lambda y: y['bitscore'], reverse=True)):
                    x['status'] = 'duplicated'
                    x['length'] = falengths[x['hit']]
                    if i > 0:
                        name = '{}_{}'.format(k, i)
                    else:
                        name = k
                        stats['duplicated'] += 1
                    b_final[name] = x
            else:
                x = v[0]
                x['length'] = falengths[x['hit']]
                if 'status' in x and x['status'] == 'fragmented':
                    stats['fragmented'] += 1
                    x['status'] = 'fragmented'
                else:
                    x['status'] = 'complete'
                    stats['single-copy'] += 1
                b_final[k] = x

        return b_final, missing, stats, Config
