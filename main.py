#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess as sp
import requests
import re

def count(novel_splicesite_outfile='splicesite.tab', stringtie_file='stringtie_file.gtf',
          abundance='abundance.tab', multi_map_frac='.95', outdir=''):
    """Parse arguments for running Stringtie
    """
    novel_splicesite_outfile = os.path.join(outdir, novel_splicesite_outfile)
    stringtie_file = os.path.join(outdir, stringtie_file)
    abundance = os.path.join(outdir, abundance)
    print('Running Stringtie\n', [p, ref, stringtie_file, abundance, multi_map_frac, bam])
    print('\n\n')
    cmd = ['stringtie']
    if p: cmd.extend(['-p', p])
    if ref: cmd.extend(['-G', ref])
    if abundance: cmd.extend(['-A', abundance])
    cmd += ['-o', stringtie_file, '-M', multi_map_frac, bam]
    sp.run(cmd)

def align(genome, sra_acc, ref, p='4', outdir='', bam='hisat.sorted.bam'):
    """Parse arguments for running Hisat2
    """
    bam = os.path.join(outdir, bam)
    print('Running Hisat2\n', [p, genome, sra_acc, novel_splicesite_outfile, bam])
    cmd = ['./hisat.sh'] + [p, genome, sra_acc, novel_splicesite_outfile, bam]
    sp.run(cmd)

def download_project(query):
    """Download all [SED]RR numbers for a given SRA project

    Returns
    -------
    srr_set : (str)
        Full set of data file identification numbers for a given project
    """
    query_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmax=5000'
    r = requests.get(query_url, params = {'term':query})
    id_list = list(map(int, re.findall("<Id>([0-9]+)</Id>", r.text)))
    srr_set = {}
    for proj in id_list:
        proj_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&format=runinfo'
        s = requests.get(proj_url, params = {'id':proj})
        srr_set.update(set(re.findall("<RUN.*accession=\"([SED]RR[0-9]+)\"", proj.text)))

    with open('run_accessions.txt', 'w') as f:
        for item in srr_set:
            f.write(item+'\n')

    return srr_set

def run_all(genome, srr_set, ref, p='4', outdir='', bam='hisat.sorted.bam',
            novel_splicesite_outfile='splicesite.tab', stringtie_file='stringtie_file.gtf',
            abundance='abundance.tab', multi_map_frac='.95'):
    """Align and count each of a set of SRR numbers
    """
    if isinstance(srr_set, str):
        with open(srr_set) as srr_set:
            srr_set = set([l.strip() for l in srr_set.readlines()])
    for sra_acc in srr_set:
        sra_acc = sra_acc.strip()
        align(genome,
              sra_acc,
              ref,
              p,
              outdir,
              '{}_{}'.format(sra_acc, bam)
              )
        count('{}_{}'.format(sra_acc, novel_splicesite_outfile),
              '{}_{}'.format(sra_acc, stringtie_file),
              '{}_{}'.format(sra_acc, abundance),
              multi_map_frac,
              outdir
              )

def run(args):
    error = not args.sra_acc and not args.file
    prjna = re.fullmatch("PRJNA[0-9]+", args.sra_acc)
    single_srr = args.sra_acc and not prjna
    do_prjna = args.sra_acc and prjna

    if error:
        print('You must use either --sra_acc or --file')
        sys.exit()
    elif single_srr:
        align(args.genome,
              args.sra_acc,
              args.reference,
              args.processes,
              args.outdir,
              args.bam
              )
        count(args.novel_splicesite_outfile,
              args.stringtie_file,
              args.abundance,
              args.multi_map_frac
              )
        return

    elif args.sra_acc: #project id (PRJNA)
        sra_acc = download_project(args.sra_acc)
    else: #A file of SRR ids
        sra_acc = args.sra_acc

    run_all(args.genome,
            sra_acc,
            args.reference,
            args.processes,
            args.outdir,
            args.bam,
            args.novel_splicesite_outfile,
            args.stringtie_file,
            args.abundance,
            args.multi_map_frac)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome', help='path of genome file')
    parser.add_argument('--sra_acc', default='', help='SRR or PRJNA number')
    parser.add_argument('--file', default='', help='path to a file with newline separated SRR numbers')
    parser.add_argument('-r', '--reference', default='', help='path to reference .gtf file')
    parser.add_argument('-p', '--processes', default='4', help='number of cores to use in run')
    parser.add_argument('-o', '--outdir', default='', help='name of directory to save everything to')
    parser.add_argument('-b', '--bam', default='hisat.sorted.bam', help='name of hisat2 output bam file')
    parser.add_argument('-nso', '--novel_splicesite_outfile', default='splicesite.tab', help='Set stringtie novel_splicesite_outfile parameter')
    parser.add_argument('-sf', '--stringtie_file', default='stringtie_file.gtf', help='name of stringtie output gtf file')
    parser.add_argument('-a', '--abundance', default='abundance.tab', help='Set stringtie -A parameter')
    parser.add_argument('-m', '--multi_map_frac', default='.95', help='Set stringtie -M parameter')
    #parser.add_argument('-all', '--all_sra_acc', action='store_true', help='Set to run a file of SRA numbers instead of just one')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    run(args)
