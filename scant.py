#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess as sp
import requests
import re

def report_params(params, values):
    print()
    for p, v in zip(params, values):
        print('{}: {}'.format(p, v))
    print()

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

def align(genome, sra_acc, p='4', outdir='', bam='hisat.sorted.bam', novel_splicesite_outfile='splicesite.tab'):
    """Parse arguments for running Hisat2
    """
    novel_splicesite_outfile = os.path.join(outdir, novel_splicesite_outfile)
    bam = os.path.join(outdir, bam)
    print('Running Hisat2')
    params = ['processes', 'genome', 'sra_acc', 'novel_splicesite_outfile', 'bam']
    values = [p, genome, sra_acc, novel_splicesite_outfile, bam]
    report_params(params, values)

    cmd = ['./hisat.sh'] + values
    sp.run(cmd)

def count(ref, bam='hisat.sorted.bam', stringtie_file='stringtie_file.gtf', abundance='abundance.tab', multi_map_frac='.95', outdir='', p='4'):
    """Parse arguments for running Stringtie
    """
    stringtie_file = os.path.join(outdir, stringtie_file)
    abundance = os.path.join(outdir, abundance)
    bam = os.path.join(outdir, bam)
    print('Running Stringtie')
    params = ['processes', 'reference', 'stringtie_file', 'abundance', 'multi_map_frac', 'bam']
    values = [p, ref, stringtie_file, abundance, multi_map_frac, bam]
    report_params(params, values)

    cmd = ['stringtie']
    if p: cmd.extend(['-p', p])
    if abundance: cmd.extend(['-A', abundance])
    cmd += ['-G', ref, '-o', stringtie_file, '-M', multi_map_frac, bam]
    sp.run(cmd)
    return stringtie_file

def merge(ref, gtfs_str, outdir=''):
    """Merges all stringtie output files
    """
    stringtie_merge_outfile = os.path.join(outdir, 'stringtie_merge.gtf')

    print('Running Stringtie --merge')
    #cmd = ['./stringtie_merge.sh'] + [ref, stringtie_merge_outfile, gtfs_str]
    cmd = 'stringtie --merge -G {} -o {} {}'.format(ref, stringtie_merge_outfile, gtfs_str)
    sp.run(cmd.split())

def compare(ref, stringtie_merge_outfile, outdir=''):
    """
    """
    stringtie_merge_outfile = os.path.join(outdir, 'stringtie_merge.gtf')
    gffcomp_annot_gtf = os.path.join(outdir, 'annotated.gtf')

    print('Running Gffcompare')
    cmd = ['./compare.sh'] + [ref, stringtie_merge_outfile]
    sp.run(cmd)

def grab_unique_exons(gffcomp_annot_gtf, outdir=''):
    """
    """
    gffcomp_annot_gtf = os.path.join(outdir, 'annotated.gtf')
    unique_exons_outfile = os.path.join(outdir, 'unique_exons.txt')

    print('Grabbing unique exons')
    with open(gffcomp_annot_gtf) as gtf:
        gtf=gtf.readlines()

    cc = 'class_code'
    cc_u = 'class_code "u"'

    gtf = []
    for i, line in enumerate(gtf):
        if cc in line:
            cc_index = line.find(cc)
            class_code = line[cc_index:cc_index+len(cc_u)]
        else:
            gtf[i]=line+class_code

    gtf = [line for line in gtf if cc_u in line]
    gtf = [line.replace(cc_u, '') for line in gtf if 'exon' in line]

    with open(unique_exons_outfile, 'w') as output:
        for line in gtf:
            output.write('{}'.format(line))

def run_all(genome, ref, srr_set, p='4', outdir='', bam='hisat.sorted.bam',
            novel_splicesite_outfile='splicesite.tab', stringtie_file='stringtie_file.gtf',
            abundance='abundance.tab', multi_map_frac='.95'):
    """Align and count each of a set of SRR numbers
    """
    if isinstance(srr_set, str):
        with open(srr_set) as srr_set:
            srr_set = set([l.strip() for l in srr_set.readlines()])
    stringtie_out_ls = []
    for sra_acc in srr_set:
        sra_acc = sra_acc.strip()
        bam_sra = '{}_{}'.format(sra_acc, bam)
        align(genome,
              sra_acc,
              p,
              outdir,
              bam_sra,
              '{}_{}'.format(sra_acc, novel_splicesite_outfile)
              )
        sto = count(ref,
                    bam_sra,
                    '{}_{}'.format(sra_acc, stringtie_file),
                    '{}_{}'.format(sra_acc, abundance),
                    multi_map_frac,
                    outdir,
                    p
                    )
        stringtie_merge_outfile.append(sto)
    merge(ref, ' '.join(stringtie_out_ls), outdir)
    compare(ref, stringtie_merge_outfile, outdir)
    grab_unique_exons(gffcomp_annot_gtf, outdir)

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
              args.processes,
              args.outdir,
              args.bam,
              args.novel_splicesite_outfile)
        count(args.reference,
              args.bam,
              args.stringtie_file,
              args.abundance,
              args.multi_map_frac,
              args.outdir,
              args.processes)
        return

    elif do_prjna:
        sra_acc = download_project(args.sra_acc)
    else:
        sra_acc = args.file

    run_all(args.genome,
            args.reference,
            sra_acc,
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
    parser.add_argument('reference', help='path to reference .gtf file')
    parser.add_argument('--sra_acc', default='', help='SRR or PRJNA number')
    parser.add_argument('--file', default='', help='path to a file with newline separated SRR numbers')
    parser.add_argument('-p', '--processes', default='4', help='number of cores to use in run')
    parser.add_argument('-o', '--outdir', default='', help='name of directory to save everything to')
    parser.add_argument('-b', '--bam', default='hisat.sorted.bam', help='name of hisat2 output bam file')
    parser.add_argument('-nso', '--novel_splicesite_outfile', default='splicesite.tab', help='Set stringtie novel_splicesite_outfile parameter')
    parser.add_argument('-sf', '--stringtie_file', default='stringtie_file.gtf', help='name of stringtie output gtf file')
    parser.add_argument('-a', '--abundance', default='abundance.tab', help='Set stringtie -A parameter')
    parser.add_argument('-m', '--multi_map_frac', default='.95', help='Set stringtie -M parameter')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    run(args)
