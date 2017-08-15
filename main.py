import sys
import os
import argparse
import subprocess as sp


def run(genome=None, sra_acc=None, ref=None, p='4', outdir='', bam='hisat.bam',
        novel_splicesite_outfile='splicesite.tab', stringtie_file='stringtie_file.gtf',
        abundance='abundance.tab', multi_map_frac='.95'):
    bam = os.path.join(outdir, bam)
    novel_splicesite_outfile = os.path.join(outdir, novel_splicesite_outfile)
    stringtie_file = os.path.join(outdir, stringtie_file)
    abundance = os.path.join(outdir, abundance)

    cmd = ['./hisat.sh'] + [p, genome, sra_acc, novel_splicesite_outfile, bam]
    sp.run(cmd)
    cmd = ['./stringtie.sh'] + [p, ref, stringtie_file, abundance, multi_map_frac, bam]
    sp.run(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome', help='full path of genome file')
    parser.add_argument('sra_acc', help='SRA Run number')
    parser.add_argument('reference', help='Path to reference .gtf file')
    parser.add_argument('-p', '--processes', default='4', help='number of cores to use in run')
    parser.add_argument('-o', '--outdir', default='', help='name of directory to save everything to')
    parser.add_argument('-b', '--bam', default='', help='name of hisat2 output bam file')
    parser.add_argument('-nso', '--novel_splicesite_outfile', default='', help='Set stringtie novel_splicesite_outfile parameter')
    parser.add_argument('-sf', '--stringtie_file', default='', help='name of stringtie output gtf file')
    parser.add_argument('-a', '--abundance', default='', help='Set stringtie -A parameter')
    parser.add_argument('-m', '--multi_map_frac', default='.95', help='Set stringtie -M parameter')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    run(args.genome,
        args.sra_acc,
        args.reference,
        args.processes,
        args.outdir,
        args.bam,
        args.novel_splicesite_outfile,
        args.stringtie_file,
        args.abundance,
        args.multi_map_frac)
