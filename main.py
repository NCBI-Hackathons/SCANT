import sys
import os
import argparse
import subprocess as sp

def run(genome, sra_acc, ref, p='4', outdir='', bam='hisat.sorted.bam',
        novel_splicesite_outfile='splicesite.tab', stringtie_file='stringtie_file.gtf',
        abundance='abundance.tab', multi_map_frac='.95'):
    bam = os.path.join(outdir, bam)
    novel_splicesite_outfile = os.path.join(outdir, novel_splicesite_outfile)
    stringtie_file = os.path.join(outdir, stringtie_file)
    abundance = os.path.join(outdir, abundance)

    print([p, genome, sra_acc, novel_splicesite_outfile, bam])
    print([p, ref, stringtie_file, abundance, multi_map_frac, bam])
    print('\n\n')

    cmd = ['./hisat.sh'] + [p, genome, sra_acc, novel_splicesite_outfile, bam]
    sp.run(cmd)
    cmd = ['./stringtie.sh'] + [p, ref, stringtie_file, abundance, multi_map_frac, bam]
    sp.run(cmd)

def run_all(genome, sra_file, ref, p='4', outdir='', bam='hisat.sorted.bam',
            novel_splicesite_outfile='splicesite.tab', stringtie_file='stringtie_file.gtf',
            abundance='abundance.tab', multi_map_frac='.95'):
    with open(srr_file) as srr:
        for sra_acc in srr:
            sra_acc = sra_acc.strip()
            run(genome,
                sra_acc,
                ref,
                p,
                outdir,
                '{}_{}'.format(sra_acc, bam),
                '{}_{}'.format(sra_acc, novel_splicesite_outfile),
                '{}_{}'.format(sra_acc, stringtie_file),
                '{}_{}'.format(sra_acc, abundance),
                '{}_{}'.format(sra_acc, multi_map_frac))

1
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome', help='full path of genome file')
    parser.add_argument('sra_acc', help='SRA Run number (or a file of numbers if --all is passed)')
    parser.add_argument('reference', help='Path to reference .gtf file')
    parser.add_argument('-p', '--processes', default='4', help='number of cores to use in run')
    parser.add_argument('-o', '--outdir', default='', help='name of directory to save everything to')
    parser.add_argument('-b', '--bam', default='hisat.sorted.bam', help='name of hisat2 output bam file')
    parser.add_argument('-nso', '--novel_splicesite_outfile', default='splicesite.tab', help='Set stringtie novel_splicesite_outfile parameter')
    parser.add_argument('-sf', '--stringtie_file', default='stringtie_file.gtf', help='name of stringtie output gtf file')
    parser.add_argument('-a', '--abundance', default='abundance.tab', help='Set stringtie -A parameter')
    parser.add_argument('-m', '--multi_map_frac', default='.95', help='Set stringtie -M parameter')
    parser.add_argument('-all', '--all_sra_acc', action='store_true', help='Set to run a file of SRA numbers instead of just one')


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    if args.all_sra_acc:
        run_all(args.genome,
                args.sra_acc,
                args.reference,
                args.processes,
                args.outdir,
                args.bam,
                args.novel_splicesite_outfile,
                args.stringtie_file,
                args.abundance,
                args.multi_map_frac)
    else:
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
