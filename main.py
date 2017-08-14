import sys
import argparse
import subprocess as sp

class HiSat2GFFComp:

    def __init__(self, genome=None, sra_acc=None, p=4, outfile='output.bam'):
        self.genome = genome
        self.sra_acc = sra_acc
        self.p = p
        self.outfile = outfile

    def run(self):
        cmd = 'hisat2 -p {} -x {} --sra-acc {} | samtools view -Sb - > {}'.format(self.p, self.genome, self.sra_acc, self.outfile)
        sp.run(cmd.split())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome', help='full path of genome file')
    parser.add_argument('sra_acc', help='SRA Run number')
    parser.add_argument('-p', '--processes', default=4, help='number of cores to use in run')
    parser.add_argument('-o', '--outfile', default=None, help='name of file to save counts to')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    counter = BasicCounter(args.fasta, args.outfile, int(args.kmer),
                           args.nonbinary, args.uncentered,
                           args.unstandardized, label=args.label)
    counter.make_count_file()
