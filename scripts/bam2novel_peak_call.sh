# it takes a sorted bam file as input
# and outputs a bed file with novel exons

# read bam input file
if [ "$#" -lt "4" ]; then
  echo "error - inconsistent number of arguments."
  exit 1
fi

bam=$1
gtfbed=$2
outdir=$3
prefix=$4

logfn="$outdir/${prefix}.log"

if [ ! -f $bam ]; then
  echo "error - bam file does not exits."
  exit 1
fi

# check paired or single-ended bam file
n_pair_ended=$(samtools view -c -f 1 $bam)

# peak call with Macs2
if [ $n_pair_ended -gt 0 ]; then
  echo "Paired end .bam"
  macs2 callpeak --keep-dup all --nomodel -f BAMPE -B --SPMR -t $bam  --name $prefix --outdir $outdir 2>&1 | tee $logfn
else
  echo "Single end .bam"
  macs2 callpeak --keep-dup all --nomodel -f BAM -B --SPMR -t $bam  --name $prefix --outdir $outdir  2>&1 | tee $logfn

fi

peakbed="$outdir/${prefix}_summits.bed"
novelbed="$outdir/${prefix}_novel.bed"

# find novel regions
bedtools intersect -v -a $peakbed -b $gtfbed > $novelbed
