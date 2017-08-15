# it takes a sorted bam file as input
# and outputs a bed file with novel exons

# read bam input file
if [ "$#" -lt "3" ]; then
  echo "error - inconsistent number of arguments."
  exit 1
fi

bam=$1
gtfbed=$2
outdir=$3

if [ ! -f $bam ]; then
  echo "error - bam file does not exits."
  exit 1
fi

# check paired or single-ended bam file
n_pair_ended=$(samtools view -c -f 1 $bam)

# peak call with Macs2
peakbed=""

# find novel regions
bedtools intersect -v -a $peakbed -b $gtfbed 
