# it takes a directory with novel bed files,
# suffix to filter out bed files,
# and outputs #sampels contributed to each merged region.

# read bam input file
if [ "$#" -lt "3" ]; then
  echo "error - inconsistent number of arguments."
  exit 1
fi

novel_bed_dir=$1
bed_suffix=$2
outdir=$3


