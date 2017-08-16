# it takes a directory with novel bed files,
# suffix to filter out bed files,
# and outputs #sampels contributed to each merged region.

# read bam input file
if [ "$#" -lt "4" ]; then
  echo "error - inconsistent number of arguments."
  exit 1
fi

bed_dir=$1
bed_suffix=$2
outdir=$3
sort_bed=$4

concatenated_bed="$outdir/concatenated.bed"
sorted_bed="$outdir/concatenated_sorted.bed"
concensus_bed="$outdir/concensus.bed"

# merge bed files and count #cells for each 
printf  "" > $concatenated_bed
counter=0
for fn in "$bed_dir/"*"$bed_suffix"; 
do
  counter=$((counter+1))
  awk -v counter=$counter '{print $0 "\t" counter}' $fn >> $concatenated_bed
done

sort -k1,1 -k2,2n -k3,3n -o $sorted_bed  $concatenated_bed

bedtools merge -i $sorted_bed -c 4 -o count_distinct > $concensus_bed    

# clean files
rm $concatenated_bed
rm $sorted_bed

