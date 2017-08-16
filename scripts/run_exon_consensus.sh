scriptfn=/home/ubuntu/alorchhota/SingleCellTranscriptEvidence/scripts/exon_consensus.sh

### toy example
# bed_dir=/home/ubuntu/data/peakCallingTest/toy_bed
# bed_suffix=".bed"
# outdir=/home/ubuntu/data/peakCallingTest/toy_bed/out

### real example
bed_dir=/home/ubuntu/data/novel_peakcall
bed_suffix="_novel.bed"
outdir=/home/ubuntu/data/novel_peakcall

sh $scriptfn $bed_dir $bed_suffix $outdir

