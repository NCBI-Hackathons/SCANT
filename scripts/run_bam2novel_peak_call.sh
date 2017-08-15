srrfn=$1
bam2novel_peak_call_script="/home/ubuntu/greenkidneybean/SingleCellTranscriptEvidence/scripts/bam2novel_peak_call.sh"
peak_call_dir=/home/ubuntu/data/novel_peakcall
bamdir=/home/ubuntu/data/alignments
gtf_bed_fn=/home/ubuntu/data/annot/ref_GRCm38.p4_top_level.ens.bed
cat $srrfn | while read acc; 
do 
  echo $acc; 
  bamfn="$bamdir/${acc}_hisat.sorted.bam"
  sh $bam2novel_peak_call_script $bamfn $gtf_bed_fn $peak_call_dir $acc
done

