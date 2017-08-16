script_fn='/home/ubuntu/alorchhota/SingleCellTranscriptEvidence/scripts/fasta_grab.R'
genome_fn='/home/ubuntu/data/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'
gtf_fn='/home/ubuntu/alorchhota/data_test/unique_test.gtf'
out_fa_fn='/home/ubuntu/alorchhota/data_test/unique_test.fa'

Rscript $script_fn -g $genome_fn -gtf $gtf_fn -o $out_fa_fn
