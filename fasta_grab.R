### given a gtf file and the genome fasta file, 
### this script spits out transcript sequences in a fasta file.

library(Biostrings)
library(data.table)
library(argparser)

args <- arg_parser("program");
args <- add_argument(args, '-g',
                     help='genome fasta file path',
                     default='/home/ubuntu/data/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz')
args <- add_argument(args, '-gtf',
                     help='gtf file path',
                     default='/home/ubuntu/alorchhota/data_test/unique_test.gtf')
args <- add_argument(args, '-o',
                     help='output fasta file path',
                     default='/home/ubuntu/alorchhota/data_test/unique_test.fa')

argv = parse_args(args)
genome_fn = argv$g
gtf_fn = argv$gtf
out_fa_fn = argv$o

# i/o utils
read_data <- function(fn, sep = '\t',  header = T, row.names=T, stringsAsFactors = F, check.names = F){
  data_df = fread(fn,
                  sep = sep,
                  header = header,
                  stringsAsFactors = stringsAsFactors,
                  check.names = check.names,
                  data.table = FALSE)
  if (row.names == TRUE){
    rownames(data_df) = data_df[,1]
    data_df = data_df[,-1]
  }
  return(data_df)
}

# gtf util
get_gtf_attr <- function(str, attr){
    pattern = paste0(attr, " \"[^\\\"]+\"")
    m = regexec(pattern, str)[[1]]
    if(m>0){
      attr_val = substr(str, m+nchar(attr)+2, m+attr(m, 'match.length')-2)
      attr_val
    } else {
      attr_val = NA
    }
    return(attr_val)
}


# read genome 
genome <- readDNAStringSet(genome_fn)

# find chr to genome mapping
genome_names <- names(genome)
chromosomes <- sapply(genome_names, function(s) strsplit(s, split = ' ')[[1]][1])
chr2genome_name <- names(chromosomes)
names(chr2genome_name) <- as.character(chromosomes)

# read gtf and parse gene_id, transcript_id, range and strand
gtf_data = read_data(gtf_fn, header = F, row.names = F)
gtf_data$gene_id = as.character(sapply(gtf_data$V9, get_gtf_attr, attr = 'gene_id'))
gtf_data$transcript_id = as.character(sapply(gtf_data$V9, get_gtf_attr, attr = 'transcript_id'))
gtf_data$title = apply(gtf_data, 1, function(row){
  paste0('>transcript_id=', row['transcript_id'], ' gene_id=', row['gene_id'], ' range=chr', row['V1'], ':', row['V4'], '-', row['V5'], ' strand=', row['V7'])
})

# for each transcript, spit out fasta sequences
get_transcript_fa_from_gtf <- function(df){
    if(nrow(df) <= 0)
       return(NA)
    chr = as.vector(df[,1])[1]
    tr_title=as.character(as.vector(df[1,'title']))
    indexes1 = as.vector(df[,4])+1   # 0-base to 1-base
    indexes2 = as.vector(df[,5])     # end index is not included
    genome_name = chr2genome_name[[chr]]
    
    sequences = sapply(1:nrow(df), function(idx){
        as.character(subseq(genome[[genome_name]], start=indexes1[idx], end=indexes2[idx]))
    })
    combined_seq = paste0(sequences, collapse='')
    tr_entry = paste(tr_title, combined_seq, sep='\n')
    return(tr_entry)
}

fa_entries = sapply(split(gtf_data, gtf_data$transcript_id),  get_transcript_fa_from_gtf)
all_entries = paste(fa_entries, collapse="\n")
cat(all_entries, file=out_fa_fn)

