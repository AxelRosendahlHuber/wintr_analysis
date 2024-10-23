# goal: The aim is to include non-coding protein data in the genomic coding file
library(dndscv) # dNdScv is required as the current package extends on dNdScv 
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(tidyverse)

data("refcds_hg19", package = "dndscv") # load the RefCDS object (and load the gr_genes object)
data("submod_192r_3w", package = "dndscv")

# load the transcripts from the biomart downloaded table. Filter for transcripts which do not generate a protein 
transcripts = fread("../data/biomart_transcripts.txt") |> 
  filter(`Protein stable ID` != "")

# get the protein ids from dNdScv refCDS object
protein_ids = sapply(RefCDS, \(x) x$protein_id)

# reformat the transcript sequences
transcript_sequences = transcripts |> 
  dplyr::filter(`Protein stable ID` %in% protein_ids) |> 
  mutate(length = `Transcript start (bp)` - `Transcript end (bp)`) |> 
  rename(protein_id = `Protein stable ID`, 
         seqnames = `Chromosome/scaffold name`,
         start = `Transcript start (bp)`, 
         end = `Transcript end (bp)`, 
         strand = Strand, 
         gene_name = `Gene name`) |> 
  select(seqnames, start, end, strand, protein_id, gene_name) |> 
  arrange(match(protein_id, protein_ids))

# edge case CDKN2A codes for tho different proteins using different transcripts
# replace CDKN2A in the gene list with the CDKN2A.p16ARF and CDKN2A.p16INK4A transcripts (keep same transcript table)
df_gene_protein = data.frame(gene_name = sapply(RefCDS, \(x) x$gene_name), 
                             protein_id = sapply(RefCDS, \(x) x$protein_id))
cdkn2a_protein_sequences = df_gene_protein |> 
  filter(grepl("CDKN2A\\.", gene_name))
transcript_sequences =  transcript_sequences |> 
  mutate(gene_name = case_when(protein_id == "ENSP00000462950" ~ "CDKN2A.p14arf", 
                               protein_id == "ENSP00000307101" ~ "CDKN2A.p16INK4a", 
                               .default =  gene_name))

all(transcript_sequences$protein_id == protein_ids) # quick check if the order of the protein ids is indeed the same
gr_transcripts = plyranges::as_granges(transcript_sequences |> select(-protein_id))


# combine the gr_transcripts and the gr_genes sets: 
transcript_genes = unique(gr_transcripts$gene_name) # can possibly be removed - the number of transcripts and genes is similar
gr_transcripts_ns = gr_transcripts
strand(gr_transcripts_ns) = "*"
gr_genes$gene_name = gr_genes$names
gr_total = c(gr_transcripts_ns, gr_genes)
gr_transcripts_all = GRangesList(split(gr_total, gr_total$gene_name)) |> 
  GenomicRanges::reduce() |> 
  unlist()

# assign the strand to the genes: 
unique_gr_names = as.data.frame(gr_transcripts) |> 
  select(gene_name, strand) |> unique() |> 
  column_to_rownames("gene_name")
strand(gr_transcripts_all) = unique_gr_names[names(gr_transcripts_all), "strand"]

seqlevelsStyle(gr_transcripts_all) = "UCSC"
sequences = getSeq(Hsapiens, gr_transcripts_all + 1) # calculate the trinucleotide frequencies for the individual reads

counts = trinucleotideFrequency(sequences)

# summarize counts as the names of the RefCDS 
gene_counts = counts |> as.data.frame() |> 
  dplyr::mutate(gene_name = names(gr_transcripts_all)) |> 
  group_by(gene_name) |> 
  summarize(across(everything(), sum)) |> 
  column_to_rownames("gene_name") |> 
  as.matrix()

trinucs = substr(rownames(substmodel), 1,3)
counts_match = gene_counts[,rep(1:64, each = 3)]

RefCDS_5 = RefCDS
for (i in 1:length(RefCDS)) { 
  print(i)
  
  gene_name = RefCDS[[i]]$gene_name
  print(gene_name)
  counts = counts_match[gene_name, ] - rowSums(RefCDS[[i]]$L) 
  
  if (min(counts) < 0) {stop("not possible scenario - exons larger than transcript")}
  
  RefCDS_5[[i]]$L = cbind(RefCDS_5[[i]]$L, counts)
  RefCDS_5[[i]]$intervals_transcript = as.numeric(transcript_sequences[i, c("start", "end")])
}

saveRDS(RefCDS_5, "../data/RefCDS_wgs.rds")
seqlevelsStyle(gr_transcripts_all) = "NCBI"
saveRDS(gr_transcripts_all, "../data/gr_transcripts.rds")