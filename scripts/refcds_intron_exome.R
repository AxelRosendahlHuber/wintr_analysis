# goal: The aim is to include non-coding protein data in the genomic coding file
library(dndscv)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(tidyverse)

data("refcds_hg19", package = "dndscv") # load the RefCDS object 
data("submod_192r_3w", package = "dndscv")

# load the exonic capture panels 
exonic_capture = import.bed("~/Nextcloud/Documents/wintr/data/exome_probesets/sorted-SeqCap_EZ_ExomeV3_Plus_UTR_hg19_primary_annotated.bed.gz")
exonic_capture$gene_name = gsub("gene_symbol\\|", "", exonic_capture$name)
exonic_capture = exonic_capture[!is.na(exonic_capture$name),]
seqlevelsStyle(exonic_capture) = "NCBI"

# get the intersect of the exonic areas and the genes 
gene_names = sapply(RefCDS, \(x) x$gene_name)
setdiff(unique(exonic_capture$name), gene_names)
setdiff(gene_names, unique(exonic_capture$name))
common_genes = intersect(gene_names, unique(exonic_capture$gene_name))

# filter for the common genes
RefCDS_select = RefCDS[gene_names %in% common_genes]
exonic_capture_gr = exonic_capture[exonic_capture$gene_name %in% common_genes]

# get the right transcripts
transcripts = fread("data/biomart_transcripts.txt") |> 
  filter(`Protein stable ID` != "")
protein_ids = sapply(RefCDS_select, \(x) x$protein_id)
transcript_sequences = transcripts |> 
  dplyr::filter(`Protein stable ID` %in% protein_ids) |> 
  rename(strand = Strand, gene_name = `Gene name`, protein_id = `Protein stable ID`) |> 
  select(gene_name, strand, protein_id)

# edge case CDKN2A codes for tho different proteins using different transcripts
# replace CDKN2A in the gene list with the CDKN2A.p16ARF and CDKN2A.p16INK4A transcripts (keep same transcript table)
df_gene_protein = data.frame(gene_name = sapply(RefCDS_select, \(x) x$gene_name), 
                             protein_id = sapply(RefCDS_select, \(x) x$protein_id))
cdkn2a_protein_sequences = df_gene_protein |> 
  filter(grepl("CDKN2A\\.", gene_name))
transcript_sequences =  transcript_sequences |> 
  mutate(gene_name = case_when(protein_id == "ENSP00000462950" ~ "CDKN2A.p14arf", 
                               protein_id == "ENSP00000307101" ~ "CDKN2A.p16INK4a", 
                               .default =  gene_name)) |> 
  as.data.frame()
rownames(transcript_sequences) = transcript_sequences$gene_name

# combine the exonic_capture and the gr_genes sets: 
gr_genes$gene_name = gr_genes$names
gr_total = c(gr_genes, exonic_capture_gr)
gr_capture_all = GRangesList(split(gr_total, gr_total$gene_name)) |> 
  GenomicRanges::reduce() |> unlist()

# assign the strand to the genes: 
gr_capture_all$gene_name = names(gr_capture_all)
transcript_names = transcript_sequences$gene_name
gr_capture_all = gr_capture_all[gr_capture_all$gene_name %in% transcript_names]
strand(gr_capture_all) = transcript_sequences[names(gr_capture_all), "strand"]
seqlevelsStyle(gr_capture_all) = "UCSC"

sequences = getSeq(Hsapiens, gr_capture_all + 1) # calculate the trinucleotide frequencies for the individual reads
counts = trinucleotideFrequency(sequences)

# summarize counts as the names of the RefCDS 
gene_counts = counts |> as.data.table() |> 
  dplyr::mutate(gene_name = names(gr_capture_all)) 
gene_counts = gene_counts[, lapply(.SD, sum), by = "gene_name"] |> 
  column_to_rownames("gene_name") |> 
  as.matrix()

trinucs = substr(rownames(substmodel), 1,3)
counts_match = gene_counts[,rep(1:64, each = 3)]

RefCDS_select = RefCDS[match(rownames(counts_match), gene_names)]
for (i in 1:length(RefCDS_select)) { 
  print(i)
  
  gene_name = RefCDS_select[[i]]$gene_name
  print(gene_name)
  counts = counts_match[gene_name, ] - rowSums(RefCDS_select[[i]]$L) 
  
  warnings = character()
  if (min(counts) < 0) {
    counts[counts < 0] = 0
    warnings = c(warnings, gene_name)
    }
  
  L = as.matrix(cbind(RefCDS_select[[i]]$L, counts))
  dimnames(L) = NULL
  RefCDS_select[[i]]$L = L
  
  if (length(warnings > 0)) {
    warning(paste0("not possible scenario - exons larger than transcript for genes: ", warnings, collapse = " "))
    # for some reason, the gene SYNM has one value more. As this is the only single case across all regions in all genes, we will set the value to 0.
  }
}

saveRDS(RefCDS_select, "data/RefCDS_wes.rds")
seqlevelsStyle(gr_capture_all) = "NCBI"
saveRDS(gr_capture_all, "data/gr_capture.rds")