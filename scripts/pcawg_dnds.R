library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(dndscv)
library(data.table)
library(patchwork)
library(wintr)
source("R/dndscv_5.R")

# loop over PCAWG cohort and save the data:
pcawg_files = list.files("/workspace/projects/mut_risk/raw_data/mutation_data/PCAWG_unfiltered/", pattern = "tsv.gz", full.names = TRUE)
names(pcawg_files) = gsub(".tsv.gz", "", basename(pcawg_files))
name = "Cervix-AdenoCA"

day = as.Date(Sys.time(), "%Y%b%e") |> as.character()
outdir = paste0("~/Nextcloud/Documents/benchmark/", day, "/")
if (!dir.exists(outdir)) {dir.create(outdir)}

for (name in names(pcawg_files)) {
  print(name)
  pcawg_file = pcawg_files[name]

  mutations = data.table::fread(pcawg_file)

  mutations = mutations |>
    select(Donor_ID, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2)
  mutations = as.data.frame(mutations)

  dndscv_vanilla = dndscv(mutations)

  # read the dndscv-intron files
  RefCDS_wgs = readRDS("data/RefCDS_wgs.rds")
  gr_genes = readRDS("data/gr_genes.rds")
  gr_transcripts = readRDS("data/gr_transcripts.rds")
  sm = read.delim("data/substmodel_introns.tsv") |>
    as.matrix()

  dndscv_intron_results = dndscv_intron(as.data.frame(mutations), refdb = RefCDS_wgs, gr_transcripts = gr_transcripts, sm = sm, max_muts_per_gene_per_sample = 5)

  list_results = list(loc = dndscv_vanilla$sel_cv, cv = dndscv_vanilla$sel_cv,
                              loc_intron = dndscv_intron_results$sel_locin, cv_intron = dndscv_intron_results$sel_cv_intr_neut,
                              globaldnds = dndscv_vanilla$globaldnds,
                              globaldnds_intron = dndscv_intron_results$globaldnds)

  saveRDS(list_results, paste0(outdir, "/", name, ".rds"))
}
