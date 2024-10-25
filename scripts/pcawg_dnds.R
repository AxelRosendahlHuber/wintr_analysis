library(GenomicRanges)
library(dndscv)
library(data.table)
library(patchwork)
library(wintr)
library(tidyverse)

# loop over PCAWG cohort and save the data:
pcawg_files = list.files("/workspace/projects/mut_risk/raw_data/mutation_data/PCAWG_unfiltered/", pattern = "tsv.gz", full.names = TRUE)
names(pcawg_files) = gsub(".tsv.gz", "", basename(pcawg_files))
name = "Cervix-AdenoCA"

day = as.Date(Sys.time(), "%Y%b%e") |> as.character()
outdir = paste0("~/Nextcloud/Documents/benchmark/", day, "/")
if (!dir.exists(outdir)) {dir.create(outdir)}

name = names(pcawg_files)[3]

for (name in names(pcawg_files)) {
  print(name)
  pcawg_file = pcawg_files[name]

  mutations = data.table::fread(pcawg_file)

  mutations = mutations |>
    select(Donor_ID, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2)
  
  mutations = as.data.frame(mutations)

  dndscv_vanilla = dndscv(mutations)
  dndscv_intron_results = dndscv_intron(as.data.frame(mutations), max_muts_per_gene_per_sample = 5)

  list_results = list(loc = dndscv_vanilla$sel_cv, cv = dndscv_vanilla$sel_cv,
                              loc_intron = dndscv_intron_results$sel_locin, cv_intron = dndscv_intron_results$sel_cv_intr_neut,
                              globaldnds = dndscv_vanilla$globaldnds,
                              globaldnds_intron = dndscv_intron_results$globaldnds, 
                              muts_per_sample = muts_per_sample)
  
  saveRDS(list_results, paste0(outdir, "/", name, ".rds"))
}