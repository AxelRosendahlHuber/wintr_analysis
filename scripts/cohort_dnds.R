# run dnds for cohort 
library(GenomicRanges)
library(dndscv)
library(data.table)
library(patchwork)
library(wintr)
library(tidyverse)

# loop over cohort cohort and save the data:
args = commandArgs(trailingOnly = TRUE)
print(args)
cohort_folder = as.character(args[1])
cohort_files = list.files(cohort_folder, full.names = TRUE)
cohort_type = as.character(args[3])

names(cohort_files) = gsub("\\..*", "", basename(cohort_files))
day = as.Date(Sys.time(), "%Y%b%e") |> as.character()
outdir = paste0("output/", cohort_type, "_", day,  "/")
if (!dir.exists(outdir)) {dir.create(outdir)}

# analyze the file from the cohort  
cohort_file = cohort_files[as.numeric(args[2])]
name = names(cohort_files)[as.numeric(args[2])]
print(cohort_file)
mutations = data.table::fread(cohort_file)

if (cohort_type == "PCAWG") {
  mutations = mutations |>
    select(Donor_ID, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2)
}

mutations = as.data.frame(mutations)
muts_per_sample = mutations |>
  group_by(sampleID) |> 
  count("n_muts")

dndscv_vanilla = dndscv(mutations)
dndscv_intron_results = dndscv_intron(as.data.frame(mutations), max_muts_per_gene_per_sample = 5)

list_results = list(loc = dndscv_vanilla$sel_cv, cv = dndscv_vanilla$sel_cv,
                    loc_intron = dndscv_intron_results$sel_locin, cv_intron = dndscv_intron_results$sel_cv_intr_neut,
                    globaldnds = dndscv_vanilla$globaldnds,
                    globaldnds_intron = dndscv_intron_results$globaldnds, 
                    muts_per_sample = muts_per_sample)

saveRDS(list_results, paste0(outdir, "/", name, ".rds"))



