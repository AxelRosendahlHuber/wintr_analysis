library(data.table)
library(tidyverse)

setwd("~/Nextcloud/Documents/wintr_analysis")
get_overlap_cgc = function(list, ngenes = 50) {
  
  gene_names = list$gene_name
  overlaps = data.frame(matrix(0, nrow = ngenes, ncol = 3))
  colnames(overlaps) = c("in CGC", "not in CGC", "pval")
  for (i in 1:nrow(overlaps)) {
    nms = gene_names[1:i]
    overlaps[i,"in CGC"] = length(intersect(cgc_genes, nms))
    overlaps[i,"not in CGC"] = length(nms) -  overlaps[i,1]
    overlaps[i, "pval"] = ifelse(list[i,2] < 0.05, "qval < 0.05", ifelse(list[i,3] < 0.05, "pval < 0.05", "non-signif"))
  }
  overlaps$pval = factor(overlaps$pval, levels = c("pval < 0.05", "qval < 0.05", "non-signif"))
  return(as.data.frame(overlaps))
}

plot_CGC_plots = function(result_list, threshold_samples = 0, threshold_muts = 0) {
  
  # get the number of muts, and the number of samples for each of the cohorts: 
  muts_per_sample = lapply(result_list, \(x) x$muts_per_sample) |> rbindlist(idcol = "cohort")
  
  cohort_stats = muts_per_sample |> 
    group_by(cohort) |> 
    summarize(n_samples = dplyr::n(), 
              n_muts = sum(n))
  
  cohort_stats = cohort_stats |> 
    filter(n_samples >= threshold_samples, n_muts >= threshold_muts)
  result_list = result_list[cohort_stats$cohort]
  if (length(result_list) == 0) {return(NULL)}
  
  
  dnds_intron_results = lapply(result_list, \(x) x$globaldnds_intron) |> 
    rbindlist(idcol = "cohort") |> 
    mutate(model =  "dndscv_intron")
  
  dnds_results = lapply(result_list, \(x) x$globaldnds) |> 
    rbindlist(idcol = "cohort") |> 
    mutate(model = "dndscv")
  
  df_dnds = rbind(dnds_results, dnds_intron_results)
  
  df_dnds |> 
    ggplot(aes(y = mle, x = cohort, , fill = name, shape = model)) + 
    geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5) + 
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = cilow, ymax = cihigh), position = position_dodge(width = 0.5), width = 0.1) + 
    facet_grid(name ~ . ) +
    theme_classic() + 
    cowplot::panel_border() + 
    scale_shape_manual(values = c(21, 22)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    labs(y = "dNdS (omega)",x = NULL, title = "PCAWG") + 
    scale_y_continuous(limits =  c(0, 4))
  
  # missense
  CGC_intron_loc = lapply(result_list, \(x) get_overlap_cgc(x$loc_intron |> arrange(qmis_loc, pmis_loc) |> select(gene_name, qmis_loc, pmis_loc))) |> 
    rbindlist(idcol = "cohort")
  CGC_cv = lapply(result_list, \(x) get_overlap_cgc(x$cv |> arrange(qmis_cv, pmis_cv) |> select(gene_name, qmis_cv, pmis_cv))) |> 
    rbindlist(idcol = "cohort")
  CGC_cv_intron = lapply(result_list, \(x) get_overlap_cgc(x$cv_intron |> arrange(qmis_cv, pmis_cv) |> select(gene_name, qmis_cv, pmis_cv))) |> 
    rbindlist(idcol = "cohort")
  CGC_models = rbindlist(list(intron_loc = CGC_intron_loc, cv = CGC_cv, cv_intron = CGC_cv_intron), idcol = "model")
  
  missense_plot = CGC_models |> 
    ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model)) + 
    facet_wrap(cohort ~ .) + 
    ggpp::geom_text_npc(cohort_stats, mapping = aes(npcx = 0.01, npcy = 0.9, 
                                                    label = paste0("N cases: ", format(n_samples, big.mark = ","), "\nN muts: ", format(n_muts, big.mark = ","))), 
                        size = 2) + 
    geom_line() + 
    theme_bw()
  
  # truncating
  CGC_intron_loc = lapply(result_list, \(x) get_overlap_cgc(x$loc_intron |> arrange(qtrunc_loc, ptrunc_loc) |> select(gene_name, qtrunc_loc, ptrunc_loc))) |> 
    rbindlist(idcol = "cohort")
  CGC_cv = lapply(result_list, \(x) get_overlap_cgc(x$cv |> arrange(qtrunc_cv, ptrunc_cv) |> select(gene_name, qtrunc_cv, ptrunc_cv))) |> 
    rbindlist(idcol = "cohort")
  CGC_cv_intron = lapply(result_list, \(x) get_overlap_cgc(x$cv_intron |> arrange(qtrunc_cv, ptrunc_cv) |> select(gene_name, qtrunc_cv, ptrunc_cv))) |> 
    rbindlist(idcol = "cohort")
  CGC_models = rbindlist(list(intron_loc = CGC_intron_loc, cv = CGC_cv, cv_intron = CGC_cv_intron), idcol = "model")
  
  truncating_plot = CGC_models |> 
    ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model)) + 
    ggpp::geom_text_npc(cohort_stats, mapping = aes(npcx = 0.01, npcy = 0.9, 
                                                    label = paste0("N cases: ", format(n_samples, big.mark = ","), "\nN muts: ", format(n_muts, big.mark = ","))), 
                        size = 2) + 
    facet_wrap(cohort ~ .) + 
    geom_line() + 
    theme_bw()
  
  # all substititions
  CGC_intron_loc = lapply(result_list, \(x) get_overlap_cgc(x$loc_intron |> arrange(qall_loc, pall_loc) |> select(gene_name, qall_loc, pall_loc))) |> 
    rbindlist(idcol = "cohort")
  CGC_cv = lapply(result_list, \(x) get_overlap_cgc(x$cv |> arrange(qallsubs_cv, pallsubs_cv) |> select(gene_name, qallsubs_cv, pallsubs_cv))) |> 
    rbindlist(idcol = "cohort")
  CGC_cv_intron = lapply(result_list, \(x) get_overlap_cgc(x$cv_intron |> arrange(qallsubs_cv, pallsubs_cv) |> select(gene_name, qallsubs_cv, pallsubs_cv))) |> 
    rbindlist(idcol = "cohort")
  CGC_models = rbindlist(list(intron_loc = CGC_intron_loc, cv = CGC_cv, cv_intron = CGC_cv_intron), idcol = "model")
  
  all_plot = CGC_models |> 
    ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model)) + 
    ggpp::geom_text_npc(cohort_stats, mapping = aes(npcx = 0.01, npcy = 0.9, 
                                                    label = paste0("N cases: ", format(n_samples, big.mark = ","), "\nN muts: ", format(n_muts, big.mark = ","))), 
                        size = 2) + 
    facet_wrap(cohort ~ .) + 
    geom_line() + 
    theme_bw()
  
  plot_list = list(missense_plot = missense_plot, 
                   truncating_plot = truncating_plot, 
                   all_plot = all_plot)
  
  
  # indels 
  indel_cohorts_idx = sapply(result_list, \(x) "qind_cv" %in% colnames(x$cv))
  indel_cohorts_idx2 = sapply(result_list, \(x) "pind_intr_cv" %in% colnames(x$cv_intron))
  indel_cohorts_idx == indel_cohorts_idx2
  
  if (sum(indel_cohorts_idx) > 0) {
    result_list_indels = result_list[indel_cohorts_idx]
    CGC_cv = lapply(result_list_indels, \(x) get_overlap_cgc(x$cv |> arrange(qind_cv, pind_cv) |> select(gene_name, qind_cv, pind_cv))) |> 
      rbindlist(idcol = "cohort")
    CGC_cv_intron = lapply(result_list_indels, \(x) get_overlap_cgc(x$cv_intron |> arrange(qind_cv, pind_cv) |> select(gene_name, qind_cv, pind_cv))) |> 
      rbindlist(idcol = "cohort")
    CGC_cv_intron_neut = lapply(result_list_indels, \(x) get_overlap_cgc(x$cv_intron |> arrange(qind_intr_cv, pind_intr_cv) |> select(gene_name, qind_intr_cv, pind_intr_cv))) |> 
      rbindlist(idcol = "cohort")
    CGC_models = rbindlist(list(cv = CGC_cv, cv_intron = CGC_cv_intron, cv_intron_neut = CGC_cv_intron_neut), idcol = "model")
    
    indel_plot = CGC_models |> 
      ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model)) + 
      ggpp::geom_text_npc(cohort_stats, mapping = aes(npcx = 0.01, npcy = 0.9, 
                                                      label = paste0("N cases: ", format(n_samples, big.mark = ","), "\nN muts: ", format(n_muts, big.mark = ","))), 
                          size = 2) + 
      facet_wrap(cohort ~ .) + 
      geom_line() + 
      theme_bw()
    plot_list$indel_plot = indel_plot
  }
    
  return(plot_list)
}

# plot the number of significant genes: 
plot_sig_genes = function(result_list) {
  list_drivers = list()
  y = lapply(result_list, \(x) x$loc_intron |> filter(qall_loc < 0.05) |> pull(gene_name) %in% cgc_genes)
  list_drivers[["loc_intron"]] = t(as.matrix(table(stack(y)))) |> 
    as.data.frame()
  y = lapply(result_list, \(x) x$cv_intron |> filter(qallsubs_cv < 0.05) |> pull(gene_name) %in% cgc_genes)
  list_drivers[["cv_intron"]] = t(as.matrix(table(stack((y))))) |> 
    as.data.frame()
  y = lapply(result_list, \(x) x$cv |> filter(qallsubs_cv < 0.05) |> pull(gene_name) %in% cgc_genes)
  list_drivers[["cv"]] = t(as.matrix(table(stack((y))))) |> 
    as.data.frame()
  all_mut_list = rbindlist(list_drivers, idcol = "model") |> 
    mutate(values = as.character(values), 
      values = ifelse(values == "TRUE", "in_CGC", "not_in_CGC")) |> 
    pivot_wider(values_from = "Freq", names_from = "values" )
  
  all_mut_list |> 
    pivot_longer(-c(model, ind)) |> 
    ggplot(aes(x = model, y = value, fill = name)) + 
    geom_col(position = position_dodge()) + 
    facet_wrap(ind ~ .  , scales = "free") + 
    coord_flip()
}

# functions for plotting and measuring in the overlap of the new driver gene method:
cgc_genes = read.delim("~/Downloads/Census_allTue Oct 15 09_47_41 2024.tsv")$Gene.Symbol
cgc_genes =  c(cgc_genes, "CDKN2A.p16INK4a", "CDKN2A.p14arf")

basename(folder)
folder = "/home/arosendahl/Nextcloud/Documents/wintr_analysis/output/CBIOP_2024-11-20/"
for (folder in dir("~/workspace/projects/mut_risk/wintr_analysis/output/", pattern = "-20", recursive = FALSE, full.names = TRUE)) {
  
  print(folder)
  name = basename(folder)
  # read in PCAWG RDS files: 
  files = list.files(folder, full.names = TRUE)
  names(files) = gsub(".rds", "", basename(files))
  result_list = lapply(files, readRDS)
  plot_list = plot_CGC_plots(result_list, threshold_samples = 50, threshold_muts = 5e4)
  
  for (i in names(plot_list)) {
    print(i)
    ggsave(paste0("plots/",name, "_", i, ".png"), plot_list[[i]],  width = 10, height = 8)
  }
  
  # plot number of significant genes: 
  sig_genes = plot_sig_genes(result_list)
  ggsave(paste0("plots/", name, "_n_sig_genes_all.png"), , width = 20, height = 8)
}

# summarize the cohorts 

# next steps: 
# find out: What are the genes which are more significant across cohorts? Can we make lists of names for each of the different cohorts?