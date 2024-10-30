library(data.table)
library(tidyverse)

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

plot_CGC_plots = function(output_files) {
  
  
  result_list = lapply(output_files, readRDS)
  names(output_files)
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
    rbindlist(idcol = "name")
  CGC_cv = lapply(result_list, \(x) get_overlap_cgc(x$cv |> arrange(qmis_cv, pmis_cv) |> select(gene_name, qmis_cv, pmis_cv))) |> 
    rbindlist(idcol = "name")
  CGC_cv_intron = lapply(result_list, \(x) get_overlap_cgc(x$cv_intron |> arrange(qmis_cv, pmis_cv) |> select(gene_name, qmis_cv, pmis_cv))) |> 
    rbindlist(idcol = "name")
  CGC_models = rbindlist(list(intron_loc = CGC_intron_loc, cv = CGC_cv, cv_intron = CGC_cv_intron), idcol = "model")
  
  missense_plot = CGC_models |> 
    ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model)) + 
    facet_wrap(name ~ .) + 
    geom_line() + 
    theme_bw()
  
  # truncating
  CGC_intron_loc = lapply(result_list, \(x) get_overlap_cgc(x$loc_intron |> arrange(qtrunc_loc, ptrunc_loc) |> select(gene_name, qtrunc_loc, ptrunc_loc))) |> 
    rbindlist(idcol = "name")
  CGC_cv = lapply(result_list, \(x) get_overlap_cgc(x$cv |> arrange(qtrunc_cv, ptrunc_cv) |> select(gene_name, qtrunc_cv, ptrunc_cv))) |> 
    rbindlist(idcol = "name")
  CGC_cv_intron = lapply(result_list, \(x) get_overlap_cgc(x$cv_intron |> arrange(qtrunc_cv, ptrunc_cv) |> select(gene_name, qtrunc_cv, ptrunc_cv))) |> 
    rbindlist(idcol = "name")
  CGC_models = rbindlist(list(intron_loc = CGC_intron_loc, cv = CGC_cv, cv_intron = CGC_cv_intron), idcol = "model")
  
  truncating_plot = CGC_models |> 
    ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model)) + 
    facet_wrap(name ~ .) + 
    geom_line() + 
    theme_bw()
  
  # all
  CGC_intron_loc = lapply(result_list, \(x) get_overlap_cgc(x$loc_intron |> arrange(qall_loc, pall_loc) |> select(gene_name, qall_loc, pall_loc))) |> 
    rbindlist(idcol = "name")
  CGC_cv = lapply(result_list, \(x) get_overlap_cgc(x$cv |> arrange(qallsubs_cv, pallsubs_cv) |> select(gene_name, qallsubs_cv, pallsubs_cv))) |> 
    rbindlist(idcol = "name")
  CGC_cv_intron = lapply(result_list, \(x) get_overlap_cgc(x$cv_intron |> arrange(qallsubs_cv, pallsubs_cv) |> select(gene_name, qallsubs_cv, pallsubs_cv))) |> 
    rbindlist(idcol = "name")
  CGC_models = rbindlist(list(intron_loc = CGC_intron_loc, cv = CGC_cv, cv_intron = CGC_cv_intron), idcol = "model")
  
  all_plot = CGC_models |> 
    ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model)) + 
    facet_wrap(name ~ .) + 
    geom_line() + 
    theme_bw()
  
  return(list(missense_plot = missense_plot, 
              truncating_plot = truncating_plot, 
              all_plot = all_plot))
}


# functions for plotting and measuring in the overlap of the new driver gene method:
cgc_genes = read.delim("~/Downloads/Census_allTue Oct 15 09_47_41 2024.tsv")$Gene.Symbol
cgc_genes =  c(cgc_genes, "CDKN2A.p16INK4a", "CDKN2A.p14arf")

# read in PCAWG RDS files: 
PCAWG_files = list.files("~/Nextcloud/Documents/benchmark/2024-10-23/", full.names = TRUE)
names(PCAWG_files) = gsub(".rds", "", basename(PCAWG_files))
PCAWG_plot_list = plot_CGC_plots(PCAWG_files)

ggsave("plots/PCAWG_missense_plot.png", PCAWG_plot_list$missense_plot, width = 10, height = 8)
ggsave("plots/PCAWG_truncating_plot.png", PCAWG_plot_list$truncating_plot, width = 10, height = 8)
ggsave("plots/PCAWG_all_plot.png", PCAWG_plot_list$all_plot, width = 10, height = 8)


hartwig_files = list.files("~/Nextcloud/Documents/benchmark/hartwig_2024-10-22/")


## development below
# check for significant results of the analysis. I think we know enough about the method at the moment to warrant furhter search. 
# Conclusions: The method is not much better, but could have specific use cases in wgs-sequenced samples
# also, it could be potentially of use to better measure selection. 
filter_significant_genes =  lapply(result_list,\(x) x$cv_intron) |> rbindlist(idcol = "cohort") |> 
  mutate(across(starts_with("q"), \(x) ifelse(x < 0.05, 1,0))) |> 
  group_by(cohort) |> 
  filter(qmis_cv == 1 | qtrunc_cv == 1, qglobal_cv ==1, qallsubs_cv ==1 ) |> 
  mutate(CGC = gene_name %in% cgc_genes) |> 
  pivot_longer(c(starts_with("q"))) |> 
  filter(value == 1)

ggplot(filter_significant_genes, aes(x = name, fill = CGC)) + 
  geom_bar() + 
  facet_wrap(cohort ~ ., scale = "free_y")

