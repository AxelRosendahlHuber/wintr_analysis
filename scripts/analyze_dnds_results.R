library(wintr)

# functions for plotting and measuring in the overlap of the new driver gene method:
cgc_genes = read.delim("~/Downloads/Census_allTue Oct 15 09_47_41 2024.tsv")$Gene.Symbol
cgc_genes =  c(cgc_genes, "CDKN2A.p16INK4a", "CDKN2A.p14arf")


get_overlap_cgc = function(list, qval_column) {

  gene_names = list$gene_name
  overlaps = data.frame(matrix(0, nrow = 75, ncol = 3))
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

# get the overlap for the CGC genes
list = list(cv = dndscv_vanilla$sel_cv |> arrange(qmis_cv, pmis_cv, gene_name) |> select(gene_name, qmis_cv, pmis_cv),
            loc = dndscv_vanilla$sel_loc  |>  arrange(qmis_loc, pmis_loc, gene_name) |> select(gene_name, qmis_loc, pmis_loc),
            loc_intron = dndscv_intron_results$sel_locin |>  arrange(qmis_loc, pmis_loc, gene_name) |> select(gene_name, qmis_loc, pmis_loc),
            #cv_mutrate_intron = dndscv_intron_results$sel_cv |>  arrange(qmis_cv, pmis_cv, gene_name) |> select(gene_name, qmis_cv, pmis_cv),
            #cv_covariate_intron = dndscv_intron_results$sel_cv_intr |>  arrange(qmis_cv, pmis_cv, gene_name) |> select(gene_name, qmis_cv, pmis_cv),
            cv_intron_neut = dndscv_intron_results$sel_cv_intr_neut |>  arrange(qmis_cv, pmis_cv, gene_name) |> select(gene_name, qmis_cv, pmis_cv))

missense = lapply(list, get_overlap_cgc) |>
  rbindlist(idcol = "model") |>
  ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model, linewidth = pval, group = model)) +
  geom_line(alpha = 0.7) + theme_classic()  +
  scale_linewidth_manual(values = c(1,2, 0.5)) +
  labs(title = name, subtitle = "missense")

# get the overlap for the CGC genes
list = list(cv = dndscv_vanilla$sel_cv |> arrange(qallsubs_cv, pallsubs_cv, gene_name) |>  select(gene_name, qallsubs_cv, pallsubs_cv),
            loc = dndscv_vanilla$sel_loc  |>  arrange(qall_loc, pall_loc, gene_name) |>  select(gene_name, qall_loc, pall_loc),
            loc_intron = dndscv_intron_results$sel_locin |>  arrange(qall_loc, pall_loc, gene_name) |>  select(gene_name, qall_loc, pall_loc),
            #cv_mutrate_intron = dndscv_intron_results$sel_cv |>  arrange(qallsubs_cv, pallsubs_cv, gene_name) |>  select(gene_name, qallsubs_cv, pallsubs_cv),
            #cv_covariate_intron = dndscv_intron_results$sel_cv_intr |>  arrange(qallsubs_cv, pallsubs_cv, gene_name) |> select(gene_name, qallsubs_cv, pallsubs_cv),
            cv_intron_neut = dndscv_intron_results$sel_cv_intr_neut |>  arrange(qallsubs_cv, pallsubs_cv, gene_name) |> select(gene_name, qallsubs_cv, pallsubs_cv))
all_subs = lapply(list, get_overlap_cgc) |>
  rbindlist(idcol = "model") |>
  ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model, linewidth = pval, group = model)) +
  geom_line(alpha = 0.7) + theme_classic() +
  scale_linewidth_manual(values = c(1,2, 0.5)) +
  labs(subtitle = 'all subs')

list = list(cv = dndscv_vanilla$sel_cv |> arrange(qtrunc_cv, ptrunc_cv, gene_name) |>  select(gene_name, qtrunc_cv, ptrunc_cv),
            loc = dndscv_vanilla$sel_loc  |>  arrange(qtrunc_loc, ptrunc_loc, gene_name) |>  select(gene_name, qtrunc_loc, ptrunc_loc),
            loc_intron = dndscv_intron_results$sel_locin |>  arrange(qtrunc_loc, ptrunc_loc, gene_name) |>  select(gene_name, qtrunc_loc, ptrunc_loc),
            #cv_mutrate_intron = dndscv_intron_results$sel_cv |>  arrange(qtrunc_cv, ptrunc_cv, gene_name) |>  select(gene_name, qtrunc_cv, ptrunc_cv),
            #cv_covariate_intron = dndscv_intron_results$sel_cv_intr |>  arrange(qtrunc_cv, ptrunc_cv, gene_name) |> select(gene_name, qtrunc_cv, ptrunc_cv),
            cv_intron_neut = dndscv_intron_results$sel_cv_intr_neut |>  arrange(qtrunc_cv, ptrunc_cv, gene_name) |> select(gene_name, qtrunc_cv, ptrunc_cv))
trunc = lapply(list, get_overlap_cgc) |>
  rbindlist(idcol = "model") |>
  ggplot(aes(x = `not in CGC`, y = `in CGC`, color = model, linewidth = pval, group = model)) +
  geom_line(alpha = 0.7) + theme_classic() +
  scale_linewidth_manual(values = c(1,2, 0.5)) +
  labs(subtitle = 'trunc')

line_plot = missense +  trunc  + all_subs + plot_layout(guides = "collect")

# line plot all muts (including indels)
if ("qglobal_cv" %in% colnames(dndscv_vanilla$sel_cv)) {
  list  = list(cv = dndscv_vanilla$sel_cv |> arrange(qglobal_cv, pglobal_cv, gene_name) |>  select(gene_name, qglobal_cv, pglobal_cv),
               #cv_mutrate_intron = dndscv_intron_results$sel_cv |> arrange(qglobal_cv, pglobal_cv, gene_name) |>  select(gene_name, qglobal_cv, pglobal_cv),
               #cv_covariate_intron = dndscv_intron_results$sel_cv_intr |> arrange(qglobal_cv, pglobal_cv, gene_name) |>  select(gene_name, qglobal_cv, pglobal_cv),
               cv_intron_neut = dndscv_intron_results$sel_cv_intr_neut |> arrange(qglobal_cv, pglobal_cv, gene_name) |>  select(gene_name, qglobal_cv, pglobal_cv))
  global = lapply(list, get_overlap_cgc) |>
    rbindlist(idcol = "model") |>
    ggplot(aes(x = `not in CGC` +`not in CGC`, y = `in CGC`, color = model, linewidth = pval, group = model)) +
    geom_line(alpha = 0.7) + theme_classic() +
    scale_linewidth_manual(values = c(1,2, 0.5)) +
    labs(subtitle = 'global', title = name)

  global = lapply(list, get_overlap_cgc) |>
    rbindlist(idcol = "model") |>
    ggplot(aes(x = `not in CGC` + `in CGC`, y = `in CGC`, color = model, linewidth = pval, group = model)) +
    geom_line(alpha = 0.7) + theme_classic() +
    scale_linewidth_manual(values = c(1,2, 0.5)) +
    labs(subtitle = 'global', title = name)


  df = data.frame(total = sapply(list, \(x) filter(x, qglobal_cv< 0.05) |> nrow()),
                  overlap_CGC = sapply(list, \(x) filter(x, qglobal_cv< 0.05) |> pull(gene_name) |> intersect(cgc_genes) |> length())) |>
    mutate(not_inCGC = total - overlap_CGC)
  barplot_global = df |> rownames_to_column("model") |>
    pivot_longer(-c(model, total)) |>
    ggplot(aes(x = model, y = value, fill = name)) + geom_col() +
    labs(x = NULL, y = "# of genes", title = name)


  ggsave(paste0(plotdir, name, "_barplot_global.png"), barplot_global, width = 5, height = 3)
  ggsave(paste0(plotdir, name, "_lineplot_global.png"), global, width = 5, height = 4)
}

#### bar plot #####
# missense
list = list(cv = dndscv_vanilla$sel_cv |> filter(qmis_cv < 0.05),
            loc = dndscv_vanilla$sel_loc  |>  filter(qmis_loc  < 0.05),
            loc_intron = dndscv_intron_results$sel_locin |>  filter(qmis_loc  < 0.05),
            #cv_mutrate_intron = dndscv_intron_results$sel_cv |>  filter(qmis_cv  < 0.05),
            #cv_covariate_intron = dndscv_intron_results$sel_cv_intr |>  filter(qmis_cv  < 0.05),
            cv_intron_neut = dndscv_intron_results$sel_cv_intr_neut |>  filter(qmis_cv  < 0.05))
missense = data.frame(overlap_cgc = sapply(list, \(x) length(intersect(x$gene_name, cgc_genes))),
                      total = sapply(list, nrow)) |> rownames_to_column("model")

# truncating
list = list(cv = dndscv_vanilla$sel_cv |> filter(qtrunc_cv < 0.05),
            loc = dndscv_vanilla$sel_loc  |>  filter(qtrunc_loc  < 0.05),
            loc_intron = dndscv_intron_results$sel_locin |>  filter(qtrunc_loc  < 0.05),
            #cv_mutrate_intron = dndscv_intron_results$sel_cv |>  filter(qtrunc_cv  < 0.05),
            #cv_covariate_intron = dndscv_intron_results$sel_cv_intr |>  filter(qtrunc_cv  < 0.05),
            cv_intron_neut = dndscv_intron_results$sel_cv_intr_neut |>  filter(qtrunc_cv  < 0.05))
trunc = data.frame(overlap_cgc = sapply(list, \(x) length(intersect(x$gene_name, cgc_genes))),
                   total = sapply(list, nrow)) |> rownames_to_column("model")

## all subs
list = list(cv = dndscv_vanilla$sel_cv |> filter(qallsubs_cv < 0.05),
            loc = dndscv_vanilla$sel_loc  |>  filter(qall_loc  < 0.05),
            loc_intron = dndscv_intron_results$sel_locin |>  filter(qall_loc  < 0.05),
            #cv_mutrate_intron = dndscv_intron_results$sel_cv |>  filter(qallsubs_cv  < 0.05),
            #cv_covariate_intron = dndscv_intron_results$sel_cv_intr |>  filter(qallsubs_cv  < 0.05),
            cv_intron_neut = dndscv_intron_results$sel_cv_intr_neut |>  filter(qallsubs_cv  < 0.05))
allsubs = data.frame(overlap_cgc = sapply(list, \(x) length(intersect(x$gene_name, cgc_genes))),
                     total = sapply(list, nrow)) |> rownames_to_column("model")

barplot = rbindlist(list(all_subs = allsubs, missense = missense, truncating = trunc), idcol = "type") |>
  mutate(non_cgc = total - overlap_cgc,
         type = factor(type, levels = c("missense", "truncating", "all_subs")),
         model = factor(model, levels = c("loc", "cv", "cv_mutrate_intron", "loc_intron", "cv_covariate_intron", "cv_intron_neut"))) |>
  pivot_longer(-c(model, type, total)) |>
  ggplot(aes(x = model, y = value, fill = name)) +
  facet_grid(. ~ type) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = NULL, y = "# genes", title = name)

ggsave(paste0(plotdir, name, "_barplot.png"), barplot, width = 5, height = 3)
ggsave(paste0(plotdir, name, "_lineplot.png"), line_plot, width = 10, height = 4)
rm(mutations)
gc()
