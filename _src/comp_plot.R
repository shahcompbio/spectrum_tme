
remove_xaxis <- theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.line.x = element_blank())

remove_guides <- guides(color = F, fill = F, shape = F, alpha = F)

## rank a comp tbl by a certain rank_column and rank_value
rank_by <- function(comp_tbl, rank_column = cell_type, rank_value = "T cell") {
  rank_column <- enquo(rank_column)
  comp_tbl %>% 
    arrange(!!rank_column != rank_value, desc(nrel)) %>% 
    mutate(sample_id_lvl = ordered(sample_id, levels = rev(unique(sample_id)))) %>%
    group_by(!!rank_column) %>% 
    mutate(sample_id_rank = row_number(sample_id_lvl)) %>%
    mutate(sample_id_rank = scales::rescale(sample_id_rank, c(-1, 1))) %>% 
    ungroup %>% 
    mutate(!!rank_column := ordered(!!rank_column, levels = rev(unique(c(rank_value, names(clrs[[as_label(rank_column)]]))))))
}

wilcoxon_wrapper <- function(comp_tbl_rank, rank_column, rank_value, 
                             test_var, test_value) {
  test_var <- enquo(test_var)
  rank_column <- enquo(rank_column)
  comp_tbl_rank <- filter(comp_tbl_rank, !!rank_column == rank_value)
  result <- tryCatch({
    x <- filter(comp_tbl_rank, !!test_var != test_value)$sample_id_rank
    y <- filter(comp_tbl_rank, !!test_var == test_value)$sample_id_rank
    wilcox.test(x = x, y = y, paired = F)$p.value
  }, error = function(e) NA)
  return(result)
}

wilcoxon_test <- function(comp_tbl_rank, rank_column, rank_value, test_var) {
  test_var <- enquo(test_var)
  rank_column <- enquo(rank_column)
  all_values <- unique(pull(comp_tbl_rank, !!test_var))
  lapply(all_values, function(x) wilcoxon_wrapper(comp_tbl_rank, !!rank_column, rank_value, !!test_var, x)) %>% 
    setNames(all_values) %>% 
    unlist %>% 
    enframe(as_label(test_var), "pval")
}

## composition bar plots
plot_comp_bar <- function(comp_tbl_rank, x, y, fill, nmax = 10000, facet = F, super_type = NULL) {
  x <- enquo(x)
  y <- enquo(y)
  fill <- enquo(fill)
  facet <- enquo(facet)
  if (as_label(y) == "nrel") ylab <- paste0("% cells")
  if (as_label(y) == "n") ylab <- paste0("# cells")
  p <- ggplot(comp_tbl_rank, aes(!!x, !!y, fill = !!fill)) +
    geom_bar(stat = "identity", position = position_stack(), width = 1) +
    theme(axis.title.y = element_text(margin = margin(0, -1, 0, 1, unit = "npc")),
          strip.text.y = element_blank(),
          strip.background.y = element_blank(),
          plot.margin = grid::unit(c(0.04, 0.01, 0.02, 0.01), "npc")) +
    remove_xaxis + 
    labs(x = "",
         y = ylab)
    if (!is.na(as_label(facet)) != F) p <- p + facet_grid(cols = vars(!!facet), space = "free_x", scales = "free_x")
  if (as_label(fill) != "cluster_label") p <- p + scale_fill_manual(values = clrs[[as_label(fill)]])
  if (as_label(fill) == "cluster_label") p <- p + scale_fill_manual(values = clrs[[as_label(fill)]][[super_type]])
  if (as_label(y) == "nrel") p <- p +
    scale_y_continuous(expand = c(0, 0),
                       breaks = c(0, 50, 100),
                       labels = c("0", "50", "100")) +
    theme(strip.text.x = element_blank(),
          strip.background.x = element_blank())
  if (as_label(y) == "n") p <- p +
    scale_y_continuous(expand = c(0, 0),
                       breaks = c(0, nmax/2, nmax),
                       limits = c(0, nmax),
                       labels = c("", as.character(c(nmax/2, nmax)))) +
    expand_limits(y = c(0, nmax)) +
    theme(panel.grid.major.y = element_line(linetype = 1, color = "grey90", size = 0.5),
          plot.title = element_text(face = "plain", size = 14,
                                    hjust = 0.5, vjust = 0.5))
  return(p)
  
}


## composition box rank plots
plot_comp_box <- function(comp_tbl_rank, x, y, color, rank_column, rank_value, pcut = 0.01) {
  x <- enquo(x)
  y <- enquo(y)
  color <- enquo(color)
  rank_column <- enquo(rank_column)
  comp_tbl_rank <- filter(comp_tbl_rank, !!rank_column == rank_value)
  
  wilcoxon_tbl <- wilcoxon_test(comp_tbl_rank, !!rank_column, rank_value, !!y) %>% 
    mutate(fdr = p.adjust(pval, method = "BH"),
           pstar = ifelse(pval < pcut, "*", ""))
  
  p <- ggplot(comp_tbl_rank) +
    geom_tile(aes(!!y, !!x, fill = !!color),
              alpha = 0.8) +
    geom_boxplot(aes(!!y, !!x, fill = !!color),
                 width = 0.35, size = 2.5, color = "white", outlier.shape = NA) +
    geom_boxplot(aes(!!y, !!x, color = !!color),
                 width = 0.35, size = 1, fill = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "black") +
    scale_fill_manual(values = clrs[[as_label(color)]]) +
    scale_color_manual(values = clrs[[as_label(color)]]) +
    coord_flip(clip = "off", ylim = c(-1.01, 1.01)) +
    theme(axis.title.y = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          plot.margin = ggplot2::margin(0.04, 0.05, 0.02, 0.01, "npc")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-1, 0, 1)) +
    labs(y = "Scaled rank") +
    geom_text(aes(x = !!y, y = 1.05, label = pstar), 
              data = wilcoxon_tbl, hjust = 0.5, vjust = 1, size = 5, angle = 90)
  return(p)
}

plot_comp_vector <- function(comp_tbl_rank, x, y, shape,
                             vector_column, vector_value,
                             rank_column, rank_value) {
  x <- enquo(x)
  y <- enquo(y)
  shape <- enquo(shape)
  vector_column <- enquo(vector_column)
  rank_column <- enquo(rank_column)
  
  comp_tbl_rank <- filter(comp_tbl_rank, !!rank_column == rank_value) %>% 
    mutate(vector_group = !!vector_column == vector_value) %>% 
    mutate(vector_group = ifelse(vector_group, "xstart", "xend"))
  
  comp_tbl_vector <- group_by(comp_tbl_rank, !!y) %>% 
    mutate(median_rank = median(!!x)) %>% 
    group_by(vector_group, median_rank, !!y) %>% 
    summarise(median_group_rank = median(!!x)) %>% 
    ungroup %>% 
    spread(vector_group, median_group_rank) %>% 
    mutate(vector_color = ifelse(xend > xstart, vector_value, "Other")) %>% 
    arrange(median_rank) %>% 
    mutate(!!y := ordered(!!y, levels = unique(!!y)))
  
  comp_tbl_rank <- comp_tbl_rank %>% 
    mutate(!!y := ordered(!!y, levels = levels(pull(comp_tbl_vector, !!y))))
  
  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_tile(aes(median_rank, !!y), color = "black",
              data = comp_tbl_vector) +
    geom_point(aes(!!x, !!y, shape = !!shape), color = "black",
               data = comp_tbl_rank) +
    geom_segment(aes(x = xstart, xend = xend, y = !!y, yend = !!y, 
                     color = vector_color),
                 arrow = arrow(type = "open", length = unit(0.04, "npc")),
                 data = comp_tbl_vector) +
    scale_shape_manual(values = shps[[as_label(shape)]]) +
    scale_color_manual(values = clrs[[as_label(vector_column)]]) +
    scale_x_continuous(expand = c(0, 0), breaks = c(-1, 0, 1), labels = c(-1, 0, 1)) +
    labs(x = "Scaled rank", y = "Patient", shape = "Sample", color = "Fraction\nin non-adnexa") +
    theme(axis.title.y = element_blank(),
          strip.text = element_blank(),
          plot.margin = ggplot2::margin(0.02, 0.01, 0.02, 0.01, "npc")) +
    # facet_grid(cols = vars(facet),
    #            space = "free_x", scales = "free_x") +
    expand_limits(x = c(-1, 1))
  return(p)
}


default_comp_grid_list <- function(comp_tbl, rank_column, rank_value, 
                                   n_bar = T, nrel_bar = T, 
                                   mutsig_box = T, site_box = T, vec_plot = T,
                                   super_type = NULL) {
  rank_column <- enquo(rank_column)
  comp_tbl_rank <- rank_by(comp_tbl, !!rank_column, rank_value)
  plist <- list()
  if (n_bar) {
    plist$pbar1 <- plot_comp_bar(comp_tbl_rank, sample_id_lvl, n, 
                                 !!rank_column, facet = sort_short_x, 
                                 super_type = super_type) +
      remove_guides
  }
  if (nrel_bar) {
    plist$pbar2 <- plot_comp_bar(comp_tbl_rank, sample_id_lvl, nrel, 
                                 !!rank_column, facet = sort_short_x, 
                                 super_type = super_type) +
      remove_guides
  }
  if (mutsig_box) {
    plist$pbox1 <- plot_comp_box(comp_tbl_rank, sample_id_rank, consensus_signature, 
                                 consensus_signature, !!rank_column, rank_value) + 
      remove_xaxis + remove_guides
  }
  if (site_box) {
    plist$pbox2 <- plot_comp_box(comp_tbl_rank, sample_id_rank, tumor_supersite, 
                                 tumor_supersite, !!rank_column, rank_value) + 
      remove_xaxis + remove_guides
  }
  if (vec_plot) {
    plist$pvec <- plot_comp_vector(comp_tbl_rank, sample_id_rank, patient_id_short,
                                   tumor_megasite, tumor_megasite, "Adnexa", 
                                   !!rank_column, rank_value) +
      remove_guides
  }
  return(plist)
}

# comp_tbl_sample %>% 
#   filter(sort_short_x == "CD45+") %>% 
#   default_comp_grid(cell_type, "T cell")
