## plotting themes --------------------------------

theme_cowplot2 <- function(...) {
  theme_cowplot(font_size = 16, font_family = "sans", ...) %+replace%
    theme(strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.border = element_blank())
}
theme_set(theme_cowplot2())

remove_xaxis <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.line.x = element_blank())

remove_guides <- guides(color = F, fill = F, shape = F, alpha = F)


## cohort marker genes ----------------------------

markers_v6 <- yaml::read_yaml("/home/uhlitzf/spectrum_tme/_data/small/signatures/hgsc_v6_major.yaml")

helper_markers <- function(x) select(unnest(enframe(x, "subtype", "gene"), cols = gene), gene, subtype)
markers_v6_super <- lapply(yaml::read_yaml("/home/uhlitzf/spectrum_tme/_data/small/signatures/hgsc_v6_super.yaml"), helper_markers)

## load color code --------------------------------

clrs <- yaml::read_yaml("/home/uhlitzf/spectrum_tme/_data/small/signatures/hgsc_v6_colors.yaml") %>% 
  lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

clrs$patient_id_short <- clrs$patient_id
names(clrs$patient_id_short) <- str_remove_all(names(clrs$patient_id), "SPECTRUM-OV-")

## load meta data ----------------------------------

meta_tbl <- read_excel("/home/uhlitzf/spectrum_tme/_data/small/MSK SPECTRUM - Single cell RNA-seq_v6.xlsx", sheet = 2) %>% 
  filter(!(patient_id %in% c("SPECTRUM-OV-100", "SPECTRUM-OV-099")),
         therapy == "pre-Rx") %>% 
  rename(sample = isabl_id) %>% 
  distinct(sample, .keep_all = T) %>% 
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sort_short = str_remove_all(sort_parameters, "singlet, live, "),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite)))

## load mutational signatures ----------------------

signature_tbl <- read_tsv("/home/uhlitzf/spectrum_tme/_data/small/mutational_signatures_summary.tsv") %>%
  select(patient_id, consensus_signature) %>% 
  bind_rows(tibble(patient_id = unique(sort(meta_tbl$patient_id[!(meta_tbl$patient_id %in% .$patient_id)])), consensus_signature = "NA")) %>% 
  mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature))) %>% 
  arrange(consensus_signature) %>% 
  distinct(patient_id, .keep_all = T)

meta_tbl <- left_join(meta_tbl, signature_tbl, by = "patient_id")

## cell type sort fraction -------------------------

cell_type_super_lookup <- c(
  B.cell = "Immune",
  Plasma.cell = "Immune",
  T.cell = "Immune", 
  Monocyte = "Immune", 
  Myeloid.cell = "Immune", 
  Dendritic.cell = "Immune", 
  Endothelial.cell = "Stromal",
  Fibroblast = "Stromal", 
  Ovarian.cancer.cell = "Stromal", 
  Other = "Other"
)
