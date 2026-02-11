## Selection plots overview ##


#### Heatmaps per Ampliconic cluster ####

# ---- Packages ---- #
library(dplyr)
library(ggplot2)
library(gridExtra)


# ---- load in the data ---- # 
# X chromosome pairwise dN/dS
selection_df_x <- read.csv("~/Downloads/Bootstrap_all_families_dN_dS_between_within_species_x_updated2.csv")
# Y chromosome pairwise dN/dS
#selection_df_x <- read.csv("~/Downloads/Bootstrap_all_families_dN_dS_between_within_species_y_updated2.csv")

# ---- Plot function ---- #
# define species order
species_order <- c(
  "PanPan", "PanTro", "HomSap", "GorGor", "PonAbe", "PonPyg", "SymSyn", "MacFas")

make_heatmap <- function(df, fam, clu) {
  d <- df %>%
    filter(Family == fam, Cluster == clu) %>%
    mutate(
      Selection = tolower(selection),
      a = pmin(Species1, Species2),
      b = pmax(Species1, Species2)
    ) %>%
    distinct(a, b, .keep_all = TRUE)
  
  #species <- sort(unique(c(d$Species1, d$Species2)))
  species <- species_order[species_order %in% unique(c(d$Species1, d$Species2))] # order species according to the phylogeny
  
  grid <- expand.grid(Species1 = species, Species2 = species, stringsAsFactors = FALSE) %>%
    mutate(
      i = match(Species1, species),
      j = match(Species2, species),
      keep_lower = i >= j,
      a = pmin(Species1, Species2),
      b = pmax(Species1, Species2)
    ) %>%
    filter(keep_lower) %>%
    left_join(d %>% select(a, b, dNdS, Selection), by = c("a","b")) %>%
    mutate(
      Selection_fill = ifelse(Selection %in% "nonsignificant", NA, Selection),
      dNdS_label     = ifelse(is.na(dNdS), "", sprintf("%.2f", dNdS)),
      Species1 = factor(Species1, levels = species),
      Species2 = factor(Species2, levels = rev(species))
    )
  
  ggplot(grid, aes(Species1, Species2)) +
    geom_tile(aes(fill = Selection_fill), color = "grey85", linewidth = 0.3) +
    geom_text(aes(label = dNdS_label), size = 2) +
    scale_fill_manual(values = c(positive="#E64B35", purifying="#4DBBD5"),
                      na.value = NA, breaks=c("positive","purifying"),
                      labels=c("Positive","Purifying"), name="Selection") +
    coord_equal() +
    labs(title = paste("dN/dS —", fam, "/", clu), x = "Species", y = "Species") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")
}

## PLOT
combos <- selection_df_x %>%
  distinct(Family, Cluster) %>%
  arrange(Family, Cluster)



# combos and make_heatmap() assumed from earlier
plots <- lapply(seq_len(nrow(combos)), function(k) {
  make_heatmap(selection_df_x, combos$Family[k], combos$Cluster[k])
})

# 12 per page: 3 rows × 4 cols
g <- marrangeGrob(grobs = plots, nrow = 4, ncol = 3, top = NULL)
g

# save multipage PDF
ggsave("all_heatmaps_12perpage_Y_updated2_try10.pdf", g, width = 16, height = 12)  


## Zoomed in per species 
# Define species name mapping
species_names <- c(
  "PanPan" = "Bonobo",
  "PanTro" = "Chimpanzee",
  "HomSap" = "Human",
  "GorGor" = "Gorilla",
  "PonAbe" = "S. orangutan",
  "PonPyg" = "B. orangutan",
  "SymSyn" = "Siamang",
  "MacFas" = "Macaque"
)

# Modified plot function with species name mapping
make_heatmap <- function(df, fam, clu) {
  d <- df %>%
    filter(Family == fam, Cluster == clu) %>%
    mutate(
      Selection = tolower(selection),
      a = pmin(Species1, Species2),
      b = pmax(Species1, Species2)
    ) %>%
    distinct(a, b, .keep_all = TRUE)
  
  # Order species according to phylogeny
  species <- species_order[species_order %in% unique(c(d$Species1, d$Species2))]
  
  grid <- expand.grid(Species1 = species, Species2 = species, stringsAsFactors = FALSE) %>%
    mutate(
      i = match(Species1, species),
      j = match(Species2, species),
      keep_lower = i >= j,
      a = pmin(Species1, Species2),
      b = pmax(Species1, Species2)
    ) %>%
    filter(keep_lower) %>%
    left_join(d %>% select(a, b, dNdS, Selection), by = c("a","b")) %>%
    mutate(
      Selection_fill = ifelse(Selection %in% "nonsignificant", NA, Selection),
      dNdS_label     = ifelse(is.na(dNdS), "", sprintf("%.2f", dNdS)),
      # Replace species codes with actual names
      Species1_name = species_names[Species1],
      Species2_name = species_names[Species2],
      Species1_name = factor(Species1_name, levels = species_names[species]),
      Species2_name = factor(Species2_name, levels = rev(species_names[species]))
    )
  
  ggplot(grid, aes(Species1_name, Species2_name)) +
    geom_tile(aes(fill = Selection_fill), color = "grey85", linewidth = 0.3) +
    geom_text(aes(label = dNdS_label), size = 5) +
    scale_fill_manual(values = c(positive="#E64B35", purifying="#4DBBD5"),
                      na.value = NA, breaks=c("positive","purifying"),
                      labels=c("Positive","Purifying"), name="Selection") +
    coord_equal() +
    labs(title = paste("dN/dS —", fam, "/", clu), x = "Species", y = "Species") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold"))
}

# For a single gene family-cluster combination
# Specify  gene family and cluster
my_family <- "VCX"  # Family name
my_cluster <- "VCX"  # Cluster name

# Create single plot
single_plot <- make_heatmap(selection_df_x, my_family, my_cluster)
print(single_plot)

# Save single plot
ggsave(paste0("heatmap_", my_family, "_", my_cluster, ".pdf"), 
       single_plot, width = 8, height = 6)











#### Branch test ####
# Trees based on branch tests #
library(ape)
library(ggtree)
library(dplyr)
library(tibble)
library(phytools)


# 1) Read and clean the omega tree
#CT45 # Hominini as forward branch
#omega_txt <- "(((((PanTro #3.78843 , PanPan #3.78843 ) #7.38139 , HomSap #3.78843 ) #7.38139 , GorGor #1.10141 ) #1.10742 , (PonPyg #1.10141 , PonAbe #1.10141 ) #1.10742 ) #1.10742 , SymSyn #1.10141 , MacFas #1.10141 );"

# MAGEB6 # African great ape as forward branch
omega_txt <- "(((((PanTro #1.36542 , PanPan #1.36542 ) #1.35846 , HomSap #1.36542 ) #1.35846 , GorGor #1.36542 ) #1.35846 , (PonPyg #0.84137 , PonAbe #0.84137 ) #0.837795 ) #0.837795 , SymSyn #0.84137 , MacFas #0.84137 );"

# HSFX4 # Hominini as forward branch
## Gorilla is NOT PRESENT HERE
#omega_txt <- "((((PanTro #5.76397 , PanPan #5.76397 ) #4.11815 , HomSap #4.11814 ) #4.11815 , (PonPyg #1.23987 , PonAbe #1.23987 ) #1.23986 ) #1.23986 , SymSyn #1.23987 , MacFas #1.23987 );"

# DAZ1 #Hominini as forward branch
#omega_txt <- "(((((PanTro #2.90829 , PanPan #2.90829 ) #2.01915 , HomSap #2.90829 ) #2.01915 , GorGor #0.70517 ) #0.632031 , (PonPyg #0.70517 , PonAbe #0.70517 ) #0.632031 ) #0.632031 , SymSyn #0.70517 , MacFas #0.70517 );"

# clean the tree with omega values
omega_tree <- read.tree(text = gsub("\\s+", "", omega_txt))

# tip labels look like "PanTro#3.78843" -> split into name and omega
tip_omega <- as.numeric(sub(".*#", "", omega_tree$tip.label))
omega_tree$tip.label <- sub("#.*$", "", omega_tree$tip.label)

# internal node omegas are in node.label like "#7.38139"
omega_tree$node.label <- as.numeric(sub("^#", "", omega_tree$node.label))


# 2) Read the big timetree and prune
big_tree <- read.tree("~/Downloads/Craig_Kumar_Hedges_final_timetree.nwk")

keep_these <- c("Pan_troglodytes", "Pan_paniscus", "Homo_sapiens",
                "Gorilla_gorilla", "Pongo_pygmaeus", "Pongo_abelii",
                "Symphalangus_syndactylus", "Macaca_fascicularis")

time_tree <- keep.tip(big_tree, keep_these)

# root with Macaca outgroup to match your intended orientation
time_tree <- root(time_tree, outgroup="Macaca_fascicularis", resolve.root=TRUE)


# 3) Rename omega-tree tips to match timetree
name_map <- c(
  PanTro="Pan_troglodytes",
  PanPan="Pan_paniscus",
  HomSap="Homo_sapiens",
  GorGor="Gorilla_gorilla",
  PonPyg="Pongo_pygmaeus",
  PonAbe="Pongo_abelii",
  SymSyn="Symphalangus_syndactylus",
  MacFas="Macaca_fascicularis"
)
omega_tree$tip.label <- unname(name_map[omega_tree$tip.label])

# make a named vector for tip omegas in full names
tip_omega_full <- tip_omega
names(tip_omega_full) <- omega_tree$tip.label


##  3b)  CHANGE TIP LABELS  
display_labels <- c(
  "Pan_troglodytes" = "Chimpanzee",
  "Pan_paniscus" = "Bonobo",
  "Homo_sapiens" = "Human",
  "Gorilla_gorilla" = "Gorilla",
  "Pongo_pygmaeus" = "B. orangutan",
  "Pongo_abelii" = "S. orangutan",
  "Symphalangus_syndactylus" = "Siamang",
  "Macaca_fascicularis" = "Macaque"
)

# Store original labels for matching omega values
original_labels <- time_tree$tip.label

# Apply new display labels
time_tree$tip.label <- display_labels[time_tree$tip.label]

# Update tip_omega_full to use new labels
names(tip_omega_full) <- display_labels[original_labels]


# 4)Transfer internal omegas by MRCA mapping 
# Need to use original labels for MRCA matching
temp_tree <- time_tree
temp_tree$tip.label <- original_labels

internal_nodes_omega <- (Ntip(omega_tree)+1):(Ntip(omega_tree)+omega_tree$Nnode)
omega_internal_vals  <- omega_tree$node.label

# start with all NA internal labels on timetree
time_tree$node.label <- rep(NA_real_, time_tree$Nnode)

for(i in seq_along(internal_nodes_omega)){
  n_om <- internal_nodes_omega[i]
  desc <- getDescendants(omega_tree, n_om)
  desc_tips <- omega_tree$tip.label[desc[desc <= Ntip(omega_tree)]]
  
  # MRCA in time tree (using original labels)
  mrca_t <- getMRCA(temp_tree, desc_tips)
  
  # assign (if multiple omegas land on same MRCA, keep the more specific/latest one)
  idx_t <- mrca_t - Ntip(time_tree)
  time_tree$node.label[idx_t] <- omega_internal_vals[i]
}


## 5) Build plotting dataframe for BOTH tips + internal nodes
node_df <- tibble(
  node  = (Ntip(time_tree)+1):(Ntip(time_tree)+time_tree$Nnode),
  omega = time_tree$node.label
)

tip_df <- tibble(
  node  = 1:Ntip(time_tree),
  omega = tip_omega_full[time_tree$tip.label]
)

all_omega_df <- bind_rows(tip_df, node_df) %>%
  mutate(omega_str = ifelse(is.na(omega), NA, sprintf("%.2f", omega)))


## 6) Find Forward clade for coloring 
# Use original labels for finding MRCA
tips_for_color_original <- c("Pan_troglodytes", "Pan_paniscus", "Homo_sapiens", "Gorilla_gorilla") # change according to forward branch labeling
mrca_time <- getMRCA(temp_tree, tips_for_color_original)
desc_time <- sort(getDescendants(time_tree, mrca_time))

# parent node of the Homo–Pan clade
parent_node <- time_tree$edge[ time_tree$edge[,2] == mrca_time , 1 ]


##  7) FINAL PLOT WITHOUT AXIS + WITH TITLE
max_age <- max(node.depth.edgelength(time_tree))

p_final <- ggtree(time_tree) %<+% all_omega_df +
  geom_tree(aes(
    color = case_when(
      node %in% desc_time ~ "clade",
      node == mrca_time   ~ "clade",
      TRUE                ~ "other"
    )
  ), size=1.5) +
  scale_color_manual(values=c(
    clade = "red",
    other = "black"
  ), guide="none") +
  geom_text2(aes(x = branch, label=omega_str), vjust=-0.6, size=3.5, na.rm=TRUE) +
  geom_tiplab(aes(label = label),
              size=5, align=TRUE, linetype=0, offset=0.6, fontface="italic") +
  theme_tree2() +
  theme(
    plot.margin = margin(5.5, 90, 5.5, 5.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  xlim(0, max_age + 6) +
  ggtitle("MAGEB6") +
  theme(plot.title = element_text(hjust = 0.2, vjust=-15.5, size = 14, face = "bold"))

p_final


# save the plot
ggsave("~/Downloads/MAGEB6_branchtest.pdf", p_final, width = 8, height = 5, units = "in")

