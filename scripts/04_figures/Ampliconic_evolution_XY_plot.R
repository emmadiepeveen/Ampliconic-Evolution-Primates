# All panels for Figure 1 Evolution of ampliconic gene families on the X and Y chromosome


#### Ampliconic / Multicopy/ Single copy/ Absent ####

# ------------------ Packages ------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## FILES NEEDED for Both X and Y ###
# 1. Ampliconic cluster per Gene family 
# 2. Gene position data

#### Load and process data ####
counts <- read.csv("~/Downloads/cluster_counts_perspecies_x_updated.tsv", sep = "\t") %>%
  mutate(Family = case_when(
    Family == "ARMCX1" ~ "ARMCX",
    Family == "CENPVL2" ~ "CENPVL2",
    Family == "CSAG2" ~ "CSAG",
    Family == "CSF2RA" ~ "CSF2RA",
    Family == "CT45A2" ~ "CT45",
    Family == "CT47C1" ~ "CT47",
    Family == "CXorf49" ~ "CXorf49",
    Family == "CXorf51A" ~ "CXorf51",
    Family == "DMRTC1B" ~ "DMRTC1",
    Family == "EOLA1" ~ "EOLA",
    Family == "ETDB" ~ "ETD",
    Family == "F8A1" ~ "F8",
    Family == "FAM156A" ~ "FAM156",
    Family == "FAM236C" ~ "FAM236",
    Family == "GAGE12F" ~ "GAGE",
    Family == "GPRASP2" ~ "GPRASP1",
    Family == "H2AB2" ~ "H2AB",
    Family == "H2BW1" ~ "H2BW",
    Family == "HSFX2" ~ "HSFX1",
    Family == "HSFX4" ~ "HSFX4",
    Family == "IKBKG" ~ "IKBKG",
    Family == "INTS6L" ~ "SAGE",
    Family == "LAGE3" ~ "CTAG",
    Family == "MAGEB16" ~ "MAGEB",
    Family == "MAGED" ~ "MAGED",
    Family == "NUDT10" ~ "NUDT",
    Family == "NXF3" ~ "NXF",
    Family == "PABPC1L2B" ~ "PABPC",
    Family == "OPN1LW" ~ "OPN1LW",
    Family == "PNMA6A" ~ "PNMA",
    Family == "PWWP4" ~ "PWWP4",
    Family == "RAB40A" ~ "RAB40A",
    Family == "RHOXF2" ~ "RHOXF2B",
    Family == "RPL36A" ~ "RPL36A-HNRNPH2",
    Family == "SMIM10" ~ "SMIM10",
    Family == "SPACA5" ~ "SPACA",
    Family == "SPANXA1" ~ "SPANX",
    Family == "SPIN4" ~ "SPIN",
    Family == "SSX4" ~ "SSX",
    Family == "TBL1X" ~ "TBL1X",
    Family == "TCEAL8" ~ "TCEAL8",
    Family == "TCP11X2" ~ "TCP11X",
    Family == "TEX28" ~ "TEX28",
    Family == "TMEM185A" ~ "TMEM185",
    Family == "TMSB15A" ~ "TMSB15",
    Family == "VCX3A" ~ "VCX",
    Family == "XAGE1A" ~ "XAGE",
    Family == "ZXDB" ~ "ZXD",
    Family == "collagen alpha-4(IV) chain-like" ~ "collagen alpha-4(IV) chain-like",
    Family == "endogenous retrovirus group K member 6 Env polyprotein-like" ~ "retrovirus K6",
    Family == "putative uncharacterized protein FLJ39060" ~ "FLJ39060",
    Family == "uncharacterized LOC115932372" ~ "LOC129138873",
    Family == "uncharacterized LOC129475109" ~ "LOC129475109",
    TRUE ~ Family
  ))

counts_y <- read.csv("~/Downloads/cluster_counts_perspecies_y_updated.tsv", sep = "\t") %>%
  mutate(Family = case_when(
    Family == "BPY2" ~ "BPY2",
    Family == "CDY1" ~ "CDY",
    Family == "DAZ1" ~ "DAZ",
    Family == "HSFY1" ~ "HSFY",
    Family == "MTRNR2-like 17" ~ "MTRNR2-like 17",
    Family == "RBMY1B" ~ "RBMY",
    Family == "TATA-box binding protein associated factor 11 like protein 2-like" ~ "TAF11L2",
    Family == "TSPY8" ~ "TSPY",
    Family == "VCY1B" ~ "VCY",
    Family == "adenylate kinase isoenzyme 6-like" ~ "AKI6",
    Family == "centriole and centriolar satellite protein OFD1-like" ~ "CCS-OFD1",
    Family == "endogenous retrovirus group K member 19 Env polyprotein-like" ~ "retrovirus K19",
    Family == "glutamate dehydrogenase 1, mitochondrial-like" ~ "GLUD1Y",
    Family == "keratin, type I cytoskeletal 18-like" ~ "KRT18Y",
    Family == "proline-rich protein, Y-linked" ~ "PRPY",
    Family == "protein FAM47A-like" ~ "FAM47AY",
    Family == "protein FRG1-like" ~ "FRG1Y",
    Family == "zinc finger protein 285-like" ~ "ZNF",
    TRUE ~ Family
  ))

#### Process status classification ####
species_cols <- setdiff(names(counts), c("Family", "Cluster"))

process_status <- function(data, species_cols) {
  data %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "Species",
      values_to = "count"
    ) %>%
    group_by(Family, Species) %>%
    summarise(
      any_ge2 = any(count >= 2, na.rm = TRUE),
      ones    = sum(count == 1, na.rm = TRUE),
      zeros   = sum(count == 0, na.rm = TRUE),
      nonzero = sum(count > 0, na.rm = TRUE),
      total   = n(),
      .groups = "drop"
    ) %>%
    mutate(Status = case_when(
      any_ge2 ~ "Ampliconic",
      ones >= 2 ~ "Multi-copy",
      ones == 1 & zeros == total - 1 ~ "Single-copy",
      nonzero == 0 ~ "Absent",
      TRUE ~ "Mixed"
    )) %>%
    select(Family, Species, Status) %>%
    pivot_wider(names_from = Species, values_from = Status)
}

result <- process_status(counts, species_cols) %>% mutate(Chromosome = "X")
result_y <- process_status(counts_y, setdiff(names(counts_y), c("Family", "Cluster"))) %>% mutate(Chromosome = "Y")

#### Load gene position data ####
xgenes <- read.csv("~/Downloads/all_families_gene_details_with_clusters_X.tsv", sep = "\t") %>%
  mutate(gene_family_symbol = case_when(
    gene_family_symbol == "ARMCX1" ~ "ARMCX",
    gene_family_symbol == "CENPVL2" ~ "CENPVL2",
    gene_family_symbol == "CSAG2" ~ "CSAG",
    gene_family_symbol == "CSF2RA" ~ "CSF2RA",
    gene_family_symbol == "CT45A2" ~ "CT45",
    gene_family_symbol == "CT47C1" ~ "CT47",
    gene_family_symbol == "CXorf49" ~ "CXorf49",
    gene_family_symbol == "TRO" ~ "MAGED",
    gene_family_symbol == "CXorf51A" ~ "CXorf51",
    gene_family_symbol == "DMRTC1B" ~ "DMRTC1",
    gene_family_symbol == "EOLA1" ~ "EOLA",
    gene_family_symbol == "ETDB" ~ "ETD",
    gene_family_symbol == "F8A1" ~ "F8",
    gene_family_symbol == "FAM156A" ~ "FAM156",
    gene_family_symbol == "FAM236C" ~ "FAM236",
    gene_family_symbol == "GAGE12F" ~ "GAGE",
    gene_family_symbol == "GPRASP2" ~ "GPRASP1",
    gene_family_symbol == "H2AB2" ~ "H2AB",
    gene_family_symbol == "H2BW1" ~ "H2BW",
    gene_family_symbol == "HSFX2" ~ "HSFX1",
    gene_family_symbol == "HSFX4" ~ "HSFX4",
    gene_family_symbol == "IKBKG" ~ "IKBKG",
    gene_family_symbol == "INTS6L" ~ "SAGE",
    gene_family_symbol == "LAGE3" ~ "CTAG",
    gene_family_symbol == "MAGEB16" ~ "MAGEB",
    gene_family_symbol == "TRO" ~ "MAGED",
    gene_family_symbol == "NUDT10" ~ "NUDT",
    gene_family_symbol == "NXF3" ~ "NXF",
    gene_family_symbol == "PABPC1L2B" ~ "PABPC",
    gene_family_symbol == "OPN1LW" ~ "OPN1LW",
    gene_family_symbol == "PNMA6A" ~ "PNMA",
    gene_family_symbol == "PWWP4" ~ "PWWP4",
    gene_family_symbol == "RAB40A" ~ "RAB40A",
    gene_family_symbol == "RHOXF2" ~ "RHOXF2B",
    gene_family_symbol == "RPL36A" ~ "RPL36A-HNRNPH2",
    gene_family_symbol == "SMIM10" ~ "SMIM10",
    gene_family_symbol == "SPACA5" ~ "SPACA",
    gene_family_symbol == "SPANXA1" ~ "SPANX",
    gene_family_symbol == "SPIN4" ~ "SPIN",
    gene_family_symbol == "SSX4" ~ "SSX",
    gene_family_symbol == "TBL1X" ~ "TBL1X",
    gene_family_symbol == "TCEAL8" ~ "TCEAL8",
    gene_family_symbol == "TCP11X2" ~ "TCP11X",
    gene_family_symbol == "TEX28" ~ "TEX28",
    gene_family_symbol == "TMEM185A" ~ "TMEM185",
    gene_family_symbol == "TMSB15A" ~ "TMSB15",
    gene_family_symbol == "VCX3A" ~ "VCX",
    gene_family_symbol == "XAGE1A" ~ "XAGE",
    gene_family_symbol == "ZXDB" ~ "ZXD",
    gene_family_symbol == "collagen alpha-4(IV) chain-like" ~ "collagen alpha-4(IV) chain-like",
    gene_family_symbol == "endogenous retrovirus group K member 6 Env polyprotein-like" ~ "retrovirus K6",
    gene_family_symbol == "putative uncharacterized protein FLJ39060" ~ "FLJ39060",
    gene_family_symbol == "uncharacterized LOC115932372" ~ "LOC129138873",
    gene_family_symbol == "uncharacterized LOC129475109" ~ "LOC129475109",
    TRUE ~ gene_family_symbol
  ))

ygenes <- read.csv("~/Downloads/all_families_gene_details_with_clusters_Y.tsv", sep = "\t") %>%
  mutate(gene_family_symbol = case_when(
    gene_family_symbol == "BPY2" ~ "BPY2",
    gene_family_symbol == "CDY1" ~ "CDY",
    gene_family_symbol == "DAZ1" ~ "DAZ",
    gene_family_symbol == "HSFY1" ~ "HSFY",
    gene_family_symbol == "MTRNR2-like 17" ~ "MTRNR2-like 17",
    gene_family_symbol == "RBMY1B" ~ "RBMY",
    gene_family_symbol == "TATA-box binding protein associated factor 11 like protein 2-like" ~ "TAF11L2",
    gene_family_symbol == "TSPY8" ~ "TSPY",
    gene_family_symbol == "VCY1B" ~ "VCY",
    gene_family_symbol == "adenylate kinase isoenzyme 6-like" ~ "AKI6",
    gene_family_symbol == "centriole and centriolar satellite protein OFD1-like" ~ "CCS-OFD1",
    gene_family_symbol == "endogenous retrovirus group K member 19 Env polyprotein-like" ~ "retrovirus K19",
    gene_family_symbol == "glutamate dehydrogenase 1, mitochondrial-like" ~ "GLUD1Y",
    gene_family_symbol == "keratin, type I cytoskeletal 18-like" ~ "KRT18Y",
    gene_family_symbol == "proline-rich protein, Y-linked" ~ "PRPY",
    gene_family_symbol == "protein FAM47A-like" ~ "FAM47AY",
    gene_family_symbol == "protein FRG1-like" ~ "FRG1Y",
    gene_family_symbol == "zinc finger protein 285-like" ~ "ZNF",
    TRUE ~ gene_family_symbol
  ))

#### Calculate X chromosome strata ordering ####
breaks <- c(0, 2394410, 3333048, 10875912, 47005769, 57347714, 61248490, Inf)
labels <- c("PAR", "Strata 5", "Strata 4", "Strata 3", "Strata 2",
            "Centromere Satellite", "Strata 1")

x_rep <- xgenes %>%
  mutate(midpoint = (Start + End) / 2) %>%
  group_by(Species, gene_family_symbol) %>%
  summarise(pos = max(midpoint, na.rm = TRUE), .groups = "drop")

species_scale <- x_rep %>%
  group_by(Species) %>%
  summarise(min_pos = min(pos, na.rm = TRUE),
            max_pos = max(pos, na.rm = TRUE), .groups = "drop")

x_rel <- x_rep %>%
  left_join(species_scale, by = "Species") %>%
  mutate(rel_pos = (pos - min_pos) / (max_pos - min_pos))

human_scale <- species_scale %>% filter(Species == "HomSap")
human_min <- human_scale$min_pos
human_max <- human_scale$max_pos

x_rel <- x_rel %>%
  mutate(human_coord = human_min + rel_pos * (human_max - human_min),
         Stratum = cut(human_coord, breaks = breaks, labels = labels,
                       right = FALSE, include.lowest = TRUE))

family_strata_X <- x_rel %>%
  group_by(gene_family_symbol) %>%
  summarise(
    n_species_present = n_distinct(Species),
    first_human_coord = min(human_coord, na.rm = TRUE),
    representative_stratum = cut(first_human_coord, breaks = breaks,
                                 labels = labels, right = FALSE,
                                 include.lowest = TRUE),
    .groups = "drop"
  )

strata_order <- c("PAR", "Strata 5", "Strata 4", "Strata 3", "Strata 2",
                  "Centromere Satellite", "Strata 1")

family_order_X <- family_strata_X %>%
  mutate(representative_stratum = factor(representative_stratum, levels = strata_order)) %>%
  arrange(representative_stratum, first_human_coord) %>%
  pull(gene_family_symbol)

#### Calculate Y chromosome ordering ####
y_rep <- ygenes %>%
  mutate(midpoint = (Start + End) / 2) %>%
  group_by(Species, gene_family_symbol) %>%
  summarise(pos = max(midpoint, na.rm = TRUE), .groups = "drop")

y_bounds <- y_rep %>%
  group_by(Species) %>%
  summarise(min_pos = min(pos, na.rm = TRUE),
            max_pos = max(pos, na.rm = TRUE), .groups = "drop")

y_rep_rel <- y_rep %>%
  left_join(y_bounds, by = "Species") %>%
  mutate(rel_pos = if_else(max_pos > min_pos, (pos - min_pos) / (max_pos - min_pos), 0.5)) %>%
  select(Species, gene_family_symbol, rel_pos)

y_hom <- y_rep_rel %>%
  filter(Species == "HomSap") %>%
  transmute(gene_family_symbol, hom_rel_pos = rel_pos)

y_cons <- y_rep_rel %>%
  filter(Species != "HomSap") %>%
  group_by(gene_family_symbol) %>%
  summarise(consensus_rel_pos = median(rel_pos, na.rm = TRUE),
            present_in_others = n(), .groups = "drop")

y_anchors <- full_join(y_hom, y_cons, by = "gene_family_symbol") %>%
  mutate(anchor_rel_pos = if_else(!is.na(hom_rel_pos), hom_rel_pos, consensus_rel_pos),
         present_in = if_else(!is.na(hom_rel_pos), Inf, present_in_others))

family_order_Y <- y_anchors %>%
  arrange(anchor_rel_pos, desc(present_in), gene_family_symbol) %>%
  pull(gene_family_symbol)

#### Create plotting data with correct ordering ####
combined <- bind_rows(result, result_y)

plot_data <- combined %>%
  pivot_longer(cols = -c(Family, Chromosome), names_to = "Species", values_to = "Status") %>%
  mutate(Family_plot = ifelse(Chromosome == "X",
                              paste0("X::", Family),
                              paste0("Y::", Family)))

# Build complete factor levels
levels_X <- paste0("X::", family_order_X)
levels_Y <- paste0("Y::", family_order_Y)

plot_data <- plot_data %>%
  mutate(Family_plot = factor(Family_plot, levels = c(levels_X, levels_Y)))

# Create label map
lvl <- levels(plot_data$Family_plot)
label_map <- setNames(gsub("^(X|Y)::", "", lvl, perl = TRUE), lvl)

# Species order
species_order <- c("MacFas", "SymSyn", "PonPyg", "PonAbe", "GorGor", "HomSap", "PanTro", "PanPan")
plot_data <- plot_data %>%
  mutate(Species = factor(Species, levels = species_order))

#### Create main tile plot ####
status_cols <- c(
  "Ampliconic"  = "darkblue",
  "Multi-copy"  = "#4682B4",#3B73B9
  "Single-copy" = "#808080", #E6A141
  "Absent"      = "white"
)

p <- ggplot(plot_data, aes(x = Family_plot, y = Species, fill = Status)) +
  geom_tile(width = 0.92, height = 0.92) +
  scale_fill_manual(values = status_cols, na.value = "white", drop = FALSE) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_discrete(position = "top", labels = label_map) +
  labs(x = "Gene family", y = "Species", fill = "Class") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1, size = 8, face = "bold"),
    strip.placement = "inside",
    strip.text.x = element_text(face = "bold", size = 7),
    strip.background.x = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.spacing.x = unit(8, "pt")
  )

#### Create strata annotation ####
x_levels <- levels(plot_data$Family_plot)[grepl("^X::", levels(plot_data$Family_plot))]
y_levels <- levels(plot_data$Family_plot)[grepl("^Y::", levels(plot_data$Family_plot))]
nX <- length(x_levels)
nY <- length(y_levels)

strata_tiles_X <- tibble(Family_plot = x_levels) %>%
  mutate(Family = sub("^X::", "", Family_plot), Chromosome = "X") %>%
  left_join(
    family_strata_X %>%
      transmute(Family = gene_family_symbol,
                Stratum = as.character(representative_stratum)),
    by = "Family"
  )

r <- rle(strata_tiles_X$Stratum)
seg_df <- tibble(
  Stratum = r$values,
  len = r$lengths,
  start = cumsum(c(1, head(r$lengths, -1))),
  end = cumsum(r$lengths)
) %>%
  filter(!is.na(Stratum)) %>%
  mutate(xmin_raw = start - 0.5, xmax_raw = end + 0.5)

gap <- 0.18
seg_df <- seg_df %>%
  mutate(
    xmin = xmin_raw + gap / 2,
    xmax = xmax_raw - gap / 2,
    xmid = (start + end) / 2,
    label = gsub("^Strata\\s*", "S", Stratum) |>
      gsub("^Centromere Satellite$", "Cent.", x = _),
    Chromosome = "X"
  )

extent_df <- bind_rows(
  tibble(Chromosome = "X", x = c(0.5, nX + 0.5), y = 1),
  tibble(Chromosome = "Y", x = c(0.5, nY + 0.5), y = 1)
)

strata_cols <- c(
  "PAR"                  = "#636363",
  "Strata 5"             = "#2166ac",
  "Strata 4"             = "#92c5de",
  "Strata 3"             = "#80cdc1",
  "Strata 2"             = "#0571b0",
  "Strata 1"             = "#a6611a",
  "Centromere Satellite" = "#bdbdbd"
)

p_strata <- ggplot() +
  geom_blank(data = extent_df, aes(x = x, y = y)) +
  geom_segment(
    data = seg_df,
    aes(x = xmin, xend = xmax, y = 1, yend = 1, color = Stratum),
    linewidth = 2, lineend = "butt"
  ) +
  geom_text(
    data = seg_df,
    aes(x = xmid, y = 0.7, label = label),
    size = 3.2, fontface = "bold", vjust = 1
  ) +
  scale_color_manual(values = strata_cols, guide = "none") +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(expand = expansion(0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 11) +
  theme(
    panel.spacing.x = unit(8, "pt"),
    plot.margin = margin(t = 0, r = 5, b = 6, l = 5),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

#### Combine plots ####
p_final <- p / p_strata + patchwork::plot_layout(heights = c(10, 1))

#### Save plot ####
n_species <- length(unique(plot_data$Species))
n_families_X <- sum(combined$Chromosome == "X")
n_families_Y <- sum(combined$Chromosome == "Y")
total_families <- n_families_X + n_families_Y

tile_size <- 0.4
plot_width <- total_families * tile_size + 1.7
plot_height <- n_species * tile_size + 2

ggsave("plot_withstrata_final4.0.jpg", p_final, width = plot_width, height = plot_height, dpi = 400)

p_final


### X and Y chromosome separately ####
#### Create species name mapping ####
species_name_map <- c(
  "MacFas" = "Macaque",
  "SymSyn" = "Siamang",
  "PonPyg" = "B. orangutan",
  "PonAbe" = "S. orangutan",
  "GorGor" = "Gorilla",
  "HomSap" = "Human",
  "PanTro" = "Chimpanzee",
  "PanPan" = "Bonobo"
)

#### Create separate X and Y plots ####

# Filter data for X chromosome only
plot_data_X <- plot_data %>%
  filter(Chromosome == "X")

# Filter data for Y chromosome only
plot_data_Y <- plot_data %>%
  filter(Chromosome == "Y")

#### X chromosome plot ####
p_X <- ggplot(plot_data_X, aes(x = Family_plot, y = Species, fill = Status)) +
  geom_tile(width = 0.92, height = 0.92, color = "black", size = 0.3) +
  scale_fill_manual(values = status_cols, na.value = "white", drop = FALSE) +
  scale_x_discrete(position = "top", labels = label_map) +
  scale_y_discrete(labels = species_name_map) + 
  labs(x = NULL, y = NULL, fill = "Class") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 15), #face = "bold"
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 15),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 15),      # Bigger title
    legend.text = element_text(size = 14),                       # Bigger text
    legend.key.size = unit(0.8, "cm"),                          # Bigger color boxes
    #plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    
  )

#### Create X-specific strata annotation ####
extent_df_X <- tibble(Chromosome = "X", x = c(0.5, nX + 0.5), y = 1)

p_strata_X <- ggplot() +
  geom_blank(data = extent_df_X, aes(x = x, y = y)) +
  geom_segment(
    data = seg_df,
    aes(x = xmin, xend = xmax, y = 1, yend = 1, color = Stratum),
    linewidth = 2, lineend = "butt"
  ) +
  geom_text(
    data = seg_df,
    aes(x = xmid, y = 0.7, label = label),
    size = 4.2, fontface = "bold", vjust = 1
  ) +
  scale_color_manual(values = strata_cols, guide = "none") +
  scale_x_continuous(expand = expansion(0), limits = c(0.5, nX + 0.5)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(t = 0, r = 5, b = 6, l = 5),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Combine X plot with X-specific strata annotation
p_X_final <- p_X / p_strata_X + patchwork::plot_layout(heights = c(10, 1))
p_X_final

#### Y chromosome plot ####
width_ratio <- n_families_X / n_families_Y

p_Y <- ggplot(plot_data_Y, aes(x = Family_plot, y = Species, fill = Status)) +
  geom_tile(width = 0.92, height = 0.92, color = "black", size = 0.3) +
  scale_fill_manual(values = status_cols, na.value = "white", drop = FALSE) +
  scale_x_discrete(position = "top", labels = label_map) +
  scale_y_discrete(labels = species_name_map) + 
  labs(x = NULL, y = NULL, fill = "Class") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 15), #face = "italic"
    legend.position = "none",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# Y chromosome doesn't need strata annotation, but you can add a spacer if you want alignment
p_Y_final <- p_Y

#### Save separate plots ####
# X chromosome
n_families_X <- sum(plot_data_X$Chromosome == "X" & !duplicated(plot_data_X$Family_plot))
plot_width_X <- n_families_X * tile_size + 1.7
ggsave("~/Downloads/plot_X_chromosome_strata_final.jpg", p_X_final, width = plot_width_X, height = plot_height, dpi = 400)

# Y chromosome
n_families_Y <- sum(plot_data_Y$Chromosome == "Y" & !duplicated(plot_data_Y$Family_plot))
plot_width_Y <- n_families_Y * tile_size + 1.7
ggsave("~/Downloads/plot_Y_chromosome_tiles.jpg", p_Y_final, width = plot_width_Y, height = plot_height, dpi = 400)

# View the plots
p_X_final
p_Y_final

#### Save plots with identical tile sizes (both width AND height) ####


# Save -  adjust margins to account for legend space
plot_width_X <- (n_families_X * fixed_tile_width) + margin_width + 2.5  # Extra for legend
plot_width_Y <- (n_families_Y * fixed_tile_width) + margin_width


ggsave("~/Downloads/plot_X_chromosome.pdf", p_X_final, width = plot_width_X, height = plot_height, units = "in")
ggsave("~/Downloads/plot_Y_chromosome.pdf", p_Y_final, width = plot_width_Y, height = plot_height, units = "in")




#### STACKED bar plot of proportions ####

library(dplyr)
library(tidyr)
library(ggplot2)

#### Prepare data for stacked barplots ####
# For X chromosome
bar_data_X <- plot_data %>%
  filter(Chromosome == "X") %>%
  group_by(Species, Status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Species) %>%
  mutate(
    total = sum(count),
    proportion = count / total  # Changed from * 100 to get proportion instead of percentage
  ) %>%
  ungroup() %>%
  mutate(Chromosome = "X")

# For Y chromosome
bar_data_Y <- plot_data %>%
  filter(Chromosome == "Y") %>%
  group_by(Species, Status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Species) %>%
  mutate(
    total = sum(count),
    proportion = count / total  # Changed from * 100 to get proportion instead of percentage
  ) %>%
  ungroup() %>%
  mutate(Chromosome = "Y")

# Define species order
species_order <- c("MacFas", "SymSyn", "PonPyg", "PonAbe", "GorGor", "HomSap", "PanTro", "PanPan")

# Make sure Status is a factor with the right order
bar_data_X <- bar_data_X %>%
  mutate(Status = factor(Status, levels = c("Absent", "Single-copy", "Multi-copy", "Ampliconic")))

bar_data_Y <- bar_data_Y %>%
  mutate(Status = factor(Status, levels = c("Absent", "Single-copy", "Multi-copy", "Ampliconic")))

# Apply species order and create ordered common names
bar_data_X <- bar_data_X %>%
  mutate(
    Species = factor(Species, levels = species_order),
    Species_name = species_name_map[as.character(Species)],
    Species_name = factor(Species_name, levels = species_name_map[species_order])
  )

bar_data_Y <- bar_data_Y %>%
  mutate(
    Species = factor(Species, levels = species_order),
    Species_name = species_name_map[as.character(Species)],
    Species_name = factor(Species_name, levels = species_name_map[species_order])
  )

# Create plots with proportion (0-1) instead of percentage
p_bar_X <- ggplot(bar_data_X, aes(x = Species_name, y = proportion, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = status_cols, drop = FALSE) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +  # Changed to 0-1
  labs(x = NULL, y = "Proportion of gene family classes", fill = "Class", title = "X Chromosome") +  # Changed label
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, hjust = 1, vjust = 1),
    axis.text.y = element_text(size=15),
    #axis.title.y = element_text(face = "bold", size = 12),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  )+
  coord_flip()

p_bar_Y <- ggplot(bar_data_Y, aes(x = Species_name, y = proportion, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = status_cols, drop = FALSE) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +  # Changed to 0-1
  labs(x = NULL, y = "Proportion of gene family classes", fill = "Class", title = "Y Chromosome") +  # Changed label
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, hjust = 1, vjust = 1),
    axis.text.y = element_text(size=15),
    #axis.title.y = element_text(face = "bold", size = 12),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  ) +
  coord_flip()

p_bar_X
p_bar_Y

#### Save barplots ####
ggsave("~/Downloads/barplot_X_chromosome.pdf", p_bar_X, width = 10, height = 6, units = "in")
ggsave("~/Downloads/barplot_Y_chromosome.pdf", p_bar_Y, width = 10, height = 6, units = "in")



# Plot for Overview ampliconic clusters

### BUBBLE PLOTS ###
# ------------------ Packages ------------------
library(readxl)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(readr)

# common species names
species_name_map <- c(
  "MacFas" = "Macaque",
  "SymSyn" = "Siamang",
  "PonPyg" = "B.orangutan",
  "PonAbe" = "S.orangutan",
  "GorGor" = "Gorilla",
  "HomSap" = "Human",
  "PanTro" = "Chimpanzee",
  "PanPan" = "Bonobo"
)

## FILES NEEDED for Both X and Y ###
# 1. Ampliconic Gene family summary, containing expression in Sertoli, leydig, Spermatogenesis, and copy number counts per species
# 2. Ampliconic cluster per Gene family 
# 3. dNdS data per Gene Family, on both dN/dS of phylogeny, as site-specific test


#### X chromosome ####

# ------------------ Inputs ------------------
# 1) Family-level summary (wide)
ampl_summary <- read_excel("~/Documents/PhD Aarhus/Amplicons/Final scripts/Amplicon_summary_count_expression_genefamily.xlsx") %>%
  filter(Chromosome == "X") %>%
  mutate(Gene_family = `Gene family`)

# 2) Cluster-level counts
cluster_counts_perspecies_y <- read_delim("~/Downloads/cluster_counts_perspecies_x_updated.tsv", 
                                          delim = "\t") %>%
  mutate(Family = case_when(
    Family == "ARMCX1" ~ "ARMCX",
    Family == "CENPVL2" ~ "CENPVL2",
    Family == "CSAG2" ~ "CSAG",
    Family == "CSF2RA" ~ "CSF2RA",
    Family == "CT45A2" ~ "CT45",
    Family == "CT47C1" ~ "CT47",
    Family == "CXorf49" ~ "CXorf49",
    Family == "MAGED" ~ "MAGED",
    Family == "CXorf51A" ~ "CXorf51",
    Family == "DMRTC1B" ~ "DMRTC1",
    Family == "EOLA1" ~ "EOLA",
    Family == "ETDB" ~ "ETD",
    Family == "F8A1" ~ "F8",
    Family == "FAM156A" ~ "FAM156",
    Family == "FAM236C" ~ "FAM236",
    Family == "GAGE12F" ~ "GAGE",
    Family == "GPRASP2" ~ "GPRASP1",
    Family == "H2AB2" ~ "H2AB",
    Family == "H2BW1" ~ "H2BW",
    Family == "HSFX2" ~ "HSFX1",
    Family == "HSFX4" ~ "HSFX4",
    Family == "IKBKG" ~ "IKBKG",
    Family == "INTS6L" ~ "SAGE",
    Family == "LAGE3" ~ "CTAG",
    Family == "MAGEB16" ~ "MAGEB",
    Family == "NUDT10" ~ "NUDT",
    Family == "NXF3" ~ "NXF",
    Family == "PABPC1L2B" ~ "PABPC",
    Family == "OPN1LW" ~ "OPN1LW",
    Family == "PNMA6A" ~ "PNMA",
    Family == "PWWP4" ~ "PWWP4",
    Family == "RAB40A" ~ "RAB40A",
    Family == "RHOXF2" ~ "RHOXF2B",
    Family == "RPL36A" ~ "RPL36A-HNRNPH2",
    Family == "SMIM10" ~ "SMIM10",
    Family == "SPACA5" ~ "SPACA",
    Family == "SPANXA1" ~ "SPANX",
    Family == "SPIN4" ~ "SPIN",
    Family == "SSX4" ~ "SSX",
    Family == "TBL1X" ~ "TBL1X",
    Family == "TCEAL8" ~ "TCEAL8",
    Family == "TCP11X2" ~ "TCP11X",
    Family == "TEX28" ~ "TEX28",
    Family == "TMEM185A" ~ "TMEM185A",
    Family == "TMSB15A" ~ "TMSB15B",
    Family == "MAGED" ~ "MAGED",
    Family == "VCX3A" ~ "VCX",
    Family == "XAGE1A" ~ "XAGE",
    Family == "ZXDB" ~ "ZXDB",
    Family == "collagen alpha-4(IV) chain-like" ~ "collagen alpha-4(IV) chain-like",
    Family == "endogenous retrovirus group K member 6 Env polyprotein-like" ~ "endogenous retrovirus group K member 113 Env polyprotein-like",
    Family == "putative uncharacterized protein FLJ39060" ~ "FLJ39060",
    Family == "uncharacterized LOC115932372" ~ "LOC129138873",
    Family == "uncharacterized LOC129475109" ~ "LOC129475109"))

# 3) --- dN/dS ---
dnds <- read_csv("~/Downloads/codeml_site_tests_median_omega_updated_X copy.csv")%>%
  dplyr::rename(Gene_family = Family) %>%
  mutate(Cluster = sub(".*_cluster_", "", Cluster)) %>%
  mutate(
    omega = coalesce(`Median M0 omega`),
    omega = ifelse(is.infinite(omega), NA, omega),
    site_test = factor(`site-test`, levels = c("No","Yes"))
  ) %>%
  select(Gene_family, Cluster, omega, site_test) %>%
  mutate(Gene_family = recode(Gene_family,
                              "TMSB"   = "TMSB15B",
                              "H2A"   = "H2AB",
                              "endogenous"   = "endogenous retrovirus group K member 113 Env polyprotein-like",
                              "SPACA5"   = "SPACA",
                              "XAGE1" ="XAGE",
                              "NUDT10"   = "NUDT",
                              "CENPVL"   = "CENPVL2",
                              "FLJ39060"   = "FLJ39060",
                              "RHOXF2" ="RHOXF2B",
                              "ZXD"   = "ZXDB",
                              "RPL36A"   = "RPL36A-HNRNPH2",
                              "TCP11X2"   = "TCP11X",
                              "GPRASP"   = "GPRASP1",
                              "INTS6L"   = "SAGE",
                              "CT45A"   = "CT45",
                              "LAGE3"   = "CTAG",
                              "F8A1"   = "F8",
                              "MAGED1" = "MAGED",
                              "LOC115932372"   = "LOC129138873",
                              "collagen"   = "collagen alpha-4(IV) chain-like",
                              "LOC129475109"   = "LOC129475109",
                              "HSFX" ="HSFX4"))

# Species columns + display order
species_cols  <- c("GorGor", "HomSap", "PanPan", "PanTro", "MacFas", "PonAbe", "PonPyg", "SymSyn")
species_order <- c("PanPan", "PanTro", "HomSap", "GorGor", "PonAbe", "PonPyg", "SymSyn", "MacFas")
column_order  <- c(species_order, "Sertoli", "Leydig", "Spermatogenesis", "dN/dS", "Site-test")

# ------------------ Family annotations ------------------
ampl_sl_long <- ampl_summary %>%
  select(Gene_family, Sertoli, Leydig) %>%
  pivot_longer(c(Sertoli, Leydig), names_to = "CellType", values_to = "Present") %>%
  mutate(CellType = factor(CellType, levels = column_order))

ampl_sp <- ampl_summary %>%
  select(Gene_family, Spermatogenesis) %>%
  mutate(Spermatogenesis = factor(Spermatogenesis,
                                  levels = c("Pre MSCI","Pre + post MSCI","Post MSCI","Unknown")))

gene_levels <- ampl_sp %>%
  distinct(Gene_family, Spermatogenesis) %>%
  arrange(Spermatogenesis, Gene_family) %>%
  pull(Gene_family)

ampl_sl_long$Gene_family <- factor(ampl_sl_long$Gene_family, levels = gene_levels)
ampl_sp$Gene_family      <- factor(ampl_sp$Gene_family,      levels = gene_levels)

competition_df <- tibble::tibble(
  Species = factor(species_order, levels = column_order),
  Competition = c("high","high","low","low","intermediate","intermediate","low","intermediate"),
  Gene_family = factor("Sperm competition", levels = c("Sperm competition", gene_levels))
)

# ------------------ Filter clusters ------------------
drops_all_singlecopy <- cluster_counts_perspecies_y %>%
  mutate(
    species_present     = rowSums(across(all_of(species_cols), ~ (.x > 0) & !is.na(.x))),
    nonzero_eq_one      = rowSums(across(all_of(species_cols), ~ (.x == 1) & !is.na(.x))),
    all_nonzero_are_one = (species_present > 0) & (nonzero_eq_one == species_present)
  ) %>%
  filter(all_nonzero_are_one) %>%
  transmute(Gene_family = Family, Cluster)

# ------------------ Long cluster table ------------------
cluster_long <- cluster_counts_perspecies_y %>%
  dplyr::rename(Gene_family = Family) %>%
  dplyr::anti_join(drops_all_singlecopy, by = c("Gene_family", "Cluster")) %>%
  pivot_longer(all_of(species_cols), names_to = "Species", values_to = "Count") %>%
  dplyr::filter(!is.na(Count), Count != 0) %>%
  dplyr::mutate(
    Species = factor(Species, levels = column_order),
    UniqueID = paste0(Gene_family, "|||", Cluster),
    ClusterDisplay = Cluster
  )

dnds <- dnds %>%
  mutate(UniqueID = paste0(Gene_family, "|||", Cluster),
         ClusterDisplay = Cluster)

# ------------------ Build y-axis ------------------
fam_order <- gene_levels[gene_levels %in% unique(cluster_long$Gene_family)]

make_block <- function(f) {
  labs <- cluster_long %>%
    filter(Gene_family == f) %>%
    distinct(UniqueID, ClusterDisplay) %>%
    arrange(ClusterDisplay) %>%
    pull(UniqueID)
  c(labs, paste0(f, "___SPACER___"))
}
y_levels <- c("Sperm competition", unlist(lapply(fam_order, make_block), use.names = FALSE))

cluster_display_map <- cluster_long %>%
  distinct(UniqueID, ClusterDisplay) %>%
  deframe()

family_brackets <- cluster_long %>%
  distinct(Gene_family, UniqueID, ClusterDisplay) %>%
  mutate(
    y_pos = match(UniqueID, y_levels),
    Gene_family = as.character(Gene_family)
  ) %>%
  group_by(Gene_family) %>%
  summarise(
    y_min = min(y_pos),
    y_max = max(y_pos),
    y_mid = (min(y_pos) + max(y_pos)) / 2,
    .groups = "drop"
  )

label_vec <- setNames(
  sapply(y_levels, function(x) {
    if (grepl("___SPACER___", x, fixed = TRUE)) {
      return("")
    } else if (x == "Sperm competition") {
      return("Sperm competition")
    } else if (x %in% names(cluster_display_map)) {
      return(cluster_display_map[x])
    } else {
      return(x)
    }
  }),
  y_levels
)

# numeric x position for each real column
x_levels <- column_order
x_map    <- setNames(seq_along(x_levels), x_levels)

cluster_bubbles <- cluster_long %>%
  transmute(
    x     = x_map[as.character(Species)],
    y_lab = factor(UniqueID, levels = y_levels),
    Count, Species
  )

comp_row <- competition_df %>%
  mutate(
    x     = x_map[as.character(Species)],
    y_lab = factor("Sperm competition", levels = y_levels)
  )

fam2clusters <- cluster_long %>% distinct(Gene_family, UniqueID)

sl_tiles <- ampl_sl_long %>%
  inner_join(fam2clusters, by = "Gene_family") %>%
  transmute(
    x     = x_map[as.character(CellType)],
    y_lab = factor(UniqueID, levels = y_levels),
    fill_val = Present
  )

sp_tiles <- ampl_sp %>%
  inner_join(fam2clusters, by = "Gene_family") %>%
  transmute(
    x     = x_map[["Spermatogenesis"]],
    y_lab = factor(UniqueID, levels = y_levels),
    fill_val = Spermatogenesis
  )

dnds_tiles <- dnds %>%
  mutate(
    x     = x_map[["dN/dS"]],
    y_lab = factor(UniqueID, levels = y_levels),
    label = ifelse(is.na(omega), "NA", sprintf("%.2f", omega)),
    group = ifelse(is.na(omega), "NA", ifelse(omega > 1, "dN/dS>1", "dN/dS≤1"))
  ) %>% 
  filter(!is.na(y_lab))

site_tiles <- dnds %>%
  mutate(
    x       = x_map[["Site-test"]],
    y_lab   = factor(UniqueID, levels = y_levels),
    fill_val = site_test
  ) %>% 
  filter(!is.na(y_lab))

cluster_labels_df <- tibble::tibble(
  y_lab = factor(names(label_vec), levels = y_levels),
  label = unname(label_vec),
  x     = 0.5  # Position just inside the left edge
) %>%
  dplyr::filter(label != "")

# ------------------ Plot ------------------
species_cols_palette <- c(
  "GorGor"="#66C2A5","HomSap"="#FC8D62","PanPan"="#8DA0CB","PanTro"="#E78AC3",
  "MacFas"="#A6D854","PonAbe"="#FFD92F","PonPyg"="#E5C494","SymSyn"="#B3B3B3"
)
fill_palette <- c(
  "high"="red","intermediate"="gold","low"="skyblue",
  "Pre MSCI"="skyblue","Post MSCI"="orange","Pre + post MSCI"="purple",
  "Unknown"="grey80","Yes"="forestgreen","No"="white",
  "dN/dS>1"="red", "dN/dS≤1"="skyblue", "NA"="grey80"
)

p <- ggplot() +
  # bubbles
  geom_point(data = cluster_bubbles,
             aes(x = x, y = y_lab, size = Count, color = Species),
             alpha = 0.85) +
  geom_text(data = cluster_bubbles,
            aes(x = x, y = y_lab, label = Count),
            size = 4, color = "black") +
  scale_size(range = c(2, 15), guide = "none") +
  scale_color_manual(values = species_cols_palette, guide = "none") +
  
  # competition header
  geom_tile(data = comp_row,
            aes(x = x, y = y_lab, fill = Competition),
            width = 0.9, height = 0.9) +
  geom_text(data = comp_row,
            aes(x = x, y = y_lab, label = Competition),
            size =4, color = "black") +
  
  # family annotation tiles
  geom_tile(data = sl_tiles,
            aes(x = x, y = y_lab, fill = fill_val),
            width = 0.4, height = 0.7,
            color = "black", size = 0.2) +  # Add black border +
  geom_tile(data = sp_tiles,
            aes(x = x, y = y_lab, fill = fill_val),
            width = 0.4, height = 0.7) +
  
  # dN/dS + site-test
  geom_tile(data = dnds_tiles,
            aes(x = x, y = y_lab, fill = group),
            width = 0.4, height = 0.7) +
  geom_text(data = dnds_tiles,
            aes(x = x, y = y_lab, label = label),
            size = 4, color = "black") +
  geom_tile(data = site_tiles,
            aes(x = x, y = y_lab, fill = fill_val),
            width = 0.4, height = 0.7,
            color = "black", size = 0.2) +  # Add black border +
  
  #---- cluster names at left edge of panel ----
geom_text(data = cluster_labels_df,
          aes(x = 0.5, y = y_lab, label = label),
          hjust = 1, size = 3.7) +
  
  # ---- family brackets and names OUTSIDE panel (using annotation_custom) ----

scale_fill_manual(values = fill_palette,
                  breaks = c("Pre MSCI","Pre + post MSCI","Post MSCI","Unknown","Yes","No"),
                  name = "Expression") +
  
  scale_y_discrete(limits = y_levels, labels = NULL, drop = FALSE) +
  
  scale_x_continuous(
    breaks = seq_along(x_levels),
    labels = c(species_name_map[species_order], "Sertoli", "Leydig", "Gametogenesis", "dN/dS", "Site-test"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  # clip="off" allows drawing outside the panel
  coord_cartesian(xlim = c(0.75, length(column_order) + 0.5), clip = "off") +
  
  theme_minimal() + 
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin  = margin(10, 10, 10, 320, "pt"),
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "cm")
  ) +
  labs(x = NULL, y = NULL)

#  add family labels and brackets as annotations OUTSIDE the plot area
for (i in 1:nrow(family_brackets)) {
  # Vertical line
  p <- p + annotate("segment",
                    x = -1.08, xend = -1.08,
                    y = family_brackets$y_min[i] - 0.4,
                    yend = family_brackets$y_max[i] + 0.4,
                    size = 0.5)
  # Top horizontal
  p <- p + annotate("segment",
                    x = -1.08, xend = -1.08,
                    y = family_brackets$y_min[i] - 0.4,
                    yend = family_brackets$y_min[i] - 0.4,
                    size = 0.5)
  # Bottom horizontal
  p <- p + annotate("segment",
                    x = -1.08, xend = -1.08,
                    y = family_brackets$y_max[i] + 0.4,
                    yend = family_brackets$y_max[i] + 0.4,
                    size = 0.5)
  # Family name
  p <- p + annotate("text",
                    x = -1.12,
                    y = family_brackets$y_mid[i],
                    label = family_brackets$Gene_family[i],
                    hjust = 1, size = 4, fontface = "bold")
}

#print(p)


ggsave("~/Downloads/X_chrom_clusters_filtered_withdNdS_updated_try23.pdf", p, width = 20, height = 25, units = "in")



#### Y Chromosome ####
# ------------------ Packages ------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)

# ------------------ Inputs ------------------
# 1) Family-level summary (wide)
ampl_summary <- read_excel("~/Documents/PhD Aarhus/Amplicons/Amplicon_summary_count_expression_genefamily.xlsx") %>%
  filter(Chromosome == "Y") %>%
  mutate(Gene_family = `Gene family`)

# 2) Cluster-level counts
cluster_counts_perspecies_y <- read_delim("~/Downloads/cluster_counts_perspecies_y_updated.tsv", 
                                          delim = "\t") %>%
  mutate(Family = case_when(
    Family == "BPY2" ~ "BPY",
    Family == "CDY1" ~ "CDY",
    Family == "DAZ1" ~ "DAZ",
    Family == "HSFY1" ~ "HSFY",
    Family == "MTRNR2-like 17" ~ "MTRNR2-like 17",
    Family == "RBMY1B" ~ "RBMY",
    Family == "TATA-box binding protein associated factor 11 like protein 2-like" ~ "TAF11L2",
    Family == "TSPY8" ~ "TSPY",
    Family == "VCY1B" ~ "VCY",
    Family == "adenylate kinase isoenzyme 6-like" ~ "AKI6",#adenylate kinase isoenzyme 6-like" 
    Family == "centriole and centriolar satellite protein OFD1-like" ~ "CCS-OFD1", # centriole and centriolar satellite protein OFD1-like
    Family == "endogenous retrovirus group K member 19 Env polyprotein-like" ~ "retrovirus K19", #endogenous retrovirus group K member 19 Env polyprotein-like
    Family == "glutamate dehydrogenase 1, mitochondrial-like" ~ "GLUD1Y",
    Family == "keratin, type I cytoskeletal 18-like" ~ "KRT18Y",
    Family == "proline-rich protein, Y-linked" ~ "PRPY", #proline-rich protein, Y-linked
    Family == "protein FAM47A-like" ~ "FAM47AY",
    Family == "protein FRG1-like" ~ "FRG1Y",
    Family == "zinc finger protein 285-like" ~ "ZNF"))

# 3) --- dN/dS ---
library(readr)

dnds <- read_csv("~/Downloads/codeml_site_tests_median_omega_updated_Y_copy.csv") %>%
  dplyr::rename(Gene_family = Family) %>%
  mutate(Cluster = sub(".*_cluster_", "", Cluster)) %>%
  mutate(
    omega = coalesce(`Median M0 omega`),
    omega = ifelse(is.infinite(omega), NA, omega), # Turn Inf values into NA because no power. 
    site_test = factor(`site-test`, levels = c("No","Yes"))
  ) %>%
  select(Gene_family, Cluster, omega, site_test) %>%
  mutate(Gene_family = recode(Gene_family,
                              "BPY2"   = "BPY",
                              "CDY1"   = "CDY",
                              "DAZ1"   = "DAZ",
                              "FAM47A"   = "FAM47AY",
                              "HSFY1" ="HSFY",
                              "RBMY1B"   = "RBMY",
                              "TSPY8"   = "TSPY",
                              "glutamate"   = "GLUD1Y",
                              "isoenzyme"   = "AKI6",
                              "retrovirus" ="retrovirus K19",
                              "zinc"   = "ZNF"))

# Species columns + display order
species_cols  <- c("GorGor", "HomSap", "PanPan", "PanTro", "MacFas", "PonAbe", "PonPyg", "SymSyn")
species_order <- c("PanPan", "PanTro", "HomSap", "GorGor", "PonAbe", "PonPyg", "SymSyn", "MacFas")
column_order  <- c(species_order, "Sertoli", "Leydig", "Spermatogenesis", "dN/dS", "Site-test")

# ------------------ Family annotations ------------------
ampl_sl_long <- ampl_summary %>%
  select(Gene_family, Sertoli, Leydig) %>%
  pivot_longer(c(Sertoli, Leydig), names_to = "CellType", values_to = "Present") %>%
  mutate(CellType = factor(CellType, levels = column_order))

ampl_sp <- ampl_summary %>%
  select(Gene_family, Spermatogenesis) %>%
  mutate(Spermatogenesis = factor(Spermatogenesis,
                                  levels = c("Pre MSCI","Pre + post MSCI","Post MSCI","Unknown")))

gene_levels <- ampl_sp %>%
  distinct(Gene_family, Spermatogenesis) %>%
  arrange(Spermatogenesis, Gene_family) %>%
  pull(Gene_family)

ampl_sl_long$Gene_family <- factor(ampl_sl_long$Gene_family, levels = gene_levels)
ampl_sp$Gene_family      <- factor(ampl_sp$Gene_family,      levels = gene_levels)

competition_df <- tibble::tibble(
  Species = factor(species_order, levels = column_order),
  Competition = c("high","high","low","low","intermediate","intermediate","low","intermediate"),
  Gene_family = factor("Sperm competition", levels = c("Sperm competition", gene_levels))
)

# ------------------ Filter clusters ------------------
drops_all_singlecopy <- cluster_counts_perspecies_y %>%
  mutate(
    species_present     = rowSums(across(all_of(species_cols), ~ (.x > 0) & !is.na(.x))),
    nonzero_eq_one      = rowSums(across(all_of(species_cols), ~ (.x == 1) & !is.na(.x))),
    all_nonzero_are_one = (species_present > 0) & (nonzero_eq_one == species_present)
  ) %>%
  filter(all_nonzero_are_one) %>%
  transmute(Gene_family = Family, Cluster)

# ------------------ Long cluster table ------------------
cluster_long <- cluster_counts_perspecies_y %>%
  dplyr::rename(Gene_family = Family) %>%
  dplyr::anti_join(drops_all_singlecopy, by = c("Gene_family", "Cluster")) %>%
  pivot_longer(all_of(species_cols), names_to = "Species", values_to = "Count") %>%
  dplyr::filter(!is.na(Count), Count != 0) %>%
  dplyr::mutate(
    Species = factor(Species, levels = column_order),
    UniqueID = paste0(Gene_family, "|||", Cluster),
    ClusterDisplay = Cluster
  )

dnds <- dnds %>%
  mutate(UniqueID = paste0(Gene_family, "|||", Cluster),
         ClusterDisplay = Cluster)

# ------------------ Build y-axis ------------------
fam_order <- gene_levels[gene_levels %in% unique(cluster_long$Gene_family)]

make_block <- function(f) {
  labs <- cluster_long %>%
    filter(Gene_family == f) %>%
    distinct(UniqueID, ClusterDisplay) %>%
    arrange(ClusterDisplay) %>%
    pull(UniqueID)
  c(labs, paste0(f, "___SPACER___"))
}
y_levels <- c("Sperm competition", unlist(lapply(fam_order, make_block), use.names = FALSE))

cluster_display_map <- cluster_long %>%
  distinct(UniqueID, ClusterDisplay) %>%
  deframe()

family_brackets <- cluster_long %>%
  distinct(Gene_family, UniqueID, ClusterDisplay) %>%
  mutate(
    y_pos = match(UniqueID, y_levels),
    Gene_family = as.character(Gene_family)
  ) %>%
  group_by(Gene_family) %>%
  summarise(
    y_min = min(y_pos),
    y_max = max(y_pos),
    y_mid = (min(y_pos) + max(y_pos)) / 2,
    .groups = "drop"
  )

label_vec <- setNames(
  sapply(y_levels, function(x) {
    if (grepl("___SPACER___", x, fixed = TRUE)) {
      return("")
    } else if (x == "Sperm competition") {
      return("Sperm competition")
    } else if (x %in% names(cluster_display_map)) {
      return(cluster_display_map[x])
    } else {
      return(x)
    }
  }),
  y_levels
)

# numeric x position for each real column
x_levels <- column_order
x_map    <- setNames(seq_along(x_levels), x_levels)

cluster_bubbles <- cluster_long %>%
  transmute(
    x     = x_map[as.character(Species)],
    y_lab = factor(UniqueID, levels = y_levels),
    Count, Species
  )

comp_row <- competition_df %>%
  mutate(
    x     = x_map[as.character(Species)],
    y_lab = factor("Sperm competition", levels = y_levels)
  )

fam2clusters <- cluster_long %>% distinct(Gene_family, UniqueID)

sl_tiles <- ampl_sl_long %>%
  inner_join(fam2clusters, by = "Gene_family") %>%
  transmute(
    x     = x_map[as.character(CellType)],
    y_lab = factor(UniqueID, levels = y_levels),
    fill_val = Present
  )

sp_tiles <- ampl_sp %>%
  inner_join(fam2clusters, by = "Gene_family") %>%
  transmute(
    x     = x_map[["Spermatogenesis"]],
    y_lab = factor(UniqueID, levels = y_levels),
    fill_val = Spermatogenesis
  )

dnds_tiles <- dnds %>%
  mutate(
    x     = x_map[["dN/dS"]],
    y_lab = factor(UniqueID, levels = y_levels),
    label = ifelse(is.na(omega), "NA", sprintf("%.2f", omega)),
    group = ifelse(is.na(omega), "NA", ifelse(omega > 1, "dN/dS>1", "dN/dS≤1"))
  ) %>% 
  filter(!is.na(y_lab))

site_tiles <- dnds %>%
  mutate(
    x       = x_map[["Site-test"]],
    y_lab   = factor(UniqueID, levels = y_levels),
    fill_val = site_test
  ) %>% 
  filter(!is.na(y_lab))

cluster_labels_df <- tibble::tibble(
  y_lab = factor(names(label_vec), levels = y_levels),
  label = unname(label_vec),
  x     = 0.5
) %>%
  dplyr::filter(label != "")

# ------------------ Plot ------------------
species_cols_palette <- c(
  "GorGor"="#66C2A5","HomSap"="#FC8D62","PanPan"="#8DA0CB","PanTro"="#E78AC3",
  "MacFas"="#A6D854","PonAbe"="#FFD92F","PonPyg"="#E5C494","SymSyn"="#B3B3B3"
)
fill_palette <- c(
  "high"="red","intermediate"="gold","low"="skyblue",
  "Pre MSCI"="skyblue","Post MSCI"="orange","Pre + post MSCI"="purple",
  "Unknown"="grey80","Yes"="forestgreen","No"="white",
  "dN/dS>1"="red", "dN/dS≤1"="skyblue", "NA"="grey80"
)

p <- ggplot() +
  geom_point(data = cluster_bubbles,
             aes(x = x, y = y_lab, size = Count, color = Species),
             alpha = 0.85) +
  geom_text(data = cluster_bubbles,
            aes(x = x, y = y_lab, label = Count),
            size = 4, color = "black") +
  scale_size(range = c(2, 15), guide = "none") +
  scale_color_manual(values = species_cols_palette, guide = "none") +
  
  geom_tile(data = comp_row,
            aes(x = x, y = y_lab, fill = Competition),
            width = 0.9, height = 0.9) +
  geom_text(data = comp_row,
            aes(x = x, y = y_lab, label = Competition),
            size = 4, color = "black") +
  
  geom_tile(data = sl_tiles,
            aes(x = x, y = y_lab, fill = fill_val),
            width = 0.4, height = 0.7, 
            color = "black", size = 0.2) +  # Add black border +
  geom_tile(data = sp_tiles,
            aes(x = x, y = y_lab, fill = fill_val),
            width = 0.4, height = 0.7) +
  
  geom_tile(data = dnds_tiles,
            aes(x = x, y = y_lab, fill = group),
            width = 0.4, height = 0.7) +
  geom_text(data = dnds_tiles,
            aes(x = x, y = y_lab, label = label),
            size = 3, color = "black") +
  geom_tile(data = site_tiles,
            aes(x = x, y = y_lab, fill = fill_val),
            width = 0.4, height = 0.7, 
            color = "black", size = 0.2) +  # Add black border +
  
  geom_text(data = cluster_labels_df,
            aes(x = 0.5, y = y_lab, label = label),
            hjust = 1, size = 3.7) +
  
  scale_fill_manual(values = fill_palette,
                    breaks = c("Pre MSCI","Pre + post MSCI","Post MSCI","Unknown","Yes","No"), 
                    name = "Expression") +
  
  scale_y_discrete(limits = y_levels, labels = NULL, drop = FALSE) +
  
  scale_x_continuous(
    breaks = seq_along(x_levels),
    labels = c(species_name_map[species_order], "Sertoli", "Leydig", "Gametogenesis", "dN/dS", "Site-test"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  coord_cartesian(xlim = c(0.75, length(column_order) + 0.5), clip = "off") +
  
  theme_minimal() + 
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin  = margin(10, 10, 10, 320, "pt"),
    legend.position = "None",
    #legend.title = element_text(size = 14, face = "bold"),
    #legend.text = element_text(size = 12),
    #legend.key.size = unit(1.2, "cm")
  ) +
  labs(x = NULL, y = NULL)

# Add family labels and brackets
for (i in 1:nrow(family_brackets)) {
  # Vertical line
  p <- p + annotate("segment",
                    x = -1.08, xend = -1.08,
                    y = family_brackets$y_min[i] - 0.4,
                    yend = family_brackets$y_max[i] + 0.4,
                    size = 0.5)
  # Top horizontal
  p <- p + annotate("segment",
                    x = -1.08, xend = -1.08,
                    y = family_brackets$y_min[i] - 0.4,
                    yend = family_brackets$y_min[i] - 0.4,
                    size = 0.5)
  # Bottom horizontal
  p <- p + annotate("segment",
                    x = -1.08, xend = -1.08,
                    y = family_brackets$y_max[i] + 0.4,
                    yend = family_brackets$y_max[i] + 0.4,
                    size = 0.5)
  # Family name
  p <- p + annotate("text",
                    x = -1.12,
                    y = family_brackets$y_mid[i],
                    label = family_brackets$Gene_family[i],
                    hjust = 1, size = 4, fontface = "bold")
}

#print(p)

ggsave("~/Downloads/Y_chrom_clusters_filtered_withdNdS_updated_try21.jpg", plot = p, width = 20, height = 10, dpi = 1000)
