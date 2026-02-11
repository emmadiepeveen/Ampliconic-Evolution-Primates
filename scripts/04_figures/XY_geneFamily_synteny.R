###  SCRIPT for GENE family Synteny visualization ###

# --- packages ---
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(jsonlite)

## Files needed ##
# 1. All gene family data with clusters
# 2. Sequence classes information (Makova et al., 2024)
# 3. Palindrome information data (Makova et al., 2024)

#### Prepare input data ####
gene_details_X <- read_delim("~/Downloads/all_families_gene_details_with_clusters-3.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(Start = col_number(), 
                                              End = col_number()), trim_ws = TRUE)
# mutate that if gene lays in a palindrome, the class is called palindrome
gene_details_X <- gene_details_X %>%
  mutate(Class = ifelse(in_palindrome == "yes", "PALINDROME", Class))

# filter out rows with genes that are not classified as ampliconic genes (so have an NA in the gene_family_symbol)
gene_details_X<-gene_details_X %>% drop_na(gene_family_symbol)


gene_details_Y <- read_delim("~/Downloads/gene_details_updated_with_palindromes_coordinates_y.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(Start = col_number(), 
                                              End = col_number()), trim_ws = TRUE)
#mutate that if gene lays in a palindrome, the class is called palindrome
gene_details_Y <- gene_details_Y %>%
  mutate(Class = ifelse(in_palindrome == "yes", "PALINDROME", Class))

# filter out rows with genes that are not classified as ampliconic genes (so have an NA in the gene_family_symbol)
gene_details_Y<-gene_details_Y %>% drop_na(gene_family_symbol)

mutate the scientific names within the seqClass dataframe to the same as the species names
#Load the data
seqClasses <- read.table('~/Downloads/2023-11-20873B-s5/2023-11-20873B-s5/2023-11-20873B-AdditionalFile2-SeqClasses.txt', header = FALSE, sep = '\t')
colnames(seqClasses) <- c('species', 'scientific_name', 'chromosome', 'start', 'end', 'classification')

unique(seqClasses$scientific_name)

# change the names accordingly. ADD MORE THEN LOOKING AT MORE SPECIES
seqClasses <- seqClasses %>%
  mutate(scientific_name = case_when(
    scientific_name == "mPanTro3" ~ "PanTro",
    scientific_name == "CHM13" ~ "HomSap",
    scientific_name == "HG002" ~ "HomSap",
    scientific_name == "mPanPan1" ~ "PanPan", 
    scientific_name == "mGorGor1" ~ "GorGor", 
    scientific_name == "mPonAbe1" ~ "PonAbe", 
    scientific_name == "mPonPyg2" ~ "PonPyg", 
    scientific_name == "mSymSyn1" ~ "SymSyn",
    TRUE ~ scientific_name  # Keep all other values unchanged
  ))


## read in the palindrome file
palindrome <- read_delim("~/Downloads/PalindromeLocations.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         col_types = cols(start = col_number(), 
                                          end = col_number()), trim_ws = TRUE)

#add an extra column based on chromosome so that X and Y can be filtered 
palindrome <- palindrome %>%
  mutate(chromosome = case_when(
    endsWith(assembly.chr, ".chrX") ~ "chrX",
    endsWith(assembly.chr, ".chrY") ~ "chrY"))

unique(palindrome$assembly.chr)

# change the names accordingly. ADD MORE THEN LOOKING AT MORE SPECIES
palindrome <- palindrome %>%
  mutate(assembly.chr = case_when(
    assembly.chr == "mPanTro3.chrX" ~ "PanTro",
    assembly.chr == "mPanTro3.chrY" ~ "PanTro",
    assembly.chr == "CHM13.chrX" ~ "HomSap",
    assembly.chr == "HG002.chrY" ~ "HomSap",
    assembly.chr == "mGorGor1.chrX" ~ "GorGor", 
    assembly.chr == "mGorGor1.chrY" ~ "GorGor", 
    assembly.chr == "mPanPan1.chrX" ~ "PanPan", 
    assembly.chr == "mPanPan1.chrY" ~ "PanPan", 
    assembly.chr == "mPonAbe1.chrX" ~ "PonAbe", 
    assembly.chr == "mPonAbe1.chrY" ~ "PonAbe", 
    assembly.chr == "mPonPyg2.chrX" ~ "PonPyg", 
    assembly.chr == "mPonPyg2.chrY" ~ "PonPyg", 
    assembly.chr == "mSymSyn1.chrX" ~ "SymSyn",
    assembly.chr == "mSymSyn1.chrY" ~ "SymSyn",
    TRUE ~ assembly.chr  # Keep all other values unchanged
  ))
head(palindrome)

### Function to integrate palindrome data into seqClasses data -> so that all are within the same dataframe ###
update_seq_classes <- function(seqClasses, palindrome) {
  new_seqClasses <- data.frame(
    species = character(),
    scientific_name = character(),
    chromosome = character(),
    start = numeric(),
    end = numeric(),
    classification = character(),
    palindrome_name = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(seqClasses))) {
    row <- seqClasses[i, ]
    palindromes_in_region <- palindrome %>%
      filter(assembly.chr == row$scientific_name,
             chromosome == row$chromosome,
             start >= row$start,
             end <= row$end)
    
    current_start <- row$start
    
    if (nrow(palindromes_in_region) == 0) {
      new_seqClasses <- rbind(new_seqClasses, c(row, palindrome_name = NA))
    } else {
      for (j in seq_len(nrow(palindromes_in_region))) {
        pal <- palindromes_in_region[j, ]
        
        # Add the region before the palindrome
        if (current_start < pal$start) {
          new_seqClasses <- rbind(new_seqClasses, data.frame(
            species = row$species,
            scientific_name = row$scientific_name,
            chromosome = row$chromosome,
            start = current_start,
            end = pal$start - 1,
            classification = row$classification,
            palindrome_name = NA,
            stringsAsFactors = FALSE
          ))
        }
        
        # Add the palindrome region
        new_seqClasses <- rbind(new_seqClasses, data.frame(
          species = row$species,
          scientific_name = row$scientific_name,
          chromosome = row$chromosome,
          start = pal$start,
          end = pal$end,
          classification = "PALINDROME",
          palindrome_name = pal$name,
          stringsAsFactors = FALSE
        ))
        
        current_start <- pal$end + 1
      }
      
      # Add the remaining part after the last palindrome
      if (current_start <= row$end) {
        new_seqClasses <- rbind(new_seqClasses, data.frame(
          species = row$species,
          scientific_name = row$scientific_name,
          chromosome = row$chromosome,
          start = current_start,
          end = row$end,
          classification = row$classification,
          palindrome_name = NA,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(new_seqClasses)
}

# Run function
updated_seqClasses <- update_seq_classes(seqClasses, palindrome)

# Filter seqClasses for chrX
seqClasses_filtered_X <- seqClasses %>%
  filter(chromosome == "chrX")
# Filter updated seqClasses (with palindrome) for chrX
update_seq_classes_filtered_X <- updated_seqClasses %>%
  filter(chromosome == "chrX")

# Filter seqClasses for chrY
seqClasses_filtered_Y <- seqClasses %>%
  filter(chromosome == "chrY")
# Filter updated seqClasses (with palindrome) for chrY
update_seq_classes_filtered_Y <- updated_seqClasses %>%
  filter(chromosome == "chrY")


# Prepare df for arrows (as before), ensuring it is filtered for chrX
gene_details_X_arrows <- gene_details_X_sorted %>%
  filter(Species %in% unique(gene_details_X_sorted$Species)) %>%
  mutate(
    direction = ifelse(Strand == "+", 1, -1),  # Direction based on Strand
    arrow_start = ifelse(direction == -1, End, Start),
    arrow_end = ifelse(direction == -1, Start, End)
  )

# Prepare df for arrows (as before), ensuring it is filtered for chrY
gene_details_Y_arrows <- gene_details_Y_sorted %>%
  filter(Species %in% unique(gene_details_Y_sorted$Species)) %>%
  mutate(
    direction = ifelse(Strand == "+", 1, -1),  # Direction based on Strand
    arrow_start = ifelse(direction == -1, End, Start),
    arrow_end = ifelse(direction == -1, Start, End)
  )






#### SYNTENY PLOT X CHROMOSOME ####

## Update gene family names 
gene_details_X_arrows_updated <- gene_details_X_arrows %>%
  mutate(gene_family_symbol = case_when(
    gene_family_symbol == "ARMCX2" ~ "ARMCX",
    gene_family_symbol == "CENPVL2" ~ "CENPVL2",
    gene_family_symbol == "CSAG1" ~ "CSAG",
    gene_family_symbol == "CSF2RA" ~ "CSF2RA",
    gene_family_symbol == "CT45A8" ~ "CT45",
    gene_family_symbol == "CT47C1" ~ "CT47",
    gene_family_symbol == "CXorf49B" ~ "CXorf49",
    gene_family_symbol == "CXorf51B" ~ "CXorf51",
    gene_family_symbol == "DMRTC1B" ~ "DMRTC1",
    gene_family_symbol == "EOLA1" ~ "EOLA",
    gene_family_symbol == "ETDA" ~ "ETD",
    gene_family_symbol == "F8A2" ~ "F8",
    gene_family_symbol == "FAM156B" ~ "FAM156",
    gene_family_symbol == "FAM236C" ~ "FAM236",
    gene_family_symbol == "GAGE1" ~ "GAGE",
    gene_family_symbol == "GPRASP2" ~ "GPRASP",
    gene_family_symbol == "H2AB3" ~ "H2AB",
    gene_family_symbol == "H2BW2" ~ "H2BW",
    gene_family_symbol == "HSFX3" ~ "HSF31",
    gene_family_symbol == "IKBKG" ~ "IKBKG",
    gene_family_symbol == "LAGE3" ~ "CTAG",
    gene_family_symbol == "MAGEB18" ~ "MAGEB",
    gene_family_symbol == "MAGED1" ~ "MAGED",
    gene_family_symbol == "NUDT10" ~ "NUDT",
    gene_family_symbol == "NXF2" ~ "NXF",
    gene_family_symbol == "PABPC1L2B" ~ "PABPC",
    gene_family_symbol == "OPN1LW" ~ "OPN1LW",
    gene_family_symbol == "PNMA6E" ~ "PNMA",
    gene_family_symbol == "PWWP4" ~ "PWWP4",
    gene_family_symbol == "RAB40A" ~ "RAB40A",
    gene_family_symbol == "RHOXF2B" ~ "RHOXF2B",
    gene_family_symbol == "RPL36A-HNRNPH2" ~ "RPL36A-HNRNPH2",
    gene_family_symbol == "SAGE1" ~ "SAGE1",
    gene_family_symbol == "SMIM10L2A" ~ "SMIM10",
    gene_family_symbol == "SPACA5" ~ "SPACA",
    gene_family_symbol == "SPANXA1" ~ "SPANX",
    gene_family_symbol == "SPIN2A" ~ "SPIN",
    gene_family_symbol == "SSX4B" ~ "SSX",
    gene_family_symbol == "TBL1X" ~ "TBL1X",
    gene_family_symbol == "TCEAL8" ~ "TCEAL8",
    gene_family_symbol == "TCP11X2" ~ "TCP11X",
    gene_family_symbol == "TEX28" ~ "TEX28",
    gene_family_symbol == "TMEM185A" ~ "TMEM185",
    gene_family_symbol == "TMSB15B" ~ "TMSB15B",
    gene_family_symbol == "VCX2" ~ "VCX",
    gene_family_symbol == "XAGE1A" ~ "XAGE",
    gene_family_symbol == "ZXDA" ~ "ZXD",
    gene_family_symbol == "collagen alpha-4(IV) chain-like" ~ "collagen alpha-4(IV) chain-like",
    gene_family_symbol == "eendogenous retrovirus group K member 6 Env polyprotein-like" ~ "retrovirus K6",
    gene_family_symbol == "putative uncharacterized protein FLJ39060" ~ "FLJ39060",
    gene_family_symbol == "uncharacterized LOC129476473" ~ "LOC129476473",
    gene_family_symbol == "uncharacterized LOC129138873" ~ "LOC129138873",
    gene_family_symbol == "uncharacterized LOC129053094" ~ "LOC129053094",
    TRUE ~ gene_family_symbol
  ))

unique(gene_details_X_arrows_updated$gene_family_symbol)

# update the Species names to their common species names
unique(gene_details_X_arrows_updated$Species)

gene_details_X_arrows_updated <- gene_details_X_arrows_updated %>%
  mutate(Species = case_when(
    Species == "GorGor"~ "Gorilla", 
    Species == "HomSap"~ "Human",
    Species == "MacFas"~ "Macaque",
    Species == "PanPan"~ "Bonobo",
    Species == "PanTro"~ "Chimpanzee",
    Species == "PonAbe"~ "S.orangutan",
    Species == "PonPyg"~ "B.orangutan",
    Species == "SymSyn"~ "Siamang",
  ))

unique(gene_details_X_arrows_updated$Species)

# --- inputs you control ---
species_order <- c("Bonobo","Chimpanzee","Human","Gorilla","S.orangutan","B.orangutan","Siamang","Macaque")

#all the _unmoved_ families on the X:
moved_families<-c("ARMCX", "CENPVL2", "CSAG", "CSF2RA", "CT45", "CT47", "CXorf49", "CXorf51", "DMRTC1", "EOLA", "ETD", "F8", "FAM156", "FAM236", "GAGE", 
                  "GPRASP", "H2AB", "H2BW", "HSFX3", "IKBKG", "CTAG", "MAGEB", "MAGED", "NUDT", "NXF", "OPN1LW", "PABPC", "PNMA", "PWWP4", "RAB40A", 
                  "RHOXF2B", "SAGE1", "SMIM10", "SPACA", "SPIN", "SSX", "TCP11X", "TEX28", "TMEM185", 
                  "TMSB15B", "VCX", "XAGE", "ZXD", "FLJ39060", "LOC129476473", "HSFX1", 
                  "LOC129053094")

# collagen alpha-4(IV) chain-like is left out because it only occur in human
# LOC129138873 is left out because it occur in PanTro and MacFas and cannot be traced throughout the other species
# Gene families with stable multi-cluster structure - use nearest neighbor linking
# These families have persistent clusters that don't represent major rearrangements
nearest_neighbor_families <- c("H2AB", "MAGEB", "PABPC", "TMSB15B")


# all the MOVED families on the X
moved_families <- c("TBL1X","TCEAL8", "SPANX", "RPL36A-HNRNPH2", "retrovirus K6")
nearest_neighbor_families <- c("SPANX")

### FAMILIES ON THE Y! ### 
#moved_families<-c("VCY1B", "TSPY8", "RBMY1B", "HSFY1", "DAZ1", "CDY1", "BPY2", "glutamate dehydrogenase 1, mitochondrial-like", "protein FRG1-like, protein FAM47A-like")
#moved_families<-c("proline-rich protein, Y-linked", "MTRNR2-like 17", "adenylate kinase isoenzyme 6−like","endogenous retrovirus group K member 19 Env polyprotein−like" )

macfas_len_bp  <- 162126771

# centromeres
centro_tbl <- tribble(
  ~Species, ~cstart,     ~cend,
  "Bonobo",  59530707,   62384166,
  "Chimpanzee",  58108941,   60525984,
  "Human",  57347714,   61248490,
  "Gorilla",  68254530,   74519264,
  "S.orangutan",  58372707,   67812153,
  "B.orangutan",  58605191,   66367560,
  "Siamang",  70550338,   71318140,
  "Macaque",  58000000,   68000000
)

# --- helpers ---
to_num <- function(x) as.numeric(gsub("[^0-9.]", "", x))

# cluster nearby copies
cluster_positions <- function(pos_vec, tol_bp = 5e5) {
  if (length(pos_vec) == 0) return(tibble(cluster = integer(), y = numeric(), n = integer()))
  pos <- sort(pos_vec)
  grp <- cumsum(c(0, diff(pos) > tol_bp)) + 1
  tibble(pos = pos, grp = grp) |>
    group_by(grp) |>
    summarise(y = mean(pos), n = n(), .groups = "drop") |>
    mutate(cluster = row_number()) |>
    select(cluster, y, n)
}

# --- data prep ---
update_seq_classes_filtered_X_2 <- update_seq_classes_filtered_X %>%
  mutate(scientific_name = case_when(
    scientific_name == "GorGor"~ "Gorilla", 
    scientific_name == "HomSap"~ "Human",
    scientific_name == "MacFas"~ "Macaque",
    scientific_name == "PanPan"~ "Bonobo",
    scientific_name == "PanTro"~ "Chimpanzee",
    scientific_name == "PonAbe"~ "S.orangutan",
    scientific_name == "PonPyg"~ "B.orangutan",
    scientific_name == "SymSyn"~ "Siamang",
  ))

classes <- update_seq_classes_filtered_X_2 |>
  mutate(start = to_num(start), end = to_num(end)) |>
  filter(scientific_name %in% species_order, is.finite(start), is.finite(end)) |>
  mutate(species_idx = match(scientific_name, species_order))

chr_len <- classes |>
  dplyr::group_by(scientific_name) |>
  summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") |>
  dplyr::rename(Species = scientific_name)

if (!"Macaque" %in% chr_len$Species) {
  chr_len <- bind_rows(chr_len, tibble(Species = "Macaque", chr_len = macfas_len_bp))
}

chr_len <- chr_len |>
  filter(Species %in% species_order) |>
  mutate(species_idx = match(Species, species_order))

# --- gene copies ---
copies <- gene_details_X_arrows_updated |>
  filter(Species %in% species_order, gene_family_symbol %in% moved_families) |>
  mutate(arrow_start = to_num(arrow_start),
         arrow_end   = to_num(arrow_end),
         mid = (arrow_start + arrow_end)/2) |>
  inner_join(chr_len |> select(Species, chr_len, species_idx), by = "Species") |>
  filter(is.finite(mid))

# --- cluster within species ---
clusters <- copies |>
  group_by(gene_family_symbol, Species, species_idx) |>
  summarise(df = list(cluster_positions(mid, tol_bp = 5e5)), .groups = "drop") |>
  unnest(df)

# --- prepare data for D3 ---
seq_classes_df <- classes |>
  select(species = scientific_name, species_idx, start, end, classification) |>
  arrange(species_idx)

gene_clusters_df <- clusters |>
  mutate(species_idx = match(Species, species_order)) |>
  select(family = gene_family_symbol, species = Species, species_idx, position = y, copies = n, cluster)

# Create links - intelligently handle splits and merges
links_df <- map_dfr(moved_families, function(fam) {
  fam_clusters <- filter(clusters, gene_family_symbol == fam)
  links_list <- list()
  
  # Check if this family should use nearest neighbor linking
  use_nearest_neighbor <- fam %in% nearest_neighbor_families
  
  for (i in seq_len(length(species_order) - 1)) {
    sp_from <- species_order[i]
    sp_to <- species_order[i + 1]
    
    from_data <- filter(fam_clusters, Species == sp_from) %>% arrange(y)
    to_data <- filter(fam_clusters, Species == sp_to) %>% arrange(y)
    
    if (nrow(from_data) == 0 || nrow(to_data) == 0) next
    
    n_from <- nrow(from_data)
    n_to <- nrow(to_data)
    
    if (n_from == n_to) {
      # Same number - direct mapping
      for (j in seq_len(n_from)) {
        links_list[[length(links_list) + 1]] <- tibble(
          family = fam,
          source_species = sp_from,
          source_idx = i,
          source_pos = from_data$y[j],
          target_species = sp_to,
          target_idx = i + 1,
          target_pos = to_data$y[j],
          is_split = FALSE
        )
      }
    } else if (use_nearest_neighbor) {
      # Different numbers BUT use nearest neighbor matching
      # This prevents the "spaghetti" effect for multi-cluster families
      if (n_from < n_to) {
        # More clusters in target - each source connects to nearest target
        for (j in seq_len(n_from)) {
          # Find nearest target cluster
          distances <- abs(to_data$y - from_data$y[j])
          nearest_idx <- which.min(distances)
          
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[j],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[nearest_idx],
            is_split = FALSE
          )
        }
        # Connect remaining target clusters to their nearest source
        connected_targets <- sapply(seq_len(n_from), function(j) {
          which.min(abs(to_data$y - from_data$y[j]))
        })
        unconnected_targets <- setdiff(seq_len(n_to), connected_targets)
        for (k in unconnected_targets) {
          nearest_source_idx <- which.min(abs(from_data$y - to_data$y[k]))
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[nearest_source_idx],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[k],
            is_split = TRUE
          )
        }
      } else {
        # More clusters in source - each target connects to nearest source
        for (k in seq_len(n_to)) {
          # Find nearest source cluster
          distances <- abs(from_data$y - to_data$y[k])
          nearest_idx <- which.min(distances)
          
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[nearest_idx],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[k],
            is_split = FALSE
          )
        }
        # Connect remaining source clusters to their nearest target
        connected_sources <- sapply(seq_len(n_to), function(k) {
          which.min(abs(from_data$y - to_data$y[k]))
        })
        unconnected_sources <- setdiff(seq_len(n_from), connected_sources)
        for (j in unconnected_sources) {
          nearest_target_idx <- which.min(abs(to_data$y - from_data$y[j]))
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[j],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[nearest_target_idx],
            is_split = TRUE
          )
        }
      }
    } else {
      # Different numbers - all combinations (original behavior)
      for (j in seq_len(n_from)) {
        for (k in seq_len(n_to)) {
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[j],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[k],
            is_split = TRUE
          )
        }
      }
    }
  }
  
  bind_rows(links_list)
})

centromeres_df <- centro_tbl |> 
  mutate(species_idx = match(Species, species_order)) |>
  select(species = Species, species_idx, cstart, cend)

# Create chromosome lengths list for D3
chr_lengths_list <- setNames(
  as.list(chr_len$chr_len),
  chr_len$Species
)

d3_data <- list(
  species = species_order,
  max_chr_len = max(chr_len$chr_len),
  chr_lengths = chr_lengths_list,
  centromeres = centromeres_df,
  seq_classes = seq_classes_df,
  gene_clusters = gene_clusters_df,
  links = links_df,
  families = moved_families,
  colors = list(
    # Colorblind-safe palette inspired by synteny plots
    # Using Okabe-Ito base colors with variations in lightness
    
    # Orange family (highly visible)
    CSF2RA = "#E69F00",      # Bright orange
    XAGE = "#CC8800",        # Dark orange
    TCP11X = "#FF9933",      # Light orange
    EOLA = "#D4A574",        # Pale orange
    
    # Sky blue family (safe for all colorblindness types)
    VCX = "#56B4E9",         # Sky blue
    SPACA = "#4A9FD5",       # Medium blue
    NUDT = "#6AC5F5",        # Light sky blue
    TMEM185 = "#3D8CB8",    # Dark sky blue
    PNMA = "#4A9FD5",        # Medium blue
    TCEAL8 = "#3D8CB8",      # Dark sky blue
    
    # Bluish green family (distinguishable)
    GAGE = "#009E73",        # Bluish green
    ETD = "#00B589",         # Light bluish green
    FAM236 = "#008C61",      # Dark bluish green
    CTAG = "#00AA7A",        # Medium bluish green
    
    # Yellow family (high contrast)
    CENPVL2 = "#F0E442",     # Bright yellow
    PWWP4 = "#D4C534",       # Dark yellow
    OPN1LW = "#E8D84D",      # Light yellow
    
    # Blue family (strong, distinguishable)
    H2BW = "#0072B2",        # Strong blue
    ARMCX = "#005A8F",       # Dark blue
    CT45 = "#0088CC",        # Medium blue
    TMSB15B = "#005A8F",     # Dark blue
    PABPC = "#0072B2",       # Strong blue
    
    # Vermillion/red-orange family (colorblind-safe red)
    MAGED = "#D55E00",       # Vermillion
    SSX = "#B84D00",         # Dark vermillion
    SPIN = "#E66B00",        # Light vermillion
    F8 = "#C25500",          # Medium vermillion
    SPANX = "#D55E00",       # Vermillion
    
    # Purple/magenta family (distinguishable from blues)
    FAM156 = "#CC79A7",      # Reddish purple
    RAB40A = "#B3648F",      # Dark reddish purple
    MAGEB = "#D98AB8",       # Light reddish purple
    SMIM10 = "#BD6E9A",      # Medium reddish purple
    CT47 = "#9B7DB8",        # Soft purple
    TBL1X = "#9B7DB8",       # Soft purple
    
    # Brown/tan family (neutral, distinguishable)
    DMRTC1 = "#8B6914",      # Dark brown
    NXF = "#A67C00",         # Medium brown
    TEX28 = "#9A7200",       # Brown-orange
    CXorf51 = "#7A5C0F",     # Very dark brown
    LOC129053094 = "#8B6914", # Dark brown
    
    # Gray family (for less important/ambiguous genes)
    ZXD = "#6C6C6C",        # Medium gray
    CXorf49 = "#5A5A5A",     # Dark gray
    HSFX1 = "#787878",       # Light gray
    HSFX3 = "#545454",       # Very dark gray
    CSAG = "#696969",        # Medium-dark gray
    HSFX1 = "#808080",       #  Medium gray
    "LOC129476473" = "black", # Dark gray
    
    # Teal/cyan family (distinct from blues)
    SAGE1 = "#00CED1",       # Dark turquoise
    IKBKG = "#00B5B8",       # Medium turquoise
    RHOXF2B = "#00A8AB",     # Teal
    
    # Olive/green family (darker greens, safe)
    GPRASP = "#6B8E23",      # Olive drab
    "RPL36A-HNRNPH2" = "#7A9A2E",      # Yellow-green
    "FLJ39060" = "#5F7A1A",  # Dark olive
    
    # Rose/salmon family (lighter, distinguishable)
    H2AB = "#CD5C5C",        # Indian red
    
    # Additional assignments to cover all families
    INTS6L = "#4682B4",      # Steel blue (if needed)
    "retrovirus K6" = "#DAA520"  # Goldenrod
  ),
  class_colors = list(
    ANCESTRAL = "#D3D3D3",#FFD700
    PAR = "#D3D3D3",#90EE90
    AMPLICONIC = "#D3D3D3",#87CEEB
    SATELLITE = "#D3D3D3",#A0522D
    XTR = "#D3D3D3",#FFB6C1
    UNCLASSIFIED = "#D3D3D3",#D3D3D3
    PALINDROME = "#D3D3D3"#4169E1
  )
)

# Convert to JSON
json_data <- toJSON(d3_data, auto_unbox = TRUE, pretty = TRUE)

# Verify JSON was created
if (is.null(json_data) || nchar(json_data) == 0) {
  stop("Failed to create JSON data")
}

cat("JSON data created successfully, length:", nchar(json_data), "\n")

# Create standalone HTML file with updated label code
html_template <- '
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Gene Family Movement</title>
    <script src="https://d3js.org/d3.v6.min.js"></script>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 20px;
            background: white;
        }
        #chart { 
            margin: 0 auto; 
            display: block;
        }
    </style>
</head>
<body>
    <div id="chart"></div>
    <script>
    const data = DATA_PLACEHOLDER;
    
    const CHR_WIDTH = 30;
    const CENTRO_WIDTH = 17;
    
    const myWidth = 1200;
    const myHeight = 700;
    const padding = {left: 400, right: 80, top: 100, bottom: 50};'

html_content <- paste0(
  html_template,
  '
    
    const yScale = d3.scaleLinear()
      .domain([0, data.max_chr_len])
      .range([padding.top, myHeight - padding.bottom]);
    
    const xScale = d3.scalePoint()
      .domain(data.species)
      .range([padding.left + CHR_WIDTH, myWidth - padding.right - CHR_WIDTH])
      .padding(0.5);
    
    const chrInfo = {};
    
    // First pass: create basic info with actual chromosome lengths
    data.species.forEach(sp => {
      const chrData = data.centromeres.find(c => c.species === sp);
      // Find the actual chromosome length for this species
      const chrLength = data.chr_lengths ? data.chr_lengths[sp] : data.max_chr_len;
      
      chrInfo[sp] = {
        x: xScale(sp),
        chrLength: chrLength,
        yStart: yScale(0),
        yEnd: yScale(chrLength),
        centroStart: yScale(chrData.cstart),
        centroEnd: yScale(chrData.cend),
        chrWidth: CHR_WIDTH
      };
    });
    
    const svg = d3.select("#chart")
      .append("svg")
      .attr("width", myWidth)
      .attr("height", myHeight);
    
    const chrGroup = svg.append("g").attr("class", "chromosomes");
    
    data.species.forEach(sp => {
      const info = chrInfo[sp];
      const yStart = info.yStart;
      const yEnd = info.yEnd;
      const centroStart = info.centroStart;
      const centroEnd = info.centroEnd;
      const x = info.x;
      
      const topCap = CHR_WIDTH * 0.3;
      const bottomCap = CHR_WIDTH * 0.3;
      
      // Fixed centromere visual size in pixels for consistency
      const centroVisualSize = 30; // pixels
      const centroCenterY = (centroStart + centroEnd) / 2;
      const centroVisualStart = centroCenterY - centroVisualSize / 2;
      const centroVisualEnd = centroCenterY + centroVisualSize / 2;
      
      const points = [];
      const nPoints = 100;
      
      for (let i = 0; i <= nPoints; i++) {
        const t = i / nPoints;
        const y = yStart + t * (yEnd - yStart);
        
        let width;
        if (y < yStart + topCap) {
          const capT = (y - yStart) / topCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y > yEnd - bottomCap) {
          const capT = (yEnd - y) / bottomCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y >= centroVisualStart && y <= centroVisualEnd) {
          // Use fixed visual size for centromere
          const centroT = (y - centroVisualStart) / (centroVisualEnd - centroVisualStart);
          const pinchFactor = Math.sin(centroT * Math.PI);
          width = CENTRO_WIDTH/2 + (CHR_WIDTH/2 - CENTRO_WIDTH/2) * pinchFactor;
        } else {
          const chrCenter = (yStart + yEnd) / 2;
          const distFromCenter = Math.abs(y - chrCenter) / ((yEnd - yStart) / 2);
          width = CHR_WIDTH/2 * (0.95 + 0.05 * (1 - distFromCenter * 0.3));
        }
        
        points.push({y: y, width: width});
      }
      
      let pathData = "M" + (x + points[0].width) + "," + points[0].y;
      for (let i = 1; i < points.length; i++) {
        pathData += " L" + (x + points[i].width) + "," + points[i].y;
      }
      for (let i = points.length - 1; i >= 0; i--) {
        pathData += " L" + (x - points[i].width) + "," + points[i].y;
      }
      pathData += " Z";
      
      chrGroup.append("path")
        .attr("d", pathData)
        .attr("fill", "#EAEAEA")
        .attr("stroke", "#000")
        .attr("stroke-width", 2);
    });
    
    const defs = svg.append("defs");
    
    data.species.forEach(sp => {
      const info = chrInfo[sp];
      const clipPath = defs.append("clipPath")
        .attr("id", "clip-" + sp);
      
      const yStart = info.yStart;
      const yEnd = info.yEnd;
      const centroStart = info.centroStart;
      const centroEnd = info.centroEnd;
      const x = info.x;
      
      const topCap = CHR_WIDTH * 0.3;
      const bottomCap = CHR_WIDTH * 0.3;
      
      // Fixed centromere visual size in pixels for consistency
      const centroVisualSize = 30; // pixels
      const centroCenterY = (centroStart + centroEnd) / 2;
      const centroVisualStart = centroCenterY - centroVisualSize / 2;
      const centroVisualEnd = centroCenterY + centroVisualSize / 2;
      
      const points = [];
      const nPoints = 100;
      
      for (let i = 0; i <= nPoints; i++) {
        const t = i / nPoints;
        const y = yStart + t * (yEnd - yStart);
        
        let width;
        if (y < yStart + topCap) {
          const capT = (y - yStart) / topCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y > yEnd - bottomCap) {
          const capT = (yEnd - y) / bottomCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y >= centroVisualStart && y <= centroVisualEnd) {
          // Use fixed visual size for centromere
          const centroT = (y - centroVisualStart) / (centroVisualEnd - centroVisualStart);
          const pinchFactor = Math.sin(centroT * Math.PI);
          width = CENTRO_WIDTH/2 + (CHR_WIDTH/2 - CENTRO_WIDTH/2) * pinchFactor;
        } else {
          const chrCenter = (yStart + yEnd) / 2;
          const distFromCenter = Math.abs(y - chrCenter) / ((yEnd - yStart) / 2);
          width = CHR_WIDTH/2 * (0.95 + 0.05 * (1 - distFromCenter * 0.3));
        }
        
        points.push({y: y, width: width});
      }
      
      let pathData = "M" + (x + points[0].width) + "," + points[0].y;
      for (let i = 1; i < points.length; i++) {
        pathData += " L" + (x + points[i].width) + "," + points[i].y;
      }
      for (let i = points.length - 1; i >= 0; i--) {
        pathData += " L" + (x - points[i].width) + "," + points[i].y;
      }
      pathData += " Z";
      
      clipPath.append("path")
        .attr("d", pathData);
    });
    
    const classGroup = svg.append("g").attr("class", "seq-classes");
    
    data.seq_classes.forEach(cls => {
      const info = chrInfo[cls.species];
      const x = info.x;
      
      classGroup.append("rect")
        .attr("x", x - info.chrWidth/2 + 1)
        .attr("y", yScale(cls.start))
        .attr("width", info.chrWidth - 2)
        .attr("height", yScale(cls.end) - yScale(cls.start))
        .attr("fill", data.class_colors[cls.classification] || "#CCC")
        .attr("opacity", 0.5)
        .attr("clip-path", "url(#clip-" + cls.species + ")");
    });
    
    const linkGen = d3.linkHorizontal()
      .x(d => d[0])
      .y(d => d[1]);
    
    const linksGroup = svg.append("g").attr("class", "links");
    
    data.links.forEach(link => {
      const color = data.colors[link.family];
      const sourceX = xScale(link.source_species) + CHR_WIDTH/2;
      const targetX = xScale(link.target_species) - CHR_WIDTH/2;
      const sourceY = yScale(link.source_pos);
      const targetY = yScale(link.target_pos);
      const linkWidth = link.is_split ? 1 : 2;
    
    // shading !  
    //  linksGroup.append("path")
    //    .attr("d", linkGen({
    //      source: [sourceX, sourceY],
    //      target: [targetX, targetY]
    //    }))
    //    .attr("fill", "none")
    //    .attr("stroke", color)
    //    .attr("stroke-width", linkWidth * 3)
    //    .attr("opacity", 0.12);
      
      linksGroup.append("path")
        .attr("d", linkGen({
          source: [sourceX, sourceY],
          target: [targetX, targetY]
        }))
        .attr("fill", "none")
        .attr("stroke", color)
        .attr("stroke-width", linkWidth)
        .attr("opacity", link.is_split ? 0.4 : 0.7);
    });
    
    const markersGroup = svg.append("g").attr("class", "markers");
    
    data.gene_clusters.forEach(gc => {
      const info = chrInfo[gc.species];
      const x = info.x;
      const color = data.colors[gc.family];
      const y = yScale(gc.position);
      const blockWidth = info.chrWidth - 4;
      
      markersGroup.append("rect")
        .attr("x", x - blockWidth/2)
        .attr("y", y - 3)
        .attr("width", blockWidth)
        .attr("height", 6)
        .attr("fill", color)
        .attr("opacity", 0.8)
        .attr("clip-path", "url(#clip-" + gc.species + ")");
    });
    
    // Get unique families and their appearance info
    const familyLabels = [];

    // For each family, find where to label it (first or last species)
    data.families.forEach(family => {
      const familyClusters = data.gene_clusters.filter(gc => gc.family === family);
      if (familyClusters.length === 0) return;
      
      // Sort by species index
      familyClusters.sort((a, b) => a.species_idx - b.species_idx);
      const firstAppearance = familyClusters[0];
      const lastAppearance = familyClusters[familyClusters.length - 1];
      
      // Decide: label on left (first species) or right (last species)
      const firstSpeciesIdx = firstAppearance.species_idx;
      const labelOnRight = firstSpeciesIdx > 1; // If it does not appear in first 2 species, label on right
      
      const labelSpecies = labelOnRight ? lastAppearance.species : firstAppearance.species;
      const clustersInLabelSpecies = familyClusters.filter(
        gc => gc.species === labelSpecies
      );
      
      // Calculate average position
      const avgPosition = clustersInLabelSpecies.reduce((sum, gc) => sum + gc.position, 0) / 
                          clustersInLabelSpecies.length;
      
      familyLabels.push({
        family: family,
        species: labelSpecies,
        originalY: avgPosition,
        y: avgPosition,
        labelOnRight: labelOnRight
      });
    });

    // Sort by species index first, then by position within species
    familyLabels.sort((a, b) => {
      const speciesCompare = data.species.indexOf(a.species) - data.species.indexOf(b.species);
      if (speciesCompare !== 0) return speciesCompare;
      return a.y - b.y;
    });

    // Group labels by species AND position AND side
    const groupSpacingThreshold = 8;
    const MAX_GENES_PER_GROUP = 4;

    // Separate left and right labels
    const leftLabels = familyLabels.filter(l => !l.labelOnRight);
    const rightLabels = familyLabels.filter(l => l.labelOnRight);

    // Function to create groups
    function createGroups(labels) {
      const groups = [];
      if (labels.length === 0) return groups;
      
      let currentGroup = [labels[0]];
      
      for (let i = 1; i < labels.length; i++) {
        const prevLabel = labels[i - 1];
        const currLabel = labels[i];
        
        if (prevLabel.species === currLabel.species) {
          const prevYScreen = yScale(prevLabel.y);
          const currYScreen = yScale(currLabel.y);
          
          if (currYScreen - prevYScreen < groupSpacingThreshold) {
            if (currentGroup.length >= MAX_GENES_PER_GROUP) {
              groups.push(currentGroup);
              currentGroup = [currLabel];
            } else {
              currentGroup.push(currLabel);
            }
          } else {
            groups.push(currentGroup);
            currentGroup = [currLabel];
          }
        } else {
          groups.push(currentGroup);
          currentGroup = [currLabel];
        }
      }
      groups.push(currentGroup);
      return groups;
    }

    const leftGroups = createGroups(leftLabels);
    const rightGroups = createGroups(rightLabels);

    // Adjust positions to prevent overlap - LEFT SIDE
    const leftSpeciesGroups = {};
    leftGroups.forEach(group => {
      const species = group[0].species;
      if (!leftSpeciesGroups[species]) {
        leftSpeciesGroups[species] = [];
      }
      leftSpeciesGroups[species].push(group);
    });

    Object.keys(leftSpeciesGroups).forEach(species => {
      const groups = leftSpeciesGroups[species];
      for (let pass = 0; pass < 5; pass++) {
        for (let i = 1; i < groups.length; i++) {
          const prevGroup = groups[i - 1];
          const currGroup = groups[i];
          const prevGroupY = yScale(prevGroup[0].y);
          const currGroupY = yScale(currGroup[0].y);
          const requiredSpacing = 16;
          
          if (currGroupY - prevGroupY < requiredSpacing) {
            const adjustment = yScale.invert(prevGroupY + requiredSpacing) - currGroup[0].y;
            currGroup.forEach(label => {
              label.y += adjustment;
            });
          }
        }
      }
    });

    // Adjust positions to prevent overlap - RIGHT SIDE
    const rightSpeciesGroups = {};
    rightGroups.forEach(group => {
      const species = group[0].species;
      if (!rightSpeciesGroups[species]) {
        rightSpeciesGroups[species] = [];
      }
      rightSpeciesGroups[species].push(group);
    });

    Object.keys(rightSpeciesGroups).forEach(species => {
      const groups = rightSpeciesGroups[species];
      for (let pass = 0; pass < 5; pass++) {
        for (let i = 1; i < groups.length; i++) {
          const prevGroup = groups[i - 1];
          const currGroup = groups[i];
          const prevGroupY = yScale(prevGroup[0].y);
          const currGroupY = yScale(currGroup[0].y);
          const requiredSpacing = 16;
          
          if (currGroupY - prevGroupY < requiredSpacing) {
            const adjustment = yScale.invert(prevGroupY + requiredSpacing) - currGroup[0].y;
            currGroup.forEach(label => {
              label.y += adjustment;
            });
          }
        }
      }
    });

    // Draw LEFT labels
    leftGroups.forEach(group => {
      const labelSpecies = group[0].species;
      const info = chrInfo[labelSpecies];
      const x = info.x;
      
      if (group.length === 1) {
        const fam = group[0];
        const color = data.colors[fam.family] || "#999999";
        const originalYScreen = yScale(fam.originalY);
        const adjustedYScreen = yScale(fam.y);
        
        if (Math.abs(adjustedYScreen - originalYScreen) > 2) {
          markersGroup.append("line")
            .attr("x1", x - info.chrWidth/2 - 5)
            .attr("y1", originalYScreen)
            .attr("x2", x - info.chrWidth/2 - 8)
            .attr("y2", adjustedYScreen)
            .attr("stroke", color)
            .attr("stroke-width", 0.5)
            .attr("opacity", 0.5);
        }
        
        markersGroup.append("text")
          .attr("x", x - info.chrWidth/2 - 8)
          .attr("y", adjustedYScreen + 4)
          .attr("text-anchor", "end")
          .attr("font-size", "11px")
          .attr("fill", color)
          .attr("font-weight", "bold")
          .text(fam.family);
        
      } else {
        const groupY = yScale(group[0].y);
        
        group.forEach(fam => {
          const color = data.colors[fam.family] || "#999999";
          const originalYScreen = yScale(fam.originalY);
          
          if (Math.abs(groupY - originalYScreen) > 2) {
            markersGroup.append("line")
              .attr("x1", x - info.chrWidth/2 - 5)
              .attr("y1", originalYScreen)
              .attr("x2", x - info.chrWidth/2 - 8)
              .attr("y2", groupY)
              .attr("stroke", color)
              .attr("stroke-width", 0.5)
              .attr("opacity", 0.5);
          }
        });
        
        let currentX = x - info.chrWidth/2 - 8;
        
        group.forEach((fam, idx) => {
          const color = data.colors[fam.family] || "#999999";
          
          if (idx > 0) {
            const arrowText = markersGroup.append("text")
              .attr("text-anchor", "end")
              .attr("font-size", "11px")
              .attr("font-weight", "bold")
              .attr("fill", "#999")
              .text(" ~ ");
            
            const arrowBBox = arrowText.node().getBBox();
            arrowText.attr("x", currentX).attr("y", groupY + 4);
            currentX -= arrowBBox.width;
          }
          
          const geneText = markersGroup.append("text")
            .attr("text-anchor", "end")
            .attr("font-size", "11px")
            .attr("font-weight", "bold")
            .attr("fill", color)
            .text(fam.family + " ");
          
          const geneBBox = geneText.node().getBBox();
          geneText.attr("x", currentX).attr("y", groupY + 4);
          currentX -= geneBBox.width;
        });
      }
    });

    // Draw RIGHT labels (with different styling to make them stand out)
    rightGroups.forEach(group => {
      const labelSpecies = group[0].species;
      const info = chrInfo[labelSpecies];
      const x = info.x;
      
      if (group.length === 1) {
        const fam = group[0];
        const color = data.colors[fam.family] || "#999999";
        const originalYScreen = yScale(fam.originalY);
        const adjustedYScreen = yScale(fam.y);
        
        if (Math.abs(adjustedYScreen - originalYScreen) > 2) {
          markersGroup.append("line")
            .attr("x1", x + info.chrWidth/2 + 5)
            .attr("y1", originalYScreen)
            .attr("x2", x + info.chrWidth/2 + 8)
            .attr("y2", adjustedYScreen)
            .attr("stroke", color)
            .attr("stroke-width", 0.5)
            .attr("opacity", 0.5);
        }
        
        markersGroup.append("text")
          .attr("x", x + info.chrWidth/2 + 8)
          .attr("y", adjustedYScreen + 4)
          .attr("text-anchor", "start")
          .attr("font-size", "11px")
          .attr("fill", color)
          .attr("font-weight", "bold")
          .attr("font-style", "italic")
          .text(fam.family);
        
      } else {
        const groupY = yScale(group[0].y);
        
        group.forEach(fam => {
          const color = data.colors[fam.family] || "#999999";
          const originalYScreen = yScale(fam.originalY);
          
          if (Math.abs(groupY - originalYScreen) > 2) {
            markersGroup.append("line")
              .attr("x1", x + info.chrWidth/2 + 5)
              .attr("y1", originalYScreen)
              .attr("x2", x + info.chrWidth/2 + 8)
              .attr("y2", groupY)
              .attr("stroke", color)
              .attr("stroke-width", 0.5)
              .attr("opacity", 0.5);
          }
        });
        
        let currentX = x + info.chrWidth/2 + 8;
        
        group.forEach((fam, idx) => {
          const color = data.colors[fam.family] || "#999999";
          
          if (idx > 0) {
            const arrowText = markersGroup.append("text")
              .attr("text-anchor", "start")
              .attr("font-size", "11px")
              .attr("font-weight", "bold")
              .attr("fill", "#999")
              .text(" ~ ");
            
            const arrowBBox = arrowText.node().getBBox();
            arrowText.attr("x", currentX).attr("y", groupY + 4);
            currentX += arrowBBox.width;
          }
          
          const geneText = markersGroup.append("text")
            .attr("text-anchor", "start")
            .attr("font-size", "11px")
            .attr("font-weight", "bold")
            .attr("font-style", "italic")
            .attr("fill", color)
            .text(fam.family + " ");
          
          const geneBBox = geneText.node().getBBox();
          geneText.attr("x", currentX).attr("y", groupY + 4);
          currentX += geneBBox.width;
        });
      }
    });

    // SPECIES LABELS (on top)
    const labelGroup = svg.append("g").attr("class", "labels");

    data.species.forEach(sp => {
      labelGroup.append("text")
        .attr("x", xScale(sp))
        .attr("y", padding.top - 15)
        .attr("text-anchor", "middle")
        .attr("font-size", "13.4px")
        .attr("font-weight", "bold")
        .text(sp);
    });

    // Y-AXIS (chromosome length on right)
    const yAxis = d3.axisRight(yScale)
      .tickFormat(d => Math.round(d/1e6) + " Mb")
      .ticks(10);

    svg.append("g")
      .attr("transform", "translate(" + (myWidth - padding.right + 20) + ", 0)")
      .call(yAxis)
      .selectAll("text")
      .style("font-size", "14px")
      .selectAll("line")
      .attr("x2", 15);

    // TITLE
    svg.append("text")
      .attr("x", myWidth/2)
      .attr("y", 25)
      .attr("text-anchor", "middle")
      .attr("font-size", "18px")
      .attr("font-weight", "bold")
      .text("Gene Family Movement Across Primate Species");
      
      // Add a download SVG button
      const downloadButton = document.createElement("button");
      downloadButton.textContent = "Download SVG";
      downloadButton.style.cssText = "position: fixed; top: 10px; right: 10px; padding: 10px 20px; font-size: 14px; cursor: pointer; z-index: 1000;";
      document.body.appendChild(downloadButton);
      
      downloadButton.addEventListener("click", function() {
        const svgElement = document.querySelector("svg");
        const svgData = new XMLSerializer().serializeToString(svgElement);
        const svgBlob = new Blob([svgData], {type: "image/svg+xml;charset=utf-8"});
        const svgUrl = URL.createObjectURL(svgBlob);
        const downloadLink = document.createElement("a");
        downloadLink.href = svgUrl;
        downloadLink.download = "gene_movement.svg";
        document.body.appendChild(downloadLink);
        downloadLink.click();
        document.body.removeChild(downloadLink);
        URL.revokeObjectURL(svgUrl);
        });

    console.log("Visualization complete!");
    console.log("Left labels:", leftGroups.length, "Right labels:", rightGroups.length);
  </script>
    </body>
    </html>
    ')

# Replace placeholder with actual data
html_content <- sub("DATA_PLACEHOLDER", json_data, html_content, fixed = TRUE)

# Write HTML file
output_file <- "gene_movement_changingX_updated_background_commonnames_NEW3.html"
writeLines(html_content, output_file)
cat("HTML file created:", output_file, "\n")
cat("Open it in your browser to view the visualization.\n")






#### SYNTENY PLOT Y CHROMOSOME ####
# --- packages ---
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(jsonlite)

# Change the input file
## Update gene family names 
gene_details_Y_arrows_updated <- gene_details_Y_arrows %>%
  mutate(gene_family_symbol = case_when(
    gene_family_symbol == "BPY2" ~ "BPY2",
    gene_family_symbol == "CDY1" ~ "CDY",
    gene_family_symbol == "DAZ1" ~ "DAZ",
    gene_family_symbol == "HSFY1" ~ "HSFY",
    gene_family_symbol == "RBMY1B" ~ "RBMY",
    gene_family_symbol == "TSPY8" ~ "TSPY8",
    gene_family_symbol == "centriole and centriolar satellite protein OFD1-like" ~ "CCS-OFD1",
    gene_family_symbol == "proline-rich protein, Y-linked" ~ "PRPY",
    gene_family_symbol == "protein FRG1-like" ~ "FRG1Y",
    gene_family_symbol == "keratin, type I cytoskeletal 18-like" ~ "KRT18Y",
    gene_family_symbol == "VCY1B" ~ "VCY",
    gene_family_symbol == "MTRNR2-like 17" ~ "MTRNR2-like 17",
    gene_family_symbol == "glutamate dehydrogenase 1, mitochondrial-like" ~ "GLUD1Y",
    gene_family_symbol == "adenylate kinase isoenzyme 6-like" ~ "AKI6",
    gene_family_symbol == "endogenous retrovirus group K member 19 Env polyprotein-like" ~ "retrovirus K19",
    gene_family_symbol == "protein FAM47A-like" ~ "FAM47AY",
    gene_family_symbol == "zinc finger protein 285-like" ~ "ZNF",
    gene_family_symbol == "TATA-box binding protein associated factor 11 like protein 2-like" ~ "TAF11L2",
    TRUE ~ gene_family_symbol
  ))


# Add all your gene mappings here
unique(gene_details_Y_arrows_updated$gene_family_symbol)

# update the Species names to their common species names
unique(gene_details_Y_arrows_updated$Species)

gene_details_Y_arrows_updated <- gene_details_Y_arrows_updated %>%
  mutate(Species = case_when(
    Species == "GorGor"~ "Gorilla", 
    Species == "HomSap"~ "Human",
    Species == "MacFas"~ "Macaque",
    Species == "PanPan"~ "Bonobo",
    Species == "PanTro"~ "Chimpanzee",
    Species == "PonAbe"~ "S.orangutan",
    Species == "PonPyg"~ "B.orangutan",
    Species == "SymSyn"~ "Siamang",
  ))

unique(gene_details_Y_arrows_updated$Species)

# --- inputs you control ---
species_order <- c("Bonobo","Chimpanzee","Human","Gorilla","S.orangutan","B.orangutan","Siamang","Macaque")

### FAMILIES ON THE Y! ### 
# moved families 
#moved_families<-c("VCY1B", "TSPY8", "RBMY1B", "HSFY1", "DAZ1", "CDY1", "BPY2", "glutamate dehydrogenase 1, mitochondrial-like", "protein FRG1-like, protein FAM47A-like")
moved_families<-c("BPY2","CDY","DAZ","HSFY","RBMY","TSPY8","PRPY","FRG1Y","VCY","KRT18Y" ,"GLUD1Y","AKI6" , "retrovirus K19" ,"FAM47AY",
                  "ZNF")

# MTRNR2-like is left out because it only occur in PanTRO
# TAF11L2 is left out because it only occurs in SymSyn
# CSS-OFD1 is left out because it only occurs in GorGor
#stationary families 
#moved_families<-c("proline-rich protein, Y-linked", "MTRNR2-like 17", "adenylate kinase isoenzyme 6−like","endogenous retrovirus group K member 19 Env polyprotein−like" )
nearest_neighbor_families<-c("VCY", "TSPY", "RBMY", "HSFY", "DAZ", "CDY", "BPY2")


macfas_len_bp  <- 17124583

# centromeres
centro_tbl <- tribble(
  ~Species, ~cstart,     ~cend,
  "Bonobo",  31243512,   34904854,
  "Chimpanzee",  22014363,   23334108,
  "Human",  20959054,   21232750,
  "Gorilla",  20754406,   25389901,
  "S.orangutan",  18152839,   35332926,
  "B.orangutan",  18427511,   25366665,
  "Siamang",  13028709,   17328173, 
  "Macaque",  11987208,   17124583
)

# --- helpers ---
to_num <- function(x) as.numeric(gsub("[^0-9.]", "", x))

# cluster nearby copies
cluster_positions <- function(pos_vec, tol_bp = 5e5) {
  if (length(pos_vec) == 0) return(tibble(cluster = integer(), y = numeric(), n = integer()))
  pos <- sort(pos_vec)
  grp <- cumsum(c(0, diff(pos) > tol_bp)) + 1
  tibble(pos = pos, grp = grp) |>
    group_by(grp) |>
    summarise(y = mean(pos), n = n(), .groups = "drop") |>
    mutate(cluster = row_number()) |>
    select(cluster, y, n)
}

# --- data prep ---
update_seq_classes_filtered_Y_2 <- update_seq_classes_filtered_Y %>%
  mutate(scientific_name = case_when(
    scientific_name == "GorGor"~ "Gorilla", 
    scientific_name == "HomSap"~ "Human",
    scientific_name == "MacFas"~ "Macaque",
    scientific_name == "PanPan"~ "Bonobo",
    scientific_name == "PanTro"~ "Chimpanzee",
    scientific_name == "PonAbe"~ "S.orangutan",
    scientific_name == "PonPyg"~ "B.orangutan",
    scientific_name == "SymSyn"~ "Siamang",
  ))

classes <- update_seq_classes_filtered_Y_2 |>
  mutate(start = to_num(start), end = to_num(end)) |>
  filter(scientific_name %in% species_order, is.finite(start), is.finite(end)) |>
  mutate(species_idx = match(scientific_name, species_order))

chr_len <- classes |>
  dplyr::group_by(scientific_name) |>
  summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") |>
  dplyr::rename(Species = scientific_name)

if (!"Macaque" %in% chr_len$Species) {
  chr_len <- bind_rows(chr_len, tibble(Species = "Macaque", chr_len = macfas_len_bp))
}

chr_len <- chr_len |>
  filter(Species %in% species_order) |>
  mutate(species_idx = match(Species, species_order))

# --- gene copies ---
copies <- gene_details_Y_arrows_updated |>
  filter(Species %in% species_order, gene_family_symbol %in% moved_families) |>
  mutate(arrow_start = to_num(arrow_start),
         arrow_end   = to_num(arrow_end),
         mid = (arrow_start + arrow_end)/2) |>
  inner_join(chr_len |> select(Species, chr_len, species_idx), by = "Species") |>
  filter(is.finite(mid))

# --- cluster within species ---
clusters <- copies |>
  group_by(gene_family_symbol, Species, species_idx) |>
  summarise(df = list(cluster_positions(mid, tol_bp = 5e5)), .groups = "drop") |>
  unnest(df)

# --- prepare data for D3 ---
seq_classes_df <- classes |>
  select(species = scientific_name, species_idx, start, end, classification) |>
  arrange(species_idx)

gene_clusters_df <- clusters |>
  mutate(species_idx = match(Species, species_order)) |>
  select(family = gene_family_symbol, species = Species, species_idx, position = y, copies = n, cluster)

# Create links - intelligently handle splits and merges
links_df <- map_dfr(moved_families, function(fam) {
  fam_clusters <- filter(clusters, gene_family_symbol == fam)
  links_list <- list()
  
  # Check if this family should use nearest neighbor linking
  use_nearest_neighbor <- fam %in% nearest_neighbor_families
  
  for (i in seq_len(length(species_order) - 1)) {
    sp_from <- species_order[i]
    sp_to <- species_order[i + 1]
    
    from_data <- filter(fam_clusters, Species == sp_from) %>% arrange(y)
    to_data <- filter(fam_clusters, Species == sp_to) %>% arrange(y)
    
    if (nrow(from_data) == 0 || nrow(to_data) == 0) next
    
    n_from <- nrow(from_data)
    n_to <- nrow(to_data)
    
    if (n_from == n_to) {
      # Same number - direct mapping
      for (j in seq_len(n_from)) {
        links_list[[length(links_list) + 1]] <- tibble(
          family = fam,
          source_species = sp_from,
          source_idx = i,
          source_pos = from_data$y[j],
          target_species = sp_to,
          target_idx = i + 1,
          target_pos = to_data$y[j],
          is_split = FALSE
        )
      }
    } else if (use_nearest_neighbor) {
      # Different numbers BUT use nearest neighbor matching
      # This prevents the "spaghetti" effect for multi-cluster families
      if (n_from < n_to) {
        # More clusters in target - each source connects to nearest target
        for (j in seq_len(n_from)) {
          # Find nearest target cluster
          distances <- abs(to_data$y - from_data$y[j])
          nearest_idx <- which.min(distances)
          
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[j],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[nearest_idx],
            is_split = FALSE
          )
        }
        # Connect remaining target clusters to their nearest source
        connected_targets <- sapply(seq_len(n_from), function(j) {
          which.min(abs(to_data$y - from_data$y[j]))
        })
        unconnected_targets <- setdiff(seq_len(n_to), connected_targets)
        for (k in unconnected_targets) {
          nearest_source_idx <- which.min(abs(from_data$y - to_data$y[k]))
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[nearest_source_idx],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[k],
            is_split = TRUE
          )
        }
      } else {
        # More clusters in source - each target connects to nearest source
        for (k in seq_len(n_to)) {
          # Find nearest source cluster
          distances <- abs(from_data$y - to_data$y[k])
          nearest_idx <- which.min(distances)
          
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[nearest_idx],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[k],
            is_split = FALSE
          )
        }
        # Connect remaining source clusters to their nearest target
        connected_sources <- sapply(seq_len(n_to), function(k) {
          which.min(abs(from_data$y - to_data$y[k]))
        })
        unconnected_sources <- setdiff(seq_len(n_from), connected_sources)
        for (j in unconnected_sources) {
          nearest_target_idx <- which.min(abs(to_data$y - from_data$y[j]))
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[j],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[nearest_target_idx],
            is_split = TRUE
          )
        }
      }
    } else {
      # Different numbers - all combinations (original behavior)
      for (j in seq_len(n_from)) {
        for (k in seq_len(n_to)) {
          links_list[[length(links_list) + 1]] <- tibble(
            family = fam,
            source_species = sp_from,
            source_idx = i,
            source_pos = from_data$y[j],
            target_species = sp_to,
            target_idx = i + 1,
            target_pos = to_data$y[k],
            is_split = TRUE
          )
        }
      }
    }
  }
  
  bind_rows(links_list)
})

centromeres_df <- centro_tbl |> 
  mutate(species_idx = match(Species, species_order)) |>
  select(species = Species, species_idx, cstart, cend)

# Create chromosome lengths list for D3
chr_lengths_list <- setNames(
  as.list(chr_len$chr_len),
  chr_len$Species
)

d3_data <- list(
  species = species_order,
  max_chr_len = max(chr_len$chr_len),
  chr_lengths = chr_lengths_list,
  centromeres = centromeres_df,
  seq_classes = seq_classes_df,
  gene_clusters = gene_clusters_df,
  links = links_df,
  families = moved_families,
  colors = list(
    # Y chromosome gene families - colorblind-safe palette
    
    # Ampliconic/multi-copy families (bright, attention-grabbing colors)
    BPY2 = "#E69F00",        # Bright orange
    CDY = "#D55E00",         # Vermillion (colorblind-safe red)
    DAZ = "#CC79A7",         # Reddish purple
    HSFY = "#56B4E9",        # Sky blue
    RBMY = "#009E73",        # Bluish green
    TSPY8 = "#F0E442",       # Bright yellow
    VCY = "#0072B2",         # Strong blue
    
    # Single-copy/conserved genes (medium saturation)
    "CCS-OFD1" = "#8B6914",  # Dark brown
    PRPY = "#5F7A1A",        # Dark turquoise
    FRG1Y = "#B84D00",       # Dark vermillion
    KRT18Y = "#4A9FD5",      # Medium blue
    TAF11L2 = "#6B8E23",     # Olive drab
    
    # Retroviruses and LOC genes (distinct but neutral)
    "MTRNR2-like 17" = "#9B7DB8",  # Soft purple
    "retrovirus K19" = "#A67C00",  # Medium brown
    
    # Mitochondrial-related (yellow-green spectrum)
    GLUD1Y = "#7A9A2E",      # Yellow-green
    AKI6 = "#00B589",        # Light bluish green
    
    # Gene families (cooler tones)
    FAM47AY = "#BD6E9A",     # Medium reddish purple
    ZNF = "#005A8F"          # Dark blue
    # Y nonmoved genes
    #`proline-rich protein, Y-linked` ="#C8955C", 
    #`MTRNR2-like 17` ="#B85A5A",
    #`adenylate kinase isoenzyme 6−like`="#9C7BB8",
    #`endogenous retrovirus group K member 19 Env polyprotein−like`="#5B9B82"
    #ARMCX1="#3498DB")
  ),
  class_colors = list(
    ANCESTRAL = "#D3D3D3",#FFD700
    PAR = "#D3D3D3",#90EE90
    AMPLICONIC = "#D3D3D3",#87CEEB
    SATELLITE = "#D3D3D3",#A0522D
    XTR = "#D3D3D3",#FFB6C1
    UNCLASSIFIED = "#D3D3D3",#D3D3D3
    PALINDROME = "#D3D3D3"#4169E1
  )
)

# Convert to JSON
json_data <- toJSON(d3_data, auto_unbox = TRUE, pretty = TRUE)

# Verify JSON was created
if (is.null(json_data) || nchar(json_data) == 0) {
  stop("Failed to create JSON data")
}

cat("JSON data created successfully, length:", nchar(json_data), "\n")

# Create standalone HTML file with updated label code
html_template <- '
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Gene Family Movement</title>
    <script src="https://d3js.org/d3.v6.min.js"></script>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 20px;
            background: white;
        }
        #chart { 
            margin: 0 auto; 
            display: block;
        }
    </style>
</head>
<body>
    <div id="chart"></div>
    <script>
    const data = DATA_PLACEHOLDER;
    
    const CHR_WIDTH = 30;
    const CENTRO_WIDTH = 17;
    
    const myWidth = 1200;
    const myHeight = 700;
    const padding = {left: 400, right: 80, top: 100, bottom: 50};'

html_content <- paste0(
  html_template,
  '
    
    const yScale = d3.scaleLinear()
      .domain([0, data.max_chr_len])
      .range([padding.top, myHeight - padding.bottom]);
    
    const xScale = d3.scalePoint()
      .domain(data.species)
      .range([padding.left + CHR_WIDTH, myWidth - padding.right - CHR_WIDTH])
      .padding(0.5);
    
    const chrInfo = {};
    
    // First pass: create basic info with actual chromosome lengths
    data.species.forEach(sp => {
      const chrData = data.centromeres.find(c => c.species === sp);
      // Find the actual chromosome length for this species
      const chrLength = data.chr_lengths ? data.chr_lengths[sp] : data.max_chr_len;
      
      chrInfo[sp] = {
        x: xScale(sp),
        chrLength: chrLength,
        yStart: yScale(0),
        yEnd: yScale(chrLength),
        centroStart: yScale(chrData.cstart),
        centroEnd: yScale(chrData.cend),
        chrWidth: CHR_WIDTH
      };
    });
    
    const svg = d3.select("#chart")
      .append("svg")
      .attr("width", myWidth)
      .attr("height", myHeight);
    
    const chrGroup = svg.append("g").attr("class", "chromosomes");
    
    data.species.forEach(sp => {
      const info = chrInfo[sp];
      const yStart = info.yStart;
      const yEnd = info.yEnd;
      const centroStart = info.centroStart;
      const centroEnd = info.centroEnd;
      const x = info.x;
      
      const topCap = CHR_WIDTH * 0.3;
      const bottomCap = CHR_WIDTH * 0.3;
      
      // Fixed centromere visual size in pixels for consistency
      const centroVisualSize = 30; // pixels
      const centroCenterY = (centroStart + centroEnd) / 2;
      const centroVisualStart = centroCenterY - centroVisualSize / 2;
      const centroVisualEnd = centroCenterY + centroVisualSize / 2;
      
      const points = [];
      const nPoints = 100;
      
      for (let i = 0; i <= nPoints; i++) {
        const t = i / nPoints;
        const y = yStart + t * (yEnd - yStart);
        
        let width;
        if (y < yStart + topCap) {
          const capT = (y - yStart) / topCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y > yEnd - bottomCap) {
          const capT = (yEnd - y) / bottomCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y >= centroVisualStart && y <= centroVisualEnd) {
          // Use fixed visual size for centromere
          const centroT = (y - centroVisualStart) / (centroVisualEnd - centroVisualStart);
          const pinchFactor = Math.sin(centroT * Math.PI);
          width = CENTRO_WIDTH/2 + (CHR_WIDTH/2 - CENTRO_WIDTH/2) * pinchFactor;
        } else {
          const chrCenter = (yStart + yEnd) / 2;
          const distFromCenter = Math.abs(y - chrCenter) / ((yEnd - yStart) / 2);
          width = CHR_WIDTH/2 * (0.95 + 0.05 * (1 - distFromCenter * 0.3));
        }
        
        points.push({y: y, width: width});
      }
      
      let pathData = "M" + (x + points[0].width) + "," + points[0].y;
      for (let i = 1; i < points.length; i++) {
        pathData += " L" + (x + points[i].width) + "," + points[i].y;
      }
      for (let i = points.length - 1; i >= 0; i--) {
        pathData += " L" + (x - points[i].width) + "," + points[i].y;
      }
      pathData += " Z";
      
      chrGroup.append("path")
        .attr("d", pathData)
        .attr("fill", "#EAEAEA")
        .attr("stroke", "#000")
        .attr("stroke-width", 2);
    });
    
    const defs = svg.append("defs");
    
    data.species.forEach(sp => {
      const info = chrInfo[sp];
      const clipPath = defs.append("clipPath")
        .attr("id", "clip-" + sp);
      
      const yStart = info.yStart;
      const yEnd = info.yEnd;
      const centroStart = info.centroStart;
      const centroEnd = info.centroEnd;
      const x = info.x;
      
      const topCap = CHR_WIDTH * 0.3;
      const bottomCap = CHR_WIDTH * 0.3;
      
      // Fixed centromere visual size in pixels for consistency
      const centroVisualSize = 30; // pixels
      const centroCenterY = (centroStart + centroEnd) / 2;
      const centroVisualStart = centroCenterY - centroVisualSize / 2;
      const centroVisualEnd = centroCenterY + centroVisualSize / 2;
      
      const points = [];
      const nPoints = 100;
      
      for (let i = 0; i <= nPoints; i++) {
        const t = i / nPoints;
        const y = yStart + t * (yEnd - yStart);
        
        let width;
        if (y < yStart + topCap) {
          const capT = (y - yStart) / topCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y > yEnd - bottomCap) {
          const capT = (yEnd - y) / bottomCap;
          width = CHR_WIDTH/2 * Math.sqrt(1 - Math.pow(1 - capT, 2));
        } else if (y >= centroVisualStart && y <= centroVisualEnd) {
          // Use fixed visual size for centromere
          const centroT = (y - centroVisualStart) / (centroVisualEnd - centroVisualStart);
          const pinchFactor = Math.sin(centroT * Math.PI);
          width = CENTRO_WIDTH/2 + (CHR_WIDTH/2 - CENTRO_WIDTH/2) * pinchFactor;
        } else {
          const chrCenter = (yStart + yEnd) / 2;
          const distFromCenter = Math.abs(y - chrCenter) / ((yEnd - yStart) / 2);
          width = CHR_WIDTH/2 * (0.95 + 0.05 * (1 - distFromCenter * 0.3));
        }
        
        points.push({y: y, width: width});
      }
      
      let pathData = "M" + (x + points[0].width) + "," + points[0].y;
      for (let i = 1; i < points.length; i++) {
        pathData += " L" + (x + points[i].width) + "," + points[i].y;
      }
      for (let i = points.length - 1; i >= 0; i--) {
        pathData += " L" + (x - points[i].width) + "," + points[i].y;
      }
      pathData += " Z";
      
      clipPath.append("path")
        .attr("d", pathData);
    });
    
    const classGroup = svg.append("g").attr("class", "seq-classes");
    
    data.seq_classes.forEach(cls => {
      const info = chrInfo[cls.species];
      const x = info.x;
      
      classGroup.append("rect")
        .attr("x", x - info.chrWidth/2 + 1)
        .attr("y", yScale(cls.start))
        .attr("width", info.chrWidth - 2)
        .attr("height", yScale(cls.end) - yScale(cls.start))
        .attr("fill", data.class_colors[cls.classification] || "#CCC")
        .attr("opacity", 0.5)
        .attr("clip-path", "url(#clip-" + cls.species + ")");
    });
    
    const linkGen = d3.linkHorizontal()
      .x(d => d[0])
      .y(d => d[1]);
    
    const linksGroup = svg.append("g").attr("class", "links");
    
    data.links.forEach(link => {
      const color = data.colors[link.family];
      const sourceX = xScale(link.source_species) + CHR_WIDTH/2;
      const targetX = xScale(link.target_species) - CHR_WIDTH/2;
      const sourceY = yScale(link.source_pos);
      const targetY = yScale(link.target_pos);
      const linkWidth = link.is_split ? 1 : 2;
    
    // shading of the lines (outcommanded! )  
    //  linksGroup.append("path")
    //    .attr("d", linkGen({
    //      source: [sourceX, sourceY],
    //      target: [targetX, targetY]
    //    }))
    //    .attr("fill", "none")
    //    .attr("stroke", color)
    //    .attr("stroke-width", linkWidth * 3)
    //    .attr("opacity", 0.12);
      
      linksGroup.append("path")
        .attr("d", linkGen({
          source: [sourceX, sourceY],
          target: [targetX, targetY]
        }))
        .attr("fill", "none")
        .attr("stroke", color)
        .attr("stroke-width", linkWidth)
        .attr("opacity", link.is_split ? 0.4 : 0.7);
    });
    
    const markersGroup = svg.append("g").attr("class", "markers");
    
    data.gene_clusters.forEach(gc => {
      const info = chrInfo[gc.species];
      const x = info.x;
      const color = data.colors[gc.family];
      const y = yScale(gc.position);
      const blockWidth = info.chrWidth - 4;
      
      markersGroup.append("rect")
        .attr("x", x - blockWidth/2)
        .attr("y", y - 3)
        .attr("width", blockWidth)
        .attr("height", 6)
        .attr("fill", color)
        .attr("opacity", 0.8)
        .attr("clip-path", "url(#clip-" + gc.species + ")");
    });
    
    // Get unique families and their appearance info
    const familyLabels = [];

    // For each family, find where to label it (first or last species)
    data.families.forEach(family => {
      const familyClusters = data.gene_clusters.filter(gc => gc.family === family);
      if (familyClusters.length === 0) return;
      
      // Sort by species index
      familyClusters.sort((a, b) => a.species_idx - b.species_idx);
      const firstAppearance = familyClusters[0];
      const lastAppearance = familyClusters[familyClusters.length - 1];
      
      // Decide: label on left (first species) or right (last species)
      const firstSpeciesIdx = firstAppearance.species_idx;
      const labelOnRight = firstSpeciesIdx > 1; // If it does not appear in first 2 species, label on right
      
      const labelSpecies = labelOnRight ? lastAppearance.species : firstAppearance.species;
      const clustersInLabelSpecies = familyClusters.filter(
        gc => gc.species === labelSpecies
      );
      
      // Calculate average position
      const avgPosition = clustersInLabelSpecies.reduce((sum, gc) => sum + gc.position, 0) / 
                          clustersInLabelSpecies.length;
      
      familyLabels.push({
        family: family,
        species: labelSpecies,
        originalY: avgPosition,
        y: avgPosition,
        labelOnRight: labelOnRight
      });
    });

    // Sort by species index first, then by position within species
    familyLabels.sort((a, b) => {
      const speciesCompare = data.species.indexOf(a.species) - data.species.indexOf(b.species);
      if (speciesCompare !== 0) return speciesCompare;
      return a.y - b.y;
    });

    // Group labels by species AND position AND side
    const groupSpacingThreshold = 8;
    const MAX_GENES_PER_GROUP = 4;

    // Separate left and right labels
    const leftLabels = familyLabels.filter(l => !l.labelOnRight);
    const rightLabels = familyLabels.filter(l => l.labelOnRight);

    // Function to create groups
    function createGroups(labels) {
      const groups = [];
      if (labels.length === 0) return groups;
      
      let currentGroup = [labels[0]];
      
      for (let i = 1; i < labels.length; i++) {
        const prevLabel = labels[i - 1];
        const currLabel = labels[i];
        
        if (prevLabel.species === currLabel.species) {
          const prevYScreen = yScale(prevLabel.y);
          const currYScreen = yScale(currLabel.y);
          
          if (currYScreen - prevYScreen < groupSpacingThreshold) {
            if (currentGroup.length >= MAX_GENES_PER_GROUP) {
              groups.push(currentGroup);
              currentGroup = [currLabel];
            } else {
              currentGroup.push(currLabel);
            }
          } else {
            groups.push(currentGroup);
            currentGroup = [currLabel];
          }
        } else {
          groups.push(currentGroup);
          currentGroup = [currLabel];
        }
      }
      groups.push(currentGroup);
      return groups;
    }

    const leftGroups = createGroups(leftLabels);
    const rightGroups = createGroups(rightLabels);

    // Adjust positions to prevent overlap - LEFT SIDE
    const leftSpeciesGroups = {};
    leftGroups.forEach(group => {
      const species = group[0].species;
      if (!leftSpeciesGroups[species]) {
        leftSpeciesGroups[species] = [];
      }
      leftSpeciesGroups[species].push(group);
    });

    Object.keys(leftSpeciesGroups).forEach(species => {
      const groups = leftSpeciesGroups[species];
      for (let pass = 0; pass < 5; pass++) {
        for (let i = 1; i < groups.length; i++) {
          const prevGroup = groups[i - 1];
          const currGroup = groups[i];
          const prevGroupY = yScale(prevGroup[0].y);
          const currGroupY = yScale(currGroup[0].y);
          const requiredSpacing = 16;
          
          if (currGroupY - prevGroupY < requiredSpacing) {
            const adjustment = yScale.invert(prevGroupY + requiredSpacing) - currGroup[0].y;
            currGroup.forEach(label => {
              label.y += adjustment;
            });
          }
        }
      }
    });

    // Adjust positions to prevent overlap - RIGHT SIDE
    const rightSpeciesGroups = {};
    rightGroups.forEach(group => {
      const species = group[0].species;
      if (!rightSpeciesGroups[species]) {
        rightSpeciesGroups[species] = [];
      }
      rightSpeciesGroups[species].push(group);
    });

    Object.keys(rightSpeciesGroups).forEach(species => {
      const groups = rightSpeciesGroups[species];
      for (let pass = 0; pass < 5; pass++) {
        for (let i = 1; i < groups.length; i++) {
          const prevGroup = groups[i - 1];
          const currGroup = groups[i];
          const prevGroupY = yScale(prevGroup[0].y);
          const currGroupY = yScale(currGroup[0].y);
          const requiredSpacing = 16;
          
          if (currGroupY - prevGroupY < requiredSpacing) {
            const adjustment = yScale.invert(prevGroupY + requiredSpacing) - currGroup[0].y;
            currGroup.forEach(label => {
              label.y += adjustment;
            });
          }
        }
      }
    });

    // Draw LEFT labels
    leftGroups.forEach(group => {
      const labelSpecies = group[0].species;
      const info = chrInfo[labelSpecies];
      const x = info.x;
      
      if (group.length === 1) {
        const fam = group[0];
        const color = data.colors[fam.family] || "#999999";
        const originalYScreen = yScale(fam.originalY);
        const adjustedYScreen = yScale(fam.y);
        
        if (Math.abs(adjustedYScreen - originalYScreen) > 2) {
          markersGroup.append("line")
            .attr("x1", x - info.chrWidth/2 - 5)
            .attr("y1", originalYScreen)
            .attr("x2", x - info.chrWidth/2 - 8)
            .attr("y2", adjustedYScreen)
            .attr("stroke", color)
            .attr("stroke-width", 0.5)
            .attr("opacity", 0.5);
        }
        
        markersGroup.append("text")
          .attr("x", x - info.chrWidth/2 - 8)
          .attr("y", adjustedYScreen + 4)
          .attr("text-anchor", "end")
          .attr("font-size", "11px")
          .attr("fill", color)
          .attr("font-weight", "bold")
          .text(fam.family);
        
      } else {
        const groupY = yScale(group[0].y);
        
        group.forEach(fam => {
          const color = data.colors[fam.family] || "#999999";
          const originalYScreen = yScale(fam.originalY);
          
          if (Math.abs(groupY - originalYScreen) > 2) {
            markersGroup.append("line")
              .attr("x1", x - info.chrWidth/2 - 5)
              .attr("y1", originalYScreen)
              .attr("x2", x - info.chrWidth/2 - 8)
              .attr("y2", groupY)
              .attr("stroke", color)
              .attr("stroke-width", 0.5)
              .attr("opacity", 0.5);
          }
        });
        
        let currentX = x - info.chrWidth/2 - 8;
        
        group.forEach((fam, idx) => {
          const color = data.colors[fam.family] || "#999999";
          
          if (idx > 0) {
            const arrowText = markersGroup.append("text")
              .attr("text-anchor", "end")
              .attr("font-size", "11px")
              .attr("font-weight", "bold")
              .attr("fill", "#999")
              .text(" ~ ");
            
            const arrowBBox = arrowText.node().getBBox();
            arrowText.attr("x", currentX).attr("y", groupY + 4);
            currentX -= arrowBBox.width;
          }
          
          const geneText = markersGroup.append("text")
            .attr("text-anchor", "end")
            .attr("font-size", "11px")
            .attr("font-weight", "bold")
            .attr("fill", color)
            .text(fam.family + " ");
          
          const geneBBox = geneText.node().getBBox();
          geneText.attr("x", currentX).attr("y", groupY + 4);
          currentX -= geneBBox.width;
        });
      }
    });

    // Draw RIGHT labels (with different styling to make them stand out)
    rightGroups.forEach(group => {
      const labelSpecies = group[0].species;
      const info = chrInfo[labelSpecies];
      const x = info.x;
      
      if (group.length === 1) {
        const fam = group[0];
        const color = data.colors[fam.family] || "#999999";
        const originalYScreen = yScale(fam.originalY);
        const adjustedYScreen = yScale(fam.y);
        
        if (Math.abs(adjustedYScreen - originalYScreen) > 2) {
          markersGroup.append("line")
            .attr("x1", x + info.chrWidth/2 + 5)
            .attr("y1", originalYScreen)
            .attr("x2", x + info.chrWidth/2 + 8)
            .attr("y2", adjustedYScreen)
            .attr("stroke", color)
            .attr("stroke-width", 0.5)
            .attr("opacity", 0.5);
        }
        
        markersGroup.append("text")
          .attr("x", x + info.chrWidth/2 + 8)
          .attr("y", adjustedYScreen + 4)
          .attr("text-anchor", "start")
          .attr("font-size", "11px")
          .attr("fill", color)
          .attr("font-weight", "bold")
          .attr("font-style", "italic")
          .text(fam.family);
        
      } else {
        const groupY = yScale(group[0].y);
        
        group.forEach(fam => {
          const color = data.colors[fam.family] || "#999999";
          const originalYScreen = yScale(fam.originalY);
          
          if (Math.abs(groupY - originalYScreen) > 2) {
            markersGroup.append("line")
              .attr("x1", x + info.chrWidth/2 + 5)
              .attr("y1", originalYScreen)
              .attr("x2", x + info.chrWidth/2 + 8)
              .attr("y2", groupY)
              .attr("stroke", color)
              .attr("stroke-width", 0.5)
              .attr("opacity", 0.5);
          }
        });
        
        let currentX = x + info.chrWidth/2 + 8;
        
        group.forEach((fam, idx) => {
          const color = data.colors[fam.family] || "#999999";
          
          if (idx > 0) {
            const arrowText = markersGroup.append("text")
              .attr("text-anchor", "start")
              .attr("font-size", "11px")
              .attr("font-weight", "bold")
              .attr("fill", "#999")
              .text(" ~ ");
            
            const arrowBBox = arrowText.node().getBBox();
            arrowText.attr("x", currentX).attr("y", groupY + 4);
            currentX += arrowBBox.width;
          }
          
          const geneText = markersGroup.append("text")
            .attr("text-anchor", "start")
            .attr("font-size", "11px")
            .attr("font-weight", "bold")
            .attr("font-style", "italic")
            .attr("fill", color)
            .text(fam.family + " ");
          
          const geneBBox = geneText.node().getBBox();
          geneText.attr("x", currentX).attr("y", groupY + 4);
          currentX += geneBBox.width;
        });
      }
    });

    // SPECIES LABELS (on top)
    const labelGroup = svg.append("g").attr("class", "labels");

    data.species.forEach(sp => {
      labelGroup.append("text")
        .attr("x", xScale(sp))
        .attr("y", padding.top - 15)
        .attr("text-anchor", "middle")
        .attr("font-size", "13.4px")
        .attr("font-weight", "bold")
        .text(sp);
    });

    // Y-AXIS (chromosome length on right)
    const yAxis = d3.axisRight(yScale)
      .tickFormat(d => Math.round(d/1e6) + " Mb")
      .ticks(10);

    svg.append("g")
      .attr("transform", "translate(" + (myWidth - padding.right + 20) + ", 0)")
      .call(yAxis)
      .selectAll("text")
      .style("font-size", "14px")
      .selectAll("line")
      .attr("x2", 15);

    // TITLE
    svg.append("text")
      .attr("x", myWidth/2)
      .attr("y", 25)
      .attr("text-anchor", "middle")
      .attr("font-size", "18px")
      .attr("font-weight", "bold")
      .text("Gene Family Movement Across Primate Species");
    
   // Add a download SVG button
      const downloadButton = document.createElement("button");
      downloadButton.textContent = "Download SVG";
      downloadButton.style.cssText = "position: fixed; top: 10px; right: 10px; padding: 10px 20px; font-size: 14px; cursor: pointer; z-index: 1000;";
      document.body.appendChild(downloadButton);
      
      downloadButton.addEventListener("click", function() {
        const svgElement = document.querySelector("svg");
        const svgData = new XMLSerializer().serializeToString(svgElement);
        const svgBlob = new Blob([svgData], {type: "image/svg+xml;charset=utf-8"});
        const svgUrl = URL.createObjectURL(svgBlob);
        const downloadLink = document.createElement("a");
        downloadLink.href = svgUrl;
        downloadLink.download = "gene_movement.svg";
        document.body.appendChild(downloadLink);
        downloadLink.click();
        document.body.removeChild(downloadLink);
        URL.revokeObjectURL(svgUrl);
        });

    console.log("Visualization complete!");
    console.log("Left labels:", leftGroups.length, "Right labels:", rightGroups.length);
  </script>
    </body>
    </html>
    ')

# Replace placeholder with actual data
html_content <- sub("DATA_PLACEHOLDER", json_data, html_content, fixed = TRUE)

# Write HTML file
output_file <- "gene_movement_changingY_updated_nonbackground_commonnames_NE2W.html"
writeLines(html_content, output_file)
cat("HTML file created:", output_file, "\n")
cat("Open it in your browser to view the visualization.\n")
