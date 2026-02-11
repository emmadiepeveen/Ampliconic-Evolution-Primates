### Supplementary figures ### 
# This script is for making supplementary figures 1,4,5 and 6. 


#### Supplementary figure 1 Copy number ~ variance in copy number ####
# --- packages ---
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)

# --- Data ---
df<- read_excel("~/Documents/PhD Aarhus/Amplicons/Amplicon_summary_count_expression_genefamily.xlsx") 


# 1) Identify species columns
species_cols <- names(df)[vapply(df, is.numeric, logical(1))]

gene_col <- intersect(names(df), c("Gene_family","GeneFamily","GeneFamilyID","Family","Gene"))[1]

# 2) median and variance per gene family
df_summarized <- df %>%
  rowwise() %>%
  mutate(
    median_copy = median(c_across(all_of(species_cols)), na.rm = TRUE),
    var_copy    = var(c_across(all_of(species_cols)),    na.rm = TRUE),
    non_na      = sum(!is.na(c_across(all_of(species_cols))))
  ) %>%
  ungroup() %>%
  # need at least 2 values to compute a variance
  filter(non_na >= 2)


# 3) highlight the Cancer-Testis genes: 
CT_genes <- c(
  "VCX", "MAGEB", "SSX", "GAGE", "PAGE", "CENPVL", "XAGE", "CXorf49",
  "NXF", "CT45", "CT47", "SAGE", "SPANX", "CXorf51", "CTAG",
  "TSPY", "VCY", "RBMY"
)

# Add a column marking CT genes
df_summarized <- df_summarized %>%
  mutate(
    CT_flag = `Gene family` %in% CT_genes,
    CT_group = if_else(CT_flag, "CT gene", "Other")
  )

chrom_colors <- c(X = "#E377C2", Y = "#1F77B4")  # pink for X, blue for Y

# 4) plot Variance copy number ~ median copy number
ggplot(df_summarized, aes(x = median_copy, y = var_copy)) +
  # Dots: chromosome colours
  geom_point(aes(color = Chromosome), size = 4, alpha = 0.85) +
  # Labels: same colour as dots, CT genes bold + with a star
  ggrepel::geom_text_repel(
    aes(
      label = ifelse(CT_flag, paste0(`Gene family`, "*"), `Gene family`),
      color = Chromosome,
      fontface = ifelse(CT_flag, "bold", "plain")
    ),
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
    scale_color_manual(
    values = chrom_colors,
    labels = c(X = "X chromosome", Y = "Y chromosome")
  ) +
  geom_smooth(method = "lm", se = FALSE, col = "gray") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(Log[10]*"(Median copy number across species)"),
    y = expression(Log[10]*"(Variance in copy number across species)"),
    color = "* = Cancer-Testis gene"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

ggsave("~/Downloads/median_vs_variance_XYgenes.pdf", width = 15, height = 12, dpi = 300)





### Supplementary figure 5 & 6 Gene location and orienation per species ####


# Packages
library(readr)
library(dplyr)    
library(tidyr)
library(ggplot2)
library(patchwork)


# X chromosome data
gene_details_X <- read_delim("~/Downloads/all_families_gene_details_with_clusters-3.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(Start = col_number(), 
                                              End = col_number()), trim_ws = TRUE)
#mutate that if gene lays in a palindrome, the class is called palindrome
gene_details_X <- gene_details_X %>%
  mutate(Class = ifelse(in_palindrome == "yes", "PALINDROME", Class))

# filter out rows with genes that are not classified as ampliconic genes (so have an NA in the gene_family_symbol)
gene_details_X<-gene_details_X %>% drop_na(gene_family_symbol)


# Y chromosome 
gene_details_Y <- read_delim("~/Downloads/gene_details_updated_with_palindromes_coordinates_y.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(Start = col_number(), 
                                              End = col_number()), trim_ws = TRUE)


#mutate that if gene lays in a palindrome, the class is called palindrome
gene_details_Y <- gene_details_Y %>%
  mutate(Class = ifelse(in_palindrome == "yes", "PALINDROME", Class))

# filter out rows with genes that are not classified as ampliconic genes (so have an NA in the gene_family_symbol)
gene_details_Y<-gene_details_Y %>% drop_na(gene_family_symbol)


# Sort data by Species, gene_family_symbol, and Start coordinate
# Define arrow direcitons -> a new column for arrow directions, based on class. 
# Add an arrow direction column
gene_details_X_sorted <- gene_details_X %>%
  arrange(Species, gene_family_symbol, Start) %>%
  mutate(
    direction = ifelse(Strand == "+", 1, -1),  # Direction based on Strand
    arrow_start = ifelse(direction == -1, End, Start),  # Left-facing: End -> Start
    arrow_end = ifelse(direction == -1, Start, End)     # Right-facing: Start -> End
  )

gene_details_Y_sorted <- gene_details_Y %>%
  arrange(Species, gene_family_symbol, Start) %>%
  mutate(
    direction = ifelse(Strand == "+", 1, -1),  # Direction based on Strand
    arrow_start = ifelse(direction == -1, End, Start),  # Left-facing: End -> Start
    arrow_end = ifelse(direction == -1, Start, End)     # Right-facing: Start -> End
  )


# mutate the scientific names within the seqClass dataframe to the same as the species names
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


# Function to integrate palindrome data into seqClasses data -> so that all are within the same dataframe #
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



# FOR GENE COORDINATES AND ORIENTATION
# Define your gene name mapping BEFORE the function
unique(gene_details_X_arrows$gene_family_symbol)
unique(gene_name_mapping)
gene_name_mapping <- c(
  "ARMCX2" = "ARMCX",
  "CENPVL2" = "CENPVL",
  "CSAG1" = "CSAG",
  "CSF2RA" = "CSF2RA",
  "CT45A8" = "CT45",
  "CT47C1" = "CT47",
  "CXorf49B" = "CXorf49",
  "CXorf51B" = "CXorf51",
  "DMRTC1B" = "DMRTC1B",
  "EOLA1" = "EOLA",
  "ETDA" = "ETD",
  "F8A2" = "F8",
  "FAM156B" = "FAM156",
  "FAM236C" = "FAM236",
  "GAGE1" = "GAGE",
  "GPRASP2" = "GPRASP",
  "H2AB3" = "H2AB",
  "H2BW2" = "H2BW",
  "HSFX3" = "HSFX3",
  "IKBKG" = "IKBKG",
  "LAGE3" = "CTAG",
  "MAGEB18" = "MAGEB",
  "MAGED1" = "MAGED",
  "NUDT10" = "NUDT10",
  "NXF2" = "NXF",
  "OPN1LW" = "OPN1LW",
  "PABPC1L2B" = "PABPC1L2",
  "PNMA6E" = "PNMA",
  "PWWP4" = "PWWP",
  "RAB40A" = "RAB40",
  "RHOXF2B" = "RHOXF2",
  "RPL36A-HNRNPH2" = "RPL36A-HNRNPH2",
  "SAGE1" = "SAGE1",
  "SMIM10L2A" = "SMIM10",
  "SPACA5" = "SPACA",
  "SPANXA1" = "SPANX",
  "SPIN2A" = "SPIN",
  "SSX4B" = "SSX",
  "TBL1X" = "TBL1X",
  "TCEAL8" = "TCEAL8",
  "TCP11X2" = "TCP11X",
  "TEX28" = "TEX28",
  "TMEM185A" = "TMEM185A",
  "TMSB15B" = "TMSB15B",
  "VCX2" = "VCX",
  "XAGE1A" = "XAGE1",
  "ZXDA" = "ZXD",
  "putative uncharacterized protein FLJ39060" = "FLJ39060",
  "uncharacterized LOC129476473" = "LOC129476473",
  "HSFX1" = "HSFX1",
  "collagen alpha-4(IV) chain-like" = "Collagen alpha-IV chain-like",
  "endogenous retrovirus group K member 6 Env polyprotein-like" = "retrovirus K6",
  "uncharacterized LOC129138873" = "LOC129138873",
  "uncharacterized LOC129053094" = "LOC129053094"
  # Add all your gene mappings here
)

# Function to generate plot for each species
plot_species_with_background <- function(species_name, data_arrows, data_classes, max_coord, gene_mapping) {
  # Filter the seqClasses data for the specific species
  classes_data <- filter(data_classes, scientific_name == species_name)
  
  # Filter the arrows data for the specific species
  arrows_data <- filter(data_arrows, Species == species_name)
  
  # Apply manual gene name mapping
  arrows_data <- arrows_data %>%
    mutate(gene_family_simple = ifelse(
      gene_family_symbol %in% names(gene_mapping),
      gene_mapping[gene_family_symbol],
      gene_family_symbol  # Keep original if not in mapping
    ))
  
  ggplot() +
    # Background rectangles for classifications
    geom_rect(
      data = classes_data,
      aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf,
        fill = classification
      ),
      alpha = 0.25
    ) +
    # Overlay the arrow plot
    geom_segment(
      data = arrows_data,
      aes(
        x = arrow_start,
        xend = arrow_end,
        y = gene_family_simple,
        yend = gene_family_simple,
        color = Class,
        linetype = ifelse(arrow_start < arrow_end, "forward", "reverse")
      ),
      arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed"),
      linewidth = 1.2,
      show.legend = c(color = FALSE, linetype = TRUE)  # Hide color legend, show linetype legend
    ) +
    labs(
      title = species_name,  # Just species name
      x = "Coordinate",
      y = "Gene Family Symbol"
    ) +
    # Set x-axis limits to the maximum coordinate across all species
    coord_cartesian(xlim = c(0, max_coord)) +
    scale_fill_manual(
      name = "Classification",  # Legend title for fill
      values = c(
        "ANCESTRAL" = "yellow",
        "PAR" = "lightgreen",
        "AMPLICONIC" = "lightblue",
        "SATELLITE" = "brown", 
        "XTR" = "lightpink", 
        "UNCLASSIFIED" = "gray", 
        "PALINDROME" = "blue"
      )
    ) +
    scale_color_manual(
      name = "Classification",  # Same legend title
      values = c(
        "PALINDROME" = "blue",
        "AMPLICONIC" = "black",
        "ANCESTRAL" = "black",
        "XTR" = "black",
        "Unknown" = "black",
        "SATELLITE" = "black"
      )
    ) +
    scale_linetype_manual(
      name = "Strand",
      values = c("forward" = "solid", "reverse" = "solid"),
      labels = c("forward" = "→ Forward", "reverse" = "← Reverse"),
      guide = guide_legend(
        override.aes = list(alpha = 0)  # Make everything invisible
      )
    ) +
    # Custom gene family order
    scale_y_discrete(
      limits = custom_gene_order  # Define this before running the function
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.5, "lines")
    )
}

# Calculate the maximum coordinate across all species
max_coordinate <- max(
  max(update_seq_classes_filtered_X$end, na.rm = TRUE),
  max(gene_details_X_arrows$arrow_end, na.rm = TRUE)
)

# Define your custom gene family order (these should be the NEW names after mapping)
custom_gene_order <- c("CSF2RA", "VCX", "TBL1X", "retrovirus K6", "SPACA", "SSX", "NUDT10", "CENPVL2", "FLJ39060", "XAGE1",
                       "FAM156", "MAGED", "GAGE", "ZXD", "SPIN", "CXorf49", "DMRTC1B", "FAM236", "PABPC1L2B", "ARMCX", "TCP11X2",
                       "GPRASP2", "NXF", "TCP11X","TCEAL8", "GPRASP","RAB40A", "H2BW", "TMSB15", "RPL36A-HNRNPH2", "RHOXF2", "CT47", "collagen alfa-IV chain-like", 
                       "LOC129138873", "ETD", "SMIM10", "SAGE", "CT45", "LOC115932372","SPANX", "CXorf51","TMEM185A","LOC129053094","EOLA", "HSFX3","HSFX1", "CSAG", "MAGEB", "PWWP", "PNMA", "OPN1LW","TEX28",
                       "IKBKG", "CTAG", "F8", "H2AB")  # Add your genes in desired order
unique(custom_gene_order)


# Generate plots for each species
species_names <- c( "PanTro", "PanPan","HomSap", "GorGor", "PonPyg", "PonAbe", "SymSyn", "MacFas")

plots_with_background_X <- lapply(species_names, function(species) {
  plot_species_with_background(species, gene_details_X_arrows, update_seq_classes_filtered_X, max_coordinate, gene_name_mapping)
})

# View plots
print(plots_with_background_X[[1]])
print(plots_with_background_X[[2]])
print(plots_with_background_X[[3]])
print(plots_with_background_X[[4]])
print(plots_with_background_X[[5]])
print(plots_with_background_X[[6]])
print(plots_with_background_X[[7]])
print(plots_with_background_X[[8]])


# save individually 
ggsave("~/Downloads/X_genelocation_PanTro.jpg", plot = plots_with_background_X[[1]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/X_genelocation_PanPan.jpg", plot = plots_with_background_X[[2]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/X_genelocation_HomSap.jpg", plot = plots_with_background_X[[3]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/X_genelocation_GorGor.jpg", plot = plots_with_background_X[[4]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/X_genelocation_PonPyg.jpg", plot = plots_with_background_X[[5]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/X_genelocation_PonAbe.jpg", plot = plots_with_background_X[[6]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/X_genelocation_SymSyn.jpg", plot = plots_with_background_X[[7]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/X_genelocation_MacFas.jpg", plot = plots_with_background_X[[8]], width = 20, height = 10, dpi = 1000)


# Y chromosome
## UPDATED Y gene orientation and location
unique(gene_details_Y_arrows$gene_family_symbol)
unique(gene_name_mapping)
gene_name_mapping <- c(
  "BPY2" = "BPY2",
  "CDY1" = "CDY",
  "DAZ1" = "DAZ",
  "HSFY1" = "HSFY",
  "RBMY1B" = "RBMY",
  "TSPY8" = "TSPY",
  "centriole and centriolar satellite protein OFD1-like" = "CCS-OFD1",
  "proline-rich protein, Y-linked" = "PRPY",
  "protein FRG1-like" = "FRG1Y",
  "VCY1B" = "VCY",
  "keratin, type I cytoskeletal 18-like" = "KRT18Y",
  "MTRNR2-like 17" = "MTRNR2-like 17",
  "glutamate dehydrogenase 1, mitochondrial-like" = "GLUD1Y",
  "adenylate kinase isoenzyme 6-like" = "FAM236",
  "endogenous retrovirus group K member 19 Env polyprotein-like" = "retrovirus K19",
  "protein FAM47A-like" = "FAM47AY",
  "zinc finger protein 285-like" = "ZNF",
  "TATA-box binding protein associated factor 11 like protein 2-like" = "TAF11L2"
  # Add all your gene mappings here
)

# Function to generate plot for each species
#Use same function as above

# Calculate the maximum coordinate across all species
max_coordinate <- max(
  max(update_seq_classes_filtered_Y$end, na.rm = TRUE),
  max(gene_details_Y_arrows$arrow_end, na.rm = TRUE)
)

# Define your custom gene family order (these should be the NEW names after mapping)
custom_gene_order <- c("FAM47AY", "TSPY", "VCY", "PRPY", "KRT18Y", "HSFY", "retrovirus K19", "MTRNR2-like 17", "RBMY", "FRG1Y",
                       "CCS-OFD1", "DAZ", "ZNF", "CDY", "TAF11L2", "BPY2", "FRG1Y","GLUD1Y","FAM236")  # Add your genes in desired order

unique(custom_gene_order)


# Generate plots for each species
species_names <- c( "PanTro", "PanPan","HomSap", "GorGor", "PonPyg", "PonAbe", "SymSyn", "MacFas")

plots_with_background_Y <- lapply(species_names, function(species) {
  plot_species_with_background(species, gene_details_Y_arrows, update_seq_classes_filtered_Y, max_coordinate, gene_name_mapping)
})

# View plots
print(plots_with_background_Y[[1]])
print(plots_with_background_Y[[2]])
print(plots_with_background_Y[[3]])
print(plots_with_background_Y[[4]])
print(plots_with_background_Y[[5]])
print(plots_with_background_Y[[6]])
print(plots_with_background_Y[[7]])
print(plots_with_background_Y[[7]])


ggsave("~/Downloads/Y_genelocation_PanTro.jpg", plot = plots_with_background_Y[[1]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/Y_genelocation_PanPan.jpg", plot = plots_with_background_Y[[2]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/Y_genelocation_HomSap.jpg", plot = plots_with_background_Y[[3]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/Y_genelocation_GorGor.jpg", plot = plots_with_background_Y[[4]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/Y_genelocation_PonPyg.jpg", plot = plots_with_background_Y[[5]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/Y_genelocation_PonAbe.jpg", plot = plots_with_background_Y[[6]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/Y_genelocation_SymSyn.jpg", plot = plots_with_background_Y[[7]], width = 20, height = 10, dpi = 1000)
ggsave("~/Downloads/Y_genelocation_MacFas.jpg", plot = plots_with_background_Y[[8]], width = 20, height = 10, dpi = 1000)

# save all
combined_plot <- wrap_plots(plots_with_background_Y, ncol = 1)
ggsave("all_species_plots_Y_updated.pdf",
       combined_plot,
       device = "pdf",
       width  = 20,    # total width
       height = 8 * length(plots_with_background_X),  # total height
       limitsize = FALSE,
       units  = "in")



#### Supplementary figure 4 Distance to palindromes ####  
# This uses dataframe made in the section for Figure 5 &6  

# 1. Calculate observed distances: For each gene, find the distance to its nearest palindrome
# 2. Calculate expected distances: Randomly place genes across the chromosome many times and calculate the same distances
# 3. Compare: Test if observed < expected (genes cluster near palindromes)

library(tidyverse)

# Function to calculate distance to nearest palindrome
distance_to_nearest_palindrome <- function(gene_pos, palindrome_ranges) {
  if (nrow(palindrome_ranges) == 0) {
    return(NA)
  }
  
  distances <- apply(palindrome_ranges, 1, function(pal) {
    pal_start <- pal[1]
    pal_end <- pal[2]
    
    # Distance is 0 if gene overlaps palindrome
    if (gene_pos >= pal_start && gene_pos <= pal_end) {
      return(0)
    }
    
    # Otherwise, distance to nearest edge
    min(abs(gene_pos - pal_start), abs(gene_pos - pal_end))
  })
  
  return(min(distances))
}


# analysis function
analyze_gene_palindrome_clustering <- function(gene_details, palindrome_data, seq_classes_data, n_permutations = 10000) {
  
  # FILTER FOR ONLY PALINDROMES
  gene_details <- gene_details %>%
    filter(Species != "MacFas") 
  
  palindrome_data <- palindrome_data %>%
    filter(classification == "PALINDROME")
  
  # Get list of species
  species_list <- unique(gene_details$Species)
  n_species <- length(species_list)
  
  results_list <- list()
  
  # Loop through species with progress indicator
  for (i in seq_along(species_list)) {
    species_name <- species_list[i]
    
    cat(sprintf("\n[%d/%d] Processing %s...\n", i, n_species, species_name))
    
    # Get data for this species
    genes <- gene_details %>% filter(Species == species_name)
    palindromes <- palindrome_data %>% filter(scientific_name == species_name)
    
    # Get palindrome ranges
    pal_ranges <- palindromes %>%
      select(start, end) %>%
      as.matrix()
    
    # CORRECT: Calculate chromosome length from ALL sequence classes for this species
    chr_length <- seq_classes_data %>%
      filter(scientific_name == species_name) %>%
      pull(end) %>%
      max(na.rm = TRUE)
    
    cat(sprintf("  Chromosome length: %d bp\n", chr_length))
    
    # Observed distances - use gene midpoint
    genes <- genes %>%
      mutate(midpoint = (Start + End) / 2)
    
    observed_distances <- sapply(genes$midpoint, function(pos) {
      distance_to_nearest_palindrome(pos, pal_ranges)
    })
    
    obs_mean <- mean(observed_distances, na.rm = TRUE)
    obs_median <- median(observed_distances, na.rm = TRUE)
    
    # Calculate confidence interval for observed distances (95% CI)
    obs_se <- sd(observed_distances, na.rm = TRUE) / sqrt(sum(!is.na(observed_distances)))
    obs_ci_lower <- obs_mean - 1.96 * obs_se
    obs_ci_upper <- obs_mean + 1.96 * obs_se
    
    # Permutation test with progress
    cat(sprintf("  Running %d permutations...\n", n_permutations))
    random_means <- numeric(n_permutations)
    
    pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)
    for (j in 1:n_permutations) {
      # Random positions across FULL chromosome
      random_positions <- runif(nrow(genes), 0, chr_length)
      random_distances <- sapply(random_positions, function(pos) {
        distance_to_nearest_palindrome(pos, pal_ranges)
      })
      random_means[j] <- mean(random_distances, na.rm = TRUE)
      
      # Update progress bar every 100 iterations
      if (j %% 100 == 0) setTxtProgressBar(pb, j)
    }
    close(pb)
    
    # P-value: proportion of random means <= observed mean
    p_value <- sum(random_means <= obs_mean) / n_permutations
    
    # Store results for this species
    results_list[[i]] <- tibble(
      species = species_name,
      n_genes = nrow(genes),
      n_palindromes = nrow(palindromes),
      chr_length = chr_length,
      observed_mean_distance = obs_mean,
      observed_median_distance = obs_median,
      observed_ci_lower = obs_ci_lower,
      observed_ci_upper = obs_ci_upper,
      expected_mean_distance = mean(random_means),
      expected_std = sd(random_means),
      p_value = p_value,
      significant = p_value < 0.05
    )
    
    cat(sprintf("  Done! P-value = %.4f\n", p_value))
  }
  
  return(bind_rows(results_list))
}

# Run the analysis with the full sequence classes data for X
results <- analyze_gene_palindrome_clustering(
  gene_details_X_arrows, 
  update_seq_classes_filtered_X,
  update_seq_classes_filtered_X,  # Pass the full seq classes data
  n_permutations = 10000
)

# Display results
print(results)

# Save results
write_csv(results, "~/Downloads/X_gene_palindrome_clustering_results.csv")
# read in later
results_x<- read_csv("~/Downloads/X_gene_palindrome_clustering_results.csv")

# for the Y
# Run the analysis
results_y <- analyze_gene_palindrome_clustering(
  gene_details_Y_arrows, 
  update_seq_classes_filtered_Y,
  update_seq_classes_filtered_Y,  # Pass the full seq classes data
  n_permutations = 10000
)

# Display results
print(results_y)

# Save results
write_csv(results_y, "~/Downloads/Y_gene_palindrome_clustering_results.csv")
# read in later
results_y<- read_csv("~/Downloads/Y_gene_palindrome_clustering_results.csv")



# Plot the results 
plot_data <- results_x %>% # change to results_y for Y chromosome
  mutate(
    # Calculate expected CI before pivoting
    expected_ci_lower = expected_mean_distance - 1.96 * expected_std,
    expected_ci_upper = expected_mean_distance + 1.96 * expected_std
  ) %>%
  select(species, 
         observed_mean_distance, observed_ci_lower, observed_ci_upper,
         expected_mean_distance, expected_ci_lower, expected_ci_upper) %>%
  pivot_longer(
    cols = c(observed_mean_distance, expected_mean_distance),
    names_to = "type",
    values_to = "distance"
  ) %>%
  mutate(
    type = ifelse(type == "observed_mean_distance", "Observed", "Expected"),
    ci_lower = ifelse(type == "Observed", observed_ci_lower, expected_ci_lower),
    ci_upper = ifelse(type == "Observed", observed_ci_upper, expected_ci_upper)
  ) %>%
  mutate(species = fct_relevel(species, 
                               "PanPan", "PanTro", "HomSap", 
                               "GorGor", "PonAbe", "PonPyg", 
                               "SymSyn")) 

# Create the plot
plot_y<- ggplot(plot_data, aes(x = species, y = distance, color = type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2, 
                position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Observed" = "#E41A1C", "Expected" = "#377EB8")) +
  labs(
    x = "Species",
    y = "Distance to nearest palindrome (bp)",
    color = "Distance type",
    title = "Gene clustering near palindromes X chromosome",
    subtitle = "Points show mean with 95% confidence intervals"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  ) 

# save
ggsave("~/Downloads/gene_palindrome_clustering_X.pdf", plot_y, width = 6, height = 6, units = "in")

# compare between species
# Calculate the clustering ratio for each species
clustering_ratios <- results %>%
  mutate(
    clustering_ratio = observed_mean_distance / expected_mean_distance,
    # Also calculate fold-change
    fold_enrichment = expected_mean_distance / observed_mean_distance
  )

# The clustering ratio shows how much closer genes are to palindromes compared to random expectation (ratio < 1 means clustering).
print(clustering_ratios %>% 
        select(species, observed_mean_distance, expected_mean_distance, 
               clustering_ratio, fold_enrichment, p_value))

clustering_ratios <- results %>%
  #filter(species != "MacFas") %>%
  mutate(
    clustering_ratio = observed_mean_distance / expected_mean_distance,
    fold_enrichment = expected_mean_distance / observed_mean_distance
  )

# Summary statistics
summary_stats <- clustering_ratios %>%
  summarise(
    mean_ratio = mean(clustering_ratio),
    sd_ratio = sd(clustering_ratio),
    min_ratio = min(clustering_ratio),
    max_ratio = max(clustering_ratio),
    range_fold = max(fold_enrichment) / min(fold_enrichment)
  )

print(summary_stats)

# Test if ratios differ between species using Kruskal-Wallis
# (non-parametric test since we only have 7 species)
kruskal.test(clustering_ratio ~ species, data = clustering_ratios)

#  visualize the differences
ggplot(clustering_ratios, aes(x = reorder(species, clustering_ratio), 
                              y = clustering_ratio)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Species",
    y = "Clustering ratio (observed/expected)",
    title = "Degree of gene clustering near palindromes",
    subtitle = "Lower ratio = stronger clustering"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))









