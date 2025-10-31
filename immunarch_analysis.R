#!/usr/bin/env Rscript
################################################################################
# Immunarch TCR Repertoire Analysis
################################################################################
#
# Immunarch is a comprehensive R package for immune repertoire analysis
# Provides: diversity metrics, clonotype tracking, V/J usage, gene usage,
#          repertoire overlap, clonal space homeostasis, and more
#
# This script integrates with the Python pipeline and provides:
# 1. Comprehensive diversity analysis
# 2. Clonotype tracking and dynamics
# 3. V/J gene usage patterns
# 4. Repertoire overlap analysis
# 5. Publication-quality visualizations
#
# Installation:
#   install.packages("immunarch")
#   install.packages("tidyverse")
#
# Author: Claude Code
# Date: 2025-10-31
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  if (!require("immunarch", quietly = TRUE)) {
    message("Installing immunarch...")
    install.packages("immunarch", repos = "http://cran.us.r-project.org")
    library(immunarch)
  } else {
    library(immunarch)
  }

  if (!require("tidyverse", quietly = TRUE)) {
    message("Installing tidyverse...")
    install.packages("tidyverse", repos = "http://cran.us.r-project.org")
    library(tidyverse)
  } else {
    library(tidyverse)
  }
})

cat("================================================================================\n")
cat("IMMUNARCH TCR REPERTOIRE ANALYSIS\n")
cat("================================================================================\n\n")

cat("📚 Immunarch provides comprehensive immune repertoire analysis tools\n")
cat("   Developed by ImmunoMind for TCR/BCR sequencing data\n\n")

################################################################################
# STEP 1: Convert data to Immunarch format
################################################################################

cat("================================================================================\n")
cat("STEP 1: PREPARING DATA FOR IMMUNARCH\n")
cat("================================================================================\n\n")

# Immunarch expects data in specific format (VDJtools or AIRR)
# We'll convert our TSV files to immunarch-compatible format

convert_to_immunarch_format <- function(data_dir = "data") {
  cat("Converting TCR data to Immunarch format...\n\n")

  # Read metadata
  metadata <- read_tsv(file.path(data_dir, "metadata.tsv"),
                      show_col_types = FALSE)

  # Extract CMV status
  metadata <- metadata %>%
    mutate(cmv_status = case_when(
      grepl("CMV\\+", disease_subtype) ~ "CMV+",
      grepl("CMV-", disease_subtype) ~ "CMV-",
      TRUE ~ "Unknown"
    )) %>%
    filter(cmv_status %in% c("CMV+", "CMV-"))

  cat(sprintf("✓ Loaded metadata for %d participants\n", nrow(metadata)))
  cat(sprintf("  - CMV+: %d\n", sum(metadata$cmv_status == "CMV+")))
  cat(sprintf("  - CMV-: %d\n", sum(metadata$cmv_status == "CMV-")))

  # Find TCR files
  tcr_files <- list.files(data_dir, pattern = "part_table_.*\\.tsv",
                         full.names = TRUE)

  cat(sprintf("\n✓ Found %d TCR repertoire files\n\n", length(tcr_files)))

  # Process each file
  immdata_list <- list()

  for (tcr_file in tcr_files) {
    participant_id <- toupper(gsub("part_table_|.tsv", "", basename(tcr_file)))

    tryCatch({
      # Read TCR data
      tcr_data <- read_tsv(tcr_file, show_col_types = FALSE)

      # Filter productive sequences
      tcr_data <- tcr_data %>%
        filter(productive == "t")

      if (nrow(tcr_data) > 0) {
        # Convert to Immunarch format
        # Required columns: Clones, Proportion, CDR3.aa, V.name, J.name
        immdata_format <- tcr_data %>%
          group_by(cdr3_aa, v_call, j_call) %>%
          summarise(Clones = n(), .groups = "drop") %>%
          mutate(
            Proportion = Clones / sum(Clones),
            CDR3.aa = cdr3_aa,
            V.name = v_call,
            J.name = j_call
          ) %>%
          select(Clones, Proportion, CDR3.aa, V.name, J.name) %>%
          arrange(desc(Clones))

        immdata_list[[participant_id]] <- immdata_format

        cat(sprintf("  %s: %d productive sequences\n",
                   participant_id, nrow(tcr_data)))
      }
    }, error = function(e) {
      cat(sprintf("  ⚠️  Error loading %s: %s\n", basename(tcr_file), e$message))
    })
  }

  cat(sprintf("\n✓ Successfully converted %d repertoires\n", length(immdata_list)))

  # Create metadata for Immunarch
  immdata_meta <- metadata %>%
    filter(participant_label %in% names(immdata_list)) %>%
    select(Sample = participant_label,
           CMV = cmv_status,
           Age = age,
           Sex = sex,
           Ancestry = ancestry) %>%
    mutate(Sample = as.character(Sample))

  # Create Immunarch data structure
  immdata <- list(data = immdata_list, meta = immdata_meta)

  return(immdata)
}

# Convert data
immdata <- convert_to_immunarch_format("data")

################################################################################
# STEP 2: Diversity Analysis
################################################################################

cat("\n================================================================================\n")
cat("STEP 2: DIVERSITY ANALYSIS\n")
cat("================================================================================\n\n")

cat("📊 Calculating multiple diversity indices:\n")
cat("   • Chao1: Species richness estimator\n")
cat("   • Shannon: Entropy-based diversity\n")
cat("   • Simpson: Probability-based diversity\n")
cat("   • Inv.Simpson: Inverse Simpson index\n\n")

# Calculate diversity metrics
div_chao <- repDiversity(immdata$data, "chao1")
div_shannon <- repDiversity(immdata$data, "shannon")
div_simpson <- repDiversity(immdata$data, "inv.simp")

# Combine diversity metrics
diversity_results <- data.frame(
  Sample = names(immdata$data),
  Chao1 = div_chao,
  Shannon = div_shannon,
  InvSimpson = div_simpson
)

# Merge with metadata
if (nrow(immdata$meta) > 0) {
  diversity_results <- diversity_results %>%
    left_join(immdata$meta, by = "Sample")
}

cat("✓ Diversity metrics calculated\n\n")

# Print summary by CMV status
if ("CMV" %in% colnames(diversity_results)) {
  cat("📊 Diversity Summary by CMV Status:\n\n")

  div_summary <- diversity_results %>%
    group_by(CMV) %>%
    summarise(
      n = n(),
      Shannon_mean = mean(Shannon, na.rm = TRUE),
      Shannon_sd = sd(Shannon, na.rm = TRUE),
      InvSimpson_mean = mean(InvSimpson, na.rm = TRUE),
      InvSimpson_sd = sd(InvSimpson, na.rm = TRUE)
    )

  print(div_summary)

  # Statistical test
  if (nrow(div_summary) == 2) {
    cmv_pos <- diversity_results %>% filter(CMV == "CMV+")
    cmv_neg <- diversity_results %>% filter(CMV == "CMV-")

    if (nrow(cmv_pos) > 0 && nrow(cmv_neg) > 0) {
      cat("\n🔬 Statistical Tests (Mann-Whitney U):\n\n")

      # Shannon test
      shannon_test <- wilcox.test(cmv_pos$Shannon, cmv_neg$Shannon)
      cat(sprintf("Shannon Entropy: p = %.4f %s\n",
                 shannon_test$p.value,
                 ifelse(shannon_test$p.value < 0.05, "***", "")))

      # Simpson test
      simpson_test <- wilcox.test(cmv_pos$InvSimpson, cmv_neg$InvSimpson)
      cat(sprintf("Inverse Simpson: p = %.4f %s\n",
                 simpson_test$p.value,
                 ifelse(simpson_test$p.value < 0.05, "***", "")))
    }
  }
}

# Save diversity results
write_csv(diversity_results, "immunarch_diversity_metrics.csv")
cat("\n✓ Saved diversity metrics to: immunarch_diversity_metrics.csv\n")

################################################################################
# STEP 3: Clonotype Abundance
################################################################################

cat("\n================================================================================\n")
cat("STEP 3: CLONOTYPE ABUNDANCE ANALYSIS\n")
cat("================================================================================\n\n")

cat("📊 Analyzing clonotype frequency distributions\n\n")

# Top clonotypes
top_clones <- repClonality(immdata$data, .method = "top", .n = 10)

cat("✓ Top clonotype analysis complete\n")

# Clonality (normalized entropy)
clonality <- repClonality(immdata$data, .method = "homeo")

cat("✓ Clonality (homeostasis) calculated\n")

################################################################################
# STEP 4: Gene Usage Analysis
################################################################################

cat("\n================================================================================\n")
cat("STEP 4: V/J GENE USAGE ANALYSIS\n")
cat("================================================================================\n\n")

cat("📊 Analyzing V and J gene segment usage patterns\n\n")

# V gene usage
tryCatch({
  v_usage <- geneUsage(immdata$data, .gene = "hs.trbv", .norm = TRUE)
  cat("✓ V gene usage calculated\n")

  # Save V gene usage
  v_usage_df <- as.data.frame(v_usage)
  write_csv(v_usage_df, "immunarch_v_gene_usage.csv")
  cat("✓ Saved V gene usage to: immunarch_v_gene_usage.csv\n")
}, error = function(e) {
  cat("⚠️  V gene analysis error:", e$message, "\n")
})

# J gene usage
tryCatch({
  j_usage <- geneUsage(immdata$data, .gene = "hs.trbj", .norm = TRUE)
  cat("✓ J gene usage calculated\n")

  # Save J gene usage
  j_usage_df <- as.data.frame(j_usage)
  write_csv(j_usage_df, "immunarch_j_gene_usage.csv")
  cat("✓ Saved J gene usage to: immunarch_j_gene_usage.csv\n")
}, error = function(e) {
  cat("⚠️  J gene analysis error:", e$message, "\n")
})

################################################################################
# STEP 5: Repertoire Overlap
################################################################################

cat("\n================================================================================\n")
cat("STEP 5: REPERTOIRE OVERLAP ANALYSIS\n")
cat("================================================================================\n\n")

cat("📊 Calculating repertoire overlap between individuals\n")
cat("   Using Jaccard index (shared clonotypes / total unique)\n\n")

if (length(immdata$data) >= 2) {
  # Calculate overlap
  overlap <- repOverlap(immdata$data, .method = "jaccard", .verbose = FALSE)

  cat("✓ Repertoire overlap calculated\n")

  # Save overlap matrix
  overlap_df <- as.data.frame(as.matrix(overlap))
  write_csv(overlap_df, "immunarch_repertoire_overlap.csv")
  cat("✓ Saved overlap matrix to: immunarch_repertoire_overlap.csv\n")
} else {
  cat("⚠️  Need at least 2 repertoires for overlap analysis\n")
}

################################################################################
# STEP 6: Visualization
################################################################################

cat("\n================================================================================\n")
cat("STEP 6: GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Create comprehensive figure
pdf("immunarch_analysis.pdf", width = 16, height = 12)

# Plot 1: Diversity comparison
if ("CMV" %in% colnames(diversity_results)) {
  p1 <- ggplot(diversity_results, aes(x = CMV, y = Shannon, fill = CMV)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("CMV+" = "lightcoral", "CMV-" = "lightblue")) +
    labs(title = "Shannon Diversity by CMV Status",
         x = "CMV Status", y = "Shannon Entropy") +
    theme_minimal() +
    theme(legend.position = "none")
  print(p1)
}

# Plot 2: V gene usage heatmap (if available)
if (exists("v_usage")) {
  tryCatch({
    print(vis(v_usage))
  }, error = function(e) {
    cat("⚠️  V gene visualization error\n")
  })
}

# Plot 3: Clonotype abundance
tryCatch({
  print(vis(top_clones))
}, error = function(e) {
  cat("⚠️  Top clones visualization error\n")
})

# Plot 4: Repertoire overlap heatmap (if available)
if (exists("overlap") && length(immdata$data) >= 2) {
  tryCatch({
    print(vis(overlap))
  }, error = function(e) {
    cat("⚠️  Overlap visualization error\n")
  })
}

dev.off()

cat("\n✓ Saved visualizations to: immunarch_analysis.pdf\n")

################################################################################
# STEP 7: Summary Report
################################################################################

cat("\n================================================================================\n")
cat("FINAL SUMMARY\n")
cat("================================================================================\n\n")

cat("📊 ANALYSIS COMPLETE!\n\n")

cat("Generated outputs:\n")
cat("  • immunarch_diversity_metrics.csv - Diversity indices per sample\n")
cat("  • immunarch_v_gene_usage.csv - V gene usage patterns\n")
cat("  • immunarch_j_gene_usage.csv - J gene usage patterns\n")
cat("  • immunarch_repertoire_overlap.csv - Pairwise repertoire similarity\n")
cat("  • immunarch_analysis.pdf - Comprehensive visualizations\n\n")

cat("📚 Immunarch provides many more analyses:\n")
cat("  • Clonal tracking over time\n")
cat("  • Kmer motif analysis\n")
cat("  • Clonal space homeostasis\n")
cat("  • Advanced diversity metrics\n")
cat("  • Public clonotype detection\n\n")

cat("For full documentation, visit:\n")
cat("  https://immunarch.com/\n\n")

cat("================================================================================\n")
cat("To run with your full dataset:\n")
cat("  1. Add TCR files to data/ directory\n")
cat("  2. Run: Rscript immunarch_analysis.R\n")
cat("================================================================================\n\n")
