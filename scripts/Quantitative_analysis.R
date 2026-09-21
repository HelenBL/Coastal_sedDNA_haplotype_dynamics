# ============================================================
# QUANTITATIVE SUMMARY OF HAPLOTYPE / ASV STRUCTURE
# ============================================================
#
# Purpose
# -------
# Quantify the haplotype-network patterns used in the manuscript.
# Each ASV is treated as a PUTATIVE HAPLOTYPE, consistently with the
# haplotype-network workflow.
#
# Main output (one row per species x marker x status):
#   species
#   marker
#   status
#   n_haplotypes
#   dominant_haplotype_frequency
#   shared_haplotypes
#   site_specific_haplotypes
#   sites_per_haplotype
#
# Additional diagnostic metrics are also exported so that every value
# can be audited and interpreted conservatively.
#
# IMPORTANT INTERPRETATION
# ------------------------
# - This script quantifies ASV / putative-haplotype structure; it does
#   NOT infer population-genetic parameters.
# - A haplotype is "shared" when detected in >= 2 sites.
# - A haplotype is "site-specific" when detected in exactly 1 site.
# - dominant_haplotype_frequency is based on POSITIVE ASV-SAMPLE
#   occurrence records, not read abundance:
#
#       detections of most frequent haplotype
#       -------------------------------------
#       detections of all haplotypes
#
#   This avoids letting read counts dominate the ecological signal.
# - Because several haplotypes can occur in the same biological sample,
#   a second measure (dominant_haplotype_sample_occupancy) is exported:
#
#       samples containing dominant haplotype
#       -------------------------------------
#       positive samples for the species
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(writexl)
library(ggplot2)


# ============================================================
# 1. PATHS
# ============================================================

# Same data folder used by the timeline workflow.
data_dir <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Data"

output_dir <- file.path(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript",
  "Haplotype_structure_quantification"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

metadata_file <- file.path(
  data_dir,
  "metadata_samples.xlsx"
)

status_file <- file.path(
  data_dir,
  "Status_selected_species.xlsx"
)

file_18S <- file.path(
  data_dir,
  "Selected_species_ASVs_18S_AbRel.xlsx"
)

file_COI <- file.path(
  data_dir,
  "Selected_species_ASVs_COI_AbRel.xlsx"
)


# ============================================================
# 2. HELPERS
# ============================================================

clean_text <- function(x) {
  stringr::str_squish(as.character(x))
}

clean_id <- function(x) {
  janitor::make_clean_names(as.character(x))
}

find_first_col <- function(
    df,
    candidates,
    required = TRUE,
    label = "column"
) {
  found <- intersect(candidates, names(df))
  
  if (length(found) > 0) {
    return(found[1])
  }
  
  if (required) {
    stop(
      "Could not identify ", label, ".\n",
      "Available columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }
  
  NULL
}

normalise_status <- function(x) {
  x <- clean_text(x)
  x <- toupper(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  
  dplyr::case_when(
    x == "NAT" ~ "NAT",
    x == "NIS" ~ "NIS",
    x == "NAT/NIS" ~ "NAT/NIS",
    x == "NIS/NAT" ~ "NIS/NAT",
    x == "CRY" ~ "CRY",
    x %in% c("CRY/NAT", "NAT/CRY") ~ "CRY/NAT",
    TRUE ~ x
  )
}

normalise_kingdom <- function(x) {
  x <- clean_text(x)
  
  dplyr::case_when(
    x %in% c("Plantae", "Viridiplantae") ~ "Plantae",
    x == "Metazoa" ~ "Metazoa",
    TRUE ~ x
  )
}

collapse_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & x != ""])))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ", ")
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}


# ============================================================
# 3. METADATA: SAMPLE -> SITE
# ============================================================

metadata <- readxl::read_excel(metadata_file) %>%
  janitor::clean_names()

metadata_sample_col <- find_first_col(
  metadata,
  c(
    "sample", "sample_id", "sample_name", "sampleid",
    "sample_code", "muestra", "muestra_id", "nombre_muestra"
  ),
  required = TRUE,
  label = "sample ID column in metadata"
)

metadata_site_col <- find_first_col(
  metadata,
  c(
    "site", "site_code", "sampling_site", "location",
    "locality", "sitio"
  ),
  required = FALSE
)

metadata_year_col <- find_first_col(
  metadata,
  c(
    "year", "sampling_year", "sample_year", "year_sampling",
    "ano", "anio"
  ),
  required = FALSE
)

metadata_key <- metadata %>%
  transmute(
    sample_original = clean_text(.data[[metadata_sample_col]]),
    sample = clean_id(.data[[metadata_sample_col]]),
    site_metadata = if (!is.null(metadata_site_col)) {
      toupper(clean_text(.data[[metadata_site_col]]))
    } else {
      NA_character_
    },
    year = if (!is.null(metadata_year_col)) {
      suppressWarnings(as.integer(.data[[metadata_year_col]]))
    } else {
      NA_integer_
    }
  ) %>%
  mutate(
    site_from_sample = toupper(
      stringr::str_extract(sample_original, "^[A-Za-z]+")
    ),
    site = dplyr::coalesce(
      na_if(site_metadata, ""),
      site_from_sample
    )
  ) %>%
  select(sample, sample_original, site, year) %>%
  filter(
    !is.na(sample),
    sample != "",
    !is.na(site),
    site != ""
  ) %>%
  distinct()

metadata_conflicts <- metadata_key %>%
  distinct(sample, site) %>%
  count(sample, name = "n_sites") %>%
  filter(n_sites > 1)

if (nrow(metadata_conflicts) > 0) {
  print(metadata_conflicts)
  stop(
    "Some biological samples map to more than one site in metadata_samples.xlsx."
  )
}

metadata_key <- metadata_key %>%
  distinct(sample, .keep_all = TRUE)


# ============================================================
# 4. STATUS TABLE
# ============================================================

status_raw <- readxl::read_excel(status_file) %>%
  janitor::clean_names()

status_species_col <- find_first_col(
  status_raw,
  c("especie", "species", "selected_species"),
  required = TRUE,
  label = "species column in Status_selected_species.xlsx"
)

status_combined_col <- find_first_col(
  status_raw,
  c("status_combinado", "combined_status", "status"),
  required = TRUE,
  label = "combined status column in Status_selected_species.xlsx"
)

# NAT is intentionally INCLUDED here, unlike the invasion timeline figure.
status_levels <- c(
  "NAT",
  "NIS",
  "NAT/NIS",
  "NIS/NAT",
  "CRY",
  "CRY/NAT"
)

species_status <- status_raw %>%
  transmute(
    species = clean_text(.data[[status_species_col]]),
    status = normalise_status(.data[[status_combined_col]])
  ) %>%
  filter(
    !is.na(species),
    species != "",
    status %in% status_levels
  ) %>%
  distinct(species, status)

status_conflicts <- species_status %>%
  count(species, name = "n_status") %>%
  filter(n_status > 1)

if (nrow(status_conflicts) > 0) {
  print(
    species_status %>%
      semi_join(status_conflicts, by = "species") %>%
      arrange(species, status)
  )
  stop(
    "Some species have more than one combined status. Resolve Status_selected_species.xlsx first."
  )
}

species_status <- species_status %>%
  mutate(
    status = factor(status, levels = status_levels)
  )

cat("\nSPECIES BY STATUS INCLUDED IN THIS ANALYSIS:\n")
print(species_status %>% count(status, .drop = FALSE))


# ============================================================
# 5. READ COLLAPSED ASV ABUNDANCE TABLES
# ============================================================

read_abundance_table <- function(file) {
  message("Reading: ", basename(file))
  
  if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    df <- readxl::read_excel(file)
  } else if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    df <- readr::read_csv(file, show_col_types = FALSE)
  } else {
    stop("Unsupported file format: ", basename(file))
  }
  
  df <- df %>% janitor::clean_names()
  
  if (!"asv" %in% names(df)) {
    stop("Column ASV not found in: ", basename(file))
  }
  
  species_col <- find_first_col(
    df,
    c("selected_species", "species"),
    required = TRUE,
    label = paste0("species column in ", basename(file))
  )
  
  df %>%
    mutate(
      asv = clean_text(asv),
      species = clean_text(.data[[species_col]])
    )
}

data_18S <- read_abundance_table(file_18S)
data_COI <- read_abundance_table(file_COI)


# ============================================================
# 6. IDENTIFY BIOLOGICAL-SAMPLE COLUMNS
# ============================================================

taxonomy_cols <- c(
  "asv", "marker", "pindent", "assignment_qual", "assignmentqual",
  "tax_id", "taxid", "domain", "kingdom", "phylum", "class",
  "order", "family", "genus", "species", "selected_species",
  "status_combinado", "mediterraneo", "atlantico", "confianza",
  "justificacion", "referencia_principal", "referencia_secundaria",
  "nota_taxonomica_cautela", "assignment_level", "assignment_type",
  "best_pident", "best_qcov", "best_bitscore", "n_blast_hits",
  "n_species_hits", "blast_species_hits"
)

identify_collapsed_sample_cols <- function(df, metadata_key) {
  # First preference: sample names that match metadata after cleaning.
  matched <- intersect(metadata_key$sample, names(df))
  
  if (length(matched) > 0) {
    return(matched)
  }
  
  # Fallback: numeric columns not recognised as taxonomy/metadata.
  candidates <- setdiff(names(df), taxonomy_cols)
  
  candidates[
    vapply(df[candidates], is.numeric, logical(1))
  ]
}

sample_cols_18S <- identify_collapsed_sample_cols(data_18S, metadata_key)
sample_cols_COI <- identify_collapsed_sample_cols(data_COI, metadata_key)

if (length(sample_cols_18S) == 0) {
  stop("No biological-sample columns detected in the 18S collapsed table.")
}

if (length(sample_cols_COI) == 0) {
  stop("No biological-sample columns detected in the COI collapsed table.")
}

cat("\nDetected biological-sample columns:\n")
cat("18S: ", length(sample_cols_18S), "\n", sep = "")
cat("COI: ", length(sample_cols_COI), "\n", sep = "")


# ============================================================
# 7. TAXONOMY / STATUS BY ASV
# ============================================================

prepare_taxonomy <- function(df, marker_name) {
  if (!"kingdom" %in% names(df)) {
    stop("The ", marker_name, " table has no kingdom column.")
  }
  
  df %>%
    transmute(
      asv = clean_text(asv),
      species = clean_text(species),
      kingdom = normalise_kingdom(kingdom),
      marker = marker_name
    ) %>%
    filter(
      !is.na(asv), asv != "",
      !is.na(species), species != "",
      kingdom %in% c("Metazoa", "Plantae")
    ) %>%
    inner_join(species_status, by = "species") %>%
    distinct(marker, asv, species, kingdom, status)
}

tax_18S <- prepare_taxonomy(data_18S, "18S")
tax_COI <- prepare_taxonomy(data_COI, "COI")

taxonomy_all <- bind_rows(tax_18S, tax_COI)

taxonomy_conflicts <- taxonomy_all %>%
  distinct(marker, asv, species) %>%
  count(marker, asv, name = "n_species") %>%
  filter(n_species > 1)

if (nrow(taxonomy_conflicts) > 0) {
  print(taxonomy_conflicts)
  stop(
    "Some marker-ASV combinations map to more than one species. Resolve taxonomy first."
  )
}


# ============================================================
# 8. BUILD ASV x SAMPLE POSITIVE OCCURRENCE TABLE
# ============================================================
#
# This is the fundamental table for all metrics below.
# Each row = one positive ASV (putative haplotype) in one biological sample.
# Read abundance is retained for diagnostics but occurrence metrics use
# presence/absence so that sequencing depth does not define dominance.
# ============================================================

build_asv_occurrences <- function(
    df,
    taxonomy_df,
    sample_cols,
    metadata_key,
    marker_name
) {
  df %>%
    select(asv, all_of(sample_cols)) %>%
    inner_join(
      taxonomy_df %>%
        select(asv, species, kingdom, status),
      by = "asv"
    ) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    mutate(
      abundance = suppressWarnings(as.numeric(abundance))
    ) %>%
    filter(
      !is.na(abundance),
      abundance > 0
    ) %>%
    mutate(marker = marker_name) %>%
    left_join(
      metadata_key %>% select(sample, site, year),
      by = "sample"
    ) %>%
    filter(
      !is.na(site),
      site != ""
    ) %>%
    group_by(
      species, kingdom, status, marker,
      asv, sample, site, year
    ) %>%
    summarise(
      abundance = max(abundance, na.rm = TRUE),
      .groups = "drop"
    )
}

occ_18S <- build_asv_occurrences(
  data_18S,
  tax_18S,
  sample_cols_18S,
  metadata_key,
  "18S"
)

occ_COI <- build_asv_occurrences(
  data_COI,
  tax_COI,
  sample_cols_COI,
  metadata_key,
  "COI"
)

asv_occurrences <- bind_rows(occ_18S, occ_COI) %>%
  mutate(
    status = factor(as.character(status), levels = status_levels),
    marker = factor(marker, levels = c("COI", "18S"))
  ) %>%
  arrange(status, species, marker, asv, site, sample)

if (nrow(asv_occurrences) == 0) {
  stop("No positive ASV occurrences were recovered after filtering.")
}


# ============================================================
# 9. SPECIES-SAMPLE DENOMINATORS
# ============================================================

species_sample_summary <- asv_occurrences %>%
  distinct(species, kingdom, status, marker, sample, site) %>%
  group_by(species, kingdom, status, marker) %>%
  summarise(
    n_positive_samples_species = n_distinct(sample),
    n_sites_species = n_distinct(site),
    species_sites = collapse_values(site),
    .groups = "drop"
  )


# ============================================================
# 10. HAPLOTYPE-LEVEL METRICS
# ============================================================
#
# One row per species x marker x ASV.
# This table is the audit trail behind the species summary.
# ============================================================

haplotype_metrics <- asv_occurrences %>%
  group_by(species, kingdom, status, marker, asv) %>%
  summarise(
    n_positive_samples = n_distinct(sample),
    n_sites = n_distinct(site),
    sites = collapse_values(site),
    total_relative_abundance = sum(abundance, na.rm = TRUE),
    mean_relative_abundance_positive = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  group_by(species, kingdom, status, marker) %>%
  
  mutate(
    total_haplotype_detection_records = sum(n_positive_samples),
    
    # Frequency of each haplotype based on positive sample detections
    haplotype_detection_frequency =
      n_positive_samples / total_haplotype_detection_records,
    
    # Haplotype occurs in at least two sites
    shared_between_sites = n_sites >= 2,
    
    # Haplotype restricted to one site
    site_specific = n_sites == 1,
    
    # Rank haplotypes from most to least frequently detected
    rank_by_detections =
      dplyr::min_rank(-n_positive_samples)
  ) %>%
  
  ungroup() %>%
  
  left_join(
    species_sample_summary,
    by = c(
      "species",
      "kingdom",
      "status",
      "marker"
    )
  ) %>%
  
  mutate(
    haplotype_sample_occupancy =
      n_positive_samples / n_positive_samples_species
  ) %>%
  
  arrange(
    status,
    species,
    marker,
    dplyr::desc(n_positive_samples),
    asv
  )

# ============================================================
# 11. DOMINANT HAPLOTYPE
# ============================================================
#
# Ties are preserved in a separate diagnostic table.
# For the main summary, if two ASVs tie for first place, all tied ASV IDs
# are reported and the shared maximum frequency is used.
# ============================================================

dominant_haplotype_ties <- haplotype_metrics %>%
  group_by(species, kingdom, status, marker) %>%
  filter(n_positive_samples == max(n_positive_samples, na.rm = TRUE)) %>%
  mutate(n_tied_dominant_haplotypes = n()) %>%
  ungroup() %>%
  select(
    species, kingdom, status, marker,
    asv,
    n_tied_dominant_haplotypes,
    n_positive_samples,
    haplotype_detection_frequency,
    haplotype_sample_occupancy,
    n_sites,
    sites
  )

dominant_summary <- dominant_haplotype_ties %>%
  group_by(species, kingdom, status, marker) %>%
  summarise(
    dominant_haplotype = collapse_values(asv),
    n_tied_dominant_haplotypes = max(n_tied_dominant_haplotypes),
    dominant_haplotype_n_detections = max(n_positive_samples),
    dominant_haplotype_frequency = max(haplotype_detection_frequency),
    dominant_haplotype_sample_occupancy = max(haplotype_sample_occupancy),
    dominant_haplotype_n_sites = max(n_sites),
    dominant_haplotype_sites = collapse_values(sites),
    .groups = "drop"
  )


# ============================================================
# 12. MAIN SPECIES x MARKER x STATUS SUMMARY
# ============================================================

haplotype_summary <- haplotype_metrics %>%
  group_by(species, kingdom, status, marker) %>%
  summarise(
    n_haplotypes = n_distinct(asv),
    
    shared_haplotypes = sum(shared_between_sites),
    site_specific_haplotypes = sum(site_specific),
    
    prop_shared_haplotypes =
      shared_haplotypes / n_haplotypes,
    
    prop_site_specific_haplotypes =
      site_specific_haplotypes / n_haplotypes,
    
    # Main requested sites_per_haplotype = MEAN number of sites occupied
    # by each ASV. Median and maximum are also exported.
    sites_per_haplotype = safe_mean(n_sites),
    median_sites_per_haplotype = safe_median(n_sites),
    max_sites_per_haplotype = safe_max(n_sites),
    
    total_haplotype_detection_records = sum(n_positive_samples),
    
    .groups = "drop"
  ) %>%
  left_join(
    species_sample_summary,
    by = c("species", "kingdom", "status", "marker")
  ) %>%
  left_join(
    dominant_summary,
    by = c("species", "kingdom", "status", "marker")
  ) %>%
  mutate(
    status = factor(as.character(status), levels = status_levels)
  ) %>%
  select(
    species,
    kingdom,
    marker,
    status,
    
    # Core requested metrics
    n_haplotypes,
    dominant_haplotype_frequency,
    shared_haplotypes,
    site_specific_haplotypes,
    sites_per_haplotype,
    
    # Useful supporting metrics
    prop_shared_haplotypes,
    prop_site_specific_haplotypes,
    median_sites_per_haplotype,
    max_sites_per_haplotype,
    dominant_haplotype,
    n_tied_dominant_haplotypes,
    dominant_haplotype_n_detections,
    dominant_haplotype_sample_occupancy,
    dominant_haplotype_n_sites,
    dominant_haplotype_sites,
    n_positive_samples_species,
    n_sites_species,
    species_sites,
    total_haplotype_detection_records
  ) %>%
  arrange(status, species, marker)


# ============================================================
# 13. STATUS-LEVEL DESCRIPTIVE SUMMARY
# ============================================================
#
# Descriptive only. We deliberately do NOT run automatic significance
# tests here because species were selected for ecological relevance and
# markers have different resolutions. Formal tests should be justified
# separately if later needed.
# ============================================================

status_marker_summary <- haplotype_summary %>%
  group_by(status, marker) %>%
  summarise(
    n_species = n_distinct(species),
    
    median_n_haplotypes = median(n_haplotypes, na.rm = TRUE),
    mean_n_haplotypes = mean(n_haplotypes, na.rm = TRUE),
    
    median_dominant_haplotype_frequency =
      median(dominant_haplotype_frequency, na.rm = TRUE),
    
    mean_dominant_haplotype_frequency =
      mean(dominant_haplotype_frequency, na.rm = TRUE),
    
    median_prop_shared_haplotypes =
      median(prop_shared_haplotypes, na.rm = TRUE),
    
    mean_prop_shared_haplotypes =
      mean(prop_shared_haplotypes, na.rm = TRUE),
    
    median_prop_site_specific_haplotypes =
      median(prop_site_specific_haplotypes, na.rm = TRUE),
    
    mean_prop_site_specific_haplotypes =
      mean(prop_site_specific_haplotypes, na.rm = TRUE),
    
    median_sites_per_haplotype =
      median(sites_per_haplotype, na.rm = TRUE),
    
    mean_sites_per_haplotype =
      mean(sites_per_haplotype, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(status, marker)


# ============================================================
# 14. SITE x HAPLOTYPE PRESENCE MATRIX (LONG)
# ============================================================
#
# Especially useful later for NIS/NAT and NAT/NIS comparisons between
# Atlantic and Mediterranean native/introduced ranges.
# ============================================================

haplotype_site_presence <- asv_occurrences %>%
  group_by(species, kingdom, status, marker, asv, site) %>%
  summarise(
    n_positive_samples = n_distinct(sample),
    total_relative_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(status, species, marker, asv, site)


# ============================================================
# 15. SIMPLE SUPPORTING FIGURES
# ============================================================
#
# These are exploratory / supplementary summaries, not intended to
# replace the species-level haplotype networks.
# ============================================================

plot_status_metric <- function(data, yvar, ylab, filename) {
  p <- ggplot(
    data,
    aes(
      x = status,
      y = .data[[yvar]]
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65
    ) +
    geom_jitter(
      aes(shape = marker),
      width = 0.12,
      height = 0,
      size = 2.8,
      alpha = 0.8
    ) +
    scale_shape_manual(
      values = c("COI" = 16, "18S" = 17),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = ylab,
      shape = "Marker"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "top"
    )
  
  ggsave(
    file.path(output_dir, paste0(filename, ".svg")),
    p,
    width = 8,
    height = 5
  )
  
  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    p,
    width = 8,
    height = 5
  )
  
  p
}

p_hap_richness <- plot_status_metric(
  haplotype_summary,
  "n_haplotypes",
  "Number of ASV haplotypes",
  "Status_haplotype_richness"
)

p_dominance <- plot_status_metric(
  haplotype_summary,
  "dominant_haplotype_frequency",
  "Dominant haplotype frequency",
  "Status_dominant_haplotype_frequency"
)

p_shared <- plot_status_metric(
  haplotype_summary,
  "prop_shared_haplotypes",
  "Proportion of haplotypes shared among sites",
  "Status_shared_haplotype_proportion"
)

p_specific <- plot_status_metric(
  haplotype_summary,
  "prop_site_specific_haplotypes",
  "Proportion of site-specific haplotypes",
  "Status_site_specific_haplotype_proportion"
)

p_sites_per_hap <- plot_status_metric(
  haplotype_summary,
  "sites_per_haplotype",
  "Mean sites per haplotype",
  "Status_sites_per_haplotype"
)


# ============================================================
# 16. EXPORT CSV FILES
# ============================================================

readr::write_csv(
  haplotype_summary,
  file.path(output_dir, "Haplotype_structure_species_marker_summary.csv")
)

readr::write_csv(
  haplotype_metrics,
  file.path(output_dir, "Haplotype_structure_haplotype_level.csv")
)

readr::write_csv(
  haplotype_site_presence,
  file.path(output_dir, "Haplotype_structure_site_presence.csv")
)

readr::write_csv(
  status_marker_summary,
  file.path(output_dir, "Haplotype_structure_status_marker_summary.csv")
)


# ============================================================
# 17. EXPORT EXCEL WORKBOOK
# ============================================================

output_xlsx <- file.path(
  output_dir,
  "Haplotype_structure_quantification.xlsx"
)

writexl::write_xlsx(
  list(
    Species_marker_summary = haplotype_summary,
    Haplotype_level = haplotype_metrics,
    Haplotype_site_presence = haplotype_site_presence,
    Status_marker_summary = status_marker_summary,
    Dominant_haplotype_ties = dominant_haplotype_ties,
    Species_sample_summary = species_sample_summary,
    ASV_positive_occurrences = asv_occurrences,
    Taxonomy_status_ASVs = taxonomy_all,
    Status_species = species_status
  ),
  path = output_xlsx
)


# ============================================================
# 18. CONSOLE SUMMARY / DIAGNOSTICS
# ============================================================

cat("\n============================================================\n")
cat("HAPLOTYPE STRUCTURE QUANTIFICATION COMPLETE\n")
cat("============================================================\n\n")

cat("Main species x marker table:\n")
print(
  haplotype_summary %>%
    select(
      species,
      marker,
      status,
      n_haplotypes,
      dominant_haplotype_frequency,
      shared_haplotypes,
      site_specific_haplotypes,
      sites_per_haplotype
    ),
  n = Inf
)

cat("\nStatus x marker descriptive summary:\n")
print(status_marker_summary, n = Inf)

cat("\nDominant-haplotype ties:\n")
print(
  dominant_haplotype_ties %>%
    filter(n_tied_dominant_haplotypes > 1),
  n = Inf
)

cat("\nFiles written to:\n", output_dir, "\n", sep = "")
cat("Excel workbook:\n", output_xlsx, "\n", sep = "")

cat("\nMetric definitions:\n")
cat("- n_haplotypes: number of distinct ASVs assigned to species x marker and detected in >=1 biological sample.\n")
cat("- dominant_haplotype_frequency: detections of most frequent ASV / detections of all ASVs for that species x marker.\n")
cat("- shared_haplotypes: number of ASVs detected in >=2 sites.\n")
cat("- site_specific_haplotypes: number of ASVs detected in exactly 1 site.\n")
cat("- sites_per_haplotype: mean number of occupied sites across ASVs.\n")
cat("- prop_shared_haplotypes and prop_site_specific_haplotypes are recommended for comparisons among taxa with different haplotype richness.\n")
cat("\nNOTE: Compare COI and 18S separately when interpreting status patterns.\n")

# ============================================================
# END
# ============================================================
