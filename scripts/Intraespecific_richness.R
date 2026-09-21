# ======================================================================
# ======================================================================
#
# FINAL REPRODUCIBLE TEMPORAL HAPLOTYPE ANALYSIS
#
# 18S = MAIN ANALYSIS
#       ANALYTICAL / GRAPHICAL PIPELINE PRESERVED
#
# COI = SUPPLEMENTARY
#       SAME ANALYSES
#       GRAPHICAL SAFEGUARDS FOR SPARSE DATA
#
# ======================================================================
# ======================================================================


# ======================================================================
# 0. CLEAN WORKSPACE
# ======================================================================

rm(list = ls())
gc()


# ======================================================================
# 1. PACKAGES
# ======================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "writexl",
  "scales",
  "svglite"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(writexl)
library(scales)
library(svglite)


# ======================================================================
# 2. PATHS
# ======================================================================

data_dir <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Haplotypes_manuscript/Data"
)

metadata_file <- file.path(
  data_dir,
  "metadata_samples.xlsx"
)

status_file <- file.path(
  data_dir,
  "Status_selected_species_OK.xlsx"
)

file_18S <- file.path(
  data_dir,
  "Selected_species_ASVs_18S_AbRel.xlsx"
)

file_COI <- file.path(
  data_dir,
  "Selected_species_ASVs_COI_AbRel.xlsx"
)


# ======================================================================
# 3. OUTPUT
# ======================================================================

output_dir <- file.path(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Figures/Intraspecific_ASV_richness_NAT_NIS/FINAL_HAPLOTYPE_TEMPORAL_ANALYSIS"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# 4. SETTINGS
# ======================================================================

BIN_WIDTH <- 50
MIN_N_STATUS <- 3
MIN_N_PAIRED <- 3

status_colors <- c(
  "NAT" = "#8FC793",
  "NIS" = "#D81B60"
)


# ======================================================================
# 5. CHECK INPUT FILES
# ======================================================================

input_files <- c(
  metadata_file,
  status_file,
  file_18S,
  file_COI
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "\nMissing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ======================================================================
# 6. READ DATA
# ======================================================================

cat(
  "\n============================================================\n",
  "READING INPUT DATA\n",
  "============================================================\n"
)

metadata_raw <- read_excel(metadata_file)
status_raw   <- read_excel(status_file)
asv18_raw    <- read_excel(file_18S)
asvCOI_raw   <- read_excel(file_COI)

cat(
  "\nMetadata:", nrow(metadata_raw), "×", ncol(metadata_raw), "\n",
  "Status:", nrow(status_raw), "×", ncol(status_raw), "\n",
  "18S:", nrow(asv18_raw), "×", ncol(asv18_raw), "\n",
  "COI:", nrow(asvCOI_raw), "×", ncol(asvCOI_raw), "\n"
)


# ======================================================================
# 7. SIMPLE NAME CLEANING
# ======================================================================

clean_names_simple <- function(df) {
  
  names(df) <- names(df) %>%
    str_trim() %>%
    str_replace_all("\\s+", "_")
  
  df
}

metadata_raw <- clean_names_simple(metadata_raw)
status_raw   <- clean_names_simple(status_raw)
asv18_raw    <- clean_names_simple(asv18_raw)
asvCOI_raw   <- clean_names_simple(asvCOI_raw)


# ======================================================================
# 8. COLUMN FINDER
# ======================================================================

find_column <- function(
    df,
    candidates,
    label,
    required = TRUE
) {
  
  nm <- names(df)
  nm_lower <- tolower(nm)
  candidates_lower <- tolower(candidates)
  
  idx <- match(
    candidates_lower,
    nm_lower
  )
  
  idx <- idx[!is.na(idx)]
  
  if (length(idx) > 0) {
    
    found <- nm[idx[1]]
    
    cat(
      label,
      "column:",
      found,
      "\n"
    )
    
    return(found)
  }
  
  if (required) {
    stop(
      "\nCould not identify ",
      label,
      " column.\nAvailable columns:\n",
      paste(nm, collapse = ", ")
    )
  }
  
  NULL
}


# ======================================================================
# 9. METADATA COLUMNS
# ======================================================================

metadata_sample_col <- find_column(
  metadata_raw,
  c(
    "sample",
    "sample_id",
    "sampleid",
    "sample_name",
    "samplename"
  ),
  "metadata sample ID"
)

metadata_site_col <- find_column(
  metadata_raw,
  c(
    "site",
    "site_code",
    "sitecode",
    "location"
  ),
  "metadata site"
)

metadata_year_col <- find_column(
  metadata_raw,
  c(
    "year",
    "age",
    "calendar_year",
    "cal_year"
  ),
  "metadata year"
)


# ======================================================================
# 10. STANDARDISE METADATA
# ======================================================================

metadata <- metadata_raw %>%
  transmute(
    sample = as.character(.data[[metadata_sample_col]]),
    site   = as.character(.data[[metadata_site_col]]),
    year   = suppressWarnings(
      as.numeric(.data[[metadata_year_col]])
    )
  ) %>%
  mutate(
    sample = str_trim(sample),
    site   = str_trim(site),
    site   = recode(site, "COL" = "CLR")
  ) %>%
  filter(
    !is.na(sample),
    sample != "",
    !is.na(site),
    site != "",
    !is.na(year)
  ) %>%
  distinct(
    sample,
    .keep_all = TRUE
  ) %>%
  mutate(
    period_start =
      floor(year / BIN_WIDTH) * BIN_WIDTH,
    
    period_end =
      period_start + BIN_WIDTH - 1,
    
    period_label =
      paste0(
        period_start,
        "–",
        period_end
      )
  )


# ======================================================================
# 11. METADATA SUMMARY
# ======================================================================

metadata_summary <- metadata %>%
  group_by(site) %>%
  summarise(
    n_samples = n(),
    min_year = min(year),
    max_year = max(year),
    n_periods = n_distinct(period_start),
    .groups = "drop"
  )

cat("\nMETADATA SUMMARY:\n")
print(metadata_summary, n = Inf)


# ======================================================================
# 12. STATUS COLUMNS
# ======================================================================

status_species_col <- find_column(
  status_raw,
  c(
    "species",
    "selected_species"
  ),
  "status species"
)

status_med_col <- find_column(
  status_raw,
  c(
    "mediterraneo",
    "mediterranean"
  ),
  "Mediterranean status"
)

status_atl_col <- find_column(
  status_raw,
  c(
    "atlantico",
    "atlantic"
  ),
  "Atlantic status"
)


# ======================================================================
# 13. STATUS TABLE
# ======================================================================

status_table <- status_raw %>%
  transmute(
    species =
      as.character(
        .data[[status_species_col]]
      ),
    
    mediterraneo =
      as.character(
        .data[[status_med_col]]
      ),
    
    atlantico =
      as.character(
        .data[[status_atl_col]]
      )
  ) %>%
  mutate(
    species =
      str_trim(species),
    
    mediterraneo =
      str_trim(
        toupper(mediterraneo)
      ),
    
    atlantico =
      str_trim(
        toupper(atlantico)
      )
  ) %>%
  filter(
    !is.na(species),
    species != ""
  ) %>%
  distinct(
    species,
    .keep_all = TRUE
  )


# ======================================================================
# 14. REGIONAL STATUS
#
# DEE = MEDITERRANEAN
# ALL OTHER SITES = ATLANTIC
# ======================================================================

assign_regional_status <- function(
    species_vector,
    site_vector
) {
  
  tmp <- tibble(
    species = species_vector,
    site = site_vector
  ) %>%
    left_join(
      status_table,
      by = "species"
    ) %>%
    mutate(
      regional_status =
        if_else(
          site == "DEE",
          mediterraneo,
          atlantico
        ),
      
      regional_status =
        str_trim(
          toupper(regional_status)
        )
    )
  
  tmp$regional_status
}


# ======================================================================
# 15. ASV / SPECIES COLUMN FINDERS
# ======================================================================

find_asv_column <- function(df) {
  
  find_column(
    df,
    c(
      "ASV",
      "asv",
      "ASV_ID",
      "asv_id",
      "sequence_id"
    ),
    "ASV"
  )
}

find_species_column <- function(df) {
  
  find_column(
    df,
    c(
      "species",
      "selected_species"
    ),
    "species"
  )
}


# ======================================================================
# 16. PREPARE ABUNDANCE TABLE
# ======================================================================

prepare_asv_long <- function(
    df,
    marker_name
) {
  
  cat(
    "\n============================================================\n",
    "PREPARING ", marker_name,
    "\n============================================================\n"
  )
  
  asv_col <- find_asv_column(df)
  species_col <- find_species_column(df)
  
  sample_cols <- intersect(
    names(df),
    metadata$sample
  )
  
  cat(
    marker_name,
    "sample columns matched to metadata:",
    length(sample_cols),
    "\n"
  )
  
  if (length(sample_cols) == 0) {
    stop(
      "\nNo sample columns for ",
      marker_name,
      " matched metadata IDs."
    )
  }
  
  long <- df %>%
    transmute(
      asv =
        as.character(
          .data[[asv_col]]
        ),
      
      species =
        as.character(
          .data[[species_col]]
        ),
      
      across(
        all_of(sample_cols)
      )
    ) %>%
    mutate(
      asv = str_trim(asv),
      species = str_trim(species)
    ) %>%
    filter(
      !is.na(asv),
      asv != "",
      !is.na(species),
      species != ""
    ) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    mutate(
      abundance =
        suppressWarnings(
          as.numeric(abundance)
        )
    ) %>%
    filter(
      !is.na(abundance),
      abundance > 0
    ) %>%
    left_join(
      metadata,
      by = "sample"
    ) %>%
    filter(
      !is.na(site),
      !is.na(year),
      !is.na(period_start)
    ) %>%
    mutate(
      marker = marker_name,
      
      regional_status =
        assign_regional_status(
          species,
          site
        )
    )
  
  cat(
    marker_name,
    "positive observations:",
    nrow(long),
    "\n"
  )
  
  cat(
    marker_name,
    "unique ASVs:",
    n_distinct(long$asv),
    "\n"
  )
  
  cat(
    marker_name,
    "unique species:",
    n_distinct(long$species),
    "\n"
  )
  
  long
}


# ======================================================================
# 17. BUILD 18S + COI
# ======================================================================

asv18_long <- prepare_asv_long(
  asv18_raw,
  "18S"
)

asvCOI_long <- prepare_asv_long(
  asvCOI_raw,
  "COI"
)

asv_samples <- bind_rows(
  asv18_long,
  asvCOI_long
)


# ======================================================================
# 18. SAMPLE EFFORT
# ======================================================================

sample_effort <- metadata %>%
  group_by(
    site,
    period_start
  ) %>%
  summarise(
    n_samples_available =
      n_distinct(sample),
    .groups = "drop"
  )


# ======================================================================
# 19. BUILD RICHNESS
# ======================================================================

richness_all <- asv_samples %>%
  filter(
    regional_status %in%
      c(
        "NAT",
        "NIS"
      )
  ) %>%
  group_by(
    marker,
    site,
    species,
    regional_status,
    period_start
  ) %>%
  summarise(
    n_asvs =
      n_distinct(asv),
    .groups = "drop"
  ) %>%
  left_join(
    sample_effort,
    by =
      c(
        "site",
        "period_start"
      )
  ) %>%
  mutate(
    period_end =
      period_start +
      BIN_WIDTH -
      1,
    
    period_label =
      paste0(
        period_start,
        "–",
        period_end
      )
  )


# ======================================================================
# 20. DIAGNOSTIC
# ======================================================================

input_summary <- richness_all %>%
  group_by(marker) %>%
  summarise(
    n_sites =
      n_distinct(site),
    
    n_species =
      n_distinct(species),
    
    n_NAT_species =
      n_distinct(
        species[
          regional_status == "NAT"
        ]
      ),
    
    n_NIS_species =
      n_distinct(
        species[
          regional_status == "NIS"
        ]
      ),
    
    n_species_period_rows =
      n(),
    
    .groups = "drop"
  )

cat(
  "\n============================================================\n",
  "INPUT SUMMARY\n",
  "============================================================\n"
)

print(input_summary)


# ======================================================================
# ======================================================================
# HELPER FUNCTIONS
# ======================================================================
# ======================================================================

safe_mean <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  mean(x)
}


safe_median <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  median(x)
}


safe_sd <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) < 2) {
    return(NA_real_)
  }
  
  sd(x)
}


significance_label <- function(p) {
  
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}


# ======================================================================
# 21. ALL TEMPORAL PAIRS
# ======================================================================

make_all_pairs <- function(periods) {
  
  periods <- sort(
    unique(
      periods[
        !is.na(periods)
      ]
    )
  )
  
  if (length(periods) < 2) {
    
    return(
      tibble(
        period_1 = numeric(),
        period_2 = numeric()
      )
    )
  }
  
  cmb <- combn(
    periods,
    2
  )
  
  tibble(
    period_1 = cmb[1, ],
    period_2 = cmb[2, ]
  )
}


# ======================================================================
# 22. Q1 — NAT vs NIS WITHIN BIN
# ======================================================================

compare_nat_nis_bin <- function(df) {
  
  nat <- df %>%
    filter(
      regional_status == "NAT"
    ) %>%
    pull(n_asvs)
  
  nis <- df %>%
    filter(
      regional_status == "NIS"
    ) %>%
    pull(n_asvs)
  
  n_nat <- length(nat)
  n_nis <- length(nis)
  
  result <- tibble(
    n_NAT = n_nat,
    n_NIS = n_nis,
    
    mean_NAT =
      safe_mean(nat),
    
    mean_NIS =
      safe_mean(nis),
    
    median_NAT =
      safe_median(nat),
    
    median_NIS =
      safe_median(nis),
    
    delta_NIS_minus_NAT =
      safe_mean(nis) -
      safe_mean(nat),
    
    W = NA_real_,
    p_raw = NA_real_
  )
  
  if (
    n_nat < MIN_N_STATUS ||
    n_nis < MIN_N_STATUS
  ) {
    return(result)
  }
  
  wt <- suppressWarnings(
    wilcox.test(
      nis,
      nat,
      paired = FALSE,
      exact = FALSE,
      alternative = "two.sided"
    )
  )
  
  result$W <-
    unname(
      wt$statistic
    )
  
  result$p_raw <-
    wt$p.value
  
  result
}


# ======================================================================
# 23. Q2 — TEMPORAL RICHNESS PAIR
# ======================================================================

compare_richness_pair <- function(
    df,
    p1,
    p2
) {
  
  d1 <- df %>%
    filter(
      period_start == p1
    ) %>%
    select(
      species,
      richness_1 = n_asvs
    )
  
  d2 <- df %>%
    filter(
      period_start == p2
    ) %>%
    select(
      species,
      richness_2 = n_asvs
    )
  
  paired <- inner_join(
    d1,
    d2,
    by = "species"
  ) %>%
    mutate(
      delta =
        richness_2 -
        richness_1
    )
  
  n_shared <- nrow(paired)
  
  if (n_shared == 0) {
    
    return(
      tibble(
        period_1 = p1,
        period_2 = p2,
        lag_years = p2 - p1,
        
        n_shared_species = 0,
        
        mean_1 = NA_real_,
        mean_2 = NA_real_,
        
        mean_change = NA_real_,
        median_change = NA_real_,
        sd_change = NA_real_,
        
        n_increase = 0,
        n_decrease = 0,
        n_equal = 0,
        
        proportion_increase = NA_real_,
        proportion_decrease = NA_real_,
        proportion_equal = NA_real_,
        
        V = NA_real_,
        p_raw = NA_real_
      )
    )
  }
  
  result <- tibble(
    period_1 = p1,
    period_2 = p2,
    lag_years = p2 - p1,
    
    n_shared_species = n_shared,
    
    mean_1 =
      mean(paired$richness_1),
    
    mean_2 =
      mean(paired$richness_2),
    
    mean_change =
      mean(paired$delta),
    
    median_change =
      median(paired$delta),
    
    sd_change =
      safe_sd(paired$delta),
    
    n_increase =
      sum(paired$delta > 0),
    
    n_decrease =
      sum(paired$delta < 0),
    
    n_equal =
      sum(paired$delta == 0),
    
    proportion_increase =
      mean(paired$delta > 0),
    
    proportion_decrease =
      mean(paired$delta < 0),
    
    proportion_equal =
      mean(paired$delta == 0),
    
    V = NA_real_,
    p_raw = NA_real_
  )
  
  if (
    n_shared <
    MIN_N_PAIRED
  ) {
    return(result)
  }
  
  if (
    all(
      paired$delta == 0
    )
  ) {
    
    result$V <- 0
    result$p_raw <- 1
    
    return(result)
  }
  
  wt <- suppressWarnings(
    wilcox.test(
      paired$richness_2,
      paired$richness_1,
      paired = TRUE,
      exact = FALSE,
      alternative = "two.sided"
    )
  )
  
  result$V <-
    unname(
      wt$statistic
    )
  
  result$p_raw <-
    wt$p.value
  
  result
}


run_all_richness_pairs <- function(df) {
  
  pairs <- make_all_pairs(
    df$period_start
  )
  
  if (nrow(pairs) == 0) {
    return(tibble())
  }
  
  map2_dfr(
    pairs$period_1,
    pairs$period_2,
    ~ compare_richness_pair(
      df,
      .x,
      .y
    )
  )
}


# ======================================================================
# 24. ASV COMPOSITION
# ======================================================================

compare_asv_sets <- function(
    df,
    p1,
    p2
) {
  
  d1 <- df %>%
    filter(
      period_start == p1
    ) %>%
    select(
      species,
      asvs_1 = asvs
    )
  
  d2 <- df %>%
    filter(
      period_start == p2
    ) %>%
    select(
      species,
      asvs_2 = asvs
    )
  
  paired <- inner_join(
    d1,
    d2,
    by = "species"
  )
  
  if (nrow(paired) == 0) {
    return(tibble())
  }
  
  paired %>%
    mutate(
      n_asvs_1 =
        map_int(
          asvs_1,
          length
        ),
      
      n_asvs_2 =
        map_int(
          asvs_2,
          length
        ),
      
      shared_asvs =
        map2(
          asvs_1,
          asvs_2,
          intersect
        ),
      
      lost_asvs =
        map2(
          asvs_1,
          asvs_2,
          setdiff
        ),
      
      gained_asvs =
        map2(
          asvs_2,
          asvs_1,
          setdiff
        ),
      
      n_shared =
        map_int(
          shared_asvs,
          length
        ),
      
      n_lost =
        map_int(
          lost_asvs,
          length
        ),
      
      n_gained =
        map_int(
          gained_asvs,
          length
        ),
      
      richness_change =
        n_asvs_2 -
        n_asvs_1,
      
      union_size =
        map2_int(
          asvs_1,
          asvs_2,
          ~ length(
            union(.x, .y)
          )
        ),
      
      jaccard_similarity =
        if_else(
          union_size > 0,
          n_shared / union_size,
          NA_real_
        ),
      
      jaccard_turnover =
        1 -
        jaccard_similarity,
      
      proportion_persisting =
        if_else(
          n_asvs_1 > 0,
          n_shared / n_asvs_1,
          NA_real_
        ),
      
      proportion_not_redetected =
        if_else(
          n_asvs_1 > 0,
          n_lost / n_asvs_1,
          NA_real_
        ),
      
      proportion_new =
        if_else(
          n_asvs_2 > 0,
          n_gained / n_asvs_2,
          NA_real_
        ),
      
      period_1 = p1,
      period_2 = p2,
      lag_years = p2 - p1
    ) %>%
    select(
      species,
      
      period_1,
      period_2,
      lag_years,
      
      n_asvs_1,
      n_asvs_2,
      richness_change,
      
      n_shared,
      n_lost,
      n_gained,
      
      proportion_persisting,
      proportion_not_redetected,
      proportion_new,
      
      jaccard_similarity,
      jaccard_turnover,
      
      shared_asvs,
      lost_asvs,
      gained_asvs
    )
}


run_all_asv_pairs <- function(df) {
  
  pairs <- make_all_pairs(
    df$period_start
  )
  
  if (nrow(pairs) == 0) {
    return(tibble())
  }
  
  map2_dfr(
    pairs$period_1,
    pairs$period_2,
    ~ compare_asv_sets(
      df,
      .x,
      .y
    )
  )
}


# ======================================================================
# 25. NAT vs NIS TURNOVER
# ======================================================================

compare_turnover_status <- function(
    df,
    response
) {
  
  nat <- df %>%
    filter(
      regional_status == "NAT"
    ) %>%
    pull(
      all_of(response)
    )
  
  nis <- df %>%
    filter(
      regional_status == "NIS"
    ) %>%
    pull(
      all_of(response)
    )
  
  nat <- nat[!is.na(nat)]
  nis <- nis[!is.na(nis)]
  
  n_nat <- length(nat)
  n_nis <- length(nis)
  
  result <- tibble(
    metric = response,
    
    n_NAT = n_nat,
    n_NIS = n_nis,
    
    mean_NAT =
      safe_mean(nat),
    
    mean_NIS =
      safe_mean(nis),
    
    median_NAT =
      safe_median(nat),
    
    median_NIS =
      safe_median(nis),
    
    delta_NIS_minus_NAT =
      safe_mean(nis) -
      safe_mean(nat),
    
    W = NA_real_,
    p_raw = NA_real_
  )
  
  if (
    n_nat < MIN_N_STATUS ||
    n_nis < MIN_N_STATUS
  ) {
    return(result)
  }
  
  wt <- suppressWarnings(
    wilcox.test(
      nis,
      nat,
      paired = FALSE,
      exact = FALSE,
      alternative = "two.sided"
    )
  )
  
  result$W <-
    unname(
      wt$statistic
    )
  
  result$p_raw <-
    wt$p.value
  
  result
}


# ======================================================================
# ======================================================================
#
# 26. MAIN ANALYSIS FUNCTION
#
# Statistical analyses are identical for 18S and COI.
#
# ONLY graphical safeguards differ for COI.
#
# ======================================================================
# ======================================================================

analyse_marker <- function(
    marker_name,
    richness_data,
    asv_data
) {
  
  cat(
    "\n\n============================================================\n",
    "ANALYSING ", marker_name,
    "\n============================================================\n"
  )
  
  marker_dir <- file.path(
    output_dir,
    marker_name
  )
  
  dir.create(
    marker_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ====================================================================
  # RICHNESS DATA
  # ====================================================================
  
  rich <- richness_data %>%
    filter(
      marker == marker_name,
      regional_status %in% c("NAT", "NIS"),
      !is.na(site),
      !is.na(species),
      !is.na(period_start),
      !is.na(n_asvs)
    )
  
  
  # ====================================================================
  # DUPLICATE CHECK
  # ====================================================================
  
  duplicate_check <- rich %>%
    count(
      site,
      species,
      regional_status,
      period_start
    ) %>%
    filter(
      n > 1
    )
  
  if (nrow(duplicate_check) > 0) {
    stop(
      "\nDuplicate species × site × period rows found for ",
      marker_name
    )
  }
  
  
  # ====================================================================
  # Q1
  # ====================================================================
  
  within_bin <- rich %>%
    group_by(
      site,
      period_start,
      period_label
    ) %>%
    group_modify(
      ~ compare_nat_nis_bin(.x)
    ) %>%
    ungroup() %>%
    mutate(
      p_FDR =
        p.adjust(
          p_raw,
          method = "BH"
        ),
      
      direction =
        case_when(
          delta_NIS_minus_NAT > 0 ~ "NIS > NAT",
          delta_NIS_minus_NAT < 0 ~ "NAT > NIS",
          delta_NIS_minus_NAT == 0 ~ "Equal",
          TRUE ~ NA_character_
        ),
      
      sig_raw =
        significance_label(p_raw),
      
      sig_FDR =
        significance_label(p_FDR),
      
      raw_significant =
        !is.na(p_raw) &
        p_raw < 0.05,
      
      FDR_significant =
        !is.na(p_FDR) &
        p_FDR < 0.05
    )
  
  
  # ====================================================================
  # Q2
  # ====================================================================
  
  all_timelags <- rich %>%
    group_by(
      site,
      regional_status
    ) %>%
    group_modify(
      ~ run_all_richness_pairs(.x)
    ) %>%
    ungroup() %>%
    group_by(
      regional_status
    ) %>%
    mutate(
      p_FDR =
        p.adjust(
          p_raw,
          method = "BH"
        )
    ) %>%
    ungroup() %>%
    mutate(
      transition =
        paste0(
          period_1,
          "→",
          period_2
        ),
      
      direction =
        case_when(
          mean_change > 0 ~ "Increase",
          mean_change < 0 ~ "Decrease",
          mean_change == 0 ~ "No change",
          TRUE ~ NA_character_
        ),
      
      sig_raw =
        significance_label(p_raw),
      
      sig_FDR =
        significance_label(p_FDR),
      
      raw_significant =
        !is.na(p_raw) &
        p_raw < 0.05,
      
      FDR_significant =
        !is.na(p_FDR) &
        p_FDR < 0.05
    )
  
  
  # ====================================================================
  # Q3 — CORRECTED LAG SUMMARY
  #
  # IMPORTANT:
  # We first rename the comparison-level variable `mean_change`
  # to `comparison_mean_change`.
  #
  # This prevents dplyr::summarise() from accidentally using a newly
  # created summary variable instead of the original vector.
  # ====================================================================
  
  lag_summary <- all_timelags %>%
    
    filter(
      !is.na(mean_change)
    ) %>%
    
    transmute(
      regional_status,
      lag_years,
      comparison_mean_change = mean_change
    ) %>%
    
    group_by(
      regional_status,
      lag_years
    ) %>%
    
    summarise(
      
      n_comparisons = n(),
      
      mean_delta =
        mean(
          comparison_mean_change,
          na.rm = TRUE
        ),
      
      median_delta =
        median(
          comparison_mean_change,
          na.rm = TRUE
        ),
      
      sd_delta =
        safe_sd(
          comparison_mean_change
        ),
      
      proportion_positive =
        mean(
          comparison_mean_change > 0,
          na.rm = TRUE
        ),
      
      proportion_negative =
        mean(
          comparison_mean_change < 0,
          na.rm = TRUE
        ),
      
      proportion_zero =
        mean(
          comparison_mean_change == 0,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  # ====================================================================
  # Q4 — DESCRIPTIVE DIFFERENTIAL CHANGE
  # ====================================================================
  
  differential_change <- all_timelags %>%
    select(
      site,
      regional_status,
      period_1,
      period_2,
      lag_years,
      n_shared_species,
      mean_change,
      median_change,
      n_increase,
      n_decrease,
      n_equal
    ) %>%
    pivot_wider(
      names_from =
        regional_status,
      
      values_from =
        c(
          n_shared_species,
          mean_change,
          median_change,
          n_increase,
          n_decrease,
          n_equal
        )
    ) %>%
    mutate(
      differential_change =
        mean_change_NIS -
        mean_change_NAT,
      
      sufficient_both =
        !is.na(n_shared_species_NAT) &
        !is.na(n_shared_species_NIS) &
        n_shared_species_NAT >= MIN_N_PAIRED &
        n_shared_species_NIS >= MIN_N_PAIRED,
      
      differential_direction =
        case_when(
          differential_change > 0 ~ "NIS more positive",
          differential_change < 0 ~ "NAT more positive",
          differential_change == 0 ~ "Equal",
          TRUE ~ NA_character_
        )
    )
  
  
  # ====================================================================
  # Q5 — ASV SETS
  # ====================================================================
  
  asv_marker <- asv_data %>%
    filter(
      marker == marker_name,
      regional_status %in% c("NAT", "NIS")
    ) %>%
    distinct(
      site,
      species,
      regional_status,
      period_start,
      asv
    )
  
  asv_sets <- asv_marker %>%
    group_by(
      site,
      species,
      regional_status,
      period_start
    ) %>%
    summarise(
      asvs =
        list(
          sort(
            unique(asv)
          )
        ),
      
      n_asvs =
        n_distinct(asv),
      
      .groups = "drop"
    )
  
  
  # ====================================================================
  # ALL ASV TEMPORAL PAIRS
  # ====================================================================
  
  asv_turnover_species <- asv_sets %>%
    group_by(
      site,
      regional_status
    ) %>%
    group_modify(
      ~ run_all_asv_pairs(.x)
    ) %>%
    ungroup()
  
  
  # ====================================================================
  # Q5 SUMMARY
  # ====================================================================
  
  asv_turnover_summary <- asv_turnover_species %>%
    group_by(
      site,
      regional_status,
      period_1,
      period_2,
      lag_years
    ) %>%
    summarise(
      n_species = n(),
      
      mean_shared =
        mean(
          n_shared,
          na.rm = TRUE
        ),
      
      mean_not_redetected =
        mean(
          n_lost,
          na.rm = TRUE
        ),
      
      mean_new =
        mean(
          n_gained,
          na.rm = TRUE
        ),
      
      mean_proportion_persisting =
        mean(
          proportion_persisting,
          na.rm = TRUE
        ),
      
      mean_proportion_not_redetected =
        mean(
          proportion_not_redetected,
          na.rm = TRUE
        ),
      
      mean_proportion_new =
        mean(
          proportion_new,
          na.rm = TRUE
        ),
      
      mean_jaccard_similarity =
        mean(
          jaccard_similarity,
          na.rm = TRUE
        ),
      
      mean_jaccard_turnover =
        mean(
          jaccard_turnover,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  # ====================================================================
  # Q6 — TURNOVER TESTS
  # ====================================================================
  
  turnover_tests <- asv_turnover_species %>%
    group_by(
      site,
      period_1,
      period_2,
      lag_years
    ) %>%
    group_modify(
      ~ bind_rows(
        compare_turnover_status(
          .x,
          "proportion_persisting"
        ),
        
        compare_turnover_status(
          .x,
          "jaccard_turnover"
        )
      )
    ) %>%
    ungroup() %>%
    group_by(metric) %>%
    mutate(
      p_FDR =
        p.adjust(
          p_raw,
          method = "BH"
        )
    ) %>%
    ungroup() %>%
    mutate(
      sig_raw =
        significance_label(p_raw),
      
      sig_FDR =
        significance_label(p_FDR),
      
      raw_significant =
        !is.na(p_raw) &
        p_raw < 0.05,
      
      FDR_significant =
        !is.na(p_FDR) &
        p_FDR < 0.05
    )
  
  
  # ====================================================================
  # GLOBAL TURNOVER
  # ====================================================================
  
  turnover_global <- asv_turnover_species %>%
    group_by(
      regional_status
    ) %>%
    summarise(
      n_species_comparisons =
        n(),
      
      mean_persistence =
        mean(
          proportion_persisting,
          na.rm = TRUE
        ),
      
      median_persistence =
        median(
          proportion_persisting,
          na.rm = TRUE
        ),
      
      mean_turnover =
        mean(
          jaccard_turnover,
          na.rm = TRUE
        ),
      
      median_turnover =
        median(
          jaccard_turnover,
          na.rm = TRUE
        ),
      
      mean_not_redetected =
        mean(
          n_lost,
          na.rm = TRUE
        ),
      
      mean_new =
        mean(
          n_gained,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  # ====================================================================
  # TURNOVER BY LAG
  # ====================================================================
  
  turnover_by_lag <- asv_turnover_species %>%
    group_by(
      regional_status,
      lag_years
    ) %>%
    summarise(
      n_species_comparisons =
        n(),
      
      mean_persistence =
        mean(
          proportion_persisting,
          na.rm = TRUE
        ),
      
      median_persistence =
        median(
          proportion_persisting,
          na.rm = TRUE
        ),
      
      mean_turnover =
        mean(
          jaccard_turnover,
          na.rm = TRUE
        ),
      
      median_turnover =
        median(
          jaccard_turnover,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  # ====================================================================
  # FIGURE 1
  #
  # NAT vs NIS WITHIN EACH 50-YEAR BIN
  #
  # Points = all species observations
  # Boxplots = only status groups with >= MIN_N_STATUS species
  #
  # Significance:
  # *   raw p < 0.05
  # **  raw p < 0.01
  # *** raw p < 0.001
  #
  # FDR-supported comparisons additionally receive "FDR"
  # ====================================================================
  
  
  # --------------------------------------------------------------------
  # PERIOD ORDER
  # --------------------------------------------------------------------
  
  period_levels <- rich %>%
    distinct(
      period_start,
      period_label
    ) %>%
    arrange(period_start) %>%
    pull(period_label)
  
  
  plot_rich <- rich %>%
    mutate(
      period_label = factor(
        period_label,
        levels = period_levels
      )
    )
  
  
  # --------------------------------------------------------------------
  # BOXPLOTS ONLY WHERE THAT STATUS HAS >= 3 SPECIES
  # --------------------------------------------------------------------
  
  box_data <- plot_rich %>%
    group_by(
      site,
      period_label,
      regional_status
    ) %>%
    mutate(
      n_group = n_distinct(species)
    ) %>%
    ungroup() %>%
    filter(
      n_group >= MIN_N_STATUS
    )
  
  
  # --------------------------------------------------------------------
  # SIGNIFICANCE ANNOTATIONS
  #
  # IMPORTANT:
  # One test = NAT vs NIS within the SAME site and SAME period.
  # --------------------------------------------------------------------
  
  annotation_q1 <- within_bin %>%
    
    # Only comparisons where a statistical test was actually possible
    filter(
      !is.na(p_raw)
    ) %>%
    
    mutate(
      
      # Raw significance
      stars_raw = case_when(
        p_raw < 0.001 ~ "***",
        p_raw < 0.01  ~ "**",
        p_raw < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      
      # Explicit FDR flag
      FDR_supported =
        !is.na(p_FDR) &
        p_FDR < 0.05,
      
      # Label shown in figure
      plot_label = case_when(
        
        FDR_supported &
          p_raw < 0.001 ~ "***\nFDR",
        
        FDR_supported &
          p_raw < 0.01 ~ "**\nFDR",
        
        FDR_supported &
          p_raw < 0.05 ~ "*\nFDR",
        
        TRUE ~ stars_raw
      )
    ) %>%
    
    # Only annotate nominally significant comparisons
    filter(
      p_raw < 0.05
    ) %>%
    
    # Find maximum richness within that site × period
    left_join(
      
      rich %>%
        group_by(
          site,
          period_start,
          period_label
        ) %>%
        summarise(
          ymax = max(
            n_asvs,
            na.rm = TRUE
          ),
          .groups = "drop"
        ),
      
      by = c(
        "site",
        "period_start",
        "period_label"
      )
    ) %>%
    
    mutate(
      
      # Put annotation clearly above the highest observation
      y_position =
        ymax +
        pmax(
          1,
          ymax * 0.20
        ),
      
      period_label = factor(
        period_label,
        levels = period_levels
      )
    )
  
  
  # --------------------------------------------------------------------
  # PRINT EXACT COMPARISONS THAT WILL RECEIVE ASTERISKS
  # --------------------------------------------------------------------
  
  cat(
    "\n============================================================\n",
    marker_name,
    " — Q1 SIGNIFICANT NAT vs NIS WITHIN-BIN COMPARISONS\n",
    "============================================================\n"
  )
  
  
  if (nrow(annotation_q1) == 0) {
    
    cat(
      "No within-bin NAT vs NIS comparisons have raw p < 0.05.\n"
    )
    
  } else {
    
    print(
      annotation_q1 %>%
        select(
          site,
          period_start,
          period_label,
          n_NAT,
          n_NIS,
          mean_NAT,
          mean_NIS,
          delta_NIS_minus_NAT,
          W,
          p_raw,
          p_FDR,
          stars_raw,
          FDR_supported
        ),
      n = Inf
    )
  }
  
  
  # --------------------------------------------------------------------
  # FIGURE
  # --------------------------------------------------------------------
  
  p1 <- ggplot(
    plot_rich,
    aes(
      x = period_label,
      y = n_asvs,
      fill = regional_status
    )
  ) +
    
    # --------------------------------------------------------------
  # Boxplots
  # --------------------------------------------------------------
  
  geom_boxplot(
    data = box_data,
    
    aes(
      group = interaction(
        period_label,
        regional_status
      )
    ),
    
    position = position_dodge(
      width = 0.72
    ),
    
    width = 0.62,
    alpha = 0.38,
    colour = "black",
    linewidth = 0.6,
    outlier.shape = NA
  ) +
    
    
    # --------------------------------------------------------------
  # Individual species observations
  # --------------------------------------------------------------
  
  geom_point(
    aes(
      colour = regional_status
    ),
    
    position = position_jitterdodge(
      jitter.width = 0.08,
      jitter.height = 0,
      dodge.width = 0.72
    ),
    
    size = 2.5,
    alpha = 0.78
  ) +
    
    
    # --------------------------------------------------------------
  # Significance stars
  # --------------------------------------------------------------
  
  geom_text(
    data = annotation_q1,
    
    aes(
      x = period_label,
      y = y_position,
      label = plot_label
    ),
    
    inherit.aes = FALSE,
    
    size = 6.5,
    fontface = "bold",
    vjust = 0
  ) +
    
    
    # --------------------------------------------------------------
  # Sites
  # --------------------------------------------------------------
  
  facet_wrap(
    ~ site,
    ncol = 1,
    scales = "free"
  ) +
    
    
    # --------------------------------------------------------------
  # Colours
  # --------------------------------------------------------------
  
  scale_fill_manual(
    values = status_colors,
    drop = FALSE
  ) +
    
    scale_colour_manual(
      values = status_colors,
      drop = FALSE
    ) +
    
    
    # --------------------------------------------------------------
  # Extra vertical space so stars are not clipped
  # --------------------------------------------------------------
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.05,
        0.25
      )
    )
  ) +
    
    
    # --------------------------------------------------------------
  # Labels
  # --------------------------------------------------------------
  
  labs(
    x = "50-year period",
    
    y = paste0(
      "Putative haplotype richness\n",
      "(number of detected ASVs)"
    ),
    
    fill = NULL,
    colour = NULL,
    
    caption = paste0(
      "NAT vs NIS within the same 50-year period: ",
      "* p < 0.05; ** p < 0.01; *** p < 0.001. ",
      "\"FDR\" indicates comparisons also supported after ",
      "Benjamini–Hochberg correction."
    )
  ) +
    
    
    # --------------------------------------------------------------
  # Theme
  # --------------------------------------------------------------
  
  theme_bw(
    base_size = 19
  ) +
    
    theme(
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 16,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = 16,
        colour = "black"
      ),
      
      axis.title = element_text(
        size = 18
      ),
      
      strip.text = element_text(
        size = 18,
        face = "bold"
      ),
      
      legend.position = "top",
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.x =
        element_blank(),
      
      # Prevent significance annotations from being clipped
      plot.margin = margin(
        15,
        15,
        15,
        15
      )
    ) +
    
    
    # Allow annotations outside normal panel range
    coord_cartesian(
      clip = "off"
    )
  
  
  # ====================================================================
  # FIGURE 2 — SAME FOR 18S AND COI
  # ====================================================================
  
  heatmap_data <- all_timelags %>%
    filter(
      !is.na(mean_change)
    ) %>%
    mutate(
      significance_plot =
        case_when(
          FDR_significant ~ "*†",
          raw_significant ~ "*",
          TRUE ~ ""
        ),
      
      period_1_f =
        factor(
          period_1,
          levels =
            sort(
              unique(
                c(
                  period_1,
                  period_2
                )
              )
            )
        ),
      
      period_2_f =
        factor(
          period_2,
          levels =
            sort(
              unique(
                c(
                  period_1,
                  period_2
                )
              )
            )
        )
    )
  
  p2 <- ggplot(
    heatmap_data,
    aes(
      x = period_1_f,
      y = period_2_f,
      fill = mean_change
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.35
    ) +
    geom_text(
      aes(
        label = significance_plot
      ),
      size = 4.5,
      fontface = "bold"
    ) +
    facet_grid(
      regional_status ~ site,
      scales = "free",
      space = "free"
    ) +
    scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      name = "Δ ASVs"
    ) +
    labs(
      x = "Earlier 50-year period",
      y = "Later 50-year period",
      caption =
        "* raw p < 0.05; † also FDR < 0.05"
    ) +
    theme_bw(
      base_size = 18
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          size = 16
        ),
      axis.text.y =
        element_text(
          size = 16
        ),
      strip.text =
        element_text(
          size = 13,
          face = "bold"
        ),
      panel.grid =
        element_blank()
    )
  
  
  # ====================================================================
  # FIGURES 3–6
  #
  # 18S:
  # ORIGINAL GRAPHICAL BEHAVIOUR PRESERVED.
  #
  # COI:
  # SPARSE-DATA SAFEGUARDS.
  # ====================================================================
  
  
  if (marker_name == "18S") {
    
    # ------------------------------------------------------------------
    # FIGURE 3 — ORIGINAL 18S
    # ------------------------------------------------------------------
    
    p3 <- ggplot(
      all_timelags %>%
        filter(
          !is.na(mean_change)
        ),
      aes(
        x = lag_years,
        y = mean_change,
        colour = regional_status
      )
    ) +
      geom_hline(
        yintercept = 0,
        linetype = 2
      ) +
      geom_point(
        size = 2.8,
        alpha = 0.60
      ) +
      geom_smooth(
        method = "loess",
        se = TRUE,
        linewidth = 1
      ) +
      facet_wrap(
        ~ site,
        ncol = 2
      ) +
      scale_colour_manual(
        values = status_colors
      ) +
      labs(
        x = "Time lag between periods (years)",
        y =
          "Mean change in putative haplotype richness",
        colour = NULL,
        caption =
          "Lines are descriptive LOESS smooths."
      ) +
      theme_bw(
        base_size = 18
      ) +
      theme(
        legend.position = "top",
        strip.text =
          element_text(
            face = "bold"
          ),
        panel.grid.minor =
          element_blank()
      )
    
    
    # ------------------------------------------------------------------
    # FIGURE 4 — ORIGINAL 18S
    # ------------------------------------------------------------------
    
    p4 <- ggplot(
      differential_change %>%
        filter(
          sufficient_both,
          !is.na(differential_change)
        ),
      aes(
        x = lag_years,
        y = differential_change
      )
    ) +
      geom_hline(
        yintercept = 0,
        linetype = 2
      ) +
      geom_point(
        size = 3,
        alpha = 0.70
      ) +
      geom_smooth(
        method = "loess",
        se = TRUE
      ) +
      facet_wrap(
        ~ site,
        ncol = 2
      ) +
      labs(
        x =
          "Time lag between periods (years)",
        y =
          "ΔNIS − ΔNAT richness",
        caption =
          paste0(
            "Positive values indicate a more positive ",
            "richness change in NIS."
          )
      ) +
      theme_bw(
        base_size = 18
      ) +
      theme(
        strip.text =
          element_text(
            face = "bold"
          ),
        panel.grid.minor =
          element_blank()
      )
    
    
    # ------------------------------------------------------------------
    # FIGURE 5 — ORIGINAL 18S
    # ------------------------------------------------------------------
    
    p5 <- ggplot(
      asv_turnover_species,
      aes(
        x = lag_years,
        y = proportion_persisting,
        colour = regional_status
      )
    ) +
      geom_point(
        alpha = 0.22,
        size = 1.8
      ) +
      geom_smooth(
        method = "loess",
        se = TRUE,
        linewidth = 1.1
      ) +
      facet_wrap(
        ~ site,
        ncol = 2
      ) +
      scale_colour_manual(
        values = status_colors
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        labels =
          percent_format()
      ) +
      labs(
        x =
          "Time lag between periods (years)",
        y =
          paste0(
            "Proportion of earlier ASVs\n",
            "detected in both periods"
          ),
        colour = NULL,
        caption =
          paste0(
            "Repeated detection does not necessarily imply ",
            "continuous biological persistence."
          )
      ) +
      theme_bw(
        base_size = 18
      ) +
      theme(
        legend.position = "top",
        strip.text =
          element_text(
            face = "bold"
          ),
        panel.grid.minor =
          element_blank()
      )
    
    
    # ------------------------------------------------------------------
    # FIGURE 6 — ORIGINAL 18S
    # ------------------------------------------------------------------
    
    p6 <- ggplot(
      asv_turnover_species,
      aes(
        x = lag_years,
        y = jaccard_turnover,
        colour = regional_status
      )
    ) +
      geom_point(
        alpha = 0.22,
        size = 1.8
      ) +
      geom_smooth(
        method = "loess",
        se = TRUE,
        linewidth = 1.1
      ) +
      facet_wrap(
        ~ site,
        ncol = 2
      ) +
      scale_colour_manual(
        values = status_colors
      ) +
      scale_y_continuous(
        limits = c(0, 1)
      ) +
      labs(
        x =
          "Time lag between periods (years)",
        y =
          paste0(
            "ASV compositional turnover\n",
            "(1 − Jaccard similarity)"
          ),
        colour = NULL,
        caption =
          paste0(
            "0 = identical detected ASV composition; ",
            "1 = no ASVs detected in both periods."
          )
      ) +
      theme_bw(
        base_size = 18
      ) +
      theme(
        legend.position = "top",
        strip.text =
          element_text(
            face = "bold"
          ),
        panel.grid.minor =
          element_blank()
      )
    
  } else {
    
    # ==================================================================
    # COI ONLY — ROBUST GRAPHICAL VERSION
    # ==================================================================
    
    
    # ------------------------------------------------------------------
    # COI FIGURE 3
    # ------------------------------------------------------------------
    
    q3_data <- all_timelags %>%
      filter(
        !is.na(site),
        !is.na(lag_years),
        !is.na(mean_change),
        !is.na(regional_status)
      )
    
    cat(
      "\nCOI — observations available for Figure 3:",
      nrow(q3_data),
      "\n"
    )
    
    if (nrow(q3_data) > 0) {
      
      q3_smooth_data <- q3_data %>%
        group_by(
          site,
          regional_status
        ) %>%
        filter(
          n() >= 5,
          n_distinct(lag_years) >= 4
        ) %>%
        ungroup()
      
      p3 <- ggplot(
        q3_data,
        aes(
          x = lag_years,
          y = mean_change,
          colour = regional_status
        )
      ) +
        geom_hline(
          yintercept = 0,
          linetype = 2
        ) +
        geom_point(
          size = 2.8,
          alpha = 0.60
        ) +
        facet_wrap(
          ~ site,
          ncol = 2
        ) +
        scale_colour_manual(
          values = status_colors
        ) +
        labs(
          x =
            "Time lag between periods (years)",
          y =
            "Mean change in putative haplotype richness",
          colour = NULL,
          caption =
            "Points represent all available temporal comparisons."
        ) +
        theme_bw(
          base_size = 18
        ) +
        theme(
          legend.position = "top",
          strip.text =
            element_text(
              face = "bold"
            ),
          panel.grid.minor =
            element_blank()
        )
      
      if (nrow(q3_smooth_data) > 0) {
        
        p3 <- p3 +
          geom_smooth(
            data = q3_smooth_data,
            method = "loess",
            formula = y ~ x,
            se = TRUE,
            linewidth = 1
          )
      }
      
    } else {
      
      p3 <- NULL
      
      message(
        "COI: Figure 3 skipped — no valid data."
      )
    }
    
    
    # ------------------------------------------------------------------
    # COI FIGURE 4
    # ------------------------------------------------------------------
    
    q4_data <- differential_change %>%
      filter(
        sufficient_both,
        !is.na(site),
        !is.na(lag_years),
        !is.na(differential_change)
      )
    
    cat(
      "COI — observations available for Figure 4:",
      nrow(q4_data),
      "\n"
    )
    
    if (nrow(q4_data) > 0) {
      
      q4_smooth_data <- q4_data %>%
        group_by(site) %>%
        filter(
          n() >= 5,
          n_distinct(lag_years) >= 4
        ) %>%
        ungroup()
      
      p4 <- ggplot(
        q4_data,
        aes(
          x = lag_years,
          y = differential_change
        )
      ) +
        geom_hline(
          yintercept = 0,
          linetype = 2
        ) +
        geom_point(
          size = 3,
          alpha = 0.70
        ) +
        facet_wrap(
          ~ site,
          ncol = 2
        ) +
        labs(
          x =
            "Time lag between periods (years)",
          y =
            "ΔNIS − ΔNAT richness",
          caption =
            paste0(
              "Positive values indicate a more positive ",
              "richness change in NIS. Descriptive effect."
            )
        ) +
        theme_bw(
          base_size = 18
        ) +
        theme(
          strip.text =
            element_text(
              face = "bold"
            ),
          panel.grid.minor =
            element_blank()
        )
      
      if (nrow(q4_smooth_data) > 0) {
        
        p4 <- p4 +
          geom_smooth(
            data = q4_smooth_data,
            method = "loess",
            formula = y ~ x,
            se = TRUE,
            linewidth = 1
          )
      }
      
    } else {
      
      p4 <- NULL
      
      message(
        "COI: Figure 4 skipped — no site × temporal pair ",
        "contains at least ",
        MIN_N_PAIRED,
        " NAT and ",
        MIN_N_PAIRED,
        " NIS shared species."
      )
    }
    
    
    # ------------------------------------------------------------------
    # COI FIGURE 5
    # ------------------------------------------------------------------
    
    q5_data <- asv_turnover_species %>%
      filter(
        !is.na(site),
        !is.na(regional_status),
        !is.na(lag_years),
        !is.na(proportion_persisting)
      )
    
    cat(
      "COI — observations available for Figure 5:",
      nrow(q5_data),
      "\n"
    )
    
    if (nrow(q5_data) > 0) {
      
      q5_smooth_data <- q5_data %>%
        group_by(
          site,
          regional_status
        ) %>%
        filter(
          n() >= 5,
          n_distinct(lag_years) >= 4
        ) %>%
        ungroup()
      
      p5 <- ggplot(
        q5_data,
        aes(
          x = lag_years,
          y = proportion_persisting,
          colour = regional_status
        )
      ) +
        geom_point(
          alpha = 0.22,
          size = 1.8
        ) +
        facet_wrap(
          ~ site,
          ncol = 2
        ) +
        scale_colour_manual(
          values = status_colors
        ) +
        scale_y_continuous(
          limits = c(0, 1),
          labels = percent_format()
        ) +
        labs(
          x =
            "Time lag between periods (years)",
          y =
            paste0(
              "Proportion of earlier ASVs\n",
              "detected in both periods"
            ),
          colour = NULL
        ) +
        theme_bw(
          base_size = 18
        ) +
        theme(
          legend.position = "top",
          strip.text =
            element_text(
              face = "bold"
            ),
          panel.grid.minor =
            element_blank()
        )
      
      if (nrow(q5_smooth_data) > 0) {
        
        p5 <- p5 +
          geom_smooth(
            data = q5_smooth_data,
            method = "loess",
            formula = y ~ x,
            se = TRUE,
            linewidth = 1.1
          )
      }
      
    } else {
      
      p5 <- NULL
      
      message(
        "COI: Figure 5 skipped — no valid persistence data."
      )
    }
    
    
    # ------------------------------------------------------------------
    # COI FIGURE 6
    # ------------------------------------------------------------------
    
    q6_data <- asv_turnover_species %>%
      filter(
        !is.na(site),
        !is.na(regional_status),
        !is.na(lag_years),
        !is.na(jaccard_turnover)
      )
    
    cat(
      "COI — observations available for Figure 6:",
      nrow(q6_data),
      "\n"
    )
    
    if (nrow(q6_data) > 0) {
      
      q6_smooth_data <- q6_data %>%
        group_by(
          site,
          regional_status
        ) %>%
        filter(
          n() >= 5,
          n_distinct(lag_years) >= 4
        ) %>%
        ungroup()
      
      p6 <- ggplot(
        q6_data,
        aes(
          x = lag_years,
          y = jaccard_turnover,
          colour = regional_status
        )
      ) +
        geom_point(
          alpha = 0.22,
          size = 1.8
        ) +
        facet_wrap(
          ~ site,
          ncol = 2
        ) +
        scale_colour_manual(
          values = status_colors
        ) +
        scale_y_continuous(
          limits = c(0, 1)
        ) +
        labs(
          x =
            "Time lag between periods (years)",
          y =
            paste0(
              "ASV compositional turnover\n",
              "(1 − Jaccard similarity)"
            ),
          colour = NULL
        ) +
        theme_bw(
          base_size = 18
        ) +
        theme(
          legend.position = "top",
          strip.text =
            element_text(
              face = "bold"
            ),
          panel.grid.minor =
            element_blank()
        )
      
      if (nrow(q6_smooth_data) > 0) {
        
        p6 <- p6 +
          geom_smooth(
            data = q6_smooth_data,
            method = "loess",
            formula = y ~ x,
            se = TRUE,
            linewidth = 1.1
          )
      }
      
    } else {
      
      p6 <- NULL
      
      message(
        "COI: Figure 6 skipped — no valid turnover data."
      )
    }
  }
  
  
  # ====================================================================
  # FIGURE 7
  #
  # Robust check used because this does not alter any 18S result.
  # ====================================================================
  
  turnover_plot_data <- asv_turnover_species %>%
    filter(
      !is.na(site),
      !is.na(regional_status),
      regional_status %in% c("NAT", "NIS"),
      !is.na(jaccard_turnover)
    )
  
  cat(
    marker_name,
    " — observations available for Figure 7:",
    nrow(turnover_plot_data),
    "\n"
  )
  
  if (nrow(turnover_plot_data) > 0) {
    
    p7 <- ggplot(
      turnover_plot_data,
      aes(
        x = regional_status,
        y = jaccard_turnover,
        fill = regional_status
      )
    ) +
      geom_boxplot(
        width = 0.60,
        alpha = 0.40,
        outlier.shape = NA
      ) +
      geom_jitter(
        aes(
          colour = regional_status
        ),
        width = 0.12,
        height = 0,
        alpha = 0.18,
        size = 1.5
      ) +
      facet_wrap(
        ~ site,
        ncol = 4
      ) +
      scale_fill_manual(
        values = status_colors,
        drop = FALSE
      ) +
      scale_colour_manual(
        values = status_colors,
        drop = FALSE
      ) +
      scale_y_continuous(
        limits = c(0, 1)
      ) +
      labs(
        x = NULL,
        y =
          paste0(
            "ASV compositional turnover\n",
            "(1 − Jaccard similarity)"
          )
      ) +
      theme_bw(
        base_size = 18
      ) +
      theme(
        legend.position = "none",
        strip.text =
          element_text(
            face = "bold"
          ),
        panel.grid.minor =
          element_blank()
      )
    
  } else {
    
    p7 <- NULL
    
    message(
      marker_name,
      ": Figure 7 skipped — no valid turnover observations."
    )
  }
  
  
  # ====================================================================
  # FIGURE 8
  # ====================================================================
  
  if (nrow(asv_turnover_species) > 0) {
    
    composition_long <- asv_turnover_species %>%
      filter(
        !is.na(site),
        !is.na(regional_status),
        regional_status %in% c("NAT", "NIS")
      ) %>%
      select(
        site,
        species,
        regional_status,
        period_1,
        period_2,
        lag_years,
        n_lost,
        n_shared,
        n_gained
      ) %>%
      pivot_longer(
        cols =
          c(
            n_lost,
            n_shared,
            n_gained
          ),
        names_to = "component",
        values_to = "n_asvs"
      ) %>%
      filter(
        !is.na(n_asvs)
      ) %>%
      mutate(
        component =
          recode(
            component,
            n_lost =
              "Not redetected",
            n_shared =
              "Detected in both",
            n_gained =
              "Newly detected"
          ),
        
        component =
          factor(
            component,
            levels =
              c(
                "Not redetected",
                "Detected in both",
                "Newly detected"
              )
          )
      )
    
  } else {
    
    composition_long <- tibble()
  }
  
  cat(
    marker_name,
    " — observations available for Figure 8:",
    nrow(composition_long),
    "\n"
  )
  
  if (nrow(composition_long) > 0) {
    
    p8 <- ggplot(
      composition_long,
      aes(
        x = regional_status,
        y = n_asvs,
        fill = regional_status
      )
    ) +
      geom_boxplot(
        alpha = 0.40,
        outlier.shape = NA
      ) +
      geom_jitter(
        aes(
          colour = regional_status
        ),
        width = 0.12,
        height = 0,
        alpha = 0.15,
        size = 1.3
      ) +
      facet_grid(
        component ~ site,
        scales = "free",
        space = "free_x"
      ) +
      scale_fill_manual(
        values = status_colors,
        drop = FALSE
      ) +
      scale_colour_manual(
        values = status_colors,
        drop = FALSE
      ) +
      labs(
        x = NULL,
        y = "Number of ASVs"
      ) +
      theme_bw(
        base_size = 18
      ) +
      theme(
        legend.position = "none",
        strip.text =
          element_text(
            face = "bold"
          ),
        panel.grid.minor =
          element_blank()
      )
    
  } else {
    
    p8 <- NULL
    
    message(
      marker_name,
      ": Figure 8 skipped — no valid composition data."
    )
  }
  
  
  # ====================================================================
  # ROBUST SAVE FUNCTION
  #
  # NULL figures are simply skipped.
  # ====================================================================
  
  save_plot <- function(
    plot_object,
    filename,
    width,
    height
  ) {
    
    if (is.null(plot_object)) {
      
      message(
        "Skipping ",
        filename,
        ": no valid data available."
      )
      
      return(
        invisible(NULL)
      )
    }
    
    for (
      extension in
      c(
        "png",
        "pdf",
        "svg"
      )
    ) {
      
      ggsave(
        filename =
          file.path(
            marker_dir,
            paste0(
              filename,
              ".",
              extension
            )
          ),
        plot = plot_object,
        width = width,
        height = height,
        dpi = 400,
        bg = "white"
      )
    }
    
    invisible(NULL)
  }
  
  
  # ====================================================================
  # SAVE FIGURES
  # ====================================================================
  
  save_plot(
    p1,
    paste0(
      marker_name,
      "_Q1_NAT_NIS_WITHIN_PERIOD"
    ),
    16,
    21
  )
  
  save_plot(
    p2,
    paste0(
      marker_name,
      "_Q2_ALL_PAIRWISE_RICHNESS_HEATMAP"
    ),
    22,
    11
  )
  
  save_plot(
    p3,
    paste0(
      marker_name,
      "_Q3_RICHNESS_CHANGE_BY_LAG"
    ),
    15,
    14
  )
  
  save_plot(
    p4,
    paste0(
      marker_name,
      "_Q4_DIFFERENTIAL_NAT_NIS_CHANGE"
    ),
    15,
    14
  )
  
  save_plot(
    p5,
    paste0(
      marker_name,
      "_Q5_ASV_PERSISTENCE_BY_LAG"
    ),
    15,
    14
  )
  
  save_plot(
    p6,
    paste0(
      marker_name,
      "_Q6_ASV_TURNOVER_BY_LAG"
    ),
    15,
    14
  )
  
  save_plot(
    p7,
    paste0(
      marker_name,
      "_Q6_NAT_NIS_TURNOVER_BY_SITE"
    ),
    16,
    9
  )
  
  save_plot(
    p8,
    paste0(
      marker_name,
      "_Q6_ASV_COMPONENTS"
    ),
    20,
    11
  )
  
  
  # ====================================================================
  # EXCEL-SAFE ASV IDENTITIES
  # ====================================================================
  
  asv_turnover_excel <- asv_turnover_species %>%
    mutate(
      shared_asvs =
        map_chr(
          shared_asvs,
          ~ paste(
            .x,
            collapse = "; "
          )
        ),
      
      lost_asvs =
        map_chr(
          lost_asvs,
          ~ paste(
            .x,
            collapse = "; "
          )
        ),
      
      gained_asvs =
        map_chr(
          gained_asvs,
          ~ paste(
            .x,
            collapse = "; "
          )
        )
    )
  
  
  # ====================================================================
  # EXPORT EXCEL
  # ====================================================================
  
  write_xlsx(
    list(
      Q1_NAT_NIS_within_bin =
        within_bin,
      
      Q2_All_richness_pairs =
        all_timelags,
      
      Q3_Lag_summary =
        lag_summary,
      
      Q4_Differential_change =
        differential_change,
      
      Q5_ASV_species_pairs =
        asv_turnover_excel,
      
      Q5_ASV_summary =
        asv_turnover_summary,
      
      Q6_Turnover_tests =
        turnover_tests,
      
      Q6_Global_turnover =
        turnover_global,
      
      Q6_Turnover_by_lag =
        turnover_by_lag
    ),
    
    file.path(
      marker_dir,
      paste0(
        "FINAL_",
        marker_name,
        "_HAPLOTYPE_TEMPORAL_ANALYSIS.xlsx"
      )
    )
  )
  
  
  # ====================================================================
  # CONSOLE SUMMARY
  # ====================================================================
  
  cat(
    "\n------------------------------------------------------------\n",
    marker_name,
    " — Q1 NAT vs NIS WITHIN BINS\n",
    "------------------------------------------------------------\n"
  )
  
  cat(
    "Valid tests:",
    sum(
      !is.na(within_bin$p_raw)
    ),
    "\n"
  )
  
  cat(
    "Raw p < 0.05:",
    sum(
      within_bin$p_raw < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "FDR < 0.05:",
    sum(
      within_bin$p_FDR < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "\n------------------------------------------------------------\n",
    marker_name,
    " — Q2 ALL TEMPORAL PAIRS\n",
    "------------------------------------------------------------\n"
  )
  
  cat(
    "All temporal comparisons:",
    nrow(all_timelags),
    "\n"
  )
  
  cat(
    "Valid paired tests:",
    sum(
      !is.na(all_timelags$p_raw)
    ),
    "\n"
  )
  
  cat(
    "Raw p < 0.05:",
    sum(
      all_timelags$p_raw < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "FDR < 0.05:",
    sum(
      all_timelags$p_FDR < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "\n------------------------------------------------------------\n",
    marker_name,
    " — GLOBAL TURNOVER\n",
    "------------------------------------------------------------\n"
  )
  
  print(
    turnover_global
  )
  
  
  cat(
    "\n------------------------------------------------------------\n",
    marker_name,
    " — RAW SIGNIFICANT TURNOVER TESTS\n",
    "------------------------------------------------------------\n"
  )
  
  print(
    turnover_tests %>%
      filter(
        raw_significant
      ) %>%
      arrange(
        p_raw
      ),
    n = Inf
  )
  
  
  cat(
    "\n------------------------------------------------------------\n",
    marker_name,
    " — FDR-SUPPORTED TURNOVER TESTS\n",
    "------------------------------------------------------------\n"
  )
  
  print(
    turnover_tests %>%
      filter(
        FDR_significant
      ) %>%
      arrange(
        p_FDR
      ),
    n = Inf
  )
  
  
  # ====================================================================
  # RETURN
  # ====================================================================
  
  list(
    richness =
      rich,
    
    within_bin =
      within_bin,
    
    all_timelags =
      all_timelags,
    
    lag_summary =
      lag_summary,
    
    differential_change =
      differential_change,
    
    asv_turnover_species =
      asv_turnover_species,
    
    asv_turnover_summary =
      asv_turnover_summary,
    
    turnover_tests =
      turnover_tests,
    
    turnover_global =
      turnover_global,
    
    turnover_by_lag =
      turnover_by_lag,
    
    plots =
      list(
        Q1 = p1,
        Q2 = p2,
        Q3 = p3,
        Q4 = p4,
        Q5 = p5,
        Q6 = p6,
        Q7 = p7,
        Q8 = p8
      )
  )
}


# ======================================================================
# ======================================================================
#
# 27. RUN 18S
#
# MAIN ANALYSIS
#
# ======================================================================
# ======================================================================

results_18S <- analyse_marker(
  marker_name = "18S",
  richness_data = richness_all,
  asv_data = asv_samples
)


# ======================================================================
# ======================================================================
#
# 28. RUN COI
#
# SUPPLEMENTARY
#
# ======================================================================
# ======================================================================

results_COI <- analyse_marker(
  marker_name = "COI",
  richness_data = richness_all,
  asv_data = asv_samples
)


# ======================================================================
# 29. FINAL MARKER SUMMARY
# ======================================================================

marker_final_summary <- bind_rows(
  
  results_18S$richness %>%
    summarise(
      marker = "18S",
      
      n_sites =
        n_distinct(site),
      
      n_species =
        n_distinct(species),
      
      n_NAT_species =
        n_distinct(
          species[
            regional_status == "NAT"
          ]
        ),
      
      n_NIS_species =
        n_distinct(
          species[
            regional_status == "NIS"
          ]
        ),
      
      n_species_period_observations =
        n()
    ),
  
  results_COI$richness %>%
    summarise(
      marker = "COI",
      
      n_sites =
        n_distinct(site),
      
      n_species =
        n_distinct(species),
      
      n_NAT_species =
        n_distinct(
          species[
            regional_status == "NAT"
          ]
        ),
      
      n_NIS_species =
        n_distinct(
          species[
            regional_status == "NIS"
          ]
        ),
      
      n_species_period_observations =
        n()
    )
)


# ======================================================================
# 30. EXPORT CROSS-MARKER SUMMARY
# ======================================================================

write_xlsx(
  list(
    Marker_representation =
      marker_final_summary,
    
    Turnover_global_18S =
      results_18S$turnover_global,
    
    Turnover_global_COI =
      results_COI$turnover_global,
    
    Turnover_by_lag_18S =
      results_18S$turnover_by_lag,
    
    Turnover_by_lag_COI =
      results_COI$turnover_by_lag
  ),
  
  file.path(
    output_dir,
    "FINAL_18S_COI_SUMMARY.xlsx"
  )
)


# ======================================================================
# 31. FINAL MESSAGE
# ======================================================================

cat(
  "\n\n============================================================\n",
  "FINAL ANALYSIS COMPLETE\n",
  "============================================================\n\n"
)

cat(
  "18S = MAIN ANALYSIS\n",
  "COI = SUPPLEMENTARY ANALYSIS\n\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n\n"
)

cat(
  "MARKER REPRESENTATION:\n"
)

print(
  marker_final_summary
)

cat(
  "\nIMPORTANT:\n",
  "- 18S analytical pipeline has not been changed.\n",
  "- COI uses the same statistical analyses.\n",
  "- Only COI plotting is protected against sparse data.\n",
  "- Missing COI figures are skipped rather than treated as errors.\n",
  "- p_raw and p_FDR are retained.\n",
  "- No comparison is removed because it fails FDR.\n",
  "- 'Not redetected' does not mean extinction.\n",
  "- 'Newly detected' does not necessarily mean biological arrival.\n",
  "\n============================================================\n"
)


# ======================================================================
# 32. WARNINGS DIAGNOSTIC
# ======================================================================

cat(
  "\nIf warnings were generated, inspect them with:\n\n",
  "warnings()\n\n"
)

# ======================================================================
# END
# ======================================================================


















# ======================================================================
# ======================================================================
# FINAL MANUSCRIPT ANALYSIS
#
# LONG-TERM TEMPORAL TURNOVER OF WITHIN-SPECIES HAPLOTYPE COMPOSITION
#
# OUTPUTS
#
# MAIN FIGURE
#   18S:
#   A = observed site-specific temporal turnover
#   B = site-specific NAT-NIS turnover contrasts with cluster-bootstrap
#       95% confidence intervals
#
# SUPPLEMENTARY FIGURE
#   Q1 haplotype richness:
#   LEFT  = COI
#   RIGHT = 18S
#
# OPTIONAL SUPPLEMENTARY FIGURE
#   COI temporal turnover (descriptive only)
#
# SUPPLEMENTARY TABLE
#   marker × site × status representation and turnover
#
# STATISTICAL TABLES
#   site-specific NAT-NIS contrasts
#   cluster-bootstrap confidence intervals
#
# IMPORTANT
#
# - Run AFTER results_18S and results_COI have been generated.
# - Does NOT modify Q1-Q6.
# - Main turnover inference is based on 18S.
# - COI is complementary because NIS representation is sparse.
# - Temporal pairs are NOT treated as independent replicates.
# - Bootstrap resampling is performed at the species level.
# ======================================================================
# ======================================================================


# ======================================================================
# 0. PACKAGES
# ======================================================================

final_packages <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "patchwork",
  "writexl",
  "scales",
  "svglite"
)

missing_final_packages <-
  final_packages[
    !final_packages %in%
      rownames(installed.packages())
  ]

if (length(missing_final_packages) > 0) {
  install.packages(missing_final_packages)
}

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(writexl)
library(scales)
library(svglite)


# ======================================================================
# 1. SETTINGS
# ======================================================================

N_BOOT <- 5000

MIN_SPECIES_STATUS <- 3

set.seed(1310)


# ======================================================================
# 2. OUTPUT DIRECTORY
# ======================================================================

final_manuscript_dir <-
  file.path(
    output_dir,
    "FINAL_MANUSCRIPT_ANALYSIS"
  )

dir.create(
  final_manuscript_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# 3. COLOURS
# ======================================================================

status_colors_final <- c(
  "NAT" = "#E69F00",
  "NIS" = "#CC3366"
)


# ======================================================================
# 4. CHECK REQUIRED OBJECTS
# ======================================================================

if (!exists("results_18S")) {
  stop("results_18S does not exist.")
}

if (!exists("results_COI")) {
  stop("results_COI does not exist.")
}


required_result_objects <- c(
  "within_bin",
  "asv_turnover_species"
)


for (nm in required_result_objects) {
  
  if (is.null(results_18S[[nm]])) {
    
    stop(
      paste0(
        "results_18S$",
        nm,
        " is missing.\n",
        "Add it to the return() list of analyse_marker()."
      )
    )
  }
  
  if (is.null(results_COI[[nm]])) {
    
    stop(
      paste0(
        "results_COI$",
        nm,
        " is missing.\n",
        "Add it to the return() list of analyse_marker()."
      )
    )
  }
}


# ======================================================================
# ======================================================================
# PART I
# PREPARE TURNOVER DATA
# ======================================================================
# ======================================================================


prepare_turnover_final <- function(
    x,
    marker_name
) {
  
  x %>%
    
    filter(
      !is.na(site),
      !is.na(species),
      !is.na(regional_status),
      !is.na(period_1),
      !is.na(period_2),
      !is.na(lag_years),
      !is.na(jaccard_turnover),
      regional_status %in% c(
        "NAT",
        "NIS"
      )
    ) %>%
    
    mutate(
      
      marker =
        marker_name,
      
      site =
        as.character(site),
      
      species =
        as.character(species),
      
      regional_status =
        factor(
          regional_status,
          levels = c(
            "NAT",
            "NIS"
          )
        ),
      
      lag_years =
        as.numeric(lag_years),
      
      lag_century =
        lag_years / 100,
      
      jaccard_turnover =
        as.numeric(
          jaccard_turnover
        )
    )
}


turnover_18S <-
  prepare_turnover_final(
    results_18S$asv_turnover_species,
    "18S"
  )


turnover_COI <-
  prepare_turnover_final(
    results_COI$asv_turnover_species,
    "COI"
  )


turnover_all <-
  bind_rows(
    turnover_18S,
    turnover_COI
  )


# ======================================================================
# 5. SANITY CHECK
# ======================================================================

if (
  any(
    turnover_all$jaccard_turnover < 0 |
    turnover_all$jaccard_turnover > 1,
    na.rm = TRUE
  )
) {
  
  stop(
    "Jaccard turnover outside [0,1]."
  )
}


# ======================================================================
# ======================================================================
# PART II
# DESCRIPTIVE REPRESENTATION
# ======================================================================
# ======================================================================


# ======================================================================
# 6. MARKER × SITE × STATUS REPRESENTATION
# ======================================================================

representation_table <-
  turnover_all %>%
  
  group_by(
    marker,
    site,
    regional_status
  ) %>%
  
  summarise(
    
    n_species =
      n_distinct(species),
    
    n_temporal_comparisons =
      n(),
    
    n_period_pairs =
      n_distinct(
        interaction(
          period_1,
          period_2
        )
      ),
    
    min_lag_years =
      min(
        lag_years,
        na.rm = TRUE
      ),
    
    max_lag_years =
      max(
        lag_years,
        na.rm = TRUE
      ),
    
    mean_turnover =
      mean(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    median_turnover =
      median(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    sd_turnover =
      if_else(
        n() > 1,
        sd(
          jaccard_turnover,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    proportion_zero =
      mean(
        jaccard_turnover == 0,
        na.rm = TRUE
      ),
    
    proportion_one =
      mean(
        jaccard_turnover == 1,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ======================================================================
# 7. SITE × STATUS × LAG
# ======================================================================

site_lag_summary <-
  turnover_all %>%
  
  group_by(
    marker,
    site,
    regional_status,
    lag_years
  ) %>%
  
  summarise(
    
    n_comparisons =
      n(),
    
    n_species =
      n_distinct(species),
    
    mean_turnover =
      mean(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    median_turnover =
      median(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    sd_turnover =
      if_else(
        n() > 1,
        sd(
          jaccard_turnover,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    .groups = "drop"
  )


# ======================================================================
# ======================================================================
# PART III
# SITE-SPECIFIC 18S NAT-NIS CONTRAST
#
# IMPORTANT:
#
# We first calculate one mean turnover per SPECIES.
#
# This prevents species contributing many temporal pairs from dominating
# the site-level NAT-NIS comparison.
#
# The estimand is therefore:
#
# mean species-level turnover in NIS
# minus
# mean species-level turnover in NAT
#
# within each site.
# ======================================================================
# ======================================================================


# ======================================================================
# 8. ONE VALUE PER SPECIES × SITE
# ======================================================================

species_site_turnover_18S <-
  turnover_18S %>%
  
  group_by(
    site,
    species,
    regional_status
  ) %>%
  
  summarise(
    
    n_temporal_pairs =
      n(),
    
    mean_turnover =
      mean(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    median_turnover =
      median(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    mean_lag =
      mean(
        lag_years,
        na.rm = TRUE
      ),
    
    max_lag =
      max(
        lag_years,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ======================================================================
# 9. OBSERVED SITE CONTRASTS
# ======================================================================

observed_site_contrasts <-
  species_site_turnover_18S %>%
  
  group_by(
    site
  ) %>%
  
  group_modify(
    ~ {
      
      d <- .x
      
      nat <-
        d %>%
        filter(
          regional_status == "NAT"
        )
      
      nis <-
        d %>%
        filter(
          regional_status == "NIS"
        )
      
      tibble(
        
        n_species_NAT =
          nrow(nat),
        
        n_species_NIS =
          nrow(nis),
        
        mean_NAT =
          ifelse(
            nrow(nat) > 0,
            mean(
              nat$mean_turnover,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
        mean_NIS =
          ifelse(
            nrow(nis) > 0,
            mean(
              nis$mean_turnover,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
        difference_NIS_minus_NAT =
          ifelse(
            nrow(nat) > 0 &
              nrow(nis) > 0,
            
            mean(
              nis$mean_turnover,
              na.rm = TRUE
            ) -
              mean(
                nat$mean_turnover,
                na.rm = TRUE
              ),
            
            NA_real_
          ),
        
        sufficient =
          nrow(nat) >=
          MIN_SPECIES_STATUS &
          nrow(nis) >=
          MIN_SPECIES_STATUS
      )
    }
  ) %>%
  
  ungroup()


# ======================================================================
# ======================================================================
# PART IV
# CLUSTER BOOTSTRAP BY SPECIES
#
# Resampling unit = species within site × status.
#
# This avoids pretending that all temporal pairs from the same species
# are independent replicates.
# ======================================================================
# ======================================================================


# ======================================================================
# 10. BOOTSTRAP FUNCTION FOR ONE SITE
# ======================================================================

bootstrap_site_contrast <- function(
    site_name,
    data,
    n_boot = 5000
) {
  
  d <-
    data %>%
    filter(
      site == site_name
    )
  
  
  nat <-
    d %>%
    filter(
      regional_status == "NAT"
    )
  
  
  nis <-
    d %>%
    filter(
      regional_status == "NIS"
    )
  
  
  n_nat <-
    nrow(nat)
  
  n_nis <-
    nrow(nis)
  
  
  if (
    n_nat < MIN_SPECIES_STATUS |
    n_nis < MIN_SPECIES_STATUS
  ) {
    
    return(
      tibble(
        
        site =
          site_name,
        
        n_species_NAT =
          n_nat,
        
        n_species_NIS =
          n_nis,
        
        estimate =
          NA_real_,
        
        CI_low =
          NA_real_,
        
        CI_high =
          NA_real_,
        
        bootstrap_probability_positive =
          NA_real_,
        
        bootstrap_probability_negative =
          NA_real_,
        
        supported_direction =
          "Insufficient species"
      )
    )
  }
  
  
  observed <-
    mean(
      nis$mean_turnover,
      na.rm = TRUE
    ) -
    mean(
      nat$mean_turnover,
      na.rm = TRUE
    )
  
  
  boot_diff <-
    replicate(
      
      n_boot,
      
      {
        
        nat_boot <-
          nat[
            sample(
              seq_len(n_nat),
              size = n_nat,
              replace = TRUE
            ),
          ]
        
        
        nis_boot <-
          nis[
            sample(
              seq_len(n_nis),
              size = n_nis,
              replace = TRUE
            ),
          ]
        
        
        mean(
          nis_boot$mean_turnover,
          na.rm = TRUE
        ) -
          mean(
            nat_boot$mean_turnover,
            na.rm = TRUE
          )
      }
    )
  
  
  CI <-
    quantile(
      boot_diff,
      probs =
        c(
          0.025,
          0.975
        ),
      na.rm = TRUE,
      names = FALSE
    )
  
  
  direction <-
    case_when(
      
      CI[1] > 0 ~
        "NIS > NAT",
      
      CI[2] < 0 ~
        "NAT > NIS",
      
      TRUE ~
        "No supported difference"
    )
  
  
  tibble(
    
    site =
      site_name,
    
    n_species_NAT =
      n_nat,
    
    n_species_NIS =
      n_nis,
    
    estimate =
      observed,
    
    CI_low =
      CI[1],
    
    CI_high =
      CI[2],
    
    bootstrap_probability_positive =
      mean(
        boot_diff > 0,
        na.rm = TRUE
      ),
    
    bootstrap_probability_negative =
      mean(
        boot_diff < 0,
        na.rm = TRUE
      ),
    
    supported_direction =
      direction
  )
}


# ======================================================================
# 11. RUN BOOTSTRAP FOR ALL SITES
# ======================================================================

sites_18S <-
  sort(
    unique(
      species_site_turnover_18S$site
    )
  )


set.seed(1310)


bootstrap_site_results <-
  map_dfr(
    sites_18S,
    bootstrap_site_contrast,
    data =
      species_site_turnover_18S,
    n_boot =
      N_BOOT
  )


cat(
  "\n============================================================\n",
  "18S SITE-SPECIFIC NAT-NIS CONTRASTS\n",
  "============================================================\n"
)

print(
  bootstrap_site_results,
  n = Inf
)


# ======================================================================
# ======================================================================
# PART V
# GLOBAL 18S CONTRAST
#
# IMPORTANT:
#
# We calculate site-specific NAT-NIS differences first and then average
# those site effects.
#
# Therefore large sites do not automatically dominate the global
# estimate merely because they contain more temporal pairs.
# ======================================================================
# ======================================================================


# ======================================================================
# 12. GLOBAL OBSERVED SITE-AVERAGED EFFECT
# ======================================================================

global_site_effect_data <-
  observed_site_contrasts %>%
  
  filter(
    sufficient,
    !is.na(
      difference_NIS_minus_NAT
    )
  )


global_observed_effect <-
  mean(
    global_site_effect_data$difference_NIS_minus_NAT,
    na.rm = TRUE
  )


# ======================================================================
# 13. GLOBAL HIERARCHICAL BOOTSTRAP
#
# Resample:
#   1. sites
#   2. species within status within selected site
#
# The final statistic is the mean NAT-NIS difference across sites.
# ======================================================================

bootstrap_global_effect <- function(
    data,
    n_boot = 5000
) {
  
  eligible_sites <-
    data %>%
    
    group_by(
      site,
      regional_status
    ) %>%
    
    summarise(
      n_species = n(),
      .groups = "drop"
    ) %>%
    
    pivot_wider(
      names_from =
        regional_status,
      values_from =
        n_species,
      values_fill =
        0
    ) %>%
    
    filter(
      NAT >= MIN_SPECIES_STATUS,
      NIS >= MIN_SPECIES_STATUS
    ) %>%
    
    pull(site)
  
  
  if (
    length(
      eligible_sites
    ) < 2
  ) {
    
    return(
      tibble(
        estimate = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
    )
  }
  
  
  boot_global <-
    replicate(
      
      n_boot,
      
      {
        
        sampled_sites <-
          sample(
            eligible_sites,
            size =
              length(
                eligible_sites
              ),
            replace = TRUE
          )
        
        
        site_differences <-
          map_dbl(
            
            sampled_sites,
            
            function(s) {
              
              d <-
                data %>%
                filter(
                  site == s
                )
              
              
              nat <-
                d %>%
                filter(
                  regional_status == "NAT"
                )
              
              
              nis <-
                d %>%
                filter(
                  regional_status == "NIS"
                )
              
              
              nat_values <-
                sample(
                  nat$mean_turnover,
                  size =
                    nrow(nat),
                  replace = TRUE
                )
              
              
              nis_values <-
                sample(
                  nis$mean_turnover,
                  size =
                    nrow(nis),
                  replace = TRUE
                )
              
              
              mean(
                nis_values,
                na.rm = TRUE
              ) -
                mean(
                  nat_values,
                  na.rm = TRUE
                )
            }
          )
        
        
        mean(
          site_differences,
          na.rm = TRUE
        )
      }
    )
  
  
  CI <-
    quantile(
      boot_global,
      probs =
        c(
          0.025,
          0.975
        ),
      na.rm = TRUE,
      names = FALSE
    )
  
  
  tibble(
    
    estimate =
      global_observed_effect,
    
    CI_low =
      CI[1],
    
    CI_high =
      CI[2],
    
    bootstrap_probability_positive =
      mean(
        boot_global > 0,
        na.rm = TRUE
      ),
    
    bootstrap_probability_negative =
      mean(
        boot_global < 0,
        na.rm = TRUE
      )
  )
}


set.seed(1310)


global_bootstrap_result <-
  bootstrap_global_effect(
    species_site_turnover_18S,
    n_boot =
      N_BOOT
  )


cat(
  "\n============================================================\n",
  "GLOBAL SITE-AVERAGED 18S NAT-NIS CONTRAST\n",
  "============================================================\n"
)

print(
  global_bootstrap_result
)


# ======================================================================
# ======================================================================
# PART VI
# TEMPORAL DISTANCE EFFECT
#
# We calculate one slope per species wherever the species has enough
# temporal information.
#
# This answers:
#
# Does compositional turnover tend to increase as periods become
# further separated?
#
# IMPORTANT:
# Species slopes are the analysis unit here.
# ======================================================================
# ======================================================================


# ======================================================================
# 14. SPECIES-SPECIFIC TEMPORAL SLOPES
# ======================================================================

species_temporal_slopes_18S <-
  turnover_18S %>%
  
  group_by(
    site,
    species,
    regional_status
  ) %>%
  
  group_modify(
    ~ {
      
      d <- .x
      
      n_lags <-
        n_distinct(
          d$lag_years
        )
      
      
      if (
        nrow(d) < 3 |
        n_lags < 2
      ) {
        
        return(
          tibble(
            n_comparisons =
              nrow(d),
            
            n_lags =
              n_lags,
            
            temporal_slope =
              NA_real_
          )
        )
      }
      
      
      fit <-
        lm(
          jaccard_turnover ~
            lag_century,
          data = d
        )
      
      
      tibble(
        
        n_comparisons =
          nrow(d),
        
        n_lags =
          n_lags,
        
        temporal_slope =
          unname(
            coef(fit)[
              "lag_century"
            ]
          )
      )
    }
  ) %>%
  
  ungroup()


# ======================================================================
# 15. SUMMARY OF SPECIES SLOPES BY SITE × STATUS
# ======================================================================

temporal_slope_summary <-
  species_temporal_slopes_18S %>%
  
  filter(
    !is.na(
      temporal_slope
    )
  ) %>%
  
  group_by(
    site,
    regional_status
  ) %>%
  
  summarise(
    
    n_species =
      n(),
    
    mean_slope =
      mean(
        temporal_slope,
        na.rm = TRUE
      ),
    
    median_slope =
      median(
        temporal_slope,
        na.rm = TRUE
      ),
    
    sd_slope =
      if_else(
        n() > 1,
        sd(
          temporal_slope,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    proportion_positive =
      mean(
        temporal_slope > 0,
        na.rm = TRUE
      ),
    
    proportion_negative =
      mean(
        temporal_slope < 0,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ======================================================================
# ======================================================================
# PART VII
# MAIN FIGURE — 18S
# ======================================================================
# ======================================================================


# ======================================================================
# 16. PANEL A
#
# SITE-SPECIFIC OBSERVED TURNOVER BY TEMPORAL LAG
# ======================================================================

main_A_data <-
  site_lag_summary %>%
  
  filter(
    marker == "18S"
  )


main_panel_A <-
  ggplot(
    main_A_data,
    aes(
      x = lag_years,
      y = mean_turnover,
      colour = regional_status,
      group = regional_status
    )
  ) +
  
  geom_line(
    linewidth = 0.8,
    alpha = 0.80
  ) +
  
  geom_point(
    aes(
      size =
        n_comparisons
    ),
    alpha = 0.90
  ) +
  
  facet_wrap(
    ~ site,
    scales = "free",
    ncol = 4
  ) +
  
  scale_colour_manual(
    values =
      status_colors_final,
    drop = FALSE
  ) +
  
  scale_size_continuous(
    name =
      "Species comparisons",
    range =
      c(
        1.5,
        6
      )
  ) +
  
  scale_y_continuous(
    limits =
      c(
        0,
        1
      ),
    breaks =
      seq(
        0,
        1,
        by = 0.25
      )
  ) +
  
  labs(
    
    x =
      "Temporal separation (years)",
    
    y =
      "Mean haplotype compositional turnover\n(Jaccard dissimilarity)",
    
    colour =
      NULL,
    
    title =
      "A  Temporal haplotype turnover within sedimentary archives"
  ) +
  
  theme_bw(
    base_size = 18
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    panel.grid.minor =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    plot.title =
      element_text(
        face = "bold"
      )
  )


# ======================================================================
# 17. PANEL B
#
# SITE-SPECIFIC NAT-NIS EFFECT SIZE
#
# Points = mean species-level difference
# Lines  = cluster-bootstrap 95% CI
#
# Positive values:
#     NIS > NAT turnover
#
# Negative values:
#     NAT > NIS turnover
#
# This is the inferential/statistical panel.
# ======================================================================

main_B_data <-
  bootstrap_site_results %>%
  
  filter(
    !is.na(
      estimate
    )
  ) %>%
  
  mutate(
    
    site =
      factor(
        site,
        levels =
          rev(
            sort(
              unique(site)
            )
          )
      ),
    
    supported =
      if_else(
        CI_low > 0 |
          CI_high < 0,
        "95% CI excludes zero",
        "95% CI overlaps zero"
      )
  )


main_panel_B <-
  ggplot(
    main_B_data,
    aes(
      x = estimate,
      y = site
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.7
  ) +
  
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.18,
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(
      shape = supported
    ),
    size = 3.2
  ) +
  
  scale_shape_manual(
    values = c(
      "95% CI excludes zero" = 16,
      "95% CI overlaps zero" = 1
    ),
    name = NULL
  ) +
  
  labs(
    
    x =
      expression(
        Delta *
          " mean turnover (NIS - NAT)"
      ),
    
    y =
      NULL,
    
    title =
      "B  Site-specific NAT–NIS contrasts"
  ) +
  
  theme_bw(
    base_size = 18
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    panel.grid.minor =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      )
  )


# ======================================================================
# 18. ADD GLOBAL EFFECT TO PANEL B AS ANNOTATION
# ======================================================================

if (
  nrow(
    global_bootstrap_result
  ) == 1 &&
  !is.na(
    global_bootstrap_result$estimate
  )
) {
  
  global_label <-
    paste0(
      
      "Site-averaged effect = ",
      
      sprintf(
        "%.2f",
        global_bootstrap_result$estimate
      ),
      
      " [95% CI ",
      
      sprintf(
        "%.2f",
        global_bootstrap_result$CI_low
      ),
      
      ", ",
      
      sprintf(
        "%.2f",
        global_bootstrap_result$CI_high
      ),
      
      "]"
    )
  
  
  main_panel_B <-
    main_panel_B +
    
    labs(
      subtitle =
        global_label
    )
}


# ======================================================================
# 19. COMBINE MAIN FIGURE
#
# A is larger because it contains eight site facets.
# ======================================================================

main_18S_figure <-
  main_panel_A /
  main_panel_B +
  
  patchwork::plot_layout(
    heights =
      c(
        2.3,
        1
      )
  )


# ======================================================================
# 20. SAVE MAIN FIGURE
# ======================================================================

ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "MAIN_18S_TEMPORAL_HAPLOTYPE_TURNOVER.png"
    ),
  
  plot =
    main_18S_figure,
  
  width = 26,
  height = 16,
  dpi = 600,
  bg = "white"
)

# SVG vector copy (DPI does not apply to vector output)
ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "MAIN_18S_TEMPORAL_HAPLOTYPE_TURNOVER.svg"
    ),
  
  plot =
    main_18S_figure,
  
  width = 26,
  height = 16,
  bg = "white"
)


ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "MAIN_18S_TEMPORAL_HAPLOTYPE_TURNOVER.pdf"
    ),
  
  plot =
    main_18S_figure,
  
  width = 26,
  height = 16,
  bg = "white"
)


# ======================================================================
# ======================================================================
# PART VIII
# SUPPLEMENTARY Q1 FIGURE
#
# COI LEFT
# 18S RIGHT
# ======================================================================
# ======================================================================


# ======================================================================
# 21. PREPARE Q1 DATA
#
# We use the already-generated Q1 data.
#
# Expected columns include:
# site
# period_start
# n_NAT
# n_NIS
# mean_NAT
# mean_NIS
# p_raw
# p_FDR
#
# However, for the boxplots we need the underlying richness observations.
#
# Therefore use dat18 / datCOI if available.
# ======================================================================
# ======================================================================
# 21. PREPARE Q1 DATA — ROBUST VERSION
#
# The final manuscript block should work even if dat18/datCOI are no
# longer present in the current R environment.
#
# Priority:
#   1. Use dat18/datCOI if they still exist.
#   2. Otherwise reconstruct the richness observations from richness_all.
#
# Q1 requires INDIVIDUAL species × site × period richness observations,
# not the already summarised within_bin table.
# ======================================================================


# ----------------------------------------------------------------------
# FUNCTION TO PREPARE Q1 DATA
# ----------------------------------------------------------------------

prepare_Q1_plot_data <- function(
    x,
    marker_name
) {
  
  x %>%
    
    filter(
      marker == marker_name,
      regional_status %in% c(
        "NAT",
        "NIS"
      ),
      !is.na(n_asvs),
      !is.na(period_start),
      !is.na(site),
      !is.na(species)
    ) %>%
    
    mutate(
      
      marker =
        marker_name,
      
      regional_status =
        factor(
          regional_status,
          levels = c(
            "NAT",
            "NIS"
          )
        ),
      
      site =
        as.character(site),
      
      species =
        as.character(species),
      
      period_start =
        as.numeric(period_start),
      
      n_asvs =
        as.numeric(n_asvs)
    )
}


# ----------------------------------------------------------------------
# 18S
# ----------------------------------------------------------------------

if (exists("dat18")) {
  
  message(
    "Q1 18S: using existing dat18 object."
  )
  
  Q1_18S_data <-
    dat18 %>%
    
    mutate(
      marker = "18S"
    ) %>%
    
    prepare_Q1_plot_data(
      marker_name = "18S"
    )
  
} else {
  
  message(
    "Q1 18S: dat18 not found. Reconstructing from richness_all."
  )
  
  if (!exists("richness_all")) {
    
    stop(
      paste0(
        "Neither dat18 nor richness_all exists.\n",
        "Q1 individual richness observations cannot be reconstructed.\n",
        "Run the first data-preparation part of the original script first."
      )
    )
  }
  
  Q1_18S_data <-
    richness_all %>%
    
    prepare_Q1_plot_data(
      marker_name = "18S"
    )
}


# ----------------------------------------------------------------------
# COI
# ----------------------------------------------------------------------

if (exists("datCOI")) {
  
  message(
    "Q1 COI: using existing datCOI object."
  )
  
  Q1_COI_data <-
    datCOI %>%
    
    mutate(
      marker = "COI"
    ) %>%
    
    prepare_Q1_plot_data(
      marker_name = "COI"
    )
  
} else {
  
  message(
    "Q1 COI: datCOI not found. Reconstructing from richness_all."
  )
  
  if (!exists("richness_all")) {
    
    stop(
      paste0(
        "Neither datCOI nor richness_all exists.\n",
        "Q1 individual richness observations cannot be reconstructed.\n",
        "Run the first data-preparation part of the original script first."
      )
    )
  }
  
  Q1_COI_data <-
    richness_all %>%
    
    prepare_Q1_plot_data(
      marker_name = "COI"
    )
}


# ----------------------------------------------------------------------
# SANITY CHECK
# ----------------------------------------------------------------------

cat(
  "\n============================================================\n",
  "Q1 DATA RECONSTRUCTED FOR FINAL SUPPLEMENTARY FIGURE\n",
  "============================================================\n"
)


Q1_check <-
  bind_rows(
    Q1_COI_data,
    Q1_18S_data
  ) %>%
  
  group_by(
    marker,
    regional_status
  ) %>%
  
  summarise(
    
    n_observations =
      n(),
    
    n_species =
      n_distinct(species),
    
    n_sites =
      n_distinct(site),
    
    n_periods =
      n_distinct(
        interaction(
          site,
          period_start
        )
      ),
    
    .groups = "drop"
  )


print(
  Q1_check,
  n = Inf
)


if (nrow(Q1_18S_data) == 0) {
  stop(
    "Q1_18S_data is empty."
  )
}


if (nrow(Q1_COI_data) == 0) {
  warning(
    "Q1_COI_data is empty."
  )
}
# ======================================================================
# 22. Q1 PLOT FUNCTION
#
# No significance stars are added to COI because no eligible
# NAT-NIS tests exist.
#
# 18S annotations use FDR-supported tests ONLY.
# ======================================================================

make_Q1_supp_plot <- function(
    plot_data,
    q1_results,
    marker_name
) {
  
  
  p <-
    ggplot(
      plot_data,
      aes(
        x =
          factor(
            period_start
          ),
        y =
          n_asvs,
        fill =
          regional_status
      )
    ) +
    
    geom_boxplot(
      position =
        position_dodge(
          width = 0.8
        ),
      width = 0.70,
      outlier.shape = NA,
      alpha = 0.80
    ) +
    
    geom_jitter(
      aes(
        colour =
          regional_status
      ),
      position =
        position_jitterdodge(
          jitter.width = 0.15,
          dodge.width = 0.8
        ),
      size = 0.8,
      alpha = 0.35
    ) +
    
    facet_wrap(
      ~ site,
      scales = "free",
      ncol = 2
    ) +
    
    scale_fill_manual(
      values =
        status_colors_final,
      drop = FALSE
    ) +
    
    scale_colour_manual(
      values =
        status_colors_final,
      drop = FALSE
    ) +
    
    labs(
      
      x =
        "50-year period",
      
      y =
        "Detected putative haplotype richness\n(number of ASVs)",
      
      fill =
        NULL,
      
      colour =
        NULL,
      
      title =
        marker_name
    ) +
    
    theme_bw(
      base_size = 18
    ) +
    
    theme(
      
      legend.position =
        "top",
      
      axis.text.x =
        element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 16
        ),
      
      strip.text =
        element_text(
          face = "bold"
        ),
      
      panel.grid.minor =
        element_blank(),
      
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        )
    )
  
  
  # ------------------------------------------------------------
  # FDR annotations only
  # ------------------------------------------------------------
  
  if (
    !is.null(
      q1_results
    ) &&
    "p_FDR" %in%
    names(
      q1_results
    )
  ) {
    
    sig <-
      q1_results %>%
      
      filter(
        !is.na(
          p_FDR
        ),
        p_FDR < 0.05
      )
    
    
    if (
      nrow(sig) > 0
    ) {
      
      y_positions <-
        plot_data %>%
        
        group_by(
          site,
          period_start
        ) %>%
        
        summarise(
          y =
            max(
              n_asvs,
              na.rm = TRUE
            ) *
            1.12,
          .groups = "drop"
        )
      
      
      sig <-
        sig %>%
        
        left_join(
          y_positions,
          by =
            c(
              "site",
              "period_start"
            )
        ) %>%
        
        mutate(
          
          label =
            case_when(
              
              p_FDR < 0.001 ~
                "***",
              
              p_FDR < 0.01 ~
                "**",
              
              TRUE ~
                "*"
            )
        )
      
      
      p <-
        p +
        
        geom_text(
          data =
            sig,
          aes(
            x =
              factor(
                period_start
              ),
            y =
              y,
            label =
              label
          ),
          inherit.aes = FALSE,
          size = 4
        )
    }
  }
  
  
  p
}


# ======================================================================
# 23. CREATE COI + 18S Q1 PANELS
# ======================================================================

supp_Q1_COI <-
  make_Q1_supp_plot(
    Q1_COI_data,
    results_COI$within_bin,
    "A  COI"
  )


supp_Q1_18S <-
  make_Q1_supp_plot(
    Q1_18S_data,
    results_18S$within_bin,
    "B  18S"
  )


# ======================================================================
# 24. COMBINE:
#
# LEFT  = COI
# RIGHT = 18S
# ======================================================================

supp_Q1_combined <-
  supp_Q1_COI +
  supp_Q1_18S +
  
  patchwork::plot_layout(
    widths =
      c(
        1,
        1
      ),
    guides =
      "collect"
  ) &
  
  theme(
    legend.position =
      "top"
  )


# ======================================================================
# 25. SAVE SUPPLEMENTARY Q1
# ======================================================================

ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "SUPP_Q1_HAPLOTYPE_RICHNESS_COI_18S.png"
    ),
  
  plot =
    supp_Q1_combined,
  
  width = 18,
  height = 10,
  dpi = 600,
  bg = "white"
)

# SVG vector copy (DPI does not apply to vector output)
ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "SUPP_Q1_HAPLOTYPE_RICHNESS_COI_18S.svg"
    ),
  
  plot =
    supp_Q1_combined,
  
  width = 18,
  height = 10,
  bg = "white"
)


ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "SUPP_Q1_HAPLOTYPE_RICHNESS_COI_18S.pdf"
    ),
  
  plot =
    supp_Q1_combined,
  
  width = 18,
  height = 10,
  bg = "white"
)


# ======================================================================
# ======================================================================
# PART IX
# OPTIONAL COI TURNOVER FIGURE
#
# Generated for transparency.
# We can decide later whether to include it.
# ======================================================================
# ======================================================================


# ======================================================================
# 26. COI SITE-SPECIFIC TURNOVER
# ======================================================================

COI_turnover_plot_data <-
  site_lag_summary %>%
  
  filter(
    marker == "COI"
  )


optional_COI_turnover <-
  ggplot(
    COI_turnover_plot_data,
    aes(
      x = lag_years,
      y = mean_turnover,
      colour = regional_status,
      group = regional_status
    )
  ) +
  
  geom_line(
    linewidth = 0.8,
    alpha = 0.75
  ) +
  
  geom_point(
    aes(
      size =
        n_comparisons
    ),
    alpha = 0.90
  ) +
  
  facet_wrap(
    ~ site,
    scales = "free",
    ncol = 4
  ) +
  
  scale_colour_manual(
    values =
      status_colors_final,
    drop = FALSE
  ) +
  
  scale_size_continuous(
    name =
      "Species comparisons",
    range =
      c(
        1.5,
        6
      )
  ) +
  
  scale_y_continuous(
    limits =
      c(
        0,
        1
      ),
    breaks =
      seq(
        0,
        1,
        by = 0.25
      )
  ) +
  
  labs(
    
    x =
      "Temporal separation (years)",
    
    y =
      "Mean haplotype compositional turnover\n(Jaccard dissimilarity)",
    
    colour =
      NULL,
    
    title =
      "COI temporal haplotype turnover"
  ) +
  
  theme_bw(
    base_size = 18
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    panel.grid.minor =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    plot.title =
      element_text(
        face = "bold"
      )
  )



ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "OPTIONAL_SUPP_COI_TEMPORAL_TURNOVER.png"
    ),
  
  plot =
    optional_COI_turnover,
  
  width = 26,
  height = 9,
  dpi = 600,
  bg = "white"
)

# SVG vector copy (DPI does not apply to vector output)
ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "OPTIONAL_SUPP_COI_TEMPORAL_TURNOVER.svg"
    ),
  
  plot =
    optional_COI_turnover,
  
  width = 26,
  height = 9,
  bg = "white"
)


ggsave(
  
  filename =
    file.path(
      final_manuscript_dir,
      "OPTIONAL_SUPP_COI_TEMPORAL_TURNOVER.pdf"
    ),
  
  plot =
    optional_COI_turnover,
  
  width = 26,
  height = 9,
  bg = "white"
)


# ======================================================================
# ======================================================================
# PART X
# SUPPLEMENTARY TABLE FOR PAPER
# ======================================================================
# ======================================================================


# ======================================================================
# 27. FORMATTED REPRESENTATION TABLE
#
# This is the table I would actually consider as Table Sx.
# ======================================================================

supplementary_table <-
  representation_table %>%
  
  arrange(
    marker,
    site,
    regional_status
  ) %>%
  
  transmute(
    
    Marker =
      marker,
    
    Site =
      site,
    
    Status =
      as.character(
        regional_status
      ),
    
    `Species (n)` =
      n_species,
    
    `Temporal comparisons (n)` =
      n_temporal_comparisons,
    
    `Period pairs (n)` =
      n_period_pairs,
    
    `Minimum temporal lag (years)` =
      min_lag_years,
    
    `Maximum temporal lag (years)` =
      max_lag_years,
    
    `Mean Jaccard turnover` =
      round(
        mean_turnover,
        3
      ),
    
    `Median Jaccard turnover` =
      round(
        median_turnover,
        3
      ),
    
    `SD` =
      round(
        sd_turnover,
        3
      ),
    
    `Proportion turnover = 0` =
      round(
        proportion_zero,
        3
      ),
    
    `Proportion turnover = 1` =
      round(
        proportion_one,
        3
      )
  )


# ======================================================================
# 28. FORMATTED SITE-CONTRAST TABLE
# ======================================================================

site_contrast_table <-
  bootstrap_site_results %>%
  
  arrange(
    site
  ) %>%
  
  transmute(
    
    Site =
      site,
    
    `NAT species (n)` =
      n_species_NAT,
    
    `NIS species (n)` =
      n_species_NIS,
    
    `NIS - NAT turnover difference` =
      round(
        estimate,
        3
      ),
    
    `95% CI lower` =
      round(
        CI_low,
        3
      ),
    
    `95% CI upper` =
      round(
        CI_high,
        3
      ),
    
    `Bootstrap Pr(NIS > NAT)` =
      round(
        bootstrap_probability_positive,
        3
      ),
    
    `Interpretation` =
      supported_direction
  )


# ======================================================================
# 29. GLOBAL EFFECT TABLE
# ======================================================================

global_effect_table <-
  global_bootstrap_result %>%
  
  transmute(
    
    Analysis =
      "Site-averaged NAT-NIS difference",
    
    `NIS - NAT turnover difference` =
      round(
        estimate,
        3
      ),
    
    `95% CI lower` =
      round(
        CI_low,
        3
      ),
    
    `95% CI upper` =
      round(
        CI_high,
        3
      ),
    
    `Bootstrap Pr(NIS > NAT)` =
      round(
        bootstrap_probability_positive,
        3
      )
  )


# ======================================================================
# 30. Q1 STATISTICS TABLE
# ======================================================================

Q1_statistics_18S <-
  results_18S$within_bin %>%
  
  mutate(
    marker = "18S"
  )


Q1_statistics_COI <-
  results_COI$within_bin %>%
  
  mutate(
    marker = "COI"
  )


Q1_statistics_all <-
  bind_rows(
    Q1_statistics_COI,
    Q1_statistics_18S
  )


# ======================================================================
# ======================================================================
# PART XI
# EXPORT FINAL WORKBOOK
# ======================================================================
# ======================================================================


writexl::write_xlsx(
  
  list(
    
    Table_S_representation =
      supplementary_table,
    
    Table_S_site_contrasts =
      site_contrast_table,
    
    Global_effect =
      global_effect_table,
    
    Q1_richness_statistics =
      Q1_statistics_all,
    
    Species_temporal_slopes =
      species_temporal_slopes_18S,
    
    Temporal_slope_summary =
      temporal_slope_summary,
    
    Site_lag_descriptive =
      site_lag_summary,
    
    Raw_turnover_18S =
      turnover_18S,
    
    Raw_turnover_COI =
      turnover_COI
    
  ),
  
  path =
    file.path(
      final_manuscript_dir,
      "FINAL_MANUSCRIPT_STATISTICS_AND_TABLES.xlsx"
    )
)


# ======================================================================
# 31. SAVE R OBJECTS FOR FULL REPRODUCIBILITY
# ======================================================================

saveRDS(
  bootstrap_site_results,
  file.path(
    final_manuscript_dir,
    "18S_SITE_BOOTSTRAP_RESULTS.rds"
  )
)


saveRDS(
  global_bootstrap_result,
  file.path(
    final_manuscript_dir,
    "18S_GLOBAL_BOOTSTRAP_RESULT.rds"
  )
)


# ======================================================================
# 32. FINAL SUMMARY
# ======================================================================

cat(
  "\n\n============================================================\n",
  "FINAL MANUSCRIPT ANALYSIS COMPLETED\n",
  "============================================================\n",
  
  "\nMAIN FIGURE:\n",
  "  MAIN_18S_TEMPORAL_HAPLOTYPE_TURNOVER.png/pdf\n",
  
  "\nSUPPLEMENTARY Q1 FIGURE:\n",
  "  SUPP_Q1_HAPLOTYPE_RICHNESS_COI_18S.png/pdf\n",
  
  "\nOPTIONAL COI TURNOVER:\n",
  "  OPTIONAL_SUPP_COI_TEMPORAL_TURNOVER.png/pdf\n",
  
  "\nFINAL TABLES + STATISTICS:\n",
  "  FINAL_MANUSCRIPT_STATISTICS_AND_TABLES.xlsx\n",
  
  "\nBOOTSTRAP:\n",
  "  ", N_BOOT, " replicates\n",
  
  "\n============================================================\n"
)










# ======================================================================
# ======================================================================
# FINAL MANUSCRIPT BLOCK
#
# TEMPORAL CHANGE IN WITHIN-SPECIES HAPLOTYPE COMPOSITION
#
# MAIN ECOLOGICAL QUESTIONS
#
# Q1. How much does detected within-species haplotype composition
#     change among sediment horizons?
#
# Q2. What constitutes that temporal compositional change?
#     - ASVs redetected in both periods
#     - ASVs detected only in the earlier period
#     - ASVs detected only in the later period
#
# Q3. Does compositional turnover tend to increase with temporal
#     separation?
#
# SECONDARY QUESTION
#
# Do these descriptive patterns differ between NAT and NIS?
#
#
# OUTPUTS
#
# MAIN:
#   MAIN_A_18S_TEMPORAL_TURNOVER
#   MAIN_B1_18S_COMPONENTS_PROPORTIONS
#   MAIN_B2_18S_COMPONENTS_COUNTS
#   MAIN_B3_18S_REDETECTION_TURNOVER
#   MAIN_AB1 / MAIN_AB2 / MAIN_AB3
#
# SUPPLEMENTARY:
#   SUPP_Q1_COI_18S_RICHNESS
#
# TABLES:
#   Table_S1_representation
#   Table_S2_composition_components
#   Table_S3_temporal_slopes
#   Q1 statistics
#
#
# TERMINOLOGY
#
# shared / redetected:
#   ASV detected in both periods
#
# not redetected later:
#   ASV detected in earlier period but not detected in later period
#
# newly detected:
#   ASV detected in later period but not detected in earlier period
#
# These terms DO NOT imply biological survival, extinction or arrival.
# ======================================================================
# ======================================================================


# ======================================================================
# 0. PACKAGES
# ======================================================================

final_packages <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "patchwork",
  "writexl",
  "scales",
  "svglite"
)

missing_final_packages <-
  final_packages[
    !final_packages %in%
      rownames(installed.packages())
  ]

if (length(missing_final_packages) > 0) {
  install.packages(missing_final_packages)
}

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(writexl)
library(scales)
library(svglite)


# ======================================================================
# 1. SETTINGS
# ======================================================================

MIN_SPECIES_DESCRIPTIVE <- 3

status_colors_final <- c(
  "NAT" = "#8FC793",
  "NIS" = "#CC3366"
)

component_colors <- c(
  "Redetected in both periods" = "#4D4D4D",
  "Not redetected later" = "#BDBDBD",
  "Newly detected later" = "#737373"
)


# ======================================================================
# 2. OUTPUT DIRECTORY
# ======================================================================

final_dir <-
  file.path(
    output_dir,
    "FINAL_MANUSCRIPT_COMPOSITION"
  )

dir.create(
  final_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# 3. CHECK OBJECTS
# ======================================================================

if (!exists("results_18S")) {
  stop("results_18S does not exist.")
}

if (!exists("results_COI")) {
  stop("results_COI does not exist.")
}

if (is.null(results_18S$asv_turnover_species)) {
  stop(
    paste0(
      "results_18S$asv_turnover_species does not exist.\n",
      "Add asv_turnover_species to the return() list of analyse_marker()."
    )
  )
}

if (is.null(results_COI$asv_turnover_species)) {
  stop(
    paste0(
      "results_COI$asv_turnover_species does not exist.\n",
      "Add asv_turnover_species to the return() list of analyse_marker()."
    )
  )
}


# ======================================================================
# 4. PREPARE TURNOVER DATA
# ======================================================================

prepare_turnover <- function(
    x,
    marker_name
) {
  
  x %>%
    
    filter(
      !is.na(site),
      !is.na(species),
      !is.na(regional_status),
      !is.na(period_1),
      !is.na(period_2),
      !is.na(lag_years),
      !is.na(jaccard_turnover),
      regional_status %in% c(
        "NAT",
        "NIS"
      )
    ) %>%
    
    mutate(
      
      marker =
        marker_name,
      
      site =
        as.character(site),
      
      species =
        as.character(species),
      
      regional_status =
        factor(
          regional_status,
          levels = c(
            "NAT",
            "NIS"
          )
        ),
      
      lag_years =
        as.numeric(lag_years),
      
      lag_century =
        lag_years / 100,
      
      jaccard_turnover =
        as.numeric(jaccard_turnover)
    )
}


turnover_18S <-
  prepare_turnover(
    results_18S$asv_turnover_species,
    "18S"
  )


turnover_COI <-
  prepare_turnover(
    results_COI$asv_turnover_species,
    "COI"
  )


turnover_all <-
  bind_rows(
    turnover_18S,
    turnover_COI
  )


# ======================================================================
# 5. CHECK Q5 COMPONENT COLUMNS
#
# Q5 terminology in the existing results:
#
# n_shared = ASVs detected in both periods
# n_lost   = ASVs detected in period 1 but not redetected in period 2
# n_gained = ASVs detected in period 2 but not detected in period 1
#
# IMPORTANT:
# "lost" and "gained" are legacy variable names only.
# In figures/tables/manuscript we use:
#
#   Redetected in both periods
#   Not redetected later
#   Newly detected later
#
# because sedaDNA detection does not demonstrate biological
# loss/extinction or gain/arrival.
# ======================================================================

candidate_shared <- c(
  "n_shared",
  "shared_asvs",
  "n_shared_asvs"
)

candidate_earlier <- c(
  "n_lost",
  "n_earlier_only",
  "n_not_redetected",
  "earlier_only_asvs"
)

candidate_later <- c(
  "n_gained",
  "n_later_only",
  "n_new",
  "later_only_asvs",
  "n_new_asvs"
)


find_column <- function(
    dat,
    candidates,
    label
) {
  
  found <-
    candidates[
      candidates %in% names(dat)
    ]
  
  if (length(found) == 0) {
    
    stop(
      paste0(
        "\nCould not find the column for: ",
        label,
        "\n\nAvailable columns are:\n",
        paste(
          names(dat),
          collapse = "\n"
        )
      )
    )
  }
  
  found[1]
}


shared_col <-
  find_column(
    turnover_18S,
    candidate_shared,
    "shared/redetected ASVs"
  )


earlier_col <-
  find_column(
    turnover_18S,
    candidate_earlier,
    "earlier-only / not-redetected-later ASVs"
  )


later_col <-
  find_column(
    turnover_18S,
    candidate_later,
    "later-only / newly-detected-later ASVs"
  )


cat(
  "\n============================================================\n",
  "Q5 COMPONENT COLUMNS DETECTED\n",
  "============================================================\n",
  "Shared / redetected:       ", shared_col, "\n",
  "Not redetected later:      ", earlier_col, "\n",
  "Newly detected later:      ", later_col, "\n",
  "============================================================\n"
)

# ======================================================================
# 6. STANDARDISE Q5 COMPONENT NAMES
# ======================================================================

standardise_components <- function(
    dat
) {
  
  dat %>%
    
    mutate(
      
      n_shared_final =
        as.numeric(
          .data[[shared_col]]
        ),
      
      n_earlier_only_final =
        as.numeric(
          .data[[earlier_col]]
        ),
      
      n_later_only_final =
        as.numeric(
          .data[[later_col]]
        )
    )
}


turnover_18S <-
  standardise_components(
    turnover_18S
  )


# ======================================================================
# 7. SANITY CHECK COMPONENTS
# ======================================================================

if (
  any(
    turnover_18S$n_shared_final < 0 |
    turnover_18S$n_earlier_only_final < 0 |
    turnover_18S$n_later_only_final < 0,
    na.rm = TRUE
  )
) {
  
  stop(
    "Negative ASV component counts detected."
  )
}


# ======================================================================
# 8. CALCULATE COMPOSITIONAL COMPONENT FRACTIONS
#
# Union:
# shared + earlier-only + later-only
#
# These three fractions sum to 1 for each species temporal comparison.
# ======================================================================

turnover_18S <-
  turnover_18S %>%
  
  mutate(
    
    n_union_final =
      n_shared_final +
      n_earlier_only_final +
      n_later_only_final,
    
    fraction_shared =
      if_else(
        n_union_final > 0,
        n_shared_final /
          n_union_final,
        NA_real_
      ),
    
    fraction_not_redetected =
      if_else(
        n_union_final > 0,
        n_earlier_only_final /
          n_union_final,
        NA_real_
      ),
    
    fraction_newly_detected =
      if_else(
        n_union_final > 0,
        n_later_only_final /
          n_union_final,
        NA_real_
      ),
    
    fraction_turnover =
      fraction_not_redetected +
      fraction_newly_detected
  )


# ======================================================================
# 9. VERIFY RELATION WITH JACCARD
#
# Jaccard dissimilarity should equal:
#
# (earlier-only + later-only) /
# (shared + earlier-only + later-only)
# ======================================================================

turnover_18S <-
  turnover_18S %>%
  
  mutate(
    
    jaccard_reconstructed =
      fraction_turnover,
    
    jaccard_difference =
      abs(
        jaccard_turnover -
          jaccard_reconstructed
      )
  )


max_jaccard_difference <-
  max(
    turnover_18S$jaccard_difference,
    na.rm = TRUE
  )


cat(
  "\nMaximum difference between stored and reconstructed Jaccard = ",
  max_jaccard_difference,
  "\n"
)


if (
  is.finite(max_jaccard_difference) &&
  max_jaccard_difference > 1e-6
) {
  
  warning(
    paste0(
      "Stored Jaccard and reconstructed Jaccard differ by up to ",
      signif(
        max_jaccard_difference,
        4
      ),
      ". Check Q5 definitions before manuscript interpretation."
    )
  )
}


# ======================================================================
# ======================================================================
# PART I
# MAIN PANEL A
# TEMPORAL TURNOVER
# ======================================================================
# ======================================================================


# ======================================================================
# 10. SITE × STATUS × LAG SUMMARY
# ======================================================================

site_lag_18S <-
  turnover_18S %>%
  
  group_by(
    site,
    regional_status,
    lag_years
  ) %>%
  
  summarise(
    
    n_comparisons =
      n(),
    
    n_species =
      n_distinct(species),
    
    mean_turnover =
      mean(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    median_turnover =
      median(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ======================================================================
# 11. PANEL A
# ======================================================================

panel_A <-
  ggplot(
    site_lag_18S,
    aes(
      x = lag_years,
      y = mean_turnover,
      colour = regional_status,
      group = regional_status
    )
  ) +
  
  geom_line(
    linewidth = 0.8,
    alpha = 0.80
  ) +
  
  geom_point(
    aes(
      size = n_comparisons
    ),
    alpha = 0.90
  ) +
  
  facet_wrap(
    ~ site,
    scales = "free",
    ncol = 4
  ) +
  
  scale_colour_manual(
    values = status_colors_final,
    drop = FALSE
  ) +
  
  scale_size_continuous(
    name = "Species comparisons",
    range = c(
      1.5,
      6
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.25
    )
  ) +
  
  labs(
    
    x =
      "Temporal separation (years)",
    
    y =
      "Mean haplotype compositional turnover\n(Jaccard dissimilarity)",
    
    colour =
      NULL,
    
    title =
      "A  Temporal haplotype turnover within sedimentary archives"
  ) +
  
  theme_bw(
    base_size = 18
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    panel.grid.minor =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    plot.title =
      element_text(
        face = "bold"
      )
  )


# ======================================================================
# ======================================================================
# PART II
# MAIN PANEL B
#
# TEMPORAL COMPONENTS OF HAPLOTYPE COMPOSITIONAL DIFFERENCES
#
# PURPOSE:
#
# Panel A asks:
#   How does compositional turnover vary with temporal separation?
#
# Panel B asks:
#   For explicit consecutive temporal comparisons within each archive,
#   what fraction of the detected ASV union was:
#
#   1. Redetected in both periods
#   2. Detected only in the earlier period
#   3. Detected only in the later period
#
# IMPORTANT:
#
# We do NOT interpret "earlier-only" as biological loss/extinction.
# We do NOT interpret "later-only" as biological gain/arrival.
#
# These are detection categories only.
#
# Biological unit:
#   species
#
# For each site × status × temporal transition:
#   each species contributes once.
#
# ======================================================================
# ======================================================================


# ======================================================================
# 12. IDENTIFY CONSECUTIVE AVAILABLE TEMPORAL COMPARISONS
#
# IMPORTANT:
#
# We do NOT simply use lag_years == 50.
#
# Instead, for each species within each site, we reconstruct the ordered
# periods in which that species was detected and retain comparisons
# between consecutive detected periods.
#
# Example:
#
# Species detected in:
#   1450, 1500, 1600
#
# Consecutive comparisons:
#   1450 -> 1500
#   1500 -> 1600
#
# We do NOT retain:
#   1450 -> 1600
#
# This avoids showing every possible overlapping temporal comparison in
# Panel B and gives the panel a direct chronological interpretation.
# ======================================================================


species_period_sequences <-
  turnover_18S %>%
  
  select(
    site,
    species,
    regional_status,
    period_1,
    period_2
  ) %>%
  
  pivot_longer(
    cols = c(
      period_1,
      period_2
    ),
    names_to = "period_position",
    values_to = "period"
  ) %>%
  
  select(
    site,
    species,
    regional_status,
    period
  ) %>%
  
  distinct() %>%
  
  arrange(
    site,
    species,
    period
  ) %>%
  
  group_by(
    site,
    species,
    regional_status
  ) %>%
  
  mutate(
    
    next_period =
      lead(period)
    
  ) %>%
  
  filter(
    !is.na(next_period)
  ) %>%
  
  transmute(
    
    site,
    species,
    regional_status,
    
    period_1 =
      period,
    
    period_2 =
      next_period
  ) %>%
  
  ungroup()


# ======================================================================
# 13. RETAIN ONLY THOSE CONSECUTIVE COMPARISONS
# ======================================================================

B_temporal_raw <-
  turnover_18S %>%
  
  inner_join(
    
    species_period_sequences,
    
    by = c(
      "site",
      "species",
      "regional_status",
      "period_1",
      "period_2"
    )
  ) %>%
  
  filter(
    
    !is.na(
      n_shared_final
    ),
    
    !is.na(
      n_earlier_only_final
    ),
    
    !is.na(
      n_later_only_final
    ),
    
    !is.na(
      fraction_shared
    ),
    
    !is.na(
      fraction_not_redetected
    ),
    
    !is.na(
      fraction_newly_detected
    )
  ) %>%
  
  mutate(
    
    transition =
      paste0(
        period_1,
        "\u2192",
        period_2
      )
  )


# ======================================================================
# 14. CHECK THAT EACH SPECIES × TRANSITION OCCURS ONLY ONCE
# ======================================================================

B_duplicates <-
  B_temporal_raw %>%
  
  count(
    site,
    species,
    regional_status,
    period_1,
    period_2,
    name = "n"
  ) %>%
  
  filter(
    n > 1
  )


if (
  nrow(
    B_duplicates
  ) > 0
) {
  
  warning(
    paste0(
      "\nPanel B: ",
      nrow(B_duplicates),
      " duplicated species × temporal transitions were detected.\n",
      "They will be averaged at species level before plotting."
    )
  )
}


# ======================================================================
# 15. SPECIES-LEVEL SUMMARY FOR EACH EXPLICIT TRANSITION
#
# This protects the biological unit if duplicates exist.
# ======================================================================

B_temporal_species <-
  B_temporal_raw %>%
  
  group_by(
    site,
    species,
    regional_status,
    period_1,
    period_2,
    transition
  ) %>%
  
  summarise(
    
    lag_years =
      mean(
        lag_years,
        na.rm = TRUE
      ),
    
    fraction_shared =
      mean(
        fraction_shared,
        na.rm = TRUE
      ),
    
    fraction_earlier_only =
      mean(
        fraction_not_redetected,
        na.rm = TRUE
      ),
    
    fraction_later_only =
      mean(
        fraction_newly_detected,
        na.rm = TRUE
      ),
    
    jaccard_turnover =
      mean(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ======================================================================
# 16. SANITY CHECK
#
# For every species × temporal comparison:
#
# shared + earlier-only + later-only = 1
#
# and:
#
# earlier-only + later-only = Jaccard turnover
# ======================================================================

B_temporal_species <-
  B_temporal_species %>%
  
  mutate(
    
    component_sum =
      fraction_shared +
      fraction_earlier_only +
      fraction_later_only,
    
    reconstructed_turnover =
      fraction_earlier_only +
      fraction_later_only,
    
    turnover_difference =
      abs(
        jaccard_turnover -
          reconstructed_turnover
      )
  )


max_component_error <-
  max(
    abs(
      B_temporal_species$component_sum - 1
    ),
    na.rm = TRUE
  )


max_turnover_error <-
  max(
    B_temporal_species$turnover_difference,
    na.rm = TRUE
  )


cat(
  "\n\n============================================================\n",
  "TEMPORAL PANEL B SANITY CHECK\n",
  "============================================================\n",
  "Maximum component-sum error: ",
  max_component_error,
  "\n",
  "Maximum Jaccard reconstruction error: ",
  max_turnover_error,
  "\n",
  "============================================================\n"
)


if (
  is.finite(max_component_error) &&
  max_component_error > 1e-6
) {
  
  warning(
    "Panel B component proportions do not sum to 1."
  )
}


if (
  is.finite(max_turnover_error) &&
  max_turnover_error > 1e-6
) {
  
  warning(
    "Panel B components do not reconstruct Jaccard turnover."
  )
}


# ======================================================================
# 17. SITE × STATUS × EXPLICIT TEMPORAL TRANSITION SUMMARY
#
# Every species receives equal weight within each transition.
#
# n_species tells us exactly how many species support each bar.
# ======================================================================
B_temporal_summary <-
  B_temporal_species %>%
  
  group_by(
    site,
    regional_status,
    period_1,
    period_2,
    transition
  ) %>%
  
  summarise(
    
    n_species =
      n_distinct(species),
    
    mean_shared =
      mean(
        fraction_shared,
        na.rm = TRUE
      ),
    
    mean_earlier_only =
      mean(
        fraction_earlier_only,
        na.rm = TRUE
      ),
    
    mean_later_only =
      mean(
        fraction_later_only,
        na.rm = TRUE
      ),
    
    mean_turnover =
      mean(
        jaccard_turnover,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    lag_years =
      period_2 - period_1
  ) %>%
  
  arrange(
    site,
    period_1,
    period_2,
    regional_status
  )

# ======================================================================
# 18. REPRESENTATION TABLE
#
# This is essential.
#
# We DO NOT automatically remove transitions with few species.
# Instead, we retain everything transparently and export the number of
# species supporting every bar.
#
# The manuscript figure can later be filtered if necessary, but only
# after examining this table.
# ======================================================================

B_temporal_representation <-
  B_temporal_summary %>%
  
  select(
    site,
    regional_status,
    period_1,
    period_2,
    transition,
    lag_years,
    n_species
  ) %>%
  
  arrange(
    site,
    period_1,
    period_2,
    regional_status
  )


cat(
  "\n\n============================================================\n",
  "SPECIES SUPPORT FOR EACH PANEL B TRANSITION\n",
  "============================================================\n"
)


print(
  B_temporal_representation,
  n = Inf
)


# ======================================================================
# 19. PREPARE LONG DATA FOR PLOTTING
#
# Terminology is deliberately symmetrical:
#
# Redetected in both periods
# Detected only in earlier period
# Detected only in later period
#
# No biological loss/gain is implied.
# ======================================================================

B_temporal_plot_data <-
  B_temporal_summary %>%
  
  select(
    site,
    regional_status,
    period_1,
    period_2,
    transition,
    lag_years,
    n_species,
    mean_shared,
    mean_earlier_only,
    mean_later_only
  ) %>%
  
  pivot_longer(
    
    cols = c(
      mean_shared,
      mean_earlier_only,
      mean_later_only
    ),
    
    names_to =
      "component",
    
    values_to =
      "proportion"
  ) %>%
  
  mutate(
    
    component =
      recode(
        
        component,
        
        mean_shared =
          "Redetected in both periods",
        
        mean_earlier_only =
          "Detected only in earlier period",
        
        mean_later_only =
          "Detected only in later period"
      ),
    
    component =
      factor(
        component,
        levels = c(
          "Redetected in both periods",
          "Detected only in earlier period",
          "Detected only in later period"
        )
      ),
    
    regional_status =
      factor(
        regional_status,
        levels = c(
          "NAT",
          "NIS"
        )
      )
  )


# ======================================================================
# 20. UPDATE COMPONENT COLOURS
#
# Keep the same neutral palette.
# ======================================================================

component_colors_temporal <- c(
  
  "Redetected in both periods" =
    "#4D4D4D",
  
  "Detected only in earlier period" =
    "#BDBDBD",
  
  "Detected only in later period" =
    "#737373"
)


# ======================================================================
# 21. CREATE A UNIQUE NUMERIC POSITION FOR EVERY SITE TRANSITION
#
# We use period_2 as the chronological x position.
#
# Conceptually:
# the bar located at 1500 describes the transition ending at 1500.
#
# The actual label explicitly gives:
# T1 -> T2
#
# Therefore there is no ambiguity about which periods are compared.
# ======================================================================

B_temporal_plot_data <-
  B_temporal_plot_data %>%
  
  mutate(
    
    transition_factor =
      factor(
        transition,
        levels =
          unique(
            transition[
              order(
                period_1,
                period_2
              )
            ]
          )
      )
  )


# ======================================================================
# 22. PANEL B — TEMPORAL COMPONENTS
#
# We facet:
#
# rows    = NAT / NIS
# columns = sites
#
# This avoids trying to encode status and component with the same fill
# aesthetic.
#
# Every stacked bar sums to 100%.
# ======================================================================

panel_B1_temporal <-
  ggplot(
    B_temporal_plot_data,
    aes(
      x = transition,
      y = proportion,
      fill = component
    )
  ) +
  
  geom_col(
    width = 0.90
  ) +
  
  facet_grid(
    regional_status ~ site,
    scales = "free",
    space = "free"
  ) +
  
  scale_y_continuous(
    
    limits = c(
      0,
      1
    ),
    
    breaks = seq(
      0,
      1,
      0.25
    ),
    
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    
    expand =
      expansion(
        mult = c(
          0,
          0.02
        )
      )
  ) +
  
  scale_fill_manual(
    
    values =
      component_colors_temporal,
    
    drop =
      FALSE
  ) +
  
  labs(
    
    x =
      "Consecutive temporal comparison",
    
    y =
      "Mean proportion of detected ASV union",
    
    fill =
      NULL,
    
    title =
      "B  Components of temporal differences in detected haplotype composition"
  ) +
  
  theme_bw(
    base_size = 18
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    legend.box =
      "horizontal",
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.x =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 16
      ),
    
    axis.title =
      element_text(
        size = 18
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 17
      )
  )


# ======================================================================
# 23. OPTIONAL VERSION WITH n SPECIES ABOVE EACH BAR
#
# IMPORTANT:
#
# This version is generated for inspection.
# I would decide whether to use it in the manuscript only after seeing
# how crowded it becomes.
# ======================================================================

B_n_labels <-
  B_temporal_summary %>%
  
  mutate(
    
    label =
      paste0(
        "n=",
        n_species
      ),
    
    regional_status =
      factor(
        regional_status,
        levels = c(
          "NAT",
          "NIS"
        )
      )
  )


panel_B1_temporal_with_n <-
  ggplot(
    B_temporal_plot_data,
    aes(
      x = transition,
      y = proportion,
      fill = component
    )
  ) +
  
  geom_col(
    width = 0.90
  ) +
  
  geom_text(
    
    data =
      B_n_labels,
    
    aes(
      x = transition,
      y = 1.025,
      label = label
    ),
    
    inherit.aes =
      FALSE,
    
    size =
      2.3,
    
    vjust =
      0
  ) +
  
  facet_grid(
    regional_status ~ site,
    scales = "free",
    space = "free"
  ) +
  
  scale_y_continuous(
    
    limits = c(
      0,
      1.10
    ),
    
    breaks = seq(
      0,
      1,
      0.25
    ),
    
    labels = function(x) {
      
      ifelse(
        x <= 1,
        scales::percent(
          x,
          accuracy = 1
        ),
        ""
      )
    },
    
    expand =
      expansion(
        mult = c(
          0,
          0.01
        )
      )
  ) +
  
  scale_fill_manual(
    
    values =
      component_colors_temporal,
    
    drop =
      FALSE
  ) +
  
  labs(
    
    x =
      "Consecutive temporal comparison",
    
    y =
      "Mean proportion of detected ASV union",
    
    fill =
      NULL,
    
    title =
      "B  Components of temporal differences in detected haplotype composition"
  ) +
  
  theme_bw(
    base_size = 18
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.x =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 16
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 17
      )
  )


# ======================================================================
# 24. BUILD FINAL MAIN FIGURE
#
# Panel A:
# turnover versus temporal separation using all available comparisons.
#
# Panel B:
# composition of explicit consecutive comparisons through the archive.
#
# Therefore the panels are complementary rather than redundant.
# ======================================================================

main_AB1_temporal <-
  panel_A /
  panel_B1_temporal +
  
  patchwork::plot_layout(
    heights = c(
      1.15,
      1
    )
  )


main_AB1_temporal_with_n <-
  panel_A /
  panel_B1_temporal_with_n +
  
  patchwork::plot_layout(
    heights = c(
      1.15,
      1
    )
  )


# ======================================================================
# 25. SAVE TEMPORAL PANEL B ALONE
# ======================================================================

ggsave(
  
  file.path(
    final_dir,
    "MAIN_B1_18S_TEMPORAL_COMPONENTS.png"
  ),
  
  panel_B1_temporal,
  
  width = 29,
  height = 9,
  dpi = 600,
  bg = "white"
)

# SVG vector copy (DPI does not apply to vector output)
ggsave(
  
  file.path(
    final_dir,
    "MAIN_B1_18S_TEMPORAL_COMPONENTS.svg"
  ),
  
  panel_B1_temporal,
  
  width = 29,
  height = 9,
  bg = "white"
)


ggsave(
  
  file.path(
    final_dir,
    "MAIN_B1_18S_TEMPORAL_COMPONENTS.pdf"
  ),
  
  panel_B1_temporal,
  
  width = 29,
  height = 9,
  bg = "white"
)


# ======================================================================
# 26. SAVE VERSION WITH SAMPLE-SIZE LABELS
# ======================================================================

ggsave(
  
  file.path(
    final_dir,
    "MAIN_B1_18S_TEMPORAL_COMPONENTS_WITH_N.png"
  ),
  
  panel_B1_temporal_with_n,
  
  width = 18,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

# SVG vector copy (DPI does not apply to vector output)
ggsave(
  
  file.path(
    final_dir,
    "MAIN_B1_18S_TEMPORAL_COMPONENTS_WITH_N.svg"
  ),
  
  panel_B1_temporal_with_n,
  
  width = 18,
  height = 8.5,
  bg = "white"
)


# ======================================================================
# 27. SAVE FINAL A + TEMPORAL B
# ======================================================================

ggsave(
  
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS.png"
  ),
  
  main_AB1_temporal,
  
  width = 18,
  height = 15,
  dpi = 600,
  bg = "white"
)

# SVG vector copy (DPI does not apply to vector output)
ggsave(
  
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS.svg"
  ),
  
  main_AB1_temporal,
  
  width = 18,
  height = 15,
  bg = "white"
)


ggsave(
  
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS.pdf"
  ),
  
  main_AB1_temporal,
  
  width = 18,
  height = 15,
  bg = "white"
)


# ======================================================================
# 28. SAVE FINAL A + TEMPORAL B WITH n
#
# Inspection version.
# ======================================================================

ggsave(
  
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS_WITH_N.png"
  ),
  
  main_AB1_temporal_with_n,
  
  width = 18,
  height = 15,
  dpi = 600,
  bg = "white"
)

# SVG vector copy (DPI does not apply to vector output)
ggsave(
  
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS_WITH_N.svg"
  ),
  
  main_AB1_temporal_with_n,
  
  width = 18,
  height = 15,
  bg = "white"
)


# ======================================================================
# 29. TEMPORAL COMPONENT TABLE FOR SUPPLEMENTARY MATERIAL
#
# This table corresponds DIRECTLY to Panel B.
# ======================================================================

Table_S2_temporal_components <-
  B_temporal_summary %>%
  
  arrange(
    site,
    period_1,
    period_2,
    regional_status
  ) %>%
  
  transmute(
    
    Marker =
      "18S",
    
    Site =
      site,
    
    Status =
      as.character(
        regional_status
      ),
    
    `Earlier period` =
      period_1,
    
    `Later period` =
      period_2,
    
    `Temporal separation (years)` =
      lag_years,
    
    `Species (n)` =
      n_species,
    
    `Mean fraction redetected in both periods` =
      round(
        mean_shared,
        3
      ),
    
    `Mean fraction detected only in earlier period` =
      round(
        mean_earlier_only,
        3
      ),
    
    `Mean fraction detected only in later period` =
      round(
        mean_later_only,
        3
      ),
    
    `Mean Jaccard turnover` =
      round(
        mean_turnover,
        3
      )
  )


# ======================================================================
# 30. PRINT PANEL B TABLE
# ======================================================================

cat(
  "\n\n============================================================\n",
  "TEMPORAL COMPONENTS USED IN MAIN PANEL B\n",
  "============================================================\n"
)


print(
  Table_S2_temporal_components,
  n = Inf
)


# ======================================================================
# 31. EXPORT FINAL TEMPORAL PANEL-B TABLES
#
# These files DO NOT replace any previous output. They add the tables
# corresponding directly to the final temporal Panel B.
# ======================================================================

writexl::write_xlsx(
  list(
    Table_S2_temporal_components = Table_S2_temporal_components,
    Panel_B_species_support = B_temporal_representation,
    Panel_B_species_level = B_temporal_species
  ),
  path = file.path(
    final_dir,
    "FINAL_MAIN_PANEL_B_TEMPORAL_TABLES.xlsx"
  )
)

write.csv(
  Table_S2_temporal_components,
  file.path(
    final_dir,
    "Table_S2_temporal_components.csv"
  ),
  row.names = FALSE
)

write.csv(
  B_temporal_representation,
  file.path(
    final_dir,
    "Panel_B_species_support.csv"
  ),
  row.names = FALSE
)

cat(
  "\n\n============================================================\n",
  "FINAL MANUSCRIPT FIGURES/TABLES READY\n",
  "============================================================\n",
  "Main figure:\n",
  "  MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS.png/pdf\n\n",
  "Inspection version with n:\n",
  "  MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS_WITH_N.png\n\n",
  "Supplementary richness figure retained:\n",
  "  SUPP_Q1_HAPLOTYPE_RICHNESS_COI_18S.png/pdf\n\n",
  "Optional COI turnover figure retained.\n\n",
  "Final temporal Panel-B tables:\n",
  "  FINAL_MAIN_PANEL_B_TEMPORAL_TABLES.xlsx\n",
  "============================================================\n"
)


# ======================================================================
# ALTERNATIVE MAIN FIGURE
# PANEL B WITH STATUS-SPECIFIC NAT / NIS COLOUR GRADIENTS
#
# Within each status:
#   dark   = redetected in both periods
#   medium = detected only in earlier period
#   light  = detected only in later period
#
# This is an additional figure. The neutral-grey main figure is retained.
# ======================================================================

B_temporal_plot_status_colours <-
  B_temporal_plot_data %>%
  mutate(
    component_status = paste(
      as.character(regional_status),
      as.character(component),
      sep = " | "
    )
  )

component_status_colors <- c(
  "NAT | Redetected in both periods"       = "#2F6B3A",
  "NAT | Detected only in earlier period"  = "#78AD7D",
  "NAT | Detected only in later period"    = "#C8E1CA",
  "NIS | Redetected in both periods"       = "#8E0F3B",
  "NIS | Detected only in earlier period"  = "#D85B86",
  "NIS | Detected only in later period"    = "#F3BDD0"
)

component_status_labels <- c(
  "NAT | Redetected in both periods"       = "NAT — redetected in both periods",
  "NAT | Detected only in earlier period"  = "NAT — detected only in earlier period",
  "NAT | Detected only in later period"    = "NAT — detected only in later period",
  "NIS | Redetected in both periods"       = "NIS — redetected in both periods",
  "NIS | Detected only in earlier period"  = "NIS — detected only in earlier period",
  "NIS | Detected only in later period"    = "NIS — detected only in later period"
)

panel_B1_temporal_status_colours <-
  ggplot(
    B_temporal_plot_status_colours,
    aes(
      x = transition,
      y = proportion,
      fill = component_status
    )
  ) +
  geom_col(
    width = 0.90
  ) +
  facet_grid(
    regional_status ~ site,
    scales = "free",
    space = "free"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = component_status_colors,
    labels = component_status_labels,
    drop = FALSE
  ) +
  labs(
    x = "Consecutive temporal comparison",
    y = "Mean proportion of detected ASV union",
    fill = NULL,
    title = "B  Components of temporal differences in detected haplotype composition"
  ) +
  theme_bw(
    base_size = 20
  ) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 17
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 15,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 16,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 20,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 21
    )
  )

main_AB1_temporal_status_colours <-
  panel_A /
  panel_B1_temporal_status_colours +
  patchwork::plot_layout(
    heights = c(1.15, 1)
  )

ggsave(
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS_STATUS_COLOURS.png"
  ),
  main_AB1_temporal_status_colours,
  width = 29,
  height = 16,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS_STATUS_COLOURS.pdf"
  ),
  main_AB1_temporal_status_colours,
  width = 29,
  height = 16,
  bg = "white"
)

ggsave(
  file.path(
    final_dir,
    "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS_STATUS_COLOURS.svg"
  ),
  main_AB1_temporal_status_colours,
  width = 29,
  height = 16,
  bg = "white"
)

cat(
  "\nAdditional manuscript figure created:\n",
  "MAIN_FINAL_18S_TURNOVER_TEMPORAL_COMPONENTS_STATUS_COLOURS",
  " (.png, .pdf, .svg)\n"
)

