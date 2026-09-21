# ============================================================
# NIS / NAT-NIS / CRY TIMELINES BY SITE, KINGDOM AND MARKER
# ============================================================

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(janitor)
library(writexl)
library(grid)


# ============================================================
# 1. PATHS
# ============================================================

data_dir <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Data"

output_dir <- file.path(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Figures",
  "NIS_timelines"
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

file_18S_reps <- file.path(
  data_dir,
  "Selected_species_ASVs_18S_Reps_AbRel.xlsx"
)

file_COI_reps <- file.path(
  data_dir,
  "Selected_species_ASVs_COI_Reps_AbRel.xlsx"
)


# ============================================================
# 2. HELPERS
# ============================================================

clean_text <- function(x) {
  stringr::str_squish(
    as.character(x)
  )
}


clean_id <- function(x) {
  janitor::make_clean_names(
    as.character(x)
  )
}


find_first_col <- function(
    df,
    candidates,
    required = TRUE,
    label = "column"
) {
  
  found <- intersect(
    candidates,
    names(df)
  )
  
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


# ------------------------------------------------------------
# STATUS NORMALISATION
# ------------------------------------------------------------

normalise_status <- function(x) {
  
  x <- clean_text(x)
  x <- toupper(x)
  
  x <- stringr::str_replace_all(
    x,
    "\\s+",
    ""
  )
  
  dplyr::case_when(
    
    x == "NIS" ~ "NIS",
    
    x == "NAT/NIS" ~ "NAT/NIS",
    
    x == "NIS/NAT" ~ "NIS/NAT",
    
    x == "CRY" ~ "CRY",
    
    x %in% c(
      "CRY/NAT",
      "NAT/CRY"
    ) ~ "CRY/NAT",
    
    TRUE ~ x
  )
}


# ------------------------------------------------------------
# KINGDOM NORMALISATION
# ------------------------------------------------------------

normalise_kingdom <- function(x) {
  
  x <- clean_text(x)
  
  dplyr::case_when(
    
    x %in% c(
      "Plantae",
      "Viridiplantae"
    ) ~ "Plantae",
    
    x == "Metazoa" ~ "Metazoa",
    
    TRUE ~ x
  )
}


safe_filename <- function(x) {
  
  x <- stringr::str_replace_all(
    x,
    "[^A-Za-z0-9]+",
    "_"
  )
  
  x <- stringr::str_replace_all(
    x,
    "_+",
    "_"
  )
  
  stringr::str_remove_all(
    x,
    "^_|_$"
  )
}


# ============================================================
# 3. METADATA
# ============================================================

metadata <- readxl::read_excel(
  metadata_file
) %>%
  janitor::clean_names()


metadata_sample_col <- find_first_col(
  metadata,
  c(
    "sample",
    "sample_id",
    "sample_name",
    "sampleid",
    "sample_code",
    "muestra",
    "muestra_id",
    "nombre_muestra"
  ),
  required = TRUE,
  label = "sample ID column in metadata"
)


metadata_year_col <- find_first_col(
  metadata,
  c(
    "year",
    "sampling_year",
    "sample_year",
    "year_sampling",
    "ano",
    "anio"
  ),
  required = TRUE,
  label = "year column in metadata"
)


metadata_site_col <- find_first_col(
  metadata,
  c(
    "site",
    "site_code",
    "sampling_site",
    "location",
    "locality",
    "sitio"
  ),
  required = FALSE
)


metadata_key <- metadata %>%
  
  transmute(
    
    sample_original =
      clean_text(
        .data[[metadata_sample_col]]
      ),
    
    sample =
      clean_id(
        .data[[metadata_sample_col]]
      ),
    
    year =
      suppressWarnings(
        as.integer(
          .data[[metadata_year_col]]
        )
      ),
    
    site_metadata =
      if (!is.null(metadata_site_col)) {
        
        toupper(
          clean_text(
            .data[[metadata_site_col]]
          )
        )
        
      } else {
        
        NA_character_
      }
  ) %>%
  
  mutate(
    
    site_from_sample =
      toupper(
        stringr::str_extract(
          sample_original,
          "^[A-Za-z]+"
        )
      ),
    
    site =
      dplyr::coalesce(
        na_if(
          site_metadata,
          ""
        ),
        site_from_sample
      )
  ) %>%
  
  select(
    sample,
    sample_original,
    site,
    year
  ) %>%
  
  filter(
    !is.na(sample),
    sample != "",
    !is.na(year)
  ) %>%
  
  distinct()


metadata_conflicts <- metadata_key %>%
  
  distinct(
    sample,
    site,
    year
  ) %>%
  
  count(
    sample,
    name = "n_metadata_combinations"
  ) %>%
  
  filter(
    n_metadata_combinations > 1
  )


if (nrow(metadata_conflicts) > 0) {
  
  print(metadata_conflicts)
  
  stop(
    "Some biological samples have more than one site/year combination ",
    "in metadata_samples.xlsx."
  )
}


metadata_key <- metadata_key %>%
  
  distinct(
    sample,
    .keep_all = TRUE
  )


# ============================================================
# 4. STATUS TABLE
# ============================================================

status_species <- readxl::read_excel(
  status_file
) %>%
  janitor::clean_names()


status_species_col <- find_first_col(
  status_species,
  c(
    "especie",
    "species",
    "selected_species"
  ),
  required = TRUE,
  label = "species column in Status_selected_species.xlsx"
)


status_combined_col <- find_first_col(
  status_species,
  c(
    "status_combinado",
    "combined_status",
    "status"
  ),
  required = TRUE,
  label = "Status combinado column"
)


status_levels <- c(
  "NIS",
  "NAT/NIS",
  "NIS/NAT",
  "CRY",
  "CRY/NAT"
)

target_statuses <- status_levels


species_status <- status_species %>%
  
  transmute(
    
    species =
      clean_text(
        .data[[status_species_col]]
      ),
    
    status =
      normalise_status(
        .data[[status_combined_col]]
      )
  ) %>%
  
  filter(
    !is.na(species),
    species != "",
    status %in% target_statuses
  ) %>%
  
  mutate(
    status = factor(
      status,
      levels = status_levels
    )
  ) %>%
  
  distinct(
    species,
    .keep_all = TRUE
  )


cat("\nSTATUS CATEGORIES RETAINED:\n")

print(
  species_status %>%
    count(
      status,
      .drop = FALSE
    )
)


# ============================================================
# 5. READ ABUNDANCE TABLES
# ============================================================

read_abundance_table <- function(file) {
  
  message(
    "Reading: ",
    basename(file)
  )
  
  
  if (
    grepl(
      "\\.xlsx$",
      file,
      ignore.case = TRUE
    )
  ) {
    
    df <- readxl::read_excel(
      file
    )
    
  } else if (
    grepl(
      "\\.csv$",
      file,
      ignore.case = TRUE
    )
  ) {
    
    df <- readr::read_csv(
      file,
      show_col_types = FALSE
    )
    
  } else {
    
    stop(
      "Unsupported file format: ",
      basename(file)
    )
  }
  
  
  df <- df %>%
    janitor::clean_names()
  
  
  if (!"asv" %in% names(df)) {
    
    stop(
      "Column ASV not found in: ",
      basename(file)
    )
  }
  
  
  species_col <- find_first_col(
    df,
    c(
      "selected_species",
      "species"
    ),
    required = TRUE,
    label = paste0(
      "species column in ",
      basename(file)
    )
  )
  
  
  df %>%
    
    mutate(
      
      asv =
        clean_text(asv),
      
      species =
        clean_text(
          .data[[species_col]]
        )
    )
}


data_18S <- read_abundance_table(
  file_18S
)

data_COI <- read_abundance_table(
  file_COI
)

data_18S_reps <- read_abundance_table(
  file_18S_reps
)

data_COI_reps <- read_abundance_table(
  file_COI_reps
)


# ============================================================
# 6. IDENTIFY SAMPLE COLUMNS
# ============================================================

taxonomy_cols <- c(
  "asv",
  "marker",
  "pindent",
  "assignment_qual",
  "assignmentqual",
  "tax_id",
  "taxid",
  "domain",
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species",
  "selected_species",
  "status_combinado",
  "mediterraneo",
  "atlantico",
  "confianza",
  "justificacion",
  "referencia_principal",
  "referencia_secundaria",
  "nota_taxonomica_cautela",
  "assignment_level",
  "assignment_type",
  "best_pident",
  "best_qcov",
  "best_bitscore",
  "n_blast_hits",
  "n_species_hits",
  "blast_species_hits"
)


identify_collapsed_sample_cols <- function(
    df,
    metadata_key
) {
  
  matched <- intersect(
    metadata_key$sample,
    names(df)
  )
  
  
  if (length(matched) > 0) {
    return(matched)
  }
  
  
  candidates <- setdiff(
    names(df),
    taxonomy_cols
  )
  
  
  candidates[
    vapply(
      df[candidates],
      is.numeric,
      logical(1)
    )
  ]
}


sample_cols_18S <- identify_collapsed_sample_cols(
  data_18S,
  metadata_key
)

sample_cols_COI <- identify_collapsed_sample_cols(
  data_COI,
  metadata_key
)


if (length(sample_cols_18S) == 0) {
  stop(
    "No biological sample columns detected in 18S collapsed table."
  )
}


if (length(sample_cols_COI) == 0) {
  stop(
    "No biological sample columns detected in COI collapsed table."
  )
}


# ============================================================
# 7. MAP PCR REPLICATES TO BIOLOGICAL SAMPLES
# ============================================================

map_replicate_columns <- function(
    reps_df,
    biological_samples
) {
  
  replicate_candidates <- setdiff(
    names(reps_df),
    taxonomy_cols
  )
  
  
  replicate_candidates <- replicate_candidates[
    vapply(
      reps_df[replicate_candidates],
      is.numeric,
      logical(1)
    )
  ]
  
  
  map_one_rep <- function(rep_col) {
    
    possible <- biological_samples[
      
      rep_col == biological_samples |
        
        stringr::str_starts(
          rep_col,
          paste0(
            biological_samples,
            "_"
          )
        )
    ]
    
    
    if (length(possible) == 0) {
      return(NA_character_)
    }
    
    
    possible[
      which.max(
        nchar(possible)
      )
    ]
  }
  
  
  tibble(
    
    replicate_col =
      replicate_candidates,
    
    sample =
      vapply(
        replicate_candidates,
        map_one_rep,
        character(1)
      )
    
  ) %>%
    
    filter(
      !is.na(sample)
    )
}


rep_map_18S <- map_replicate_columns(
  data_18S_reps,
  sample_cols_18S
)

rep_map_COI <- map_replicate_columns(
  data_COI_reps,
  sample_cols_COI
)


replicate_structure_18S <- rep_map_18S %>%
  count(
    sample,
    name = "n_pcr_columns"
  )


replicate_structure_COI <- rep_map_COI %>%
  count(
    sample,
    name = "n_pcr_columns"
  )


# ============================================================
# 8. TAXONOMY FOR TARGET SPECIES
# ============================================================

prepare_species_taxonomy <- function(
    df,
    marker_name
) {
  
  if (!"kingdom" %in% names(df)) {
    
    stop(
      "The ",
      marker_name,
      " abundance table has no kingdom column."
    )
  }
  
  
  df %>%
    
    select(
      asv,
      species,
      kingdom
    ) %>%
    
    mutate(
      
      kingdom =
        normalise_kingdom(
          kingdom
        ),
      
      marker =
        marker_name
    ) %>%
    
    filter(
      !is.na(species),
      species != "",
      kingdom %in% c(
        "Plantae",
        "Metazoa"
      )
    ) %>%
    
    inner_join(
      species_status,
      by = "species"
    ) %>%
    
    mutate(
      status = factor(
        status,
        levels = status_levels
      )
    ) %>%
    
    distinct(
      marker,
      asv,
      species,
      kingdom,
      status
    )
}


tax_18S <- prepare_species_taxonomy(
  data_18S,
  "18S"
)

tax_COI <- prepare_species_taxonomy(
  data_COI,
  "COI"
)


taxonomy_conflicts <- bind_rows(
  tax_18S,
  tax_COI
) %>%
  
  distinct(
    marker,
    asv,
    species
  ) %>%
  
  count(
    marker,
    asv,
    name = "n_species"
  ) %>%
  
  filter(
    n_species > 1
  )


if (nrow(taxonomy_conflicts) > 0) {
  
  print(taxonomy_conflicts)
  
  stop(
    "Some marker-ASV combinations map to more than one selected species. ",
    "Resolve taxonomy before plotting."
  )
}


# ============================================================
# 9. BIOLOGICAL-SAMPLE PRESENCE
# ============================================================

build_collapsed_presence <- function(
    df,
    taxonomy_df,
    sample_cols,
    metadata_key,
    marker_name
) {
  
  df %>%
    
    select(
      asv,
      all_of(sample_cols)
    ) %>%
    
    inner_join(
      
      taxonomy_df %>%
        select(
          asv,
          species,
          kingdom,
          status
        ),
      
      by = "asv"
      
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
    
    group_by(
      species,
      kingdom,
      status,
      sample
    ) %>%
    
    summarise(
      present = TRUE,
      n_asvs_detected = n_distinct(asv),
      .groups = "drop"
    ) %>%
    
    mutate(
      marker = marker_name
    ) %>%
    
    left_join(
      metadata_key,
      by = "sample"
    )
}


presence_18S <- build_collapsed_presence(
  data_18S,
  tax_18S,
  sample_cols_18S,
  metadata_key,
  "18S"
)

presence_COI <- build_collapsed_presence(
  data_COI,
  tax_COI,
  sample_cols_COI,
  metadata_key,
  "COI"
)


missing_metadata <- bind_rows(
  presence_18S,
  presence_COI
) %>%
  
  filter(
    is.na(site) |
      is.na(year)
  ) %>%
  
  distinct(
    marker,
    sample
  )


if (nrow(missing_metadata) > 0) {
  
  print(missing_metadata)
  
  stop(
    "Some biological samples could not be matched to site/year ",
    "in metadata_samples.xlsx."
  )
}


# ============================================================
# 10. PCR REPLICATE SUPPORT
# ============================================================

build_pcr_support <- function(
    reps_df,
    taxonomy_df,
    rep_map,
    marker_name
) {
  
  if (nrow(rep_map) == 0) {
    
    warning(
      "No PCR replicate columns could be mapped for ",
      marker_name
    )
    
    return(
      tibble(
        species = character(),
        kingdom = character(),
        status = character(),
        sample = character(),
        n_pcr_reps = integer(),
        marker = character()
      )
    )
  }
  
  
  rep_cols <- rep_map$replicate_col
  
  
  reps_df %>%
    
    select(
      asv,
      all_of(rep_cols)
    ) %>%
    
    inner_join(
      
      taxonomy_df %>%
        select(
          asv,
          species,
          kingdom,
          status
        ),
      
      by = "asv"
      
    ) %>%
    
    pivot_longer(
      cols = all_of(rep_cols),
      names_to = "replicate_col",
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
      rep_map,
      by = "replicate_col"
    ) %>%
    
    distinct(
      species,
      kingdom,
      status,
      sample,
      replicate_col
    ) %>%
    
    count(
      species,
      kingdom,
      status,
      sample,
      name = "n_pcr_reps"
    ) %>%
    
    mutate(
      marker = marker_name
    )
}


pcr_18S <- build_pcr_support(
  data_18S_reps,
  tax_18S,
  rep_map_18S,
  "18S"
)

pcr_COI <- build_pcr_support(
  data_COI_reps,
  tax_COI,
  rep_map_COI,
  "COI"
)


# ============================================================
# 11. FINAL TIMELINE DATA
# ============================================================

timeline_data <- bind_rows(
  presence_18S,
  presence_COI
) %>%
  
  left_join(
    
    bind_rows(
      pcr_18S,
      pcr_COI
    ) %>%
      
      select(
        marker,
        species,
        kingdom,
        status,
        sample,
        n_pcr_reps
      ),
    
    by = c(
      "marker",
      "species",
      "kingdom",
      "status",
      "sample"
    )
  ) %>%
  
  mutate(
    
    n_pcr_reps =
      as.integer(
        n_pcr_reps
      ),
    
    status =
      factor(
        as.character(status),
        levels = status_levels
      )
    
  ) %>%
  
  arrange(
    marker,
    kingdom,
    site,
    species,
    year,
    sample
  )


collapsed_without_rep_support <- timeline_data %>%
  
  filter(
    is.na(n_pcr_reps)
  ) %>%
  
  distinct(
    marker,
    species,
    sample,
    site,
    year
  )


if (nrow(collapsed_without_rep_support) > 0) {
  
  warning(
    nrow(collapsed_without_rep_support),
    " collapsed detections have no mapped positive PCR replicate. ",
    "They are exported in diagnostics and excluded from plots."
  )
}


timeline_plot_data <- timeline_data %>%
  
  filter(
    !is.na(n_pcr_reps),
    n_pcr_reps > 0
  ) %>%
  
  mutate(
    status = factor(
      as.character(status),
      levels = status_levels
    )
  )


# ============================================================
# 12. COLOURS
# ============================================================

status_colors <- c(
  "NIS" = "#CD3333",
  "NAT/NIS" = "#97FFFF",
  "NIS/NAT" = "#9A32CD",
  "CRY" = "#CDAD00",
  "CRY/NAT" = "#0072B2"
)


site_order <- c(
  "DEE",
  "STO",
  "CCO",
  "CLR",
  "CMB",
  "CSD",
  "HIT",
  "HIO"
)


# ============================================================
# 13. LEGEND GUIDES
# ============================================================

status_guide <- guide_legend(
  
  order = 1,
  
  override.aes = list(
    shape = 21,
    size = 5,
    colour = "black",
    alpha = 1,
    stroke = 0.6
  )
)


kingdom_guide <- guide_legend(
  
  order = 2,
  
  override.aes = list(
    size = 4.5,
    fill = "grey80",
    colour = "black"
  )
)


size_guide <- guide_legend(
  
  order = 3,
  
  override.aes = list(
    shape = 21,
    fill = "grey80",
    colour = "black"
  )
)


# ============================================================
# 14. INDIVIDUAL SITE × MARKER × KINGDOM TIMELINES
# ============================================================

plot_one_timeline <- function(
    df,
    site_name,
    marker_name,
    kingdom_name
) {
  
  plot_df <- df %>%
    
    filter(
      site == site_name,
      marker == marker_name,
      kingdom == kingdom_name
    )
  
  
  if (nrow(plot_df) == 0) {
    return(NULL)
  }
  
  
  # Alphabetical species order
  species_order <- sort(
    unique(
      as.character(
        plot_df$species
      )
    )
  )
  
  
  plot_df <- plot_df %>%
    
    mutate(
      
      species = factor(
        species,
        levels = rev(
          species_order
        )
      ),
      
      status = factor(
        as.character(status),
        levels = status_levels
      )
    )
  
  
  size_breaks <- sort(
    unique(
      plot_df$n_pcr_reps
    )
  )
  
  
  p <- ggplot(
    plot_df,
    aes(
      x = year,
      y = species
    )
  ) +
    
    geom_line(
      aes(
        group = species
      ),
      linewidth = 0.30,
      colour = "grey78"
    ) +
    
    geom_point(
      aes(
        size = n_pcr_reps,
        fill = status
      ),
      shape = 21,
      colour = "black",
      stroke = 0.45,
      alpha = 0.95,
      
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    
    scale_y_discrete(
      expand = expansion(
        add = c(
          0.7,
          0.7
        )
      )
    ) +
    
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    
    scale_size_continuous(
      range = c(
        2,
        5.5
      ),
      breaks = size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    
    scale_x_continuous(
      breaks = sort(
        unique(
          plot_df$year
        )
      )
    ) +
    
    guides(
      fill = status_guide,
      size = size_guide
    ) +
    
    labs(
      
      title = paste(
        site_name,
        "—",
        kingdom_name,
        "—",
        marker_name
      ),
      
      subtitle = paste0(
        "One point = one biological sample; ",
        "point size = number of positive PCR replicates"
      ),
      
      x = "Year (CE)",
      y = NULL
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    theme(
      
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      panel.grid.minor =
        element_blank(),
      
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      
      plot.subtitle = element_text(
        size = 10
      ),
      
      axis.text.y = element_text(
        face = "italic",
        size = 11,
        margin = margin(
          r = 8
        ),
        lineheight = 1.1
      ),
      
      axis.text.x = element_text(
        size = 12,
        angle = 45,
        hjust = 1
      ),
      
      legend.position = "right"
    )
  
  
  out_svg <- file.path(
    output_dir,
    paste0(
      safe_filename(site_name),
      "_",
      marker_name,
      "_",
      safe_filename(kingdom_name),
      "_timeline.svg"
    )
  )
  
  
  # More real vertical space per species
  individual_height <- max(
    7,
    3 +
      0.42 *
      dplyr::n_distinct(
        plot_df$species
      )
  )
  
  
  individual_height <- min(
    individual_height,
    40
  )
  
  
  ggsave(
    filename = out_svg,
    plot = p,
    width = 10,
    height = individual_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  invisible(p)
}


# ============================================================
# 15. GENERATE ALL INDIVIDUAL TIMELINES
# ============================================================

plot_combinations <- timeline_plot_data %>%
  
  distinct(
    site,
    marker,
    kingdom
  ) %>%
  
  filter(
    kingdom %in% c(
      "Plantae",
      "Metazoa"
    )
  ) %>%
  
  arrange(
    site,
    marker,
    kingdom
  )


timeline_plots <- list()


for (
  i in seq_len(
    nrow(
      plot_combinations
    )
  )
) {
  
  ss <- plot_combinations$site[i]
  mm <- plot_combinations$marker[i]
  kk <- plot_combinations$kingdom[i]
  
  
  message(
    "Plotting: ",
    ss,
    " | ",
    mm,
    " | ",
    kk
  )
  
  
  plot_name <- paste(
    ss,
    mm,
    kk,
    sep = "_"
  )
  
  
  timeline_plots[[plot_name]] <-
    plot_one_timeline(
      timeline_plot_data,
      site_name = ss,
      marker_name = mm,
      kingdom_name = kk
    )
}


# ============================================================
# 16. MAIN FIGURE DATA
# ============================================================

main_timeline_data <- timeline_plot_data %>%
  
  filter(
    
    site %in% site_order,
    
    marker %in% c(
      "COI",
      "18S"
    ),
    
    kingdom %in% c(
      "Metazoa",
      "Plantae"
    )
  ) %>%
  
  mutate(
    
    site = factor(
      site,
      levels = site_order
    ),
    
    marker = factor(
      marker,
      levels = c(
        "COI",
        "18S"
      )
    ),
    
    kingdom = factor(
      kingdom,
      levels = c(
        "Metazoa",
        "Plantae"
      )
    ),
    
    status = factor(
      as.character(status),
      levels = status_levels
    )
  )


# ============================================================
# 17. MAIN MANUSCRIPT FIGURE
#
# IMPORTANT:
#
# scales = "free_y"
# space  = "free_y"
#
# Sites with many species receive MORE vertical space.
# Sites with few species receive LESS vertical space.
#
# Species are alphabetical.
# ============================================================

main_species_order <- main_timeline_data %>%
  
  distinct(
    site,
    species
  ) %>%
  
  arrange(
    site,
    species
  )


main_timeline_data <- main_timeline_data %>%
  
  mutate(
    species_site = paste(
      species,
      as.character(site),
      sep = "___SITE___"
    )
  )


species_levels_main <- main_species_order %>%
  
  mutate(
    species_site = paste(
      species,
      as.character(site),
      sep = "___SITE___"
    )
  ) %>%
  
  pull(
    species_site
  )


main_timeline_data <- main_timeline_data %>%
  
  mutate(
    species_site = factor(
      species_site,
      levels = rev(
        unique(
          species_levels_main
        )
      )
    )
  )


species_axis_labels <- function(x) {
  
  sub(
    "___SITE___.*$",
    "",
    x
  )
}


if (
  nrow(
    main_timeline_data
  ) == 0
) {
  
  warning(
    "No data available for the main manuscript timeline figure."
  )
  
} else {
  
  
  main_size_breaks <- sort(
    unique(
      main_timeline_data$n_pcr_reps
    )
  )
  
  
  p_main <- ggplot(
    main_timeline_data,
    aes(
      x = year,
      y = species_site
    )
  ) +
    
    geom_point(
      aes(
        size = n_pcr_reps,
        fill = status,
        shape = kingdom
      ),
      
      colour = "black",
      stroke = 0.50,
      alpha = 0.95,
      
      # Horizontal jitter only
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    
    # --------------------------------------------------------
  # KEY CHANGE:
  # space = "free_y"
  #
  # This gives more physical vertical space to sites that
  # contain more species.
  # --------------------------------------------------------
  
  facet_grid(
    rows = vars(site),
    cols = vars(marker),
    scales = "free_y",
    space = "free_y",
    drop = FALSE
  ) +
    
    scale_y_discrete(
      labels = species_axis_labels,
      
      # Small extra padding at the top and bottom
      # of each site's species axis
      expand = expansion(
        add = c(
          0.7,
          0.7
        )
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    
    scale_shape_manual(
      values = c(
        "Metazoa" = 21,
        "Plantae" = 24
      ),
      drop = FALSE,
      name = "Kingdom"
    ) +
    
    # Slightly smaller points
    scale_size_continuous(
      range = c(
        2,
        5.5
      ),
      breaks = main_size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    
    guides(
      fill = status_guide,
      shape = kingdom_guide,
      size = size_guide
    ) +
    
    labs(
      x = "Year (CE)",
      y = NULL
    ) +
    
    theme_bw(
      base_size = 14
    ) +
    
    theme(
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      panel.grid.major.x = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      strip.text.x = element_text(
        size = 17,
        face = "bold"
      ),
      
      strip.text.y = element_text(
        size = 15,
        face = "bold"
      ),
      
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey40",
        linewidth = 0.5
      ),
      
      axis.title.x = element_text(
        size = 17,
        face = "bold",
        margin = margin(
          t = 12
        )
      ),
      
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = 0.5
      ),
      
      # Slightly smaller labels + more room
      axis.text.y = element_text(
        size = 11,
        face = "italic",
        margin = margin(
          r = 9
        ),
        lineheight = 1.1
      ),
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.position = "right",
      
      panel.spacing.x = unit(
        1,
        "lines"
      ),
      
      panel.spacing.y = unit(
        1.5,
        "lines"
      )
    )
  
  
  # ==========================================================
  # MAIN FIGURE HEIGHT
  #
  # We give approximately 0.19 inches per site-species row.
  # This is the main control for label separation.
  # ==========================================================
  
  n_site_species <- main_timeline_data %>%
    
    distinct(
      site,
      species
    ) %>%
    
    nrow()
  
  
  species_spacing_inches <- 0.19
  
  
  main_height <- max(
    28,
    6 +
      species_spacing_inches *
      n_site_species
  )
  
  
  # Safety limit
  main_height <- min(
    main_height,
    65
  )
  
  
  cat(
    "\nMain figure site-species combinations:",
    n_site_species,
    "\n"
  )
  
  cat(
    "Target vertical spacing per species:",
    species_spacing_inches,
    "inches\n"
  )
  
  cat(
    "Main figure height:",
    round(
      main_height,
      1
    ),
    "inches\n"
  )
  
  
  main_svg <- file.path(
    output_dir,
    "MAIN_MANUSCRIPT_ALL_SITES_COI_18S_timeline.svg"
  )
  
  main_pdf <- file.path(
    output_dir,
    "MAIN_MANUSCRIPT_ALL_SITES_COI_18S_timeline.pdf"
  )
  
  main_png <- file.path(
    output_dir,
    "MAIN_MANUSCRIPT_ALL_SITES_COI_18S_timeline.png"
  )
  
  
  ggsave(
    filename = main_svg,
    plot = p_main,
    width = 16,
    height = main_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  ggsave(
    filename = main_pdf,
    plot = p_main,
    width = 16,
    height = main_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  # 300 dpi to avoid ragg dimension problems
  ggsave(
    filename = main_png,
    plot = p_main,
    width = 16,
    height = main_height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  message(
    "\nMAIN MANUSCRIPT FIGURE SAVED:\n",
    main_svg
  )
  
  
  print(
    p_main
  )
}




# ============================================================
# 17B. MAIN FIGURES SEPARATED BY KINGDOM
#
# Creates the same MAIN layout separately for:
#   - Metazoa
#   - Plantae
#
# Rows    = sites
# Columns = COI / 18S
# Fill    = status
# Size    = positive PCR replicates
#
# Because each figure contains a single kingdom, point shape is
# fixed and the redundant Kingdom legend is removed.
# ============================================================

plot_main_one_kingdom <- function(
    main_data,
    kingdom_name,
    output_dir,
    species_spacing_inches = 0.19
) {
  
  kingdom_data <- main_data %>%
    filter(
      as.character(kingdom) == kingdom_name
    ) %>%
    droplevels()
  
  
  if (nrow(kingdom_data) == 0) {
    
    warning(
      "No data available for kingdom: ",
      kingdom_name
    )
    
    return(NULL)
  }
  
  
  # ----------------------------------------------------------
  # Alphabetical species order within each site
  # ----------------------------------------------------------
  
  kingdom_species_order <- kingdom_data %>%
    
    distinct(
      site,
      species
    ) %>%
    
    arrange(
      site,
      species
    ) %>%
    
    mutate(
      species_site = paste(
        species,
        as.character(site),
        sep = "___SITE___"
      )
    )
  
  
  species_levels_kingdom <-
    kingdom_species_order$species_site
  
  
  kingdom_data <- kingdom_data %>%
    
    mutate(
      
      species_site = paste(
        species,
        as.character(site),
        sep = "___SITE___"
      ),
      
      species_site = factor(
        species_site,
        levels = rev(
          unique(
            species_levels_kingdom
          )
        )
      ),
      
      status = factor(
        as.character(status),
        levels = status_levels
      )
    )
  
  
  # Metazoa = circle; Plantae = triangle
  kingdom_shape <- ifelse(
    kingdom_name == "Metazoa",
    21,
    24
  )
  
  
  kingdom_size_breaks <- sort(
    unique(
      kingdom_data$n_pcr_reps
    )
  )
  
  
  p_kingdom <- ggplot(
    kingdom_data,
    aes(
      x = year,
      y = species_site
    )
  ) +
    
    geom_point(
      
      aes(
        size = n_pcr_reps,
        fill = status
      ),
      
      shape = kingdom_shape,
      colour = "black",
      stroke = 0.50,
      alpha = 0.95,
      
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    
    facet_grid(
      rows = vars(site),
      cols = vars(marker),
      scales = "free_y",
      space = "free_y",
      drop = FALSE
    ) +
    
    scale_y_discrete(
      labels = species_axis_labels,
      expand = expansion(
        add = c(
          0.7,
          0.7
        )
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::breaks_pretty(
        n = 8
      ),
      expand = expansion(
        mult = c(
          0.02,
          0.02
        )
      )
    ) +
    
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    
    scale_size_continuous(
      range = c(
        2,
        5.5
      ),
      breaks = kingdom_size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    
    guides(
      fill = status_guide,
      size = size_guide
    ) +
    
    labs(
      title = kingdom_name,
      x = "Year (CE)",
      y = NULL
    ) +
    
    theme_bw(
      base_size = 14
    ) +
    
    theme(
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      panel.grid.major.x = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      strip.text.x = element_text(
        size = 17,
        face = "bold"
      ),
      
      strip.text.y = element_text(
        size = 15,
        face = "bold"
      ),
      
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey40",
        linewidth = 0.5
      ),
      
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 10
        )
      ),
      
      axis.title.x = element_text(
        size = 17,
        face = "bold",
        margin = margin(
          t = 12
        )
      ),
      
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = 0.5
      ),
      
      axis.text.y = element_text(
        size = 11,
        face = "italic",
        margin = margin(
          r = 9
        ),
        lineheight = 1.1
      ),
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.position = "right",
      
      panel.spacing.x = unit(
        1,
        "lines"
      ),
      
      panel.spacing.y = unit(
        1.5,
        "lines"
      )
    )
  
  
  # ----------------------------------------------------------
  # Dynamic height
  # ----------------------------------------------------------
  
  n_site_species_kingdom <- kingdom_data %>%
    
    distinct(
      site,
      species
    ) %>%
    
    nrow()
  
  
  kingdom_height <- max(
    18,
    6 +
      species_spacing_inches *
      n_site_species_kingdom
  )
  
  
  kingdom_height <- min(
    kingdom_height,
    60
  )
  
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "MAIN KINGDOM FIGURE:",
    kingdom_name,
    "\n"
  )
  
  cat(
    "Site-species combinations:",
    n_site_species_kingdom,
    "\n"
  )
  
  cat(
    "Figure height:",
    round(
      kingdom_height,
      1
    ),
    "inches\n"
  )
  
  cat(
    "============================================================\n"
  )
  
  
  kingdom_filename <- toupper(
    safe_filename(
      kingdom_name
    )
  )
  
  
  out_svg <- file.path(
    output_dir,
    paste0(
      "MAIN_",
      kingdom_filename,
      "_ALL_SITES_COI_18S_timeline.svg"
    )
  )
  
  
  out_pdf <- file.path(
    output_dir,
    paste0(
      "MAIN_",
      kingdom_filename,
      "_ALL_SITES_COI_18S_timeline.pdf"
    )
  )
  
  
  out_png <- file.path(
    output_dir,
    paste0(
      "MAIN_",
      kingdom_filename,
      "_ALL_SITES_COI_18S_timeline.png"
    )
  )
  
  
  ggsave(
    filename = out_svg,
    plot = p_kingdom,
    width = 16,
    height = kingdom_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  ggsave(
    filename = out_pdf,
    plot = p_kingdom,
    width = 16,
    height = kingdom_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  ggsave(
    filename = out_png,
    plot = p_kingdom,
    width = 16,
    height = kingdom_height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  message(
    "\n",
    kingdom_name,
    " MAIN FIGURE SAVED:\n",
    out_svg
  )
  
  
  print(
    p_kingdom
  )
  
  
  invisible(
    p_kingdom
  )
}


# ============================================================
# 17C. GENERATE METAZOA MAIN FIGURE
# ============================================================

p_main_metazoa <- plot_main_one_kingdom(
  main_data = main_timeline_data,
  kingdom_name = "Metazoa",
  output_dir = output_dir,
  species_spacing_inches = 0.19
)


# ============================================================
# 17D. GENERATE PLANTAE MAIN FIGURE
# ============================================================

p_main_plantae <- plot_main_one_kingdom(
  main_data = main_timeline_data,
  kingdom_name = "Plantae",
  output_dir = output_dir,
  species_spacing_inches = 0.19
)


# ============================================================
# 18. FREE-Y FIGURE
#
# Each SITE × MARKER panel has its own species axis.
# This is useful for seeing marker-specific species.
# ============================================================

free_y_data <- timeline_plot_data %>%
  
  filter(
    
    site %in% site_order,
    
    marker %in% c(
      "COI",
      "18S"
    ),
    
    kingdom %in% c(
      "Metazoa",
      "Plantae"
    )
  ) %>%
  
  mutate(
    
    site = factor(
      site,
      levels = site_order
    ),
    
    marker = factor(
      marker,
      levels = c(
        "COI",
        "18S"
      )
    ),
    
    kingdom = factor(
      kingdom,
      levels = c(
        "Metazoa",
        "Plantae"
      )
    ),
    
    status = factor(
      as.character(status),
      levels = status_levels
    ),
    
    species_panel = paste(
      species,
      as.character(site),
      as.character(marker),
      sep = "___PANEL___"
    )
  )


# ------------------------------------------------------------
# Alphabetical species order within each site × marker
# ------------------------------------------------------------

free_species_order <- free_y_data %>%
  
  distinct(
    site,
    marker,
    species,
    species_panel
  ) %>%
  
  arrange(
    site,
    marker,
    species
  )


free_species_levels <- free_species_order %>%
  
  pull(
    species_panel
  )


free_y_data <- free_y_data %>%
  
  mutate(
    species_panel = factor(
      species_panel,
      levels = rev(
        unique(
          free_species_levels
        )
      )
    )
  )


species_panel_labels <- function(x) {
  
  sub(
    "___PANEL___.*$",
    "",
    x
  )
}


# Initialise for robust export even if FREE-Y data are empty
species_per_free_panel <- tibble(
  site = factor(levels = site_order),
  marker = factor(levels = c("COI", "18S")),
  n_species = integer()
)


if (
  nrow(
    free_y_data
  ) > 0
) {
  
  
  free_size_breaks <- sort(
    unique(
      free_y_data$n_pcr_reps
    )
  )
  
  
  p_main_free_y <- ggplot(
    free_y_data,
    aes(
      x = year,
      y = species_panel
    )
  ) +
    
    geom_point(
      aes(
        size = n_pcr_reps,
        fill = status,
        shape = kingdom
      ),
      
      colour = "black",
      stroke = 0.50,
      alpha = 0.95,
      
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    
    facet_wrap(
      vars(
        site,
        marker
      ),
      ncol = 2,
      scales = "free_y",
      drop = TRUE
    ) +
    
    scale_y_discrete(
      labels = species_panel_labels,
      expand = expansion(
        add = c(
          0.7,
          0.7
        )
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    
    scale_shape_manual(
      values = c(
        "Metazoa" = 21,
        "Plantae" = 24
      ),
      drop = FALSE,
      name = "Kingdom"
    ) +
    
    scale_size_continuous(
      range = c(
        2,
        5.5
      ),
      breaks = free_size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    
    guides(
      fill = status_guide,
      shape = kingdom_guide,
      size = size_guide
    ) +
    
    labs(
      x = "Year (CE)",
      y = NULL
    ) +
    
    theme_bw(
      base_size = 14
    ) +
    
    theme(
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      panel.grid.major.x = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      strip.text = element_text(
        size = 16,
        face = "bold"
      ),
      
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey40",
        linewidth = 0.5
      ),
      
      axis.title.x = element_text(
        size = 17,
        face = "bold",
        margin = margin(
          t = 12
        )
      ),
      
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = 0.5
      ),
      
      axis.text.y = element_text(
        size = 11,
        face = "italic",
        margin = margin(
          r = 9
        ),
        lineheight = 1.1
      ),
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.position = "right",
      
      panel.spacing.x = unit(
        1.2,
        "lines"
      ),
      
      panel.spacing.y = unit(
        1.5,
        "lines"
      )
    )
  
  
  # ==========================================================
  # FREE-Y HEIGHT
  # ==========================================================
  
  species_per_free_panel <- free_y_data %>%
    
    distinct(
      site,
      marker,
      species
    ) %>%
    
    count(
      site,
      marker,
      name = "n_species"
    )
  
  
  max_species_free_panel <- max(
    species_per_free_panel$n_species,
    na.rm = TRUE
  )
  
  
  n_free_panels <- nrow(
    species_per_free_panel
  )
  
  
  n_panel_rows <- ceiling(
    n_free_panels / 2
  )
  
  
  # A little more vertical room than before
  height_per_panel_row <- max(
    7,
    3 +
      0.20 *
      max_species_free_panel
  )
  
  
  free_height <- max(
    24,
    n_panel_rows *
      height_per_panel_row
  )
  
  
  # Safety limit
  free_height <- min(
    free_height,
    65
  )
  
  
  cat(
    "\nFREE-Y FIGURE\n"
  )
  
  cat(
    "Number of panels:",
    n_free_panels,
    "\n"
  )
  
  cat(
    "Maximum species in one panel:",
    max_species_free_panel,
    "\n"
  )
  
  cat(
    "Number of panel rows:",
    n_panel_rows,
    "\n"
  )
  
  cat(
    "Calculated figure height:",
    round(
      free_height,
      1
    ),
    "inches\n"
  )
  
  
  free_svg <- file.path(
    output_dir,
    "MAIN_LIKE_ALL_SITES_COI_18S_FREE_Y_timeline.svg"
  )
  
  free_pdf <- file.path(
    output_dir,
    "MAIN_LIKE_ALL_SITES_COI_18S_FREE_Y_timeline.pdf"
  )
  
  free_png <- file.path(
    output_dir,
    "MAIN_LIKE_ALL_SITES_COI_18S_FREE_Y_timeline.png"
  )
  
  
  ggsave(
    filename = free_svg,
    plot = p_main_free_y,
    width = 15,
    height = free_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  ggsave(
    filename = free_pdf,
    plot = p_main_free_y,
    width = 15,
    height = free_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  ggsave(
    filename = free_png,
    plot = p_main_free_y,
    width = 15,
    height = free_height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  message(
    "\nFREE-Y FIGURE SAVED:\n",
    free_svg
  )
  
  
  print(
    p_main_free_y
  )
}




# ============================================================
# 18B. FREE-Y FIGURES SEPARATED BY KINGDOM
#
# Creates the same FREE-Y layout separately for:
#   - Metazoa
#   - Plantae
#
# Every SITE × MARKER panel keeps its own independent species
# axis, exactly as in the combined FREE-Y figure.
# ============================================================

plot_free_y_one_kingdom <- function(
    free_data,
    kingdom_name,
    output_dir
) {
  
  kingdom_free_data <- free_data %>%
    
    filter(
      as.character(kingdom) == kingdom_name
    ) %>%
    
    droplevels()
  
  
  if (nrow(kingdom_free_data) == 0) {
    
    warning(
      "No FREE-Y data available for kingdom: ",
      kingdom_name
    )
    
    return(NULL)
  }
  
  
  # ----------------------------------------------------------
  # Rebuild panel-specific species IDs and alphabetical order
  # ----------------------------------------------------------
  
  kingdom_free_data <- kingdom_free_data %>%
    
    mutate(
      
      species_panel = paste(
        species,
        as.character(site),
        as.character(marker),
        sep = "___PANEL___"
      ),
      
      status = factor(
        as.character(status),
        levels = status_levels
      )
    )
  
  
  kingdom_free_species_order <- kingdom_free_data %>%
    
    distinct(
      site,
      marker,
      species,
      species_panel
    ) %>%
    
    arrange(
      site,
      marker,
      species
    )
  
  
  kingdom_free_species_levels <-
    kingdom_free_species_order$species_panel
  
  
  kingdom_free_data <- kingdom_free_data %>%
    
    mutate(
      species_panel = factor(
        species_panel,
        levels = rev(
          unique(
            kingdom_free_species_levels
          )
        )
      )
    )
  
  
  kingdom_shape <- ifelse(
    kingdom_name == "Metazoa",
    21,
    24
  )
  
  
  kingdom_free_size_breaks <- sort(
    unique(
      kingdom_free_data$n_pcr_reps
    )
  )
  
  
  p_free_kingdom <- ggplot(
    kingdom_free_data,
    aes(
      x = year,
      y = species_panel
    )
  ) +
    
    geom_point(
      
      aes(
        size = n_pcr_reps,
        fill = status
      ),
      
      shape = kingdom_shape,
      colour = "black",
      stroke = 0.50,
      alpha = 0.95,
      
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    
    facet_wrap(
      vars(
        site,
        marker
      ),
      ncol = 2,
      scales = "free_y",
      drop = TRUE
    ) +
    
    scale_y_discrete(
      labels = species_panel_labels,
      expand = expansion(
        add = c(
          0.7,
          0.7
        )
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::breaks_pretty(
        n = 8
      ),
      expand = expansion(
        mult = c(
          0.02,
          0.02
        )
      )
    ) +
    
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    
    scale_size_continuous(
      range = c(
        2,
        5.5
      ),
      breaks = kingdom_free_size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    
    guides(
      fill = status_guide,
      size = size_guide
    ) +
    
    labs(
      title = paste0(
        kingdom_name,
        " — independent species axis by site × marker"
      ),
      x = "Year (CE)",
      y = NULL
    ) +
    
    theme_bw(
      base_size = 14
    ) +
    
    theme(
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      panel.grid.major.x = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      strip.text = element_text(
        size = 16,
        face = "bold"
      ),
      
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey40",
        linewidth = 0.5
      ),
      
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 10
        )
      ),
      
      axis.title.x = element_text(
        size = 17,
        face = "bold",
        margin = margin(
          t = 12
        )
      ),
      
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = 0.5
      ),
      
      axis.text.y = element_text(
        size = 11,
        face = "italic",
        margin = margin(
          r = 9
        ),
        lineheight = 1.1
      ),
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.position = "right",
      
      panel.spacing.x = unit(
        1.2,
        "lines"
      ),
      
      panel.spacing.y = unit(
        1.5,
        "lines"
      )
    )
  
  
  # ----------------------------------------------------------
  # Height based on the largest panel in this kingdom
  # ----------------------------------------------------------
  
  species_per_kingdom_free_panel <- kingdom_free_data %>%
    
    distinct(
      site,
      marker,
      species
    ) %>%
    
    count(
      site,
      marker,
      name = "n_species"
    )
  
  
  max_species_kingdom_free_panel <- max(
    species_per_kingdom_free_panel$n_species,
    na.rm = TRUE
  )
  
  
  n_kingdom_free_panels <- nrow(
    species_per_kingdom_free_panel
  )
  
  
  n_kingdom_panel_rows <- ceiling(
    n_kingdom_free_panels / 2
  )
  
  
  height_per_kingdom_panel_row <- max(
    6,
    3 +
      0.20 *
      max_species_kingdom_free_panel
  )
  
  
  kingdom_free_height <- max(
    18,
    n_kingdom_panel_rows *
      height_per_kingdom_panel_row
  )
  
  
  kingdom_free_height <- min(
    kingdom_free_height,
    60
  )
  
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "FREE-Y KINGDOM FIGURE:",
    kingdom_name,
    "\n"
  )
  
  cat(
    "Number of panels:",
    n_kingdom_free_panels,
    "\n"
  )
  
  cat(
    "Maximum species in one panel:",
    max_species_kingdom_free_panel,
    "\n"
  )
  
  cat(
    "Figure height:",
    round(
      kingdom_free_height,
      1
    ),
    "inches\n"
  )
  
  cat(
    "============================================================\n"
  )
  
  
  kingdom_filename <- toupper(
    safe_filename(
      kingdom_name
    )
  )
  
  
  out_svg <- file.path(
    output_dir,
    paste0(
      "FREE_Y_",
      kingdom_filename,
      "_ALL_SITES_COI_18S_timeline.svg"
    )
  )
  
  
  out_pdf <- file.path(
    output_dir,
    paste0(
      "FREE_Y_",
      kingdom_filename,
      "_ALL_SITES_COI_18S_timeline.pdf"
    )
  )
  
  
  out_png <- file.path(
    output_dir,
    paste0(
      "FREE_Y_",
      kingdom_filename,
      "_ALL_SITES_COI_18S_timeline.png"
    )
  )
  
  
  ggsave(
    filename = out_svg,
    plot = p_free_kingdom,
    width = 15,
    height = kingdom_free_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  ggsave(
    filename = out_pdf,
    plot = p_free_kingdom,
    width = 15,
    height = kingdom_free_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  ggsave(
    filename = out_png,
    plot = p_free_kingdom,
    width = 15,
    height = kingdom_free_height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  message(
    "\n",
    kingdom_name,
    " FREE-Y FIGURE SAVED:\n",
    out_svg
  )
  
  
  print(
    p_free_kingdom
  )
  
  
  invisible(
    list(
      plot = p_free_kingdom,
      species_per_panel = species_per_kingdom_free_panel
    )
  )
}


# ============================================================
# 18C. GENERATE METAZOA FREE-Y FIGURE
# ============================================================

free_y_metazoa_result <- plot_free_y_one_kingdom(
  free_data = free_y_data,
  kingdom_name = "Metazoa",
  output_dir = output_dir
)

p_free_y_metazoa <- if (
  is.null(
    free_y_metazoa_result
  )
) {
  NULL
} else {
  free_y_metazoa_result$plot
}


# ============================================================
# 18D. GENERATE PLANTAE FREE-Y FIGURE
# ============================================================

free_y_plantae_result <- plot_free_y_one_kingdom(
  free_data = free_y_data,
  kingdom_name = "Plantae",
  output_dir = output_dir
)

p_free_y_plantae <- if (
  is.null(
    free_y_plantae_result
  )
) {
  NULL
} else {
  free_y_plantae_result$plot
}



# ============================================================
# 18E. MAIN-TEXT SELECTION: ALL METAZOA + SELECTED PLANTAE
# ============================================================

main_plant_species <- c(
  "Populus pseudoglauca",
  "Chenopodium quinoa",
  "Ophioglossum reticulatum",
  "Avicennia marina",
  "Tamarix pentandra",
  "Hordeum vulgare",
  "Chloropicon laureae",
  "Vitis_pseudoreticulata",
)

# Check that requested plant names occur exactly in the plotted dataset
missing_main_plants <- setdiff(
  main_plant_species,
  unique(as.character(timeline_plot_data$species))
)

if (length(missing_main_plants) > 0) {
  warning(
    "These selected plant species were not found exactly in timeline_plot_data: ",
    paste(missing_main_plants, collapse = ", ")
  )
}

selected_main_data <- timeline_plot_data %>%
  filter(
    site %in% site_order,
    marker %in% c("COI", "18S"),
    kingdom == "Metazoa" |
      (kingdom == "Plantae" & species %in% main_plant_species)
  ) %>%
  mutate(
    site = factor(site, levels = site_order),
    marker = factor(marker, levels = c("COI", "18S")),
    kingdom = factor(kingdom, levels = c("Metazoa", "Plantae")),
    status = factor(as.character(status), levels = status_levels)
  )

# ------------------------------------------------------------
# 18F. MAIN layout: all Metazoa + selected Plantae
# Shared species axis between COI and 18S within each site
# ------------------------------------------------------------

selected_species_order <- selected_main_data %>%
  distinct(site, species) %>%
  arrange(site, species) %>%
  mutate(
    species_site = paste(species, as.character(site), sep = "___SITE___")
  )

selected_species_levels <- selected_species_order$species_site

selected_main_plot_data <- selected_main_data %>%
  mutate(
    species_site = paste(species, as.character(site), sep = "___SITE___"),
    species_site = factor(
      species_site,
      levels = rev(unique(selected_species_levels))
    )
  )

if (nrow(selected_main_plot_data) > 0) {
  
  selected_size_breaks <- sort(unique(selected_main_plot_data$n_pcr_reps))
  
  p_main_metazoa_selected_plants <- ggplot(
    selected_main_plot_data,
    aes(x = year, y = species_site)
  ) +
    geom_point(
      aes(
        size = n_pcr_reps,
        fill = status,
        shape = kingdom
      ),
      colour = "black",
      stroke = 0.50,
      alpha = 0.95,
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    facet_grid(
      rows = vars(site),
      cols = vars(marker),
      scales = "free_y",
      space = "free_y",
      drop = FALSE
    ) +
    scale_y_discrete(
      labels = species_axis_labels,
      expand = expansion(add = c(0.7, 0.7))
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    scale_shape_manual(
      values = c("Metazoa" = 21, "Plantae" = 24),
      limits = c("Metazoa", "Plantae"),
      drop = FALSE,
      name = "Kingdom"
    ) +
    scale_size_continuous(
      range = c(2, 5.5),
      breaks = selected_size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    guides(
      fill = status_guide,
      shape = kingdom_guide,
      size = size_guide
    ) +
    labs(
      x = "Year (CE)",
      y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey90"),
      strip.text.x = element_text(size = 17, face = "bold"),
      strip.text.y = element_text(size = 15, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.5),
      axis.title.x = element_text(size = 17, face = "bold", margin = margin(t = 12)),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 11, face = "italic", margin = margin(r = 9), lineheight = 1.1),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(1.5, "lines")
    )
  
  n_selected_site_species <- selected_main_plot_data %>%
    distinct(site, species) %>%
    nrow()
  
  selected_main_height <- max(
    24,
    6 + 0.19 * n_selected_site_species
  )
  selected_main_height <- min(selected_main_height, 65)
  
  selected_main_svg <- file.path(
    output_dir,
    "MAIN_METAZOA_PLUS_SELECTED_PLANTAE_ALL_SITES_COI_18S_timeline.svg"
  )
  selected_main_pdf <- file.path(
    output_dir,
    "MAIN_METAZOA_PLUS_SELECTED_PLANTAE_ALL_SITES_COI_18S_timeline.pdf"
  )
  selected_main_png <- file.path(
    output_dir,
    "MAIN_METAZOA_PLUS_SELECTED_PLANTAE_ALL_SITES_COI_18S_timeline.png"
  )
  
  ggsave(selected_main_svg, p_main_metazoa_selected_plants,
         width = 16, height = selected_main_height,
         limitsize = FALSE, bg = "white")
  ggsave(selected_main_pdf, p_main_metazoa_selected_plants,
         width = 16, height = selected_main_height,
         limitsize = FALSE, bg = "white")
  ggsave(selected_main_png, p_main_metazoa_selected_plants,
         width = 16, height = selected_main_height, dpi = 300,
         limitsize = FALSE, bg = "white")
  
  message("\nMAIN Metazoa + selected Plantae figure saved:\n", selected_main_svg)
  print(p_main_metazoa_selected_plants)
}


# ------------------------------------------------------------
# 18G. FREE-Y layout: all Metazoa + selected Plantae
# Independent species axis for each site × marker panel
# ------------------------------------------------------------

selected_free_y_data <- selected_main_data %>%
  mutate(
    species_panel = paste(
      species,
      as.character(site),
      as.character(marker),
      sep = "___PANEL___"
    )
  )

selected_free_species_order <- selected_free_y_data %>%
  distinct(site, marker, species, species_panel) %>%
  arrange(site, marker, species)

selected_free_species_levels <- selected_free_species_order$species_panel

selected_free_y_data <- selected_free_y_data %>%
  mutate(
    species_panel = factor(
      species_panel,
      levels = rev(unique(selected_free_species_levels))
    )
  )

species_per_selected_free_panel <- tibble(
  site = factor(levels = site_order),
  marker = factor(levels = c("COI", "18S")),
  n_species = integer()
)

if (nrow(selected_free_y_data) > 0) {
  
  selected_free_size_breaks <- sort(unique(selected_free_y_data$n_pcr_reps))
  
  p_free_y_metazoa_selected_plants <- ggplot(
    selected_free_y_data,
    aes(x = year, y = species_panel)
  ) +
    geom_point(
      aes(
        size = n_pcr_reps,
        fill = status,
        shape = kingdom
      ),
      colour = "black",
      stroke = 0.50,
      alpha = 0.95,
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    facet_wrap(
      vars(site, marker),
      ncol = 2,
      scales = "free_y",
      drop = TRUE
    ) +
    scale_y_discrete(
      labels = species_panel_labels,
      expand = expansion(add = c(0.7, 0.7))
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    scale_shape_manual(
      values = c("Metazoa" = 21, "Plantae" = 24),
      limits = c("Metazoa", "Plantae"),
      drop = FALSE,
      name = "Kingdom"
    ) +
    scale_size_continuous(
      range = c(2, 5.5),
      breaks = selected_free_size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    guides(
      fill = status_guide,
      shape = kingdom_guide,
      size = size_guide
    ) +
    labs(
      x = "Year (CE)",
      y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey90"),
      strip.text = element_text(size = 16, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.5),
      axis.title.x = element_text(size = 17, face = "bold", margin = margin(t = 12)),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 11, face = "italic", margin = margin(r = 9), lineheight = 1.1),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(1.5, "lines")
    )
  
  species_per_selected_free_panel <- selected_free_y_data %>%
    distinct(site, marker, species) %>%
    count(site, marker, name = "n_species")
  
  max_species_selected_free_panel <- max(
    species_per_selected_free_panel$n_species,
    na.rm = TRUE
  )
  n_selected_free_panels <- nrow(species_per_selected_free_panel)
  n_selected_panel_rows <- ceiling(n_selected_free_panels / 2)
  
  selected_height_per_panel_row <- max(
    7,
    3 + 0.20 * max_species_selected_free_panel
  )
  
  selected_free_height <- max(
    24,
    n_selected_panel_rows * selected_height_per_panel_row
  )
  selected_free_height <- min(selected_free_height, 65)
  
  selected_free_svg <- file.path(
    output_dir,
    "FREE_Y_METAZOA_PLUS_SELECTED_PLANTAE_ALL_SITES_COI_18S_timeline.svg"
  )
  selected_free_pdf <- file.path(
    output_dir,
    "FREE_Y_METAZOA_PLUS_SELECTED_PLANTAE_ALL_SITES_COI_18S_timeline.pdf"
  )
  selected_free_png <- file.path(
    output_dir,
    "FREE_Y_METAZOA_PLUS_SELECTED_PLANTAE_ALL_SITES_COI_18S_timeline.png"
  )
  
  ggsave(selected_free_svg, p_free_y_metazoa_selected_plants,
         width = 15, height = selected_free_height,
         limitsize = FALSE, bg = "white")
  ggsave(selected_free_pdf, p_free_y_metazoa_selected_plants,
         width = 15, height = selected_free_height,
         limitsize = FALSE, bg = "white")
  ggsave(selected_free_png, p_free_y_metazoa_selected_plants,
         width = 15, height = selected_free_height, dpi = 300,
         limitsize = FALSE, bg = "white")
  
  message("\nFREE-Y Metazoa + selected Plantae figure saved:\n", selected_free_svg)
  print(p_free_y_metazoa_selected_plants)
}


# ============================================================
# 18H. SUPPLEMENTARY TIMELINE:
# PLANTAE NOT SELECTED FOR THE MAIN FIGURE
# ============================================================

# ------------------------------------------------------------
# 18H.1 Select all Plantae EXCEPT those used in the main figure
# ------------------------------------------------------------

remaining_plant_data <- timeline_plot_data %>%
  filter(
    site %in% site_order,
    marker %in% c("COI", "18S"),
    kingdom == "Plantae",
    !species %in% main_plant_species
  ) %>%
  mutate(
    site = factor(site, levels = unique(site_order)),
    marker = factor(marker, levels = c("COI", "18S")),
    kingdom = factor(kingdom, levels = c("Metazoa", "Plantae")),
    status = factor(as.character(status), levels = status_levels)
  )


# ------------------------------------------------------------
# 18H.2 Check which plant species are included
# ------------------------------------------------------------

remaining_plant_species <- remaining_plant_data %>%
  distinct(species) %>%
  arrange(species)

message(
  "\nNumber of plant species NOT included in the main figure: ",
  nrow(remaining_plant_species)
)

print(remaining_plant_species)


# Optional check:
# confirm that none of the main-figure plants remain
remaining_main_plants <- intersect(
  main_plant_species,
  unique(as.character(remaining_plant_data$species))
)

if (length(remaining_main_plants) > 0) {
  warning(
    "Some main-figure plants are still present: ",
    paste(remaining_main_plants, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 18H.3 Create shared species axis between COI and 18S
# within each site
# ------------------------------------------------------------

remaining_plant_species_order <- remaining_plant_data %>%
  distinct(site, species) %>%
  arrange(site, species) %>%
  mutate(
    species_site = paste(
      species,
      as.character(site),
      sep = "___SITE___"
    )
  )

remaining_plant_species_levels <-
  remaining_plant_species_order$species_site


remaining_plant_plot_data <- remaining_plant_data %>%
  mutate(
    species_site = paste(
      species,
      as.character(site),
      sep = "___SITE___"
    ),
    species_site = factor(
      species_site,
      levels = rev(unique(remaining_plant_species_levels))
    )
  )


# ------------------------------------------------------------
# 18H.4 Labels for species axis
#
# Remove the site identifier added internally above.
# Species names will remain italic through theme(axis.text.y = ...)
# ------------------------------------------------------------

remaining_plant_axis_labels <- function(x) {
  sub("___SITE___.*$", "", x)
}

status_colors <- c(
  "NIS"     = "#CD3333",
  "NAT"     = "#2E8B57",
  "NAT/NIS" = "#97FFFF",
  "NIS/NAT" = "#9A32CD",
  "CRY"     = "#CDAD00",
  "CRY/NAT" = "#0072B2"
)
# ------------------------------------------------------------
# 18H.5 Plot
# ------------------------------------------------------------

if (nrow(remaining_plant_plot_data) > 0) {
  
  remaining_plant_size_breaks <- sort(
    unique(remaining_plant_plot_data$n_pcr_reps)
  )
  
  p_remaining_plants <- ggplot(
    remaining_plant_plot_data,
    aes(
      x = year,
      y = species_site
    )
  ) +
    
    geom_point(
      aes(
        size = n_pcr_reps,
        fill = status,
        shape = kingdom
      ),
      colour = "black",
      stroke = 0.50,
      alpha = 0.95,
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
    ) +
    
    facet_grid(
      rows = vars(site),
      cols = vars(marker),
      scales = "free_y",
      space = "free_y",
      drop = FALSE
    ) +
    
    scale_y_discrete(
      labels = remaining_plant_axis_labels,
      expand = expansion(add = c(0.7, 0.7))
    ) +
    
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    scale_fill_manual(
      values = status_colors,
      limits = status_levels,
      breaks = status_levels,
      drop = FALSE,
      name = "Status"
    ) +
    
    # All observations are Plantae, but retaining the same
    # triangle used in the main timeline keeps figure consistency
    scale_shape_manual(
      values = c("Metazoa" = 21, "Plantae" = 24),
      limits = c("Metazoa", "Plantae"),
      drop = FALSE,
      name = "Kingdom"
    ) +
    
    scale_size_continuous(
      range = c(2, 5.5),
      breaks = remaining_plant_size_breaks,
      name = "Positive PCR\nreplicates"
    ) +
    
    guides(
      fill = status_guide,
      shape = kingdom_guide,
      size = size_guide
    ) +
    
    labs(
      x = "Year (CE)",
      y = NULL
    ) +
    
    theme_bw(base_size = 14) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      panel.grid.major.x = element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      
      strip.text.x = element_text(
        size = 17,
        face = "bold"
      ),
      
      strip.text.y = element_text(
        size = 15,
        face = "bold"
      ),
      
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey40",
        linewidth = 0.5
      ),
      
      axis.title.x = element_text(
        size = 17,
        face = "bold",
        margin = margin(t = 12)
      ),
      
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = 0.5
      ),
      
      axis.text.y = element_text(
        size = 11,
        face = "italic",
        margin = margin(r = 9),
        lineheight = 1.1
      ),
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.position = "right",
      
      panel.spacing.x = unit(1, "lines"),
      
      panel.spacing.y = unit(1.5, "lines")
    )
  
  
  # ----------------------------------------------------------
  # 18H.6 Automatic figure height
  # ----------------------------------------------------------
  
  n_remaining_plant_site_species <- remaining_plant_plot_data %>%
    distinct(site, species) %>%
    nrow()
  
  remaining_plant_height <- max(
    18,
    6 + 0.19 * n_remaining_plant_site_species
  )
  
  remaining_plant_height <- min(
    remaining_plant_height,
    65
  )
  
  
  # ----------------------------------------------------------
  # 18H.7 Save
  # ----------------------------------------------------------
  
  remaining_plant_svg <- file.path(
    output_dir,
    "SUPPLEMENTARY_REMAINING_PLANTAE_ALL_SITES_COI_18S_timeline.svg"
  )
  
  remaining_plant_pdf <- file.path(
    output_dir,
    "SUPPLEMENTARY_REMAINING_PLANTAE_ALL_SITES_COI_18S_timeline.pdf"
  )
  
  remaining_plant_png <- file.path(
    output_dir,
    "SUPPLEMENTARY_REMAINING_PLANTAE_ALL_SITES_COI_18S_timeline.png"
  )
  
  
  ggsave(
    remaining_plant_svg,
    p_remaining_plants,
    width = 16,
    height = remaining_plant_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  ggsave(
    remaining_plant_pdf,
    p_remaining_plants,
    width = 16,
    height = remaining_plant_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  ggsave(
    remaining_plant_png,
    p_remaining_plants,
    width = 16,
    height = remaining_plant_height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  
  message(
    "\nSUPPLEMENTARY remaining Plantae figure saved:\n",
    remaining_plant_svg
  )
  
  print(p_remaining_plants)
  
} else {
  
  warning(
    "No remaining Plantae were found after excluding main_plant_species."
  )
}




# ============================================================
# 19. SITE × MARKER DIAGNOSTIC
# ============================================================

site_marker_diagnostic <- tidyr::expand_grid(
  
  site = site_order,
  
  marker = c(
    "COI",
    "18S"
  )
  
) %>%
  
  left_join(
    
    timeline_plot_data %>%
      
      filter(
        site %in% site_order,
        marker %in% c(
          "COI",
          "18S"
        )
      ) %>%
      
      group_by(
        site,
        marker
      ) %>%
      
      summarise(
        
        n_species =
          n_distinct(
            species
          ),
        
        n_detections =
          n(),
        
        .groups = "drop"
      ),
    
    by = c(
      "site",
      "marker"
    )
    
  ) %>%
  
  mutate(
    
    n_species =
      replace_na(
        n_species,
        0L
      ),
    
    n_detections =
      replace_na(
        n_detections,
        0L
      ),
    
    plot_has_data =
      n_detections > 0
  )


print(
  site_marker_diagnostic
)


# ============================================================
# 20. STATUS DIAGNOSTIC
# ============================================================

status_diagnostic <- timeline_plot_data %>%
  
  count(
    status,
    marker,
    kingdom,
    name = "n_detections",
    .drop = FALSE
  ) %>%
  
  arrange(
    status,
    marker,
    kingdom
  )


print(
  status_diagnostic
)


# ============================================================
# 21. EXPORT DATA + DIAGNOSTICS
# ============================================================

summary_species <- timeline_plot_data %>%
  
  group_by(
    marker,
    kingdom,
    status,
    species
  ) %>%
  
  summarise(
    
    n_sites =
      n_distinct(
        site
      ),
    
    n_samples_detected =
      n_distinct(
        sample
      ),
    
    first_year =
      min(
        year,
        na.rm = TRUE
      ),
    
    last_year =
      max(
        year,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(
    marker,
    kingdom,
    species
  )


summary_sites <- timeline_plot_data %>%
  
  group_by(
    marker,
    kingdom,
    site
  ) %>%
  
  summarise(
    
    n_species =
      n_distinct(
        species
      ),
    
    n_samples_with_detections =
      n_distinct(
        sample
      ),
    
    .groups = "drop"
  )


writexl::write_xlsx(
  
  list(
    
    Timeline_data =
      timeline_data,
    
    Timeline_plot_data =
      timeline_plot_data,
    
    Species_summary =
      summary_species,
    
    Site_summary =
      summary_sites,
    
    Status_diagnostic =
      status_diagnostic,
    
    Status_species_table =
      species_status,
    
    Replicate_structure_18S =
      replicate_structure_18S,
    
    Replicate_structure_COI =
      replicate_structure_COI,
    
    Collapsed_without_rep_support =
      collapsed_without_rep_support,
    
    Metadata_key =
      metadata_key,
    
    Site_marker_diagnostic =
      site_marker_diagnostic,
    
    Species_per_free_panel =
      species_per_free_panel
    
  ),
  
  file.path(
    output_dir,
    "NIS_CRY_timelines_data_and_diagnostics.xlsx"
  )
)


# ============================================================
# 22. FINAL SUMMARY
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "TIMELINE ANALYSIS FINISHED\n"
)

cat(
  "Species statuses included: ",
  paste(
    target_statuses,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Detections plotted: ",
  nrow(
    timeline_plot_data
  ),
  "\n",
  sep = ""
)

cat(
  "Species plotted: ",
  n_distinct(
    timeline_plot_data$species
  ),
  "\n",
  sep = ""
)

cat(
  "Individual SVG timelines created: ",
  nrow(
    plot_combinations
  ),
  "\n",
  sep = ""
)

cat(
  "Output directory:\n",
  normalizePath(
    output_dir
  ),
  "\n"
)

cat(
  "============================================================\n"
)




# ============================================================
# 19. SUPPLEMENTARY FIGURE:
# ALL NATIVE (NAT) TAXA
#
# IMPORTANT:
# This entire section is independent from the previous
# NIS / NAT-NIS / NIS-NAT / CRY / CRY-NAT pipeline.
#
# It DOES NOT modify:
#   species_status
#   status_levels
#   tax_18S / tax_COI
#   presence_18S / presence_COI
#   timeline_data
#   timeline_plot_data
#   Fig. 1
# ============================================================


# ============================================================
# 19A. CREATE AN INDEPENDENT NAT STATUS TABLE
# ============================================================

nat_species_status <- status_species %>%
  
  transmute(
    
    species =
      clean_text(
        .data[[status_species_col]]
      ),
    
    status =
      normalise_status(
        .data[[status_combined_col]]
      )
  ) %>%
  
  filter(
    !is.na(species),
    species != "",
    status == "NAT"
  ) %>%
  
  mutate(
    status = "NAT"
  ) %>%
  
  distinct(
    species,
    .keep_all = TRUE
  )


cat(
  "\n============================================================\n",
  "NATIVE SPECIES RETAINED FROM STATUS TABLE\n",
  "============================================================\n"
)

print(
  nat_species_status %>%
    arrange(species)
)

cat(
  "\nNumber of NAT species in status table: ",
  nrow(nat_species_status),
  "\n",
  sep = ""
)


if (nrow(nat_species_status) == 0) {
  
  stop(
    "No species classified as NAT were found in Status_selected_species.xlsx."
  )
}


# ============================================================
# 19B. PREPARE NAT TAXONOMY
#
# This reproduces prepare_species_taxonomy(),
# but joins ONLY to the NAT table created above.
# ============================================================

prepare_nat_taxonomy <- function(
    df,
    marker_name
) {
  
  if (!"kingdom" %in% names(df)) {
    
    stop(
      "The ",
      marker_name,
      " abundance table has no kingdom column."
    )
  }
  
  
  df %>%
    
    select(
      asv,
      species,
      kingdom
    ) %>%
    
    mutate(
      
      asv =
        clean_text(asv),
      
      species =
        clean_text(species),
      
      kingdom =
        normalise_kingdom(
          kingdom
        ),
      
      marker =
        marker_name
    ) %>%
    
    filter(
      !is.na(species),
      species != "",
      kingdom %in% c(
        "Plantae",
        "Metazoa"
      )
    ) %>%
    
    inner_join(
      nat_species_status,
      by = "species"
    ) %>%
    
    mutate(
      status = "NAT"
    ) %>%
    
    distinct(
      marker,
      asv,
      species,
      kingdom,
      status
    )
}


# ============================================================
# 19C. NAT TAXONOMY FOR EACH MARKER
# ============================================================

nat_tax_18S <- prepare_nat_taxonomy(
  data_18S,
  "18S"
)


nat_tax_COI <- prepare_nat_taxonomy(
  data_COI,
  "COI"
)


cat(
  "\n============================================================\n",
  "NAT TAXONOMY RECOVERED FROM ABUNDANCE TABLES\n",
  "============================================================\n"
)


cat("\n18S:\n")

print(
  nat_tax_18S %>%
    distinct(
      species,
      kingdom
    ) %>%
    arrange(
      kingdom,
      species
    )
)


cat("\nCOI:\n")

print(
  nat_tax_COI %>%
    distinct(
      species,
      kingdom
    ) %>%
    arrange(
      kingdom,
      species
    )
)


cat(
  "\nNumber of NAT species in 18S: ",
  n_distinct(nat_tax_18S$species),
  "\n",
  sep = ""
)


cat(
  "Number of NAT species in COI: ",
  n_distinct(nat_tax_COI$species),
  "\n",
  sep = ""
)


# ============================================================
# 19D. CHECK NAT TAXONOMY CONFLICTS
# ============================================================

nat_taxonomy_conflicts <- bind_rows(
  nat_tax_18S,
  nat_tax_COI
) %>%
  
  distinct(
    marker,
    asv,
    species
  ) %>%
  
  count(
    marker,
    asv,
    name = "n_species"
  ) %>%
  
  filter(
    n_species > 1
  )


if (nrow(nat_taxonomy_conflicts) > 0) {
  
  print(
    nat_taxonomy_conflicts
  )
  
  stop(
    "Some NAT marker-ASV combinations map to more than one species."
  )
}


# ============================================================
# 19E. BIOLOGICAL-SAMPLE PRESENCE FOR NAT TAXA
#
# Uses the SAME function already used by the main pipeline.
# ============================================================

nat_presence_18S <- build_collapsed_presence(
  data_18S,
  nat_tax_18S,
  sample_cols_18S,
  metadata_key,
  "18S"
)


nat_presence_COI <- build_collapsed_presence(
  data_COI,
  nat_tax_COI,
  sample_cols_COI,
  metadata_key,
  "COI"
)


# ============================================================
# 19F. CHECK METADATA
# ============================================================

nat_missing_metadata <- bind_rows(
  nat_presence_18S,
  nat_presence_COI
) %>%
  
  filter(
    is.na(site) |
      is.na(year)
  ) %>%
  
  distinct(
    marker,
    sample
  )


if (nrow(nat_missing_metadata) > 0) {
  
  print(
    nat_missing_metadata
  )
  
  stop(
    "Some NAT biological samples could not be matched ",
    "to site/year in metadata_samples.xlsx."
  )
}


# ============================================================
# 19G. PCR-REPLICATE SUPPORT FOR NAT TAXA
#
# Again uses exactly the SAME function as the main pipeline.
# ============================================================

nat_pcr_18S <- build_pcr_support(
  data_18S_reps,
  nat_tax_18S,
  rep_map_18S,
  "18S"
)


nat_pcr_COI <- build_pcr_support(
  data_COI_reps,
  nat_tax_COI,
  rep_map_COI,
  "COI"
)


# ============================================================
# 19H. BUILD NAT TIMELINE DATA
#
# Same logic as section 11 of the original script.
# ============================================================

nat_timeline_data <- bind_rows(
  nat_presence_18S,
  nat_presence_COI
) %>%
  
  left_join(
    
    bind_rows(
      nat_pcr_18S,
      nat_pcr_COI
    ) %>%
      
      select(
        marker,
        species,
        kingdom,
        status,
        sample,
        n_pcr_reps
      ),
    
    by = c(
      "marker",
      "species",
      "kingdom",
      "status",
      "sample"
    )
  ) %>%
  
  mutate(
    
    n_pcr_reps =
      as.integer(
        n_pcr_reps
      ),
    
    status =
      "NAT"
  ) %>%
  
  arrange(
    marker,
    kingdom,
    site,
    species,
    year,
    sample
  )


# ============================================================
# 19I. CHECK DETECTIONS WITHOUT PCR SUPPORT
# ============================================================

nat_without_rep_support <- nat_timeline_data %>%
  
  filter(
    is.na(n_pcr_reps)
  ) %>%
  
  distinct(
    marker,
    species,
    sample,
    site,
    year
  )


if (nrow(nat_without_rep_support) > 0) {
  
  warning(
    nrow(nat_without_rep_support),
    " NAT collapsed detections have no mapped positive PCR ",
    "replicate and will therefore be excluded from the figure."
  )
}


# ============================================================
# 19J. FINAL NAT PLOTTING DATA
# ============================================================

nat_timeline_plot_data <- nat_timeline_data %>%
  
  filter(
    !is.na(n_pcr_reps),
    n_pcr_reps > 0,
    site %in% unique(site_order),
    marker %in% c(
      "COI",
      "18S"
    ),
    kingdom %in% c(
      "Metazoa",
      "Plantae"
    )
  ) %>%
  
  mutate(
    
    site =
      factor(
        as.character(site),
        levels = unique(site_order)
      ),
    
    marker =
      factor(
        as.character(marker),
        levels = c(
          "COI",
          "18S"
        )
      ),
    
    kingdom =
      factor(
        as.character(kingdom),
        levels = c(
          "Metazoa",
          "Plantae"
        )
      ),
    
    status =
      "NAT"
  )


# ============================================================
# 19K. FINAL DIAGNOSTICS
# ============================================================

cat(
  "\n============================================================\n",
  "FINAL NAT TIMELINE DATA\n",
  "============================================================\n"
)


cat(
  "\nTotal NAT detections plotted: ",
  nrow(nat_timeline_plot_data),
  "\n",
  sep = ""
)


cat(
  "Total NAT species detected: ",
  n_distinct(nat_timeline_plot_data$species),
  "\n",
  sep = ""
)


cat("\nNAT species detected:\n")

print(
  nat_timeline_plot_data %>%
    
    distinct(
      species,
      kingdom
    ) %>%
    
    arrange(
      kingdom,
      species
    ),
  n = Inf
)


cat("\nNAT species by marker:\n")

print(
  nat_timeline_plot_data %>%
    
    distinct(
      marker,
      species,
      kingdom
    ) %>%
    
    count(
      marker,
      kingdom,
      name = "n_species"
    )
)


cat("\nNAT records by site and marker:\n")

print(
  nat_timeline_plot_data %>%
    
    count(
      site,
      marker,
      name = "n_records"
    )
)


if (nrow(nat_timeline_plot_data) == 0) {
  
  stop(
    "NAT taxa were found in the status table, but no NAT ",
    "detections with positive PCR support were recovered."
  )
}


# ============================================================
# 19L. SHARED SPECIES AXIS
#
# Same species row is shared between COI and 18S
# within each site, matching the main figure.
# ============================================================

nat_species_order <- nat_timeline_plot_data %>%
  
  distinct(
    site,
    species
  ) %>%
  
  arrange(
    site,
    species
  ) %>%
  
  mutate(
    
    species_site =
      paste(
        species,
        as.character(site),
        sep = "___SITE___"
      )
  )


nat_species_levels <-
  nat_species_order$species_site


nat_timeline_plot_data <- nat_timeline_plot_data %>%
  
  mutate(
    
    species_site =
      paste(
        species,
        as.character(site),
        sep = "___SITE___"
      ),
    
    species_site =
      factor(
        species_site,
        levels = rev(
          unique(
            nat_species_levels
          )
        )
      )
  )


nat_species_axis_labels <- function(x) {
  
  sub(
    "___SITE___.*$",
    "",
    x
  )
}


# ============================================================
# 19M. NAT COLOUR
#
# Change this one value if NAT already has a fixed colour
# elsewhere in the manuscript.
# ============================================================

nat_colour <- "#2E8B57"


# ============================================================
# 19N. PLOT ALL NAT TAXA
# ============================================================

nat_size_breaks <- sort(
  unique(
    nat_timeline_plot_data$n_pcr_reps
  )
)


p_all_nat <- ggplot(
  nat_timeline_plot_data,
  aes(
    x = year,
    y = species_site
  )
) +
  
  geom_point(
    
    aes(
      size = n_pcr_reps,
      shape = kingdom
    ),
    
    fill = nat_colour,
    colour = "black",
    stroke = 0.50,
    alpha = 0.95,
    
    position =
      position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      )
  ) +
  
  facet_grid(
    rows = vars(site),
    cols = vars(marker),
    scales = "free_y",
    space = "free_y",
    drop = FALSE
  ) +
  
  scale_y_discrete(
    labels = nat_species_axis_labels,
    expand = expansion(
      add = c(
        0.7,
        0.7
      )
    )
  ) +
  
  scale_x_continuous(
    breaks =
      scales::breaks_pretty(
        n = 8
      ),
    expand =
      expansion(
        mult = c(
          0.02,
          0.02
        )
      )
  ) +
  
  scale_shape_manual(
    values = c(
      "Metazoa" = 21,
      "Plantae" = 24
    ),
    limits = c(
      "Metazoa",
      "Plantae"
    ),
    drop = TRUE,
    name = "Kingdom"
  ) +
  
  scale_size_continuous(
    range = c(
      2,
      5.5
    ),
    breaks = nat_size_breaks,
    name = "Positive PCR\nreplicates"
  ) +
  
  guides(
    
    shape =
      guide_legend(
        order = 1,
        override.aes = list(
          size = 4.5,
          fill = nat_colour,
          colour = "black"
        )
      ),
    
    size =
      guide_legend(
        order = 2,
        override.aes = list(
          shape = 21,
          fill = nat_colour,
          colour = "black"
        )
      )
  ) +
  
  labs(
    x = "Year (CE)",
    y = NULL
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.y =
      element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
    
    panel.grid.major.x =
      element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
    
    strip.text.x =
      element_text(
        size = 17,
        face = "bold"
      ),
    
    strip.text.y =
      element_text(
        size = 15,
        face = "bold"
      ),
    
    strip.background =
      element_rect(
        fill = "grey95",
        colour = "grey40",
        linewidth = 0.5
      ),
    
    axis.title.x =
      element_text(
        size = 17,
        face = "bold",
        margin = margin(
          t = 12
        )
      ),
    
    axis.text.x =
      element_text(
        size = 12,
        angle = 0,
        hjust = 0.5
      ),
    
    axis.text.y =
      element_text(
        size = 11,
        face = "italic",
        margin = margin(
          r = 9
        ),
        lineheight = 1.1
      ),
    
    legend.title =
      element_text(
        size = 13,
        face = "bold"
      ),
    
    legend.text =
      element_text(
        size = 12
      ),
    
    legend.position =
      "right",
    
    panel.spacing.x =
      unit(
        1,
        "lines"
      ),
    
    panel.spacing.y =
      unit(
        1.5,
        "lines"
      )
  )


# ============================================================
# 19O. FIGURE HEIGHT
# ============================================================

n_nat_site_species <- nat_timeline_plot_data %>%
  
  distinct(
    site,
    species
  ) %>%
  
  nrow()


nat_species_spacing_inches <- 0.19


nat_figure_height <- max(
  24,
  6 +
    nat_species_spacing_inches *
    n_nat_site_species
)


nat_figure_height <- min(
  nat_figure_height,
  65
)


cat(
  "\nNAT figure site-species combinations: ",
  n_nat_site_species,
  "\n",
  sep = ""
)


cat(
  "NAT figure height: ",
  round(
    nat_figure_height,
    1
  ),
  " inches\n",
  sep = ""
)


# ============================================================
# 19P. SAVE NAT FIGURE
# ============================================================

nat_svg <- file.path(
  output_dir,
  "SUPPLEMENTARY_ALL_NATIVE_NAT_TAXA_COI_18S_timeline.svg"
)


nat_pdf <- file.path(
  output_dir,
  "SUPPLEMENTARY_ALL_NATIVE_NAT_TAXA_COI_18S_timeline.pdf"
)


nat_png <- file.path(
  output_dir,
  "SUPPLEMENTARY_ALL_NATIVE_NAT_TAXA_COI_18S_timeline.png"
)


ggsave(
  filename = nat_svg,
  plot = p_all_nat,
  width = 16,
  height = nat_figure_height,
  limitsize = FALSE,
  bg = "white"
)


ggsave(
  filename = nat_pdf,
  plot = p_all_nat,
  width = 16,
  height = nat_figure_height,
  limitsize = FALSE,
  bg = "white"
)


ggsave(
  filename = nat_png,
  plot = p_all_nat,
  width = 16,
  height = nat_figure_height,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)


message(
  "\n============================================================\n",
  "ALL NAT FIGURE SAVED\n",
  nat_svg,
  "\n============================================================"
)


print(
  p_all_nat
)