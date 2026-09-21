# ============================================================
# TEMPORAL HAPLOTYPE NETWORKS
# GENERAL 50-YEAR BINS
# 4 NETWORKS PER ROW
# NO LOWER TIMELINES
# ============================================================
#
# QUESTION:
#
# How does detected haplotype composition and spatial sharing
# change through time?
#
#
# DESIGN:
#
# - One figure per species
# - One network panel per 50-year period
# - All sites are combined within each temporal period
# - Site contribution is represented by pie charts
# - Genetic topology is FIXED across all temporal panels
# - Node coordinates are FIXED across all temporal panels
# - Undetected haplotypes remain visible in grey
# - Detected haplotypes are represented by site-coloured pies
# - Node size = number of positive biological-sample detections
# - Four temporal networks per row
# - ONE site legend at the bottom
# - NO haplotype richness / dominance / sharing timelines below
#
#
# REQUIRED OBJECTS ALREADY CREATED IN R:
#
#   data_18S_updated
#   data_COI_updated
#   all_haplotype_results
#
# ============================================================



# ============================================================
# 1. PACKAGES
# ============================================================

library(readxl)
library(readr)

library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(tibble)

library(ggplot2)

library(igraph)
library(tidygraph)
library(ggraph)
library(scatterpie)

library(scales)
library(patchwork)
library(grid)



# ============================================================
# 2. PATHS
# ============================================================

data_dir <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Haplotypes_manuscript/Data"
)


metadata_file <- file.path(
  data_dir,
  "metadata_samples.xlsx"
)


temporal_output_dir <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Haplotypes_manuscript/Figures/",
  "Haplotype_networks/",
  "TEMPORAL_HAPLOTYPE_NETWORKS_50YR_4_PER_ROW"
)


dir.create(
  temporal_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# ============================================================
# 3. SETTINGS
# ============================================================

BIN_WIDTH <- 50

NETWORKS_PER_ROW <- 4

PNG_DPI <- 500

SHOW_ASV_LABELS <- FALSE


# ------------------------------------------------------------
# TRUE:
# Keep 50-year periods in which sediment samples exist,
# even if the focal species was not detected.
#
# This is scientifically useful because:
#
# sampled + no focal species detection
#
# is different from:
#
# no sediment samples available.
#
# ------------------------------------------------------------

KEEP_SAMPLED_ZERO_DETECTION_PERIODS <- FALSE



# ============================================================
# 4. SITE ORDER
# ============================================================

site_order <- c(
  "DEE",
  "STO",
  "CCO",
  "CLR",
  "CLR",
  "CSD",
  "CMB",
  "HIT",
  "HIO"
)



# ============================================================
# 5. SITE COLOURS
# ============================================================
#
# Same colours as the original haplotype-network workflow.
#
# COL and CLR are temporarily assigned the same colour so
# either code can be plotted if it occurs in the data.
#
# ============================================================

site_colors <- c(
  
  "DEE" = "#CDCD00",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  
  "CLR" = "burlywood1",
  "CLR" = "burlywood1",
  
  "CSD" = "#76EEC6",
  "CMB" = "#458B74",
  
  "HIT" = "#CD96CD",
  "HIO" = "#8B0A50"
  
)



# ============================================================
# 6. CHECK REQUIRED OBJECTS
# ============================================================

required_objects <- c(
  "data_18S_updated",
  "data_COI_updated",
  "all_haplotype_results"
)


missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]


if (length(missing_objects) > 0) {
  
  stop(
    "\nMissing required object(s):\n",
    paste(
      missing_objects,
      collapse = "\n"
    ),
    "\n\nRun the main haplotype-network workflow first."
  )
  
}



# ============================================================
# 7. READ METADATA
# ============================================================

metadata_raw <- readxl::read_excel(
  metadata_file
) %>%
  janitor::clean_names()


cat("\n")
cat("=============================================\n")
cat("METADATA COLUMNS\n")
cat("=============================================\n\n")

print(
  names(metadata_raw)
)



# ============================================================
# 8. FIND METADATA COLUMN
# ============================================================

find_first_column <- function(
    df,
    candidates,
    required = TRUE
) {
  
  found <- candidates[
    candidates %in%
      names(df)
  ]
  
  
  if (length(found) > 0) {
    
    return(
      found[1]
    )
    
  }
  
  
  if (required) {
    
    stop(
      "Could not identify metadata column among: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
    
  }
  
  
  return(NULL)
  
}



# ============================================================
# 9. IDENTIFY SAMPLE / YEAR / SITE
# ============================================================

sample_column <- find_first_column(
  
  metadata_raw,
  
  c(
    "sample",
    "sample_id",
    "sample_name",
    "sampleid",
    "samples",
    "biological_sample",
    "biological_sample_id"
  )
  
)


year_column <- find_first_column(
  
  metadata_raw,
  
  c(
    "year",
    "sample_year",
    "sampling_year",
    "calendar_year",
    "year_ad",
    "age_year",
    "date"
  )
  
)


site_column <- find_first_column(
  
  metadata_raw,
  
  c(
    "site",
    "site_code",
    "location",
    "locality"
  ),
  
  required = FALSE
  
)



# ============================================================
# 10. YEAR CONVERSION
# ============================================================

convert_to_year <- function(x) {
  
  if (inherits(x, "Date")) {
    
    return(
      as.numeric(
        format(
          x,
          "%Y"
        )
      )
    )
    
  }
  
  
  if (
    inherits(x, "POSIXct") ||
    inherits(x, "POSIXlt")
  ) {
    
    return(
      as.numeric(
        format(
          x,
          "%Y"
        )
      )
    )
    
  }
  
  
  suppressWarnings(
    as.numeric(x)
  )
  
}



# ============================================================
# 11. STANDARDISE METADATA
# ============================================================

metadata_temporal <- tibble::tibble(
  
  sample =
    as.character(
      metadata_raw[[sample_column]]
    ),
  
  year =
    convert_to_year(
      metadata_raw[[year_column]]
    )
  
)


if (!is.null(site_column)) {
  
  metadata_temporal$site <-
    as.character(
      metadata_raw[[site_column]]
    )
  
} else {
  
  metadata_temporal$site <-
    NA_character_
  
}



# ------------------------------------------------------------
# Infer site from sample ID when necessary
# ------------------------------------------------------------

metadata_temporal <- metadata_temporal %>%
  
  mutate(
    
    site = ifelse(
      
      is.na(site) |
        site == "",
      
      stringr::str_extract(
        sample,
        "^[A-Za-z]+"
      ),
      
      site
      
    ),
    
    site =
      toupper(
        as.character(site)
      ),
    
    sample =
      as.character(sample)
    
  ) %>%
  
  filter(
    
    !is.na(sample),
    
    sample != "",
    
    !is.na(year),
    
    is.finite(year),
    
    !is.na(site)
    
  )



# ============================================================
# 12. ASSIGN GENERAL 50-YEAR BINS
# ============================================================
#
# IMPORTANT:
#
# These bins are GENERAL.
#
# They are NOT defined separately by site.
#
# Example:
#
# 1750-1799
# 1800-1849
# 1850-1899
# 1900-1949
# 1950-1999
# 2000-2049
#
# All sites falling within the same 50-year period are pooled
# into that temporal network.
#
# ============================================================

metadata_temporal <- metadata_temporal %>%
  
  mutate(
    
    bin_start =
      floor(
        year / BIN_WIDTH
      ) *
      BIN_WIDTH,
    
    bin_end =
      bin_start +
      BIN_WIDTH -
      1,
    
    bin_mid =
      bin_start +
      (
        BIN_WIDTH - 1
      ) / 2,
    
    period =
      paste0(
        bin_start,
        "\u2013",
        bin_end
      )
    
  )



# ============================================================
# 13. CHECK TEMPORAL BINS
# ============================================================

bin_check <- metadata_temporal %>%
  
  arrange(
    year,
    site
  ) %>%
  
  select(
    sample,
    site,
    year,
    bin_start,
    bin_end,
    period
  )


cat("\n")
cat("=============================================\n")
cat("50-YEAR BIN CHECK\n")
cat("=============================================\n\n")

print(
  bin_check,
  n = 100
)


readr::write_csv(
  
  bin_check,
  
  file.path(
    temporal_output_dir,
    "CHECK_50yr_temporal_bins.csv"
  )
  
)



# ============================================================
# 14. GET MARKER DATASET
# ============================================================

get_marker_dataset <- function(marker_name) {
  
  marker_name <-
    toupper(marker_name)
  
  
  if (marker_name == "18S") {
    
    return(
      data_18S_updated
    )
    
  }
  
  
  if (marker_name == "COI") {
    
    return(
      data_COI_updated
    )
    
  }
  
  
  stop(
    "marker_name must be '18S' or 'COI'."
  )
  
}



# ============================================================
# 15. GET STORED HAPLOTYPE NETWORK
# ============================================================

get_temporal_network_result <- function(
    species_name,
    marker_name
) {
  
  marker_name <-
    toupper(marker_name)
  
  
  if (
    !species_name %in%
    names(
      all_haplotype_results
    )
  ) {
    
    stop(
      "Species absent from all_haplotype_results: ",
      species_name
    )
    
  }
  
  
  if (
    !marker_name %in%
    names(
      all_haplotype_results[[species_name]]
    )
  ) {
    
    stop(
      marker_name,
      " unavailable for ",
      species_name
    )
    
  }
  
  
  result <-
    all_haplotype_results[[species_name]][[marker_name]]
  
  
  if (is.null(result)) {
    
    stop(
      "NULL result for ",
      species_name,
      " | ",
      marker_name
    )
    
  }
  
  
  return(result)
  
}
  
test_peringia <- get_temporal_network_result(
  species_name = "Peringia ulvae",
  marker_name = "COI"
)

test_peringia
  
  # ============================================================
  # 16. EXTRACT ASV NAMES
  # ============================================================
  
  extract_network_asvs <- function(result) {
    
    if (is.null(result$ASVs)) {
      
      stop(
        "result$ASVs is missing."
      )
      
    }
    
    
    if (is.data.frame(result$ASVs)) {
      
      possible_columns <- c(
        "ASV",
        "asv",
        "name",
        "haplotype"
      )
      
      
      asv_column <- possible_columns[
        possible_columns %in%
          names(
            result$ASVs
          )
      ]
      
      
      if (length(asv_column) == 0) {
        
        stop(
          "Could not identify ASV column in result$ASVs."
        )
        
      }
      
      
      asvs <-
        as.character(
          result$ASVs[[asv_column[1]]]
        )
      
    } else {
      
      asvs <-
        as.character(
          result$ASVs
        )
      
    }
    
    
    asvs <- unique(
      
      asvs[
        !is.na(asvs) &
          asvs != ""
      ]
      
    )
    
    
    return(asvs)
    
  }
  
  
  
  # ============================================================
  # 17. CONVERT DISTANCES TO EDGE TABLE
  # ============================================================
  
  distance_object_to_edges <- function(
    distance_object,
    asvs
  ) {
    
    if (is.null(distance_object)) {
      
      return(NULL)
      
    }
    
    
    # ----------------------------------------------------------
    # MATRIX
    # ----------------------------------------------------------
    
    if (is.matrix(distance_object)) {
      
      m <- distance_object
      
      
      if (
        is.null(
          rownames(m)
        ) ||
        is.null(
          colnames(m)
        )
      ) {
        
        if (
          nrow(m) ==
          length(asvs) &&
          ncol(m) ==
          length(asvs)
        ) {
          
          rownames(m) <-
            asvs
          
          colnames(m) <-
            asvs
          
        }
        
      }
      
      
      idx <- which(
        upper.tri(m),
        arr.ind = TRUE
      )
      
      
      edges <- tibble::tibble(
        
        from =
          rownames(m)[
            idx[, 1]
          ],
        
        to =
          colnames(m)[
            idx[, 2]
          ],
        
        distance =
          as.numeric(
            m[idx]
          )
        
      )
      
      
      return(
        
        edges %>%
          
          filter(
            !is.na(distance),
            is.finite(distance)
          )
        
      )
      
    }
    
    
    # ----------------------------------------------------------
    # DIST OBJECT
    # ----------------------------------------------------------
    
    if (inherits(distance_object, "dist")) {
      
      return(
        
        distance_object_to_edges(
          
          as.matrix(
            distance_object
          ),
          
          asvs
          
        )
        
      )
      
    }
    
    
    # ----------------------------------------------------------
    # DATA FRAME
    # ----------------------------------------------------------
    
    if (is.data.frame(distance_object)) {
      
      df <-
        as.data.frame(
          distance_object
        )
      
      
      from_candidates <- c(
        "from",
        "ASV1",
        "asv1",
        "Var1"
      )
      
      
      to_candidates <- c(
        "to",
        "ASV2",
        "asv2",
        "Var2"
      )
      
      
      distance_candidates <- c(
        "distance",
        "dist",
        "value",
        "Freq"
      )
      
      
      from_col <- from_candidates[
        from_candidates %in%
          names(df)
      ]
      
      
      to_col <- to_candidates[
        to_candidates %in%
          names(df)
      ]
      
      
      dist_col <- distance_candidates[
        distance_candidates %in%
          names(df)
      ]
      
      
      if (
        length(from_col) == 0 ||
        length(to_col) == 0 ||
        length(dist_col) == 0
      ) {
        
        stop(
          "Could not identify from / to / distance columns."
        )
        
      }
      
      
      edges <- df %>%
        
        transmute(
          
          from =
            as.character(
              .data[[from_col[1]]]
            ),
          
          to =
            as.character(
              .data[[to_col[1]]]
            ),
          
          distance =
            suppressWarnings(
              as.numeric(
                .data[[dist_col[1]]]
              )
            )
          
        ) %>%
        
        filter(
          
          !is.na(from),
          
          !is.na(to),
          
          from != to,
          
          !is.na(distance),
          
          is.finite(distance)
          
        )
      
      
      return(edges)
      
    }
    
    
    stop(
      "Unsupported result$distances format."
    )
    
  }
  
  
  
  # ============================================================
  # 18. NORMALISE NETWORK COORDINATES
  # ============================================================
  
  normalise_coordinate <- function(
    x,
    to = c(
      0.08,
      0.92
    )
  ) {
    
    if (
      length(
        unique(x)
      ) <= 1
    ) {
      
      return(
        rep(
          mean(to),
          length(x)
        )
      )
      
    }
    
    
    scales::rescale(
      x,
      to = to
    )
    
  }
  
  
  
  # ============================================================
  # 19. BUILD FIXED NETWORK
  # ============================================================
  #
  # The network is calculated ONCE per species.
  #
  # Its topology and coordinates remain identical in every
  # temporal panel.
  #
  # ============================================================
  
  build_fixed_temporal_network <- function(result) {
    
    asvs <-
      extract_network_asvs(
        result
      )
    
    
    n_asvs <-
      length(asvs)
    
    
    # ----------------------------------------------------------
    # SINGLE HAPLOTYPE
    # ----------------------------------------------------------
    
    if (n_asvs == 1) {
      
      coordinates <- tibble::tibble(
        
        ASV =
          asvs,
        
        x =
          0.5,
        
        y =
          0.5
        
      )
      
      
      return(
        
        list(
          
          graph =
            NULL,
          
          coordinates =
            coordinates,
          
          edges =
            NULL,
          
          n_asvs =
            1
          
        )
        
      )
      
    }
    
    
    # ----------------------------------------------------------
    # COMPLETE GENETIC DISTANCE GRAPH
    # ----------------------------------------------------------
    
    distance_edges <-
      distance_object_to_edges(
        
        result$distances,
        
        asvs
        
      )
    
    
    if (
      is.null(distance_edges) ||
      nrow(distance_edges) == 0
    ) {
      
      stop(
        "No usable genetic distances found."
      )
      
    }
    
    
    distance_edges <- distance_edges %>%
      
      filter(
        from %in% asvs,
        to %in% asvs
      )
    
    
    vertices <- tibble::tibble(
      name = asvs
    )
    
    
    graph_full <-
      igraph::graph_from_data_frame(
        
        d =
          distance_edges,
        
        directed =
          FALSE,
        
        vertices =
          vertices
        
      )
    
    
    # ----------------------------------------------------------
    # Connectivity check
    # ----------------------------------------------------------
    
    connected <-
      tryCatch(
        
        igraph::is_connected(
          graph_full
        ),
        
        error = function(e) {
          
          igraph::is.connected(
            graph_full
          )
          
        }
        
      )
    
    
    if (!connected) {
      
      stop(
        "Full genetic-distance graph is not connected."
      )
      
    }
    
    
    # ----------------------------------------------------------
    # MST
    # ----------------------------------------------------------
    
    graph_mst <-
      igraph::mst(
        
        graph_full,
        
        weights =
          igraph::E(graph_full)$distance
        
      )
    
    
    # ----------------------------------------------------------
    # FIXED STRESS LAYOUT
    # ----------------------------------------------------------
    
    graph_tbl <-
      tidygraph::as_tbl_graph(
        graph_mst
      )
    
    
    layout_df <-
      ggraph::create_layout(
        
        graph_tbl,
        
        layout =
          "stress"
        
      ) %>%
      
      as.data.frame()
    
    
    layout_df$ASV <-
      igraph::V(graph_mst)$name
    
    
    coordinates <- layout_df %>%
      
      transmute(
        
        ASV =
          as.character(ASV),
        
        x =
          normalise_coordinate(x),
        
        y =
          normalise_coordinate(y)
        
      )
    
    
    # ----------------------------------------------------------
    # MST EDGES + FIXED COORDINATES
    # ----------------------------------------------------------
    
    mst_edges <-
      igraph::as_data_frame(
        
        graph_mst,
        
        what =
          "edges"
        
      ) %>%
      
      transmute(
        
        from =
          as.character(from),
        
        to =
          as.character(to),
        
        distance =
          as.numeric(distance)
        
      ) %>%
      
      left_join(
        
        coordinates %>%
          
          rename(
            
            from = ASV,
            
            x = x,
            
            y = y
            
          ),
        
        by =
          "from"
        
      ) %>%
      
      left_join(
        
        coordinates %>%
          
          rename(
            
            to = ASV,
            
            xend = x,
            
            yend = y
            
          ),
        
        by =
          "to"
        
      )
    
    
    return(
      
      list(
        
        graph =
          graph_mst,
        
        coordinates =
          coordinates,
        
        edges =
          mst_edges,
        
        n_asvs =
          n_asvs
        
      )
      
    )
    
  }
  
  
  
  # ============================================================
  # 20. BUILD TEMPORAL OCCURRENCES
  # ============================================================
  
  build_species_temporal_data <- function(
    species_name,
    marker_name
  ) {
    
    marker_name <-
      toupper(marker_name)
    
    
    df <-
      get_marker_dataset(
        marker_name
      )
    
    
    if (!"species" %in% names(df)) {
      
      stop(
        "'species' column missing from abundance table."
      )
      
    }
    
    
    if (!"ASV" %in% names(df)) {
      
      stop(
        "'ASV' column missing from abundance table."
      )
      
    }
    
    
    if (
      !species_name %in%
      unique(
        df$species
      )
    ) {
      
      stop(
        species_name,
        " absent from ",
        marker_name,
        " data."
      )
      
    }
    
    
    # ----------------------------------------------------------
    # Biological sample columns
    # ----------------------------------------------------------
    
    sample_cols <-
      intersect(
        
        names(df),
        
        metadata_temporal$sample
        
      )
    
    
    if (length(sample_cols) == 0) {
      
      stop(
        "No abundance-table samples matched metadata."
      )
      
    }
    
    
    # ----------------------------------------------------------
    # Marker metadata
    # ----------------------------------------------------------
    
    marker_metadata <- metadata_temporal %>%
      
      filter(
        sample %in%
          sample_cols
      )
    
    
    # ----------------------------------------------------------
    # Sampling effort PER GENERAL 50-YEAR PERIOD
    # ----------------------------------------------------------
    
    sample_effort <- marker_metadata %>%
      
      group_by(
        bin_start,
        bin_end,
        bin_mid,
        period
      ) %>%
      
      summarise(
        
        n_samples_available =
          n_distinct(sample),
        
        n_sites_sampled =
          n_distinct(site),
        
        .groups =
          "drop"
        
      )
    
    
    # ----------------------------------------------------------
    # Positive focal-species occurrences
    # ----------------------------------------------------------
    
    occurrences <- df %>%
      
      filter(
        species ==
          species_name
      ) %>%
      
      select(
        ASV,
        all_of(sample_cols)
      ) %>%
      
      pivot_longer(
        
        cols =
          all_of(sample_cols),
        
        names_to =
          "sample",
        
        values_to =
          "abundance"
        
      ) %>%
      
      mutate(
        
        ASV =
          as.character(ASV),
        
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
        ASV,
        sample
      ) %>%
      
      summarise(
        
        abundance =
          sum(
            abundance,
            na.rm = TRUE
          ),
        
        .groups =
          "drop"
        
      ) %>%
      
      left_join(
        
        marker_metadata %>%
          
          select(
            
            sample,
            site,
            year,
            
            bin_start,
            bin_end,
            bin_mid,
            period
            
          ),
        
        by =
          "sample"
        
      ) %>%
      
      filter(
        
        !is.na(period),
        
        !is.na(site)
        
      )
    
    
    return(
      
      list(
        
        occurrences =
          occurrences,
        
        sample_effort =
          sample_effort,
        
        marker_metadata =
          marker_metadata,
        
        sample_columns =
          sample_cols
        
      )
      
    )
    
  }
  
  
  
  # ============================================================
  # 21. CREATE PERIOD REFERENCE
  # ============================================================
  #
  # GENERAL 50-YEAR PERIODS.
  #
  # NOT site-specific.
  #
  # Periods with sediment samples but no species detection can
  # remain in the figure.
  #
  # Completely unsampled periods are removed.
  #
  # ============================================================
  
  create_period_reference <- function(
    sample_effort,
    occurrences
  ) {
    
    detected_periods <- occurrences %>%
      
      distinct(
        bin_start,
        period
      ) %>%
      
      mutate(
        focal_species_detected =
          TRUE
      )
    
    
    period_reference <- sample_effort %>%
      
      left_join(
        
        detected_periods,
        
        by = c(
          "bin_start",
          "period"
        )
        
      ) %>%
      
      mutate(
        
        focal_species_detected =
          replace_na(
            focal_species_detected,
            FALSE
          )
        
      )
    
    
    if (!KEEP_SAMPLED_ZERO_DETECTION_PERIODS) {
      
      period_reference <- period_reference %>%
        
        filter(
          focal_species_detected
        )
      
    }
    
    
    period_reference <- period_reference %>%
      
      filter(
        n_samples_available > 0
      ) %>%
      
      arrange(
        bin_start
      ) %>%
      
      mutate(
        
        period_order =
          row_number()
        
      )
    
    
    return(
      period_reference
    )
    
  }
  
  
  
  # ============================================================
  # 22. PREPARE TEMPORAL NODE DATA
  # ============================================================
  
  prepare_temporal_nodes <- function(
    occurrences,
    period_reference,
    fixed_network
  ) {
    
    all_asvs <-
      fixed_network$coordinates$ASV
    
    
    sites_present <- sort(
      
      unique(
        
        c(
          metadata_temporal$site,
          occurrences$site
        )
        
      )
      
    )
    
    
    sites_present <- sites_present[
      !is.na(sites_present) &
        sites_present != ""
    ]
    
    
    # ----------------------------------------------------------
    # Count POSITIVE BIOLOGICAL SAMPLES
    #
    # per:
    #
    # period x ASV x site
    # ----------------------------------------------------------
    
    site_counts <- occurrences %>%
      
      distinct(
        period,
        ASV,
        site,
        sample
      ) %>%
      
      count(
        
        period,
        ASV,
        site,
        
        name =
          "n_positive_samples_site"
        
      )
    
    
    # ----------------------------------------------------------
    # Wide pie data
    # ----------------------------------------------------------
    
    pie_data <- site_counts %>%
      
      pivot_wider(
        
        names_from =
          site,
        
        values_from =
          n_positive_samples_site,
        
        values_fill =
          0
        
      )
    
    
    # ----------------------------------------------------------
    # Complete period x ASV matrix
    # ----------------------------------------------------------
    
    node_data <- tidyr::crossing(
      
      period =
        period_reference$period,
      
      ASV =
        all_asvs
      
    ) %>%
      
      left_join(
        
        pie_data,
        
        by = c(
          "period",
          "ASV"
        )
        
      )
    
    
    # ----------------------------------------------------------
    # Add any missing site columns
    # ----------------------------------------------------------
    
    for (site_i in sites_present) {
      
      if (!site_i %in% names(node_data)) {
        
        node_data[[site_i]] <- 0
        
      }
      
    }
    
    
    # ----------------------------------------------------------
    # Replace NA site counts by zero
    # ----------------------------------------------------------
    
    node_data <- node_data %>%
      
      mutate(
        
        across(
          
          all_of(
            sites_present
          ),
          
          ~ replace_na(
            .x,
            0
          )
          
        )
        
      )
    
    
    # ----------------------------------------------------------
    # Total detections per haplotype/period
    # ----------------------------------------------------------
    
    node_data$n_positive_samples <-
      rowSums(
        
        node_data[
          ,
          sites_present,
          drop = FALSE
        ],
        
        na.rm = TRUE
        
      )
    
    
    node_data <- node_data %>%
      
      mutate(
        
        detected =
          n_positive_samples > 0
        
      )
    
    
    # ----------------------------------------------------------
    # Add fixed coordinates
    # ----------------------------------------------------------
    
    node_data <- node_data %>%
      
      left_join(
        
        fixed_network$coordinates,
        
        by =
          "ASV"
        
      )
    
    
    # ----------------------------------------------------------
    # Period metadata
    # ----------------------------------------------------------
    
    node_data <- node_data %>%
      
      left_join(
        
        period_reference %>%
          
          select(
            
            period,
            period_order,
            
            bin_start,
            bin_end,
            bin_mid,
            
            n_samples_available,
            n_sites_sampled,
            focal_species_detected
            
          ),
        
        by =
          "period"
        
      )
    
    
    # ----------------------------------------------------------
    # Haplotype richness per period
    # ----------------------------------------------------------
    
    period_metrics <- node_data %>%
      
      group_by(
        period
      ) %>%
      
      summarise(
        
        n_haplotypes =
          sum(detected),
        
        n_positive_haplotype_detections =
          sum(
            n_positive_samples
          ),
        
        .groups =
          "drop"
        
      )
    
    
    # ----------------------------------------------------------
    # Number of unique biological samples with focal species
    # ----------------------------------------------------------
    
    positive_samples <- occurrences %>%
      
      distinct(
        period,
        sample
      ) %>%
      
      count(
        
        period,
        
        name =
          "n_positive_samples_species"
        
      )
    
    
    period_metrics <- period_metrics %>%
      
      left_join(
        
        positive_samples,
        
        by =
          "period"
        
      ) %>%
      
      mutate(
        
        n_positive_samples_species =
          replace_na(
            n_positive_samples_species,
            0L
          )
        
      )
    
    
    node_data <- node_data %>%
      
      left_join(
        
        period_metrics,
        
        by =
          "period"
        
      )
    
    
    # ----------------------------------------------------------
    # Strip label
    # ----------------------------------------------------------
    #
    # Example:
    #
    # 1850–1899
    # nH = 12 | n = 8
    #
    # n = focal-species positive biological samples
    #
    # ----------------------------------------------------------
    
    node_data <- node_data %>%
      
      mutate(
        
        panel_label =
          paste0(
            
            period,
            
            "\n",
            
            "nH = ",
            n_haplotypes,
            
            " | n = ",
            n_positive_samples_species
            
          )
        
      )
    
    
    # ----------------------------------------------------------
    # COMMON NODE RADIUS SCALE
    # ----------------------------------------------------------
    
    positive_values <- node_data$n_positive_samples[
      node_data$n_positive_samples > 0
    ]
    
    
    n_total_nodes <-
      fixed_network$n_asvs
    
    
    max_radius <- dplyr::case_when(
      
      n_total_nodes >= 70 ~ 0.035,
      
      n_total_nodes >= 40 ~ 0.042,
      
      n_total_nodes >= 20 ~ 0.050,
      
      n_total_nodes >= 10 ~ 0.060,
      
      TRUE ~ 0.075
      
    )
    
    
    min_radius <- max_radius * 0.55
    
    absent_radius <- max_radius * 0.35
    
    
    if (length(positive_values) == 0) {
      
      node_data$radius <-
        absent_radius
      
    } else if (
      length(
        unique(
          positive_values
        )
      ) == 1
    ) {
      
      node_data <- node_data %>%
        
        mutate(
          
          radius =
            ifelse(
              
              detected,
              
              (
                min_radius +
                  max_radius
              ) / 2,
              
              absent_radius
              
            )
          
        )
      
    } else {
      
      min_positive <-
        min(
          positive_values
        )
      
      
      max_positive <-
        max(
          positive_values
        )
      
      
      node_data <- node_data %>%
        
        mutate(
          
          radius =
            ifelse(
              
              detected,
              
              scales::rescale(
                
                sqrt(
                  n_positive_samples
                ),
                
                from = c(
                  
                  sqrt(
                    min_positive
                  ),
                  
                  sqrt(
                    max_positive
                  )
                  
                ),
                
                to = c(
                  min_radius,
                  max_radius
                )
                
              ),
              
              absent_radius
              
            )
          
        )
      
    }
    
    
    return(
      
      list(
        
        node_data =
          node_data,
        
        sites_present =
          sites_present,
        
        period_metrics =
          period_metrics
        
      )
      
    )
    
  }
  
  
  
  # ============================================================
  # 23. PREPARE EDGES FOR EACH PERIOD
  # ============================================================
  
  prepare_temporal_edges <- function(
    fixed_network,
    period_reference,
    node_data
  ) {
    
    if (is.null(fixed_network$edges)) {
      
      return(NULL)
      
    }
    
    
    # ----------------------------------------------------------
    # Explicit expansion to avoid crossing/data-frame issues
    # ----------------------------------------------------------
    
    edge_base <- fixed_network$edges %>%
      
      mutate(
        edge_id =
          row_number()
      )
    
    
    edge_periods <- tidyr::expand_grid(
      
      period =
        period_reference$period,
      
      edge_id =
        edge_base$edge_id
      
    ) %>%
      
      left_join(
        
        edge_base,
        
        by =
          "edge_id"
        
      ) %>%
      
      select(
        -edge_id
      )
    
    
    # ----------------------------------------------------------
    # Endpoint detection state
    # ----------------------------------------------------------
    
    detection_lookup <- node_data %>%
      
      select(
        period,
        ASV,
        detected
      )
    
    
    edge_periods <- edge_periods %>%
      
      left_join(
        
        detection_lookup %>%
          
          rename(
            
            from =
              ASV,
            
            from_detected =
              detected
            
          ),
        
        by = c(
          "period",
          "from"
        )
        
      ) %>%
      
      left_join(
        
        detection_lookup %>%
          
          rename(
            
            to =
              ASV,
            
            to_detected =
              detected
            
          ),
        
        by = c(
          "period",
          "to"
        )
        
      ) %>%
      
      mutate(
        
        active_edge =
          from_detected &
          to_detected
        
      )
    
    
    return(
      edge_periods
    )
    
  }
  
  
  
  # ============================================================
  # 24. TEMPORAL METRICS FOR EXPORT ONLY
  # ============================================================
  #
  # These are NOT plotted below the networks.
  #
  # They are retained for quantitative analysis / Supplementary.
  #
  # ============================================================
  
  calculate_temporal_metrics <- function(
    occurrences,
    period_reference
  ) {
    
    if (nrow(occurrences) == 0) {
      
      return(
        
        period_reference %>%
          
          transmute(
            
            period,
            bin_start,
            bin_end,
            
            n_samples_available,
            n_sites_sampled,
            
            n_haplotypes = 0,
            n_positive_samples_species = 0,
            n_positive_sites_species = 0,
            shared_haplotypes = 0,
            site_specific_haplotypes = 0,
            prop_shared_haplotypes = NA_real_,
            prop_site_specific_haplotypes = NA_real_,
            mean_sites_per_haplotype = NA_real_,
            dominant_haplotype_frequency = NA_real_
            
          )
        
      )
      
    }
    
    
    # ----------------------------------------------------------
    # Basic metrics
    # ----------------------------------------------------------
    
    basic <- occurrences %>%
      
      group_by(
        period
      ) %>%
      
      summarise(
        
        n_haplotypes =
          n_distinct(ASV),
        
        n_positive_samples_species =
          n_distinct(sample),
        
        n_positive_sites_species =
          n_distinct(site),
        
        .groups =
          "drop"
        
      )
    
    
    # ----------------------------------------------------------
    # Spatial sharing
    # ----------------------------------------------------------
    
    spatial <- occurrences %>%
      
      distinct(
        period,
        ASV,
        site
      ) %>%
      
      count(
        
        period,
        ASV,
        
        name =
          "n_sites_haplotype"
        
      ) %>%
      
      group_by(
        period
      ) %>%
      
      summarise(
        
        shared_haplotypes =
          sum(
            n_sites_haplotype >= 2
          ),
        
        site_specific_haplotypes =
          sum(
            n_sites_haplotype == 1
          ),
        
        mean_sites_per_haplotype =
          mean(
            n_sites_haplotype
          ),
        
        .groups =
          "drop"
        
      )
    
    
    # ----------------------------------------------------------
    # Dominance
    # ----------------------------------------------------------
    
    hap_detection <- occurrences %>%
      
      distinct(
        period,
        ASV,
        sample
      ) %>%
      
      count(
        
        period,
        ASV,
        
        name =
          "n_detection_records"
        
      )
    
    
    dominance <- hap_detection %>%
      
      group_by(
        period
      ) %>%
      
      summarise(
        
        dominant_haplotype_frequency =
          
          max(
            n_detection_records
          ) /
          
          sum(
            n_detection_records
          ),
        
        .groups =
          "drop"
        
      )
    
    
    # ----------------------------------------------------------
    # Combine
    # ----------------------------------------------------------
    
    metrics <- period_reference %>%
      
      select(
        
        period,
        bin_start,
        bin_end,
        
        n_samples_available,
        n_sites_sampled
        
      ) %>%
      
      left_join(
        
        basic,
        
        by =
          "period"
        
      ) %>%
      
      left_join(
        
        spatial,
        
        by =
          "period"
        
      ) %>%
      
      left_join(
        
        dominance,
        
        by =
          "period"
        
      ) %>%
      
      mutate(
        
        n_haplotypes =
          replace_na(
            n_haplotypes,
            0L
          ),
        
        n_positive_samples_species =
          replace_na(
            n_positive_samples_species,
            0L
          ),
        
        n_positive_sites_species =
          replace_na(
            n_positive_sites_species,
            0L
          ),
        
        shared_haplotypes =
          replace_na(
            shared_haplotypes,
            0L
          ),
        
        site_specific_haplotypes =
          replace_na(
            site_specific_haplotypes,
            0L
          ),
        
        prop_shared_haplotypes =
          ifelse(
            
            n_haplotypes > 0,
            
            shared_haplotypes /
              n_haplotypes,
            
            NA_real_
            
          ),
        
        prop_site_specific_haplotypes =
          ifelse(
            
            n_haplotypes > 0,
            
            site_specific_haplotypes /
              n_haplotypes,
            
            NA_real_
            
          )
        
      )
    
    
    return(metrics)
    
  }
  
  
  
  # ============================================================
  # 25. MAKE TEMPORAL NETWORK FIGURE
  # ============================================================
  
  make_temporal_network_figure <- function(
    species_name,
    marker_name,
    status_value = NULL
  ) {
    
    marker_name <-
      toupper(marker_name)
    
    
    cat("\n")
    cat("=============================================\n")
    cat("TEMPORAL NETWORK\n")
    cat("=============================================\n")
    
    cat(
      "Species: ",
      species_name,
      "\n",
      sep = ""
    )
    
    cat(
      "Marker : ",
      marker_name,
      "\n",
      sep = ""
    )
    
    cat(
      "Status : ",
      status_value,
      "\n",
      sep = ""
    )
    
    cat("=============================================\n\n")
    
    
    # ==========================================================
    # STORED COMPLETE NETWORK
    # ==========================================================
    
    result <-
      get_temporal_network_result(
        
        species_name =
          species_name,
        
        marker_name =
          marker_name
        
      )
    
    
    # ==========================================================
    # FIXED NETWORK
    # ==========================================================
    
    fixed_network <-
      build_fixed_temporal_network(
        result
      )
    
    
    # ==========================================================
    # TEMPORAL OCCURRENCES
    # ==========================================================
    
    temporal_data <-
      build_species_temporal_data(
        
        species_name =
          species_name,
        
        marker_name =
          marker_name
        
      )
    
    
    occurrences <-
      temporal_data$occurrences
    
    
    sample_effort <-
      temporal_data$sample_effort
    
    
    # ==========================================================
    # GENERAL PERIODS
    # ==========================================================
    
    period_reference <-
      create_period_reference(
        
        sample_effort =
          sample_effort,
        
        occurrences =
          occurrences
        
      )
    
    
    if (nrow(period_reference) == 0) {
      
      stop(
        "No temporal periods available for plotting."
      )
      
    }
    
    
    # ==========================================================
    # NODE DATA
    # ==========================================================
    
    temporal_nodes <-
      prepare_temporal_nodes(
        
        occurrences =
          occurrences,
        
        period_reference =
          period_reference,
        
        fixed_network =
          fixed_network
        
      )
    
    
    node_data <-
      temporal_nodes$node_data
    
    
    sites_present <-
      temporal_nodes$sites_present
    
    
    # ==========================================================
    # SITE COLOUR CHECK
    # ==========================================================
    
    unknown_sites <-
      setdiff(
        
        sites_present,
        
        names(
          site_colors
        )
        
      )
    
    
    if (length(unknown_sites) > 0) {
      
      warning(
        "Undefined site colours: ",
        paste(
          unknown_sites,
          collapse = ", "
        ),
        ". Grey will be used."
      )
      
      
      extra_colours <-
        rep(
          "grey60",
          length(
            unknown_sites
          )
        )
      
      
      names(extra_colours) <-
        unknown_sites
      
      
      site_colors_local <-
        c(
          site_colors,
          extra_colours
        )
      
    } else {
      
      site_colors_local <-
        site_colors
      
    }
    
    
    # ==========================================================
    # PERIOD LABEL ORDER
    # ==========================================================
    
    panel_reference <- node_data %>%
      
      distinct(
        
        period_order,
        panel_label
        
      ) %>%
      
      arrange(
        period_order
      )
    
    
    panel_levels <-
      panel_reference$panel_label
    
    
    node_data <- node_data %>%
      
      mutate(
        
        panel_label =
          factor(
            
            panel_label,
            
            levels =
              panel_levels
            
          )
        
      )
    
    
    # ==========================================================
    # EDGES
    # ==========================================================
    
    edge_data <-
      prepare_temporal_edges(
        
        fixed_network =
          fixed_network,
        
        period_reference =
          period_reference,
        
        node_data =
          node_data
        
      )
    
    
    if (!is.null(edge_data)) {
      
      period_label_lookup <- node_data %>%
        
        distinct(
          period,
          panel_label
        )
      
      
      edge_data <- edge_data %>%
        
        left_join(
          
          period_label_lookup,
          
          by =
            "period"
          
        ) %>%
        
        mutate(
          
          panel_label =
            factor(
              panel_label,
              levels =
                panel_levels
            )
          
        )
      
    }
    
    
    # ==========================================================
    # BASE PLOT
    # ==========================================================
    
    p <- ggplot()
    
    
    # ----------------------------------------------------------
    # Complete fixed topology in light grey
    # ----------------------------------------------------------
    
    if (!is.null(edge_data)) {
      
      p <- p +
        
        geom_segment(
          
          data =
            edge_data,
          
          aes(
            
            x =
              x,
            
            y =
              y,
            
            xend =
              xend,
            
            yend =
              yend
            
          ),
          
          colour =
            "grey82",
          
          linewidth =
            0.45,
          
          lineend =
            "round"
          
        )
      
    }
    
    
    # ----------------------------------------------------------
    # Highlight edges whose two endpoints are detected
    # ----------------------------------------------------------
    
    if (!is.null(edge_data)) {
      
      p <- p +
        
        geom_segment(
          
          data =
            edge_data %>%
            
            filter(
              active_edge
            ),
          
          aes(
            
            x =
              x,
            
            y =
              y,
            
            xend =
              xend,
            
            yend =
              yend
            
          ),
          
          colour =
            "grey38",
          
          linewidth =
            0.65,
          
          lineend =
            "round"
          
        )
      
    }
    
    
    # ==========================================================
    # UNDETECTED HAPLOTYPES
    # ==========================================================
    
    inactive_nodes <- node_data %>%
      
      filter(
        !detected
      )
    
    
    p <- p +
      
      geom_point(
        
        data =
          inactive_nodes,
        
        aes(
          x = x,
          y = y
        ),
        
        shape =
          21,
        
        size =
          1.5,
        
        fill =
          "white",
        
        colour =
          "grey67",
        
        stroke =
          0.35
        
      )
    
    
    # ==========================================================
    # DETECTED HAPLOTYPES — SITE PIES
    # ==========================================================
    
    detected_nodes <- node_data %>%
      
      filter(
        detected
      )
    
    
    if (nrow(detected_nodes) > 0) {
      
      pie_columns <-
        sites_present[
          sites_present %in%
            names(
              detected_nodes
            )
        ]
      
      
      if (length(pie_columns) == 0) {
        
        stop(
          "No site columns available for scatterpie."
        )
        
      }
      
      
      p <- p +
        
        scatterpie::geom_scatterpie(
          
          data =
            detected_nodes,
          
          aes(
            
            x =
              x,
            
            y =
              y,
            
            r =
              radius
            
          ),
          
          cols =
            pie_columns,
          
          colour =
            "black",
          
          linewidth =
            0.35
          
        ) +
        
        scale_fill_manual(
          
          values =
            site_colors_local,
          
          breaks =
            site_order[
              site_order %in%
                sites_present
            ],
          
          limits =
            site_order[
              site_order %in%
                sites_present
            ],
          
          drop =
            FALSE,
          
          name =
            "Site"
          
        )
      
    }
    
    
    # ==========================================================
    # OPTIONAL ASV LABELS
    # ==========================================================
    
    if (SHOW_ASV_LABELS) {
      
      p <- p +
        
        geom_text(
          
          data =
            detected_nodes,
          
          aes(
            
            x =
              x,
            
            y =
              y,
            
            label =
              ASV
            
          ),
          
          nudge_y =
            0.045,
          
          size =
            2.0,
          
          check_overlap =
            TRUE
          
        )
      
    }
    
    
    # ==========================================================
    # FACETS — EXACTLY 4 NETWORKS PER ROW
    # ==========================================================
    
    p <- p +
      
      facet_wrap(
        
        ~ panel_label,
        
        ncol =
          NETWORKS_PER_ROW,
        
        drop =
          FALSE
        
      ) +
      
      coord_equal(
        
        xlim = c(
          -0.03,
          1.03
        ),
        
        ylim = c(
          -0.03,
          1.03
        ),
        
        expand =
          FALSE,
        
        clip =
          "off"
        
      )
    
    
    # ==========================================================
    # TITLE
    # ==========================================================
    
    subtitle_text <- paste0(
      
      status_value,
      " | ",
      marker_name,
      " | 50-year bins | fixed genetic topology"
      
    )
    
    
    p <- p +
      
      labs(
        
        title =
          species_name,
        
        subtitle =
          subtitle_text
        
      ) +
      
      guides(
        
        fill =
          guide_legend(
            
            title.position =
              "top",
            
            title.hjust =
              0.5,
            
            nrow =
              1,
            
            byrow =
              TRUE
            
          )
        
      ) +
      
      theme_void() +
      
      theme(
        
        plot.title =
          element_text(
            
            face =
              "italic",
            
            size =
              17,
            
            hjust =
              0
            
          ),
        
        plot.subtitle =
          element_text(
            
            size =
              10,
            
            colour =
              "grey30",
            
            margin =
              margin(
                b = 9
              )
            
          ),
        
        strip.background =
          element_rect(
            
            fill =
              "grey96",
            
            colour =
              "grey82",
            
            linewidth =
              0.35
            
          ),
        
        strip.text =
          element_text(
            
            face =
              "bold",
            
            size =
              9.5,
            
            lineheight =
              1.12,
            
            margin =
              margin(
                t = 5,
                r = 5,
                b = 5,
                l = 5
              )
            
          ),
        
        panel.spacing =
          grid::unit(
            0.70,
            "lines"
          ),
        
        legend.position =
          "bottom",
        
        legend.direction =
          "horizontal",
        
        legend.title =
          element_text(
            
            face =
              "bold",
            
            size =
              10
            
          ),
        
        legend.text =
          element_text(
            size = 9
          ),
        
        legend.key.width =
          grid::unit(
            1.0,
            "lines"
          ),
        
        legend.spacing.x =
          grid::unit(
            0.20,
            "cm"
          ),
        
        plot.margin =
          margin(
            10,
            12,
            10,
            12
          )
        
      )
    
    
    # ==========================================================
    # TEMPORAL METRICS — EXPORT ONLY
    # ==========================================================
    
    temporal_metrics <-
      calculate_temporal_metrics(
        
        occurrences =
          occurrences,
        
        period_reference =
          period_reference
        
      ) %>%
      
      mutate(
        
        species =
          species_name,
        
        marker =
          marker_name,
        
        status =
          status_value
        
      ) %>%
      
      select(
        
        species,
        marker,
        status,
        
        everything()
        
      )
    
    
    # ==========================================================
    # OUTPUT NAME
    # ==========================================================
    
    safe_species <- stringr::str_replace_all(
      
      species_name,
      
      "[^A-Za-z0-9]+",
      
      "_"
      
    )
    
    
    prefix <- paste0(
      
      safe_species,
      "_",
      marker_name,
      "_50yr_TEMPORAL_HAPLOTYPE_NETWORKS"
      
    )
    
    
    # ==========================================================
    # FIGURE DIMENSIONS
    # ==========================================================
    #
    # 4 panels per row.
    #
    # No timeline panel underneath.
    #
    # Even with ~20 temporal periods this gives:
    #
    # 20 periods / 4 = 5 rows
    #
    # so the figure stays manageable.
    #
    # ==========================================================
    
    n_panels <-
      length(
        panel_levels
      )
    
    
    n_rows <-
      ceiling(
        n_panels /
          NETWORKS_PER_ROW
      )
    
    
    figure_width <-
      12.8
    
    
    figure_height <-
      1.9 +
      (
        n_rows *
          2.75
      )
    
    
    # Safe cap — should normally never be reached.
    
    figure_height <-
      min(
        figure_height,
        45
      )
    
    
    cat(
      "Temporal panels: ",
      n_panels,
      "\n",
      sep = ""
    )
    
    
    cat(
      "Rows: ",
      n_rows,
      "\n",
      sep = ""
    )
    
    
    cat(
      "Figure size: ",
      figure_width,
      " x ",
      round(
        figure_height,
        1
      ),
      " inches\n",
      sep = ""
    )
    
    
    # ==========================================================
    # EXPORT SVG
    # ==========================================================
    
    ggsave(
      
      filename =
        file.path(
          
          temporal_output_dir,
          
          paste0(
            prefix,
            ".svg"
          )
          
        ),
      
      plot =
        p,
      
      width =
        figure_width,
      
      height =
        figure_height,
      
      units =
        "in",
      
      bg =
        "white"
      
    )
    
    
    # ==========================================================
    # EXPORT PDF
    # ==========================================================
    
    ggsave(
      
      filename =
        file.path(
          
          temporal_output_dir,
          
          paste0(
            prefix,
            ".pdf"
          )
          
        ),
      
      plot =
        p,
      
      width =
        figure_width,
      
      height =
        figure_height,
      
      units =
        "in",
      
      device =
        cairo_pdf,
      
      bg =
        "white"
      
    )
    
    
    # ==========================================================
    # EXPORT PNG
    # ==========================================================
    
    ggsave(
      
      filename =
        file.path(
          
          temporal_output_dir,
          
          paste0(
            prefix,
            ".png"
          )
          
        ),
      
      plot =
        p,
      
      width =
        figure_width,
      
      height =
        figure_height,
      
      units =
        "in",
      
      dpi =
        PNG_DPI,
      
      bg =
        "white"
      
    )
    
    
    # ==========================================================
    # EXPORT TEMPORAL METRICS
    # ==========================================================
    
    readr::write_csv(
      
      temporal_metrics,
      
      file.path(
        
        temporal_output_dir,
        
        paste0(
          prefix,
          "_metrics.csv"
        )
        
      )
      
    )
    
    
    # ==========================================================
    # EXPORT NODE DATA
    # ==========================================================
    
    readr::write_csv(
      
      node_data,
      
      file.path(
        
        temporal_output_dir,
        
        paste0(
          prefix,
          "_nodes.csv"
        )
        
      )
      
    )
    
    
    # ==========================================================
    # EXPORT OCCURRENCES
    # ==========================================================
    
    readr::write_csv(
      
      occurrences,
      
      file.path(
        
        temporal_output_dir,
        
        paste0(
          prefix,
          "_occurrences.csv"
        )
        
      )
      
    )
    
    
    cat(
      "Finished: ",
      species_name,
      " | ",
      marker_name,
      "\n\n",
      sep = ""
    )
    
    
    # ==========================================================
    # RETURN OBJECTS
    # ==========================================================
    
    return(
      
      list(
        
        plot =
          p,
        
        temporal_metrics =
          temporal_metrics,
        
        node_data =
          node_data,
        
        edge_data =
          edge_data,
        
        occurrences =
          occurrences,
        
        sample_effort =
          sample_effort,
        
        period_reference =
          period_reference,
        
        fixed_network =
          fixed_network,
        
        original_result =
          result
        
      )
      
    )
    
  }
  
  
  
  # ============================================================
  # 26. FOUR SELECTED SPECIES
  # ============================================================
  
  temporal_species <- tibble::tribble(
    
    ~species,
    ~marker,
    ~status,
    
    "Peringia ulvae",
    "COI",
    "CRY/NAT",
    
    "Botrylloides israeliense",
    "18S",
    "CRY",
    
    "Ascidia ahodori",
    "COI",
    "NIS",
    
    "Phaxas pellucidus",
    "18S",
    "NIS/NAT",
    
    "Cerastoderma edule",
    "COI",
    "NAT",
    
    "Cerastoderma edule",
    "18S",
    "NAT",
    
    "Ecteinascidia turbinata",
    "COI",
    "NIS",
    
    "Obelia geniculata",
    "COI",
    "CRY",
    
    "Styela plicata",
    "COI",
    "NIS",
  )
  
  
  
  # ============================================================
  # 27. CHECK THE FOUR NETWORKS EXIST
  # ============================================================
  
  network_check <- temporal_species %>%
    
    rowwise() %>%
    
    mutate(
      
      species_exists =
        species %in%
        names(
          all_haplotype_results
        ),
      
      marker_exists =
        
        ifelse(
          
          species_exists,
          
          marker %in%
            names(
              all_haplotype_results[[species]]
            ),
          
          FALSE
          
        )
      
    ) %>%
    
    ungroup()
  
  
  cat("\n")
  cat("=============================================\n")
  cat("SELECTED SPECIES CHECK\n")
  cat("=============================================\n\n")
  
  print(
    network_check
  )
  
  
  
  if (
    any(
      !network_check$species_exists |
      !network_check$marker_exists
    )
  ) {
    
    stop(
      paste0(
        "\nAt least one selected species-marker network is ",
        "missing from all_haplotype_results.\n",
        "Check the table printed above before continuing.\n"
      )
    )
    
  }
  
  
  
  # ============================================================
  # 28. RUN ALL FOUR SPECIES AUTOMATICALLY
  # ============================================================
  
  all_temporal_networks <- list()
  
  
  for (
    ii in seq_len(
      nrow(
        temporal_species
      )
    )
  ) {
    
    species_i <-
      temporal_species$species[ii]
    
    
    marker_i <-
      temporal_species$marker[ii]
    
    
    status_i <-
      temporal_species$status[ii]
    
    
    key_i <- paste(
      
      species_i,
      
      marker_i,
      
      sep =
        "___"
      
    )
    
    
    cat("\n\n")
    cat("#############################################\n")
    cat("# RUNNING ", species_i, "\n", sep = "")
    cat("# MARKER  ", marker_i, "\n", sep = "")
    cat("#############################################\n\n")
    
    
    all_temporal_networks[[key_i]] <-
      
      make_temporal_network_figure(
        
        species_name =
          species_i,
        
        marker_name =
          marker_i,
        
        status_value =
          status_i
        
      )
    
  }
  
  
  
  # ============================================================
  # 29. PRINT ALL FOUR FIGURES IN R
  # ============================================================
  
  print(
    all_temporal_networks[["Peringia ulvae___COI"]]$plot
  )
  
  print(
    all_temporal_networks[["Botrylloides israeliense___18S"]]$plot
  )
  
  print(
    all_temporal_networks[["Ascidia ahodori___COI"]]$plot
  )
  
  print(
    all_temporal_networks[["Phaxas pellucidus___18S"]]$plot
  )
  
  
  
  # ============================================================
  # 30. COMBINE ALL TEMPORAL METRICS
  # ============================================================
  
  all_temporal_metrics <- bind_rows(
    
    lapply(
      
      all_temporal_networks,
      
      function(x) {
        
        x$temporal_metrics
        
      }
      
    )
    
  )
  
  
  readr::write_csv(
    
    all_temporal_metrics,
    
    file.path(
      
      temporal_output_dir,
      
      "ALL_FOUR_SPECIES_50yr_temporal_metrics.csv"
      
    )
    
  )
  
  
  
  # ============================================================
  # 31. SUMMARY IN CONSOLE
  # ============================================================
  
  cat("\n\n")
  cat("=============================================\n")
  cat("TEMPORAL HAPLOTYPE ANALYSIS FINISHED\n")
  cat("=============================================\n\n")
  
  
  summary_table <- all_temporal_metrics %>%
    
    group_by(
      species,
      marker,
      status
    ) %>%
    
    summarise(
      
      first_detected_period =
        ifelse(
          
          any(
            n_haplotypes > 0
          ),
          
          min(
            bin_start[
              n_haplotypes > 0
            ]
          ),
          
          NA_real_
          
        ),
      
      last_detected_period =
        ifelse(
          
          any(
            n_haplotypes > 0
          ),
          
          max(
            bin_end[
              n_haplotypes > 0
            ]
          ),
          
          NA_real_
          
        ),
      
      maximum_haplotype_richness =
        max(
          n_haplotypes,
          na.rm = TRUE
        ),
      
      .groups =
        "drop"
      
    )
  
  
  print(
    summary_table
  )
  
  
  cat("\n")
  cat("Figures saved in:\n")
  cat(temporal_output_dir)
  cat("\n\n")
  
  
  
  # ============================================================
  # 32. END
  # ============================================================