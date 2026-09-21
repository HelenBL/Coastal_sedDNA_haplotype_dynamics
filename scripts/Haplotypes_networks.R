# ============================================================
# HAPLOTYPE NETWORK FOR ONE SELECTED SPECIES
# ============================================================
#
# USAGE:
#
#   species_name <- "Abra segmentum"
#   marker_name  <- "18S"
#
# Then run:
#
#   make_haplotype_network(species_name, marker_name)
#
# Output:
#   Species_18S_haplotype_network.svg
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

library(Biostrings)
library(DECIPHER)
library(ape)
library(igraph)
library(tidygraph)
library(ggraph)
library(scatterpie)
library(scales)


# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================
setwd("C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Figures")

haplo_dir <- "Haplotype_networks"

dir.create(
  haplo_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. SITE COLOURS
# ============================================================

site_colors <- c(
  "DEE" = "#CDCD00",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  "COL" = "burlywood1",
  "CSD" = "#76EEC6",
  "CMB" = "#458B74",
  "HIT" = "#CD96CD",
  "HIO" = "#8B0A50"
)



# ============================================================
# 4. READ ASV SEQUENCES
# ============================================================
#
# IMPORTANT:
#
# Change these two filenames to your actual FASTA files
# containing the original ASV sequences.
#
# Headers can be:
#
# >ASV_123
#
# or more complicated, e.g.
#
# >something|18S|ASV_123
#
# because the function extracts ASV_xxx automatically.
#
# ============================================================

fasta_18S <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Data/fastas/DADA2.ASVs_merged_18S_all_peninsula_solent.fasta"

fasta_COI_1 <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Data/fastas/DADA2.ASVs_merged_COI_all_peninsula_solent.fasta"

fasta_COI_2 <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Data/fastas/STO_CAD_unique_new_ASVs_for_blast_COI.fasta"


read_asv_fasta <- function(files) {
  
  all_sequences <- lapply(files, function(file) {
    
    message("Reading FASTA: ", file)
    
    seqs <- Biostrings::readDNAStringSet(file)
    
    headers <- names(seqs)
    
    # Extract ASV_xxx from FASTA header
    ASV <- stringr::str_extract(
      headers,
      "ASV_[0-9]+"
    )
    
    tibble::tibble(
      ASV = ASV,
      sequence = as.character(seqs),
      source_fasta = basename(file)
    )
  })
  
  seq_df <- dplyr::bind_rows(all_sequences)
  
  # Remove headers where ASV could not be identified
  if (any(is.na(seq_df$ASV))) {
    
    warning(
      sum(is.na(seq_df$ASV)),
      " sequences did not contain an ASV_xxx identifier."
    )
  }
  
  seq_df <- seq_df %>%
    dplyr::filter(!is.na(ASV))
  
  # ----------------------------------------------------------
  # Check if the same ASV occurs more than once
  # ----------------------------------------------------------
  
  duplicated_asvs <- seq_df %>%
    dplyr::count(ASV) %>%
    dplyr::filter(n > 1)
  
  if (nrow(duplicated_asvs) > 0) {
    
    message(
      "\nDuplicated ASVs found across FASTA files: ",
      nrow(duplicated_asvs)
    )
    
    print(duplicated_asvs)
  }
  
  # ----------------------------------------------------------
  # Important check:
  # same ASV should have same sequence
  # ----------------------------------------------------------
  
  conflicting_asvs <- seq_df %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(
      n_sequences = dplyr::n_distinct(sequence),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_sequences > 1)
  
  if (nrow(conflicting_asvs) > 0) {
    
    print(conflicting_asvs)
    
    stop(
      "Some ASVs occur in more than one FASTA with different sequences."
    )
  }
  
  # If duplicated but sequence is identical, keep one
  seq_df <- seq_df %>%
    dplyr::distinct(
      ASV,
      sequence,
      .keep_all = TRUE
    ) %>%
    dplyr::distinct(
      ASV,
      .keep_all = TRUE
    )
  
  message(
    "Unique ASVs loaded: ",
    nrow(seq_df)
  )
  
  return(seq_df)
}


# ============================================================
# 5. FUNCTION TO IDENTIFY SAMPLE COLUMNS
# ============================================================

get_sample_columns <- function(df) {
  
  # A sample column starts with one of the site codes:
  #
  # CMB_...
  # CSD_...
  # CCO_...
  # etc.
  
  pattern <- paste0(
    "^(",
    paste(names(site_colors), collapse = "|"),
    ")_"
  )
  
  names(df)[
    stringr::str_detect(names(df), pattern)
  ]
}



# ============================================================
# 6. MAIN FUNCTION
# ============================================================

make_haplotype_network <- function(
    species_name,
    marker_name,
    status_value = NULL
) {
  
  # ----------------------------------------------------------
  # Normalise input
  # ----------------------------------------------------------
  
  marker_name <- toupper(marker_name)
  
  species_name <- stringr::str_squish(
    species_name
  )
  
  
  # ----------------------------------------------------------
  # Select correct dataset and FASTA
  # ----------------------------------------------------------
  
  if (marker_name == "18S") {
    
    abundance_df <- data_18S_updated
    
    fasta_files <- c(
      fasta_18S
    )
    
  } else if (marker_name == "COI") {
    
    abundance_df <- data_COI_updated
    
    fasta_files <- c(
      fasta_COI_1,
      fasta_COI_2
    )
    
  } else {
    
    stop(
      "marker_name must be either '18S' or 'COI'."
    )
  }
  
  
  # ----------------------------------------------------------
  # Check required columns
  # ----------------------------------------------------------
  
  required_cols <- c(
    "ASV",
    "species"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(abundance_df)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "Missing columns in abundance table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # Species available?
  # ----------------------------------------------------------
  
  species_available <- abundance_df %>%
    dplyr::filter(
      !is.na(species),
      species != ""
    ) %>%
    dplyr::pull(species) %>%
    stringr::str_squish() %>%
    unique()
  
  
  if (!species_name %in% species_available) {
    
    genus_search <- stringr::word(
      species_name,
      1
    )
    
    similar_species <- species_available[
      stringr::str_detect(
        species_available,
        stringr::regex(
          genus_search,
          ignore_case = TRUE
        )
      )
    ]
    
    message(
      "\nSpecies not found: ",
      species_name
    )
    
    if (length(similar_species) > 0) {
      
      message("\nPossible matches:")
      
      print(similar_species)
    }
    
    stop(
      "Species not found in the selected dataset."
    )
  }
  
  
  # ----------------------------------------------------------
  # Select ASVs belonging to selected species
  # ----------------------------------------------------------
  
  df_sp <- abundance_df %>%
    dplyr::mutate(
      species_clean = stringr::str_squish(
        as.character(species)
      )
    ) %>%
    dplyr::filter(
      species_clean == species_name
    )
  
  
  message(
    "\n========================================"
  )
  
  message(
    "Species: ",
    species_name
  )
  
  message(
    "Marker: ",
    marker_name
  )
  
  message(
    "ASVs found: ",
    dplyr::n_distinct(df_sp$ASV)
  )
  
  message(
    "========================================"
  )
  
  
  if (nrow(df_sp) == 0) {
    
    stop(
      "No ASVs found for ",
      species_name,
      " in ",
      marker_name
    )
  }
  
  
  # ----------------------------------------------------------
  # Identify sample columns
  # ----------------------------------------------------------
  
  sample_cols <- get_sample_columns(
    df_sp
  )
  
  if (length(sample_cols) == 0) {
    
    stop(
      "No sample columns were detected."
    )
  }
  
  message(
    "Sample columns detected: ",
    length(sample_cols)
  )
  
  
  # ----------------------------------------------------------
  # Convert abundance data to long format
  # ----------------------------------------------------------
  
  abundance_long <- df_sp %>%
    dplyr::select(
      ASV,
      dplyr::all_of(sample_cols)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(sample_cols),
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    dplyr::mutate(
      
      abundance = suppressWarnings(
        as.numeric(abundance)
      ),
      
      Site = stringr::str_extract(
        sample,
        "^[^_]+"
      )
      
    ) %>%
    dplyr::filter(
      !is.na(abundance),
      abundance > 0
    )
  
  
  # ----------------------------------------------------------
  # Sum abundance per ASV and site
  # ----------------------------------------------------------
  
  site_abundance <- abundance_long %>%
    dplyr::group_by(
      ASV,
      Site
    ) %>%
    dplyr::summarise(
      abundance = sum(
        abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  message("\nSites detected:")
  
  print(
    site_abundance %>%
      dplyr::group_by(Site) %>%
      dplyr::summarise(
        abundance = sum(
          abundance,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
  )
  
  
  # ==========================================================
  # READ SEQUENCES
  # ==========================================================
  
  seq_df <- read_asv_fasta(
    fasta_files
  )
  
  
  # ----------------------------------------------------------
  # Keep sequences corresponding to selected ASVs
  # ----------------------------------------------------------
  
  asv_sequences <- df_sp %>%
    dplyr::distinct(ASV) %>%
    dplyr::left_join(
      seq_df,
      by = "ASV"
    )
  
  
  # ----------------------------------------------------------
  # Check missing sequences
  # ----------------------------------------------------------
  
  missing_sequences <- asv_sequences %>%
    dplyr::filter(
      is.na(sequence) |
        sequence == ""
    )
  
  
  if (nrow(missing_sequences) > 0) {
    
    warning(
      "\nNo sequence found for these ASVs:\n",
      paste(
        missing_sequences$ASV,
        collapse = ", "
      )
    )
  }
  
  
  asv_sequences <- asv_sequences %>%
    dplyr::filter(
      !is.na(sequence),
      sequence != ""
    )
  
  
  message(
    "\nASVs with sequence: ",
    nrow(asv_sequences)
  )
  
  
  if (nrow(asv_sequences) == 0) {
    
    stop(
      "No sequences available for the selected species."
    )
  }
  
  
  # ==========================================================
  # CASE: ONLY ONE ASV
  # ==========================================================
  
  if (nrow(asv_sequences) == 1) {
    
    message(
      "\nOnly one ASV found. Creating single-haplotype plot."
    )
    
    single_asv <- asv_sequences$ASV[1]
    
    
    pie_df <- site_abundance %>%
      dplyr::filter(
        ASV == single_asv
      ) %>%
      tidyr::pivot_wider(
        names_from = Site,
        values_from = abundance,
        values_fill = 0
      )
    
    
    missing_sites <- setdiff(
      names(site_colors),
      names(pie_df)
    )
    
    
    for (ss in missing_sites) {
      
      pie_df[[ss]] <- 0
      
    }
    
    
    pie_df <- pie_df %>%
      dplyr::select(
        ASV,
        dplyr::all_of(
          names(site_colors)
        )
      ) %>%
      dplyr::mutate(
        x = 0,
        y = 0,
        radius = 0.18
      )
    
    
    p <- ggplot2::ggplot() +
      
      scatterpie::geom_scatterpie(
        data = pie_df,
        ggplot2::aes(
          x = x,
          y = y,
          r = radius
        ),
        cols = names(site_colors),
        colour = "black",
        linewidth = 0.6
      ) +
      
      ggplot2::geom_text(
        data = pie_df,
        ggplot2::aes(
          x = x,
          y = y + 0.27,
          label = ASV
        ),
        size = 4.5,
        fontface = "bold"
      ) +
      
      ggplot2::scale_fill_manual(
        values = site_colors,
        name = "Site",
        drop = TRUE
      ) +
      
      ggplot2::coord_equal(
        xlim = c(-0.5, 0.8),
        ylim = c(-0.4, 0.5)
      ) +
      
      ggplot2::labs(
        title = species_name,
        subtitle = paste0(
          marker_name,
          " — single haplotype"
        )
      ) +
      
      ggplot2::theme_void() +
      
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "italic",
          size = 18,
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 12,
          hjust = 0.5
        ),
        legend.position = "right"
      )
    
    
    safe_species <- stringr::str_replace_all(
      species_name,
      "[^A-Za-z0-9]+",
      "_"
    )
    
    # Clean status for filename
    if (!is.null(status_value) &&
        !is.na(status_value) &&
        status_value != "") {
      
      safe_status <- stringr::str_replace_all(
        status_value,
        "[^A-Za-z0-9]+",
        "_"
      )
      
      outfile <- file.path(
        haplo_dir,
        paste0(
          safe_species,
          "_",
          safe_status,
          "_",
          marker_name,
          "_haplotype_network.svg"
        )
      )
      
    } else {
      
      outfile <- file.path(
        haplo_dir,
        paste0(
          safe_species,
          "_",
          marker_name,
          "_haplotype_network.svg"
        )
      )
    }
    
    
    ggplot2::ggsave(
      filename = outfile,
      plot = p,
      width = 7,
      height = 6,
      bg = "white"
    )
    
    
    message(
      "\nSVG saved to:\n",
      outfile
    )
    
    print(p)
    
    
    return(
      invisible(
        list(
          plot = p,
          species = species_name,
          marker = marker_name,
          ASVs = asv_sequences,
          site_abundance = site_abundance,
          aligned_sequences = NULL,
          distances = NULL
        )
      )
    )
  }
  
  
  # ==========================================================
  # ALIGN SEQUENCES
  # ==========================================================
  
  message(
    "\nAligning sequences..."
  )
  
  
  seqs <- Biostrings::DNAStringSet(
    asv_sequences$sequence
  )
  
  names(seqs) <- asv_sequences$ASV
  
  
  alignment <- DECIPHER::AlignSeqs(
    seqs,
    processors = NULL
  )
  
  
  aligned_df <- tibble::tibble(
    ASV = names(alignment),
    aligned_sequence = as.character(
      alignment
    )
  )
  
  
  # ==========================================================
  # GENETIC DISTANCES
  # ==========================================================
  
  dna <- ape::as.DNAbin(
    strsplit(
      aligned_df$aligned_sequence,
      split = ""
    )
  )
  
  names(dna) <- aligned_df$ASV
  
  
  dist_matrix <- ape::dist.dna(
    dna,
    model = "N",
    pairwise.deletion = TRUE,
    as.matrix = TRUE
  )
  
  
  edge_df <- as.data.frame(
    as.table(dist_matrix),
    stringsAsFactors = FALSE
  )
  
  
  names(edge_df) <- c(
    "from",
    "to",
    "distance"
  )
  
  
  edge_df <- edge_df %>%
    dplyr::mutate(
      from = as.character(from),
      to = as.character(to),
      distance = as.numeric(distance)
    ) %>%
    dplyr::filter(
      from < to,
      !is.na(distance)
    )
  
  
  # ==========================================================
  # BUILD MINIMUM SPANNING NETWORK
  # ==========================================================
  
  graph_full <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = aligned_df %>%
      dplyr::transmute(
        name = ASV
      )
  )
  
  
  graph_net <- igraph::mst(
    graph_full,
    weights = igraph::E(
      graph_full
    )$distance
  )
  
  
  # ==========================================================
  # PREPARE SITE PIE CHARTS
  # ==========================================================
  
  site_pie <- site_abundance %>%
    tidyr::pivot_wider(
      names_from = Site,
      values_from = abundance,
      values_fill = 0
    )
  
  
  missing_sites <- setdiff(
    names(site_colors),
    names(site_pie)
  )
  
  
  for (ss in missing_sites) {
    
    site_pie[[ss]] <- 0
    
  }
  
  
  site_pie <- site_pie %>%
    dplyr::select(
      ASV,
      dplyr::all_of(
        names(site_colors)
      )
    ) %>%
    dplyr::mutate(
      
      total_abundance = rowSums(
        dplyr::across(
          dplyr::all_of(
            names(site_colors)
          )
        )
      )
      
    )
  
  
  # ----------------------------------------------------------
  # Node size
  # ----------------------------------------------------------
  
  if (
    length(
      unique(
        site_pie$total_abundance
      )
    ) == 1
  ) {
    
    site_pie$radius <- 0.12
    
  } else {
    
    site_pie$radius <- scales::rescale(
      sqrt(
        site_pie$total_abundance
      ),
      to = c(
        0.07,
        0.18
      )
    )
  }
  
  
  node_attr <- site_pie %>%
    dplyr::rename(
      name = ASV
    )
  
  
  graph_tbl <- tidygraph::as_tbl_graph(
    graph_net
  ) %>%
    tidygraph::activate(nodes) %>%
    dplyr::left_join(
      node_attr,
      by = "name"
    )
  
  
  # ==========================================================
  # PLOT NETWORK
  # ==========================================================
  
  p <- ggraph::ggraph(
    graph_tbl,
    layout = "stress"
  ) +
    
    ggraph::geom_edge_link(
      ggplot2::aes(
        width = distance
      ),
      colour = "grey40",
      alpha = 0.8
    ) +
    
    scatterpie::geom_scatterpie(
      ggplot2::aes(
        x = x,
        y = y,
        r = radius
      ),
      cols = names(site_colors),
      colour = "black",
      linewidth = 0.5
    ) +
    
    ggraph::geom_node_text(
      ggplot2::aes(
        label = name
      ),
      repel = TRUE,
      size = 3.8,
      fontface = "bold"
    ) +
    
    ggplot2::scale_fill_manual(
      values = site_colors,
      name = "Site",
      drop = TRUE
    ) +
    
    ggraph::scale_edge_width(
      range = c(
        0.5,
        2.2
      ),
      name = "Sequence distance"
    ) +
    
    ggplot2::coord_equal() +
    
    ggplot2::labs(
      title = species_name,
      subtitle = paste0(
        marker_name,
        " haplotype network"
      ),
      caption = paste0(
        "Pie sectors = abundance by site; ",
        "node size = total abundance"
      )
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      
      plot.title = ggplot2::element_text(
        size = 20,
        face = "italic",
        hjust = 0.5
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 12,
        hjust = 0.5
      ),
      
      plot.caption = ggplot2::element_text(
        size = 9
      ),
      
      legend.position = "right",
      
      legend.title = ggplot2::element_text(
        face = "bold"
      )
    )
  
  
  # ==========================================================
  # SAVE SVG
  # ==========================================================
  
  safe_species <- stringr::str_replace_all(
    species_name,
    "[^A-Za-z0-9]+",
    "_"
  )
  
  
  outfile <- file.path(
    haplo_dir,
    paste0(
      safe_species,
      "_",
      marker_name,
      "_haplotype_network.svg"
    )
  )
  
  
  ggplot2::ggsave(
    filename = outfile,
    plot = p,
    width = 9,
    height = 8,
    bg = "white"
  )
  
  
  message(
    "\n========================================"
  )
  
  message(
    "SVG saved to:\n",
    outfile
  )
  
  message(
    "========================================"
  )
  
  
  print(p)
  
  
  # ==========================================================
  # RETURN RESULTS
  # ==========================================================
  
  invisible(
    list(
      plot = p,
      species = species_name,
      marker = marker_name,
      ASVs = asv_sequences,
      aligned_sequences = aligned_df,
      site_abundance = site_abundance,
      distances = edge_df
    )
  )
  
}


# ============================================================
# 15. CHOOSE SPECIES HERE
# ============================================================

species_name <- "Carcinus aestuarii"
marker_name  <- "COI"


# ============================================================
# 16. CREATE NETWORK
# ============================================================

haplo_result <- make_haplotype_network(
  species_name,
  marker_name
)


# ============================================================
# 17. AUTOMATIC HAPLOTYPE NETWORKS FOR ALL SELECTED SPECIES
# ============================================================

# ------------------------------------------------------------
# Read species status table
# ------------------------------------------------------------

status_haplo <- readxl::read_excel(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript/Data/Status_selected_species.xlsx"
) %>%
  janitor::clean_names() %>%
  dplyr::transmute(
    
    species_name = stringr::str_squish(
      as.character(especie)
    ),
    
    status_value = stringr::str_squish(
      as.character(status_combinado)
    )
    
  ) %>%
  dplyr::filter(
    !is.na(species_name),
    species_name != ""
  ) %>%
  dplyr::distinct(
    species_name,
    .keep_all = TRUE
  )


# Check
print(status_haplo)

cat(
  "\nSpecies to process:",
  nrow(status_haplo),
  "\n"
)


# ============================================================
# FUNCTION TO PROCESS ONE SPECIES AUTOMATICALLY
# ============================================================

run_species_haplotype <- function(
    species_name,
    status_value
) {
  
  cat(
    "\n\n============================================================\n"
  )
  
  cat(
    "PROCESSING: ",
    species_name,
    " | Status: ",
    status_value,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================================\n"
  )
  
  
  # ----------------------------------------------------------
  # Determine which markers are available for this species
  # from the NEW BLAST taxonomy
  # ----------------------------------------------------------
  
  markers_available <- new_taxonomy %>%
    dplyr::mutate(
      species_clean = stringr::str_squish(
        as.character(species)
      ),
      marker = toupper(
        as.character(marker)
      )
    ) %>%
    dplyr::filter(
      species_clean == species_name
    ) %>%
    dplyr::pull(marker) %>%
    unique()
  
  
  markers_available <- intersect(
    markers_available,
    c("18S", "COI")
  )
  
  
  # ----------------------------------------------------------
  # Species absent from BLAST taxonomy
  # ----------------------------------------------------------
  
  if (length(markers_available) == 0) {
    
    message(
      "Skipping ",
      species_name,
      ": no 18S/COI ASVs found in new taxonomy."
    )
    
    return(NULL)
  }
  
  
  message(
    "Markers available: ",
    paste(
      markers_available,
      collapse = ", "
    )
  )
  
  
  # ----------------------------------------------------------
  # Generate one network for each available marker
  # ----------------------------------------------------------
  
  results_species <- list()
  
  
  for (mk in markers_available) {
    
    message(
      "\nCreating network: ",
      species_name,
      " - ",
      mk
    )
    
    
    result <- tryCatch(
      
      {
        
        make_haplotype_network(
          species_name = species_name,
          marker_name = mk,
          status_value = status_value
        )
        
      },
      
      error = function(e) {
        
        message(
          "\nERROR for ",
          species_name,
          " - ",
          mk,
          ":\n",
          e$message
        )
        
        return(NULL)
      }
      
    )
    
    
    results_species[[mk]] <- result
  }
  
  
  return(
    invisible(results_species)
  )
}


# ============================================================
# RUN ALL SPECIES
# ============================================================

all_haplotype_results <- list()


for (i in seq_len(nrow(status_haplo))) {
  
  sp <- status_haplo$species_name[i]
  
  st <- status_haplo$status_value[i]
  
  
  all_haplotype_results[[sp]] <- run_species_haplotype(
    species_name = sp,
    status_value = st
  )
}


# ============================================================
# FINISHED
# ============================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "ALL SPECIES FINISHED\n"
)

cat(
  "Networks saved in:\n",
  normalizePath(haplo_dir),
  "\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# ADD MANUALLY GENERATED NETWORKS TO all_haplotype_results
# ============================================================

# Carcinus aestuarii was generated separately and therefore
# was not included in the original all_haplotype_results list.

if (!exists("all_haplotype_results")) {
  stop("all_haplotype_results does not exist.")
}

if (!exists("haplo_result")) {
  stop(
    "haplo_result does not exist. ",
    "Run make_haplotype_network() for Carcinus aestuarii first."
  )
}

# Check that haplo_result is actually Carcinus aestuarii 18S
if (
  haplo_result$species != "Carcinus aestuarii" ||
  toupper(haplo_result$marker) != "18S"
) {
  
  stop(
    "haplo_result is not Carcinus aestuarii 18S.\n",
    "Current object contains: ",
    haplo_result$species,
    " | ",
    haplo_result$marker
  )
}

# Create species entry if it does not already exist
if (!"Carcinus aestuarii" %in% names(all_haplotype_results)) {
  
  all_haplotype_results[["Carcinus aestuarii"]] <- list()
}

# Add the 18S network result
all_haplotype_results[["Carcinus aestuarii"]][["18S"]] <- haplo_result


# ============================================================
# CHECK
# ============================================================

cat(
  "\nCarcinus aestuarii added to all_haplotype_results:\n"
)

print(
  names(
    all_haplotype_results[["Carcinus aestuarii"]]
  )
)

print(
  all_haplotype_results[["Carcinus aestuarii"]][["18S"]]$species
)

print(
  all_haplotype_results[["Carcinus aestuarii"]][["18S"]]$marker
)

cat(
  "\nNumber of ASVs: ",
  nrow(
    all_haplotype_results[["Carcinus aestuarii"]][["18S"]]$ASVs
  ),
  "\n",
  sep = ""
)



# ============================================================
# 18. FINAL MAIN FIGURE — SELECTED HAPLOTYPE NETWORKS
# ============================================================
#
# This section must be run AFTER:
#
#   all_haplotype_results
#
# has been created.
#
# FINAL FIGURE:
#   - selected ecologically informative species
#   - one marker per species
#   - same genetic-distance scale across ALL networks
#   - one common Site legend
#   - one common Sequence distance legend
#   - species names italicised
#   - status shown in each panel
#   - quantitative metrics shown below species name:
#
#       COI | nH = 4 | Dom = 0.65 | Shared = 50%
#
# IMPORTANT:
#   Network topology and node composition are taken from the
#   existing make_haplotype_network() results.
#
# ============================================================


# ============================================================
# 18.1. PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scatterpie)
library(scales)
library(readr)
library(tibble)
library(patchwork)


# ============================================================
# 18.2. OUTPUT DIRECTORY
# ============================================================

main_haplo_dir <- file.path(
  haplo_dir,
  "MAIN_FINAL_HAPLOTYPE_FIGURE"
)

dir.create(
  main_haplo_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 18.3. HAPLOTYPE METRICS FILE
# ============================================================
#
# Generated by:
#   Haplotype_structure_quantification.R
#
# ============================================================

haplotype_metrics_file <- file.path(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Haplotypes_manuscript",
  "Haplotype_structure_quantification",
  "Haplotype_structure_species_marker_summary.csv"
)


if (!file.exists(haplotype_metrics_file)) {
  
  stop(
    "Haplotype metrics file not found:\n",
    haplotype_metrics_file
  )
}


haplotype_metrics_main <- readr::read_csv(
  haplotype_metrics_file,
  show_col_types = FALSE
) %>%
  
  dplyr::mutate(
    
    species = stringr::str_squish(
      as.character(species)
    ),
    
    marker = toupper(
      as.character(marker)
    ),
    
    status = stringr::str_squish(
      as.character(status)
    )
    
  )


# ============================================================
# 18.4. SELECT NETWORKS FOR THE MAIN FIGURE
# ============================================================
#
# One row = one panel.
#
# The order here is the order used in the final figure.
#
# We prioritise COI where an informative COI network exists.
# 18S is retained for selected taxa where it provides the
# available / biologically informative contrast.
#
# To remove a species, simply delete its row.
#
# ============================================================

selected_main_networks <- tibble::tribble(
  
  ~panel_order, ~status,   ~species,                     ~marker,
  
  # ----------------------------------------------------------
  # NAT
  # ----------------------------------------------------------
  
  1,           "NAT",     "Cerastoderma edule",         "COI",
  2,           "NAT",     "Cerastoderma edule",         "18S",
  3,           "NAT",     "Phallusia mammillata",       "COI",
  
  
  # ----------------------------------------------------------
  # NIS
  # ----------------------------------------------------------
  
  4,           "NIS",     "Ascidia ahodori",            "COI",
  5,           "NIS",     "Styela plicata",             "COI",
  6,           "NIS",     "Ecteinascidia turbinata",    "COI",
  7,           "NIS",     "Diadumene lineata",    "18S",
  
  # ----------------------------------------------------------
  # CRY/NAT
  # ----------------------------------------------------------
  
  8,           "CRY/NAT", "Peringia ulvae",             "COI",
  9,           "CRY/NAT", "Botryllus schlosseri",       "COI",
  
  
  # ----------------------------------------------------------
  # CRY
  # ----------------------------------------------------------
  
  10,           "CRY",     "Botrylloides israeliense",   "18S",
  11,           "CRY",     "Amathia verticillata",       "18S",
  12,           "CRY",     "Obelia geniculata",          "COI",
  
  
  # ----------------------------------------------------------
  # NIS/NAT
  # ----------------------------------------------------------
  
  13,           "NIS/NAT", "Onchidella celtica",         "COI",
  14,           "NIS/NAT", "Phaxas pellucidus",          "18S",
  
  
  # ----------------------------------------------------------
  # NAT/NIS
  # ----------------------------------------------------------
  
  15,           "NAT/NIS", "Carcinus aestuarii",         "18S"
  
) %>%
  
  dplyr::mutate(
    
    marker = toupper(marker),
    
    species = stringr::str_squish(species),
    
    status = stringr::str_squish(status)
    
  ) %>%
  
  dplyr::arrange(panel_order)


print(selected_main_networks)


# ============================================================
# 18.5. CHECK THAT METRICS EXIST
# ============================================================

selected_metrics_check <- selected_main_networks %>%
  
  dplyr::left_join(
    
    haplotype_metrics_main %>%
      dplyr::select(
        species,
        marker,
        status_metrics = status,
        n_haplotypes,
        dominant_haplotype_frequency,
        prop_shared_haplotypes,
        prop_site_specific_haplotypes,
        sites_per_haplotype
      ),
    
    by = c(
      "species",
      "marker"
    )
    
  )


missing_metrics <- selected_metrics_check %>%
  dplyr::filter(
    is.na(n_haplotypes)
  )


if (nrow(missing_metrics) > 0) {
  
  print(missing_metrics)
  
  stop(
    "Some selected species/marker combinations are missing ",
    "from Haplotype_structure_species_marker_summary.csv."
  )
}


# ------------------------------------------------------------
# Warn if status in the main selection differs from the
# quantitative status table.
#
# We do NOT silently change it.
# ------------------------------------------------------------

status_disagreement <- selected_metrics_check %>%
  dplyr::filter(
    !is.na(status_metrics),
    status != status_metrics
  )


if (nrow(status_disagreement) > 0) {
  
  message(
    "\nWARNING: status disagreement between selected figure ",
    "and quantitative metrics:"
  )
  
  print(
    status_disagreement %>%
      dplyr::select(
        species,
        marker,
        status_selected = status,
        status_metrics
      )
  )
  
}


# ============================================================
# 18.6. FORMAT METRIC LABEL
# ============================================================

selected_metrics_main <- selected_metrics_check %>%
  
  dplyr::mutate(
    
    metric_label = paste0(
      
      marker,
      
      " | nH = ",
      n_haplotypes,
      
      " | Dom = ",
      sprintf(
        "%.2f",
        dominant_haplotype_frequency
      ),
      
      " | Shared = ",
      sprintf(
        "%.0f%%",
        100 * prop_shared_haplotypes
      )
      
    )
    
  )


print(
  selected_metrics_main %>%
    dplyr::select(
      panel_order,
      status,
      species,
      marker,
      n_haplotypes,
      dominant_haplotype_frequency,
      prop_shared_haplotypes,
      metric_label
    )
)


# ============================================================
# 18.7. FUNCTION TO RETRIEVE STORED NETWORK RESULT
# ============================================================

get_main_haplotype_result <- function(
    species_name,
    marker_name
) {
  
  marker_name <- toupper(marker_name)
  
  
  if (!species_name %in% names(all_haplotype_results)) {
    
    stop(
      "Species not found in all_haplotype_results: ",
      species_name
    )
  }
  
  
  species_result <- all_haplotype_results[[species_name]]
  
  
  if (is.null(species_result)) {
    
    stop(
      "NULL result for species: ",
      species_name
    )
  }
  
  
  if (!marker_name %in% names(species_result)) {
    
    stop(
      "Marker ",
      marker_name,
      " not available in all_haplotype_results for ",
      species_name
    )
  }
  
  
  result <- species_result[[marker_name]]
  
  
  if (is.null(result)) {
    
    stop(
      "NULL result for ",
      species_name,
      " - ",
      marker_name
    )
  }
  
  
  return(result)
}


# ============================================================
# 18.8. CHECK THAT ALL SELECTED NETWORKS EXIST
# ============================================================

network_check <- vector(
  mode = "list",
  length = nrow(selected_main_networks)
)


for (i in seq_len(nrow(selected_main_networks))) {
  
  sp <- selected_main_networks$species[i]
  mk <- selected_main_networks$marker[i]
  
  
  network_check[[i]] <- tryCatch(
    
    {
      
      result <- get_main_haplotype_result(
        species_name = sp,
        marker_name = mk
      )
      
      tibble::tibble(
        species = sp,
        marker = mk,
        network_available = TRUE,
        error = NA_character_
      )
      
    },
    
    error = function(e) {
      
      tibble::tibble(
        species = sp,
        marker = mk,
        network_available = FALSE,
        error = e$message
      )
      
    }
    
  )
  
}


network_check <- dplyr::bind_rows(
  network_check
)


print(network_check)


if (any(!network_check$network_available)) {
  
  stop(
    "One or more selected networks are unavailable. ",
    "See network_check above."
  )
}


# ============================================================
# 18.9. FUNCTION TO REBUILD THE MST
# ============================================================
#
# We rebuild the MST from the pairwise distances already
# calculated by make_haplotype_network().
#
# This allows us to:
#
#   1. determine the genetic distances actually represented
#      by network edges;
#
#   2. calculate ONE common genetic-distance scale;
#
#   3. use exactly the same edge-width mapping in every panel.
#
# ============================================================

build_main_mst <- function(result) {
  
  # Single-haplotype network
  if (
    is.null(result$distances) ||
    nrow(result$ASVs) <= 1
  ) {
    
    return(NULL)
    
  }
  
  
  edge_df <- result$distances %>%
    
    dplyr::mutate(
      from = as.character(from),
      to = as.character(to),
      distance = as.numeric(distance)
    ) %>%
    
    dplyr::filter(
      !is.na(distance),
      is.finite(distance)
    )
  
  
  if (nrow(edge_df) == 0) {
    
    return(NULL)
    
  }
  
  
  vertex_df <- result$ASVs %>%
    
    dplyr::transmute(
      name = as.character(ASV)
    ) %>%
    
    dplyr::distinct()
  
  
  graph_full <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = vertex_df
  )
  
  
  graph_mst <- igraph::mst(
    graph_full,
    weights = igraph::E(graph_full)$distance
  )
  
  
  return(graph_mst)
}


# ============================================================
# 18.10. DETERMINE GLOBAL GENETIC-DISTANCE SCALE
# ============================================================

main_mst_objects <- list()

all_main_edge_distances <- numeric()


for (i in seq_len(nrow(selected_main_networks))) {
  
  sp <- selected_main_networks$species[i]
  mk <- selected_main_networks$marker[i]
  
  
  result <- get_main_haplotype_result(
    species_name = sp,
    marker_name = mk
  )
  
  
  mst_i <- build_main_mst(
    result
  )
  
  
  key_i <- paste(
    sp,
    mk,
    sep = "___"
  )
  
  
  main_mst_objects[[key_i]] <- mst_i
  
  
  if (
    !is.null(mst_i) &&
    igraph::ecount(mst_i) > 0
  ) {
    
    all_main_edge_distances <- c(
      
      all_main_edge_distances,
      
      as.numeric(
        igraph::E(mst_i)$distance
      )
      
    )
    
  }
  
}


all_main_edge_distances <- all_main_edge_distances[
  is.finite(all_main_edge_distances)
]


if (length(all_main_edge_distances) == 0) {
  
  stop(
    "No finite genetic distances found in selected networks."
  )
}


global_distance_max <- max(
  all_main_edge_distances,
  na.rm = TRUE
)


# Always start the common scale at zero.
global_distance_limits <- c(
  0,
  global_distance_max
)


# Pretty common breaks
global_distance_breaks <- scales::breaks_pretty(
  n = 4
)(
  global_distance_limits
)


# Keep only valid breaks inside range
global_distance_breaks <- global_distance_breaks[
  global_distance_breaks >= global_distance_limits[1] &
    global_distance_breaks <= global_distance_limits[2]
]


cat(
  "\n============================================================\n"
)

cat(
  "GLOBAL GENETIC DISTANCE SCALE\n"
)

cat(
  "Minimum = ",
  global_distance_limits[1],
  "\n",
  sep = ""
)

cat(
  "Maximum = ",
  global_distance_limits[2],
  "\n",
  sep = ""
)

cat(
  "Breaks  = ",
  paste(
    global_distance_breaks,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)


# ============================================================
# 18.11. PREPARE NODE PIE DATA
# ============================================================

prepare_main_node_data <- function(result) {
  
  site_abundance <- result$site_abundance %>%
    
    dplyr::filter(
      Site %in% names(site_colors)
    )
  
  
  site_pie <- site_abundance %>%
    
    tidyr::pivot_wider(
      names_from = Site,
      values_from = abundance,
      values_fill = 0
    )
  
  
  missing_sites <- setdiff(
    names(site_colors),
    names(site_pie)
  )
  
  
  for (ss in missing_sites) {
    
    site_pie[[ss]] <- 0
    
  }
  
  
  site_pie <- site_pie %>%
    
    dplyr::select(
      ASV,
      dplyr::all_of(
        names(site_colors)
      )
    ) %>%
    
    dplyr::mutate(
      
      total_abundance = rowSums(
        dplyr::across(
          dplyr::all_of(
            names(site_colors)
          )
        ),
        na.rm = TRUE
      )
      
    )
  
  
  # ----------------------------------------------------------
  # Node sizes are rescaled WITHIN species.
  #
  # This is deliberate:
  #
  # node size represents relative dominance within each
  # network and is NOT interpreted quantitatively between
  # species.
  # ----------------------------------------------------------
  
  if (
    nrow(site_pie) == 1 ||
    length(
      unique(
        site_pie$total_abundance
      )
    ) == 1
  ) {
    
    site_pie$radius <- 0.11
    
  } else {
    
    site_pie$radius <- scales::rescale(
      
      sqrt(
        site_pie$total_abundance
      ),
      
      to = c(
        0.055,
        0.14
      )
      
    )
    
  }
  
  
  return(site_pie)
}


# ============================================================
# 18.12. FINAL NETWORK-PLOT FUNCTION
# ============================================================
#
# labels_max_nodes:
#
# ASV labels become unreadable for large networks such as
# Peringia ulvae.
#
# Individual labelled networks remain available in the
# Supplementary Material.
#
# Here labels are shown only for small networks.
#
# ============================================================

make_main_network_panel <- function(
    species_name,
    marker_name,
    status_value,
    metric_label,
    labels_max_nodes = 12
) {
  
  result <- get_main_haplotype_result(
    species_name = species_name,
    marker_name = marker_name
  )
  
  
  node_data <- prepare_main_node_data(
    result
  )
  
  
  n_nodes <- nrow(
    result$ASVs
  )
  
  
  # ==========================================================
  # SINGLE-HAPLOTYPE CASE
  # ==========================================================
  
  if (
    n_nodes == 1 ||
    is.null(result$distances)
  ) {
    
    single_asv <- as.character(
      result$ASVs$ASV[1]
    )
    
    
    single_df <- node_data %>%
      
      dplyr::filter(
        ASV == single_asv
      ) %>%
      
      dplyr::mutate(
        x = 0,
        y = 0,
        radius = 0.16
      )
    
    
    p_single <- ggplot2::ggplot() +
      
      scatterpie::geom_scatterpie(
        
        data = single_df,
        
        ggplot2::aes(
          x = x,
          y = y,
          r = radius
        ),
        
        cols = names(site_colors),
        
        colour = "black",
        linewidth = 0.45
        
      ) +
      
      ggplot2::scale_fill_manual(
        
        values = site_colors,
        
        breaks = names(site_colors),
        
        limits = names(site_colors),
        
        drop = FALSE,
        
        name = "Site"
        
      ) +
      
      ggplot2::coord_equal(
        xlim = c(-0.50, 0.50),
        ylim = c(-0.45, 0.45),
        clip = "off"
      ) +
      
      ggplot2::labs(
        
        title = species_name,
        
        subtitle = metric_label,
        
        tag = status_value
        
      ) +
      
      ggplot2::theme_void() +
      
      ggplot2::theme(
        
        plot.title = ggplot2::element_text(
          face = "italic",
          size = 13,
          hjust = 0.5,
          margin = ggplot2::margin(
            b = 3
          )
        ),
        
        plot.subtitle = ggplot2::element_text(
          size = 8.5,
          hjust = 0.5,
          margin = ggplot2::margin(
            b = 4
          )
        ),
        
        plot.tag = ggplot2::element_text(
          face = "bold",
          size = 9
        ),
        
        plot.tag.position = c(
          0.02,
          0.98
        ),
        
        legend.position = "bottom",
        
        legend.title = ggplot2::element_text(
          face = "bold",
          size = 9
        ),
        
        legend.text = ggplot2::element_text(
          size = 8
        )
        
      )
    
    
    return(
      p_single
    )
    
  }
  
  
  # ==========================================================
  # MULTI-HAPLOTYPE NETWORK
  # ==========================================================
  
  key_i <- paste(
    species_name,
    marker_name,
    sep = "___"
  )
  
  
  graph_mst <- main_mst_objects[[key_i]]
  
  
  if (is.null(graph_mst)) {
    
    stop(
      "MST could not be reconstructed for ",
      species_name,
      " - ",
      marker_name
    )
    
  }
  
  
  node_attr <- node_data %>%
    
    dplyr::rename(
      name = ASV
    )
  
  
  graph_tbl <- tidygraph::as_tbl_graph(
    graph_mst
  ) %>%
    
    tidygraph::activate(nodes) %>%
    
    dplyr::left_join(
      node_attr,
      by = "name"
    )
  
  
  p <- ggraph::ggraph(
    graph_tbl,
    layout = "stress"
  ) +
    
    # --------------------------------------------------------
  # NETWORK EDGES
  #
  # IMPORTANT:
  # ALL panels use exactly the same limits.
  # --------------------------------------------------------
  
  ggraph::geom_edge_link(
    
    ggplot2::aes(
      width = distance
    ),
    
    colour = "grey35",
    alpha = 0.80,
    lineend = "round"
    
  ) +
    
    
    # --------------------------------------------------------
  # SITE COMPOSITION WITHIN EACH HAPLOTYPE
  # --------------------------------------------------------
  
  scatterpie::geom_scatterpie(
    
    ggplot2::aes(
      x = x,
      y = y,
      r = radius
    ),
    
    cols = names(site_colors),
    
    colour = "black",
    linewidth = 0.35
    
  ) +
    
    
    # --------------------------------------------------------
  # COMMON SITE SCALE
  # --------------------------------------------------------
  
  ggplot2::scale_fill_manual(
    
    values = site_colors,
    
    breaks = names(site_colors),
    
    limits = names(site_colors),
    
    drop = FALSE,
    
    name = "Site"
    
  ) +
    
    
    # --------------------------------------------------------
  # COMMON GENETIC-DISTANCE SCALE
  # --------------------------------------------------------
  
  ggraph::scale_edge_width_continuous(
    
    limits = global_distance_limits,
    
    breaks = global_distance_breaks,
    
    range = c(
      0.35,
      2.6
    ),
    
    name = "Sequence distance"
    
  ) +
    
    
    # --------------------------------------------------------
  # LEGEND ORDER
  # --------------------------------------------------------
  
  ggplot2::guides(
    
    fill = ggplot2::guide_legend(
      order = 1,
      title.position = "top"
    ),
    
    edge_width = ggplot2::guide_legend(
      order = 2,
      title.position = "top"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # ASV LABELS ONLY FOR SMALL NETWORKS
  # ----------------------------------------------------------
  
  if (n_nodes <= labels_max_nodes) {
    
    p <- p +
      
      ggraph::geom_node_text(
        
        ggplot2::aes(
          label = name
        ),
        
        repel = TRUE,
        size = 2.6,
        fontface = "bold"
        
      )
    
  }
  
  
  # ----------------------------------------------------------
  # TITLES / METRICS / STATUS
  # ----------------------------------------------------------
  
  p <- p +
    
    ggplot2::coord_equal(
      clip = "off"
    ) +
    
    ggplot2::labs(
      
      title = species_name,
      
      subtitle = metric_label,
      
      tag = status_value
      
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      
      plot.title = ggplot2::element_text(
        face = "italic",
        size = 13,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 3
        )
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 8.5,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 4
        )
      ),
      
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 9
      ),
      
      plot.tag.position = c(
        0.02,
        0.98
      ),
      
      plot.margin = ggplot2::margin(
        5,
        7,
        5,
        7
      ),
      
      legend.position = "bottom",
      
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 9
      ),
      
      legend.text = ggplot2::element_text(
        size = 8
      ),
      
      legend.key.width = grid::unit(
        1.1,
        "cm"
      )
      
    )
  
  
  return(p)
}


# ============================================================
# 18.13. BUILD ALL MAIN PANELS
# ============================================================

main_network_plots <- vector(
  mode = "list",
  length = nrow(selected_metrics_main)
)


names(main_network_plots) <- paste0(
  "panel_",
  selected_metrics_main$panel_order
)


for (i in seq_len(nrow(selected_metrics_main))) {
  
  sp <- selected_metrics_main$species[i]
  mk <- selected_metrics_main$marker[i]
  st <- selected_metrics_main$status[i]
  ml <- selected_metrics_main$metric_label[i]
  
  
  message(
    "\nBuilding MAIN panel: ",
    sp,
    " | ",
    mk,
    " | ",
    st
  )
  
  
  main_network_plots[[i]] <- make_main_network_panel(
    
    species_name = sp,
    
    marker_name = mk,
    
    status_value = st,
    
    metric_label = ml,
    
    labels_max_nodes = 12
    
  )
  
}


# ============================================================
# 18.14. FINAL COMPOSITE FIGURE
# ============================================================
#
# 4 columns works well for 13 networks:
#
# Row 1: panels 1–4
# Row 2: panels 5–8
# Row 3: panels 9–12
# Row 4: panel 13
#
# Patchwork collects identical legends.
#
# Because ALL multi-haplotype panels have the same
# scale_edge_width limits and breaks, only ONE genetic-
# distance legend is shown.
#
# ============================================================

main_haplotype_figure <- patchwork::wrap_plots(
  
  main_network_plots,
  
  ncol = 4,
  
  guides = "collect"
  
) +
  
  patchwork::plot_annotation(
    
    title = "Haplotype structure across biogeographic and invasion-status categories",
    
    subtitle = paste0(
      "Node sectors represent site composition; ",
      "edge width represents pairwise sequence distance"
    ),
    
    theme = ggplot2::theme(
      
      plot.title = ggplot2::element_text(
        size = 17,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 10.5,
        hjust = 0.5
      )
      
    )
    
  ) &
  
  ggplot2::theme(
    
    legend.position = "bottom",
    
    legend.box = "horizontal",
    
    legend.box.just = "center"
    
  )


# ============================================================
# 18.15. SHOW FIGURE
# ============================================================

print(
  main_haplotype_figure
)


# ============================================================
# 18.16. EXPORT FINAL MAIN FIGURE
# ============================================================

main_svg <- file.path(
  main_haplo_dir,
  "MAIN_SELECTED_HAPLOTYPE_NETWORKS.svg"
)

main_pdf <- file.path(
  main_haplo_dir,
  "MAIN_SELECTED_HAPLOTYPE_NETWORKS.pdf"
)

main_png <- file.path(
  main_haplo_dir,
  "MAIN_SELECTED_HAPLOTYPE_NETWORKS.png"
)


# ------------------------------------------------------------
# SVG
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename = main_svg,
  
  plot = main_haplotype_figure,
  
  width = 16,
  height = 18,
  
  units = "in",
  
  bg = "white"
  
)


# ------------------------------------------------------------
# PDF
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename = main_pdf,
  
  plot = main_haplotype_figure,
  
  width = 16,
  height = 18,
  
  units = "in",
  
  device = cairo_pdf,
  
  bg = "white"
  
)


# ------------------------------------------------------------
# PNG
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename = main_png,
  
  plot = main_haplotype_figure,
  
  width = 16,
  height = 18,
  
  units = "in",
  
  dpi = 400,
  
  bg = "white"
  
)


# ============================================================
# 18.17. EXPORT THE EXACT METRICS USED IN THE MAIN FIGURE
# ============================================================
#
# This is useful for:
#   - figure caption
#   - Results
#   - Supplementary table
#   - reproducibility
#
# ============================================================

main_metrics_export <- selected_metrics_main %>%
  
  dplyr::select(
    
    panel_order,
    status,
    species,
    marker,
    
    n_haplotypes,
    
    dominant_haplotype_frequency,
    
    shared_haplotypes,
    prop_shared_haplotypes,
    
    site_specific_haplotypes,
    prop_site_specific_haplotypes,
    
    sites_per_haplotype,
    
    n_positive_samples_species,
    n_sites_species,
    
    metric_label
    
  )


readr::write_csv(
  
  main_metrics_export,
  
  file.path(
    main_haplo_dir,
    "MAIN_SELECTED_HAPLOTYPE_NETWORKS_metrics.csv"
  )
  
)


# ============================================================
# 18.18. EXPORT GLOBAL GENETIC-DISTANCE INFORMATION
# ============================================================

distance_scale_export <- tibble::tibble(
  
  parameter = c(
    "global_min_distance",
    "global_max_distance"
  ),
  
  value = c(
    global_distance_limits[1],
    global_distance_limits[2]
  )
  
)


readr::write_csv(
  
  distance_scale_export,
  
  file.path(
    main_haplo_dir,
    "MAIN_genetic_distance_scale.csv"
  )
  
)


# ============================================================
# 18.19. FINAL MESSAGE
# ============================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "FINAL MAIN HAPLOTYPE FIGURE FINISHED\n"
)

cat(
  "\nFiles saved in:\n",
  normalizePath(main_haplo_dir),
  "\n",
  sep = ""
)

cat(
  "\nSVG:\n",
  main_svg,
  "\n",
  sep = ""
)

cat(
  "\nPDF:\n",
  main_pdf,
  "\n",
  sep = ""
)

cat(
  "\nPNG:\n",
  main_png,
  "\n",
  sep = ""
)

cat(
  "\nCommon sequence-distance range: ",
  global_distance_limits[1],
  " – ",
  global_distance_limits[2],
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)




# ============================================================
# 19. PUBLICATION-QUALITY MAIN HAPLOTYPE FIGURE
# ============================================================
#
# Requirements already existing in the workspace:
#
#   all_haplotype_results
#   selected_metrics_main
#   site_colors
#
# This version:
#
#   - DOES NOT use the previous composite figure
#   - rebuilds only the selected networks
#   - normalises network coordinates within every panel
#   - recalculates genetic distance as % sequence divergence
#   - uses ONE common genetic-distance scale
#   - uses ONE manual Site legend
#   - uses ONE manual Sequence divergence legend
#   - removes ASV labels from Main Figure
#   - keeps quantitative metrics under each species
#   - exports SVG, PDF and PNG
#
# ============================================================


# ============================================================
# 19.1. PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(Biostrings)
library(ape)
library(igraph)
library(tidygraph)
library(ggraph)
library(scatterpie)
library(scales)
library(patchwork)
library(tibble)
library(readr)
library(grid)


# ============================================================
# 19.2. OUTPUT DIRECTORY
# ============================================================

paper_haplo_dir <- file.path(
  haplo_dir,
  "MAIN_FINAL_HAPLOTYPE_FIGURE_PUBLICATION"
)

dir.create(
  paper_haplo_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 19.3. CHECK REQUIRED OBJECTS
# ============================================================

required_objects <- c(
  "all_haplotype_results",
  "selected_metrics_main",
  "site_colors"
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
    "These required objects are missing:\n",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}


# ============================================================
# 19.4. RETRIEVE ONE STORED NETWORK
# ============================================================

# ============================================================
# 19.4. RETRIEVE ONE STORED NETWORK
# ============================================================

get_publication_result <- function(
    species_name,
    marker_name
) {
  
  marker_name <- toupper(marker_name)
  
  
  # ----------------------------------------------------------
  # Check species
  # ----------------------------------------------------------
  
  if (!species_name %in% names(all_haplotype_results)) {
    
    stop(
      "Species not found in all_haplotype_results: ",
      species_name
    )
  }
  
  
  # ----------------------------------------------------------
  # Check marker
  # ----------------------------------------------------------
  
  if (
    !marker_name %in%
    names(all_haplotype_results[[species_name]])
  ) {
    
    stop(
      "Marker ",
      marker_name,
      " not found for ",
      species_name
    )
  }
  
  
  # ----------------------------------------------------------
  # Retrieve stored result
  # ----------------------------------------------------------
  
  result <- all_haplotype_results[[species_name]][[marker_name]]
  
  
  # ----------------------------------------------------------
  # Check result
  # ----------------------------------------------------------
  
  if (is.null(result)) {
    
    stop(
      "NULL result for ",
      species_name,
      " - ",
      marker_name
    )
  }
  
  
  return(result)
}


# ============================================================
# 19.5. RECOMPUTE GENETIC DISTANCE
# ============================================================
#
# IMPORTANT SCIENTIFIC CHANGE
#
# Previous networks used:
#
#   ape::dist.dna(model = "N")
#
# i.e. absolute numbers of differences.
#
# For one legend shared by COI and 18S, this is not ideal
# because marker lengths can differ.
#
# Here we calculate:
#
#   model = "raw"
#
# which gives the proportion of differing sites.
#
# We multiply by 100 so that the figure shows:
#
#   Sequence divergence (%)
#
# ============================================================

calculate_publication_distances <- function(result) {
  
  n_asvs <- nrow(
    result$ASVs
  )
  
  
  # ----------------------------------------------------------
  # Single haplotype
  # ----------------------------------------------------------
  
  if (n_asvs <= 1) {
    
    return(NULL)
  }
  
  
  # ----------------------------------------------------------
  # Prefer the alignment already stored in the network object
  # ----------------------------------------------------------
  
  aligned_df <- result$aligned_sequences
  
  
  if (
    is.null(aligned_df) ||
    nrow(aligned_df) <= 1
  ) {
    
    stop(
      "Aligned sequences unavailable for ",
      result$species,
      " - ",
      result$marker
    )
  }
  
  
  aligned_df <- aligned_df %>%
    
    dplyr::filter(
      !is.na(aligned_sequence),
      aligned_sequence != ""
    )
  
  
  # ----------------------------------------------------------
  # Convert alignment to DNAbin
  # ----------------------------------------------------------
  
  dna_list <- strsplit(
    aligned_df$aligned_sequence,
    split = ""
  )
  
  
  names(dna_list) <- aligned_df$ASV
  
  
  dna <- ape::as.DNAbin(
    dna_list
  )
  
  
  # ----------------------------------------------------------
  # Pairwise proportional sequence divergence
  # ----------------------------------------------------------
  
  dist_matrix <- ape::dist.dna(
    
    dna,
    
    model = "raw",
    
    pairwise.deletion = TRUE,
    
    as.matrix = TRUE
    
  )
  
  
  edge_df <- as.data.frame(
    as.table(dist_matrix),
    stringsAsFactors = FALSE
  )
  
  
  names(edge_df) <- c(
    "from",
    "to",
    "distance_prop"
  )
  
  
  edge_df <- edge_df %>%
    
    dplyr::mutate(
      
      from = as.character(from),
      
      to = as.character(to),
      
      distance_prop =
        as.numeric(distance_prop),
      
      distance_pct =
        distance_prop * 100
      
    ) %>%
    
    dplyr::filter(
      
      from < to,
      
      !is.na(distance_pct),
      
      is.finite(distance_pct)
      
    )
  
  
  return(edge_df)
}


# ============================================================
# 19.6. BUILD PUBLICATION MST
# ============================================================

build_publication_mst <- function(result) {
  
  edge_df <- calculate_publication_distances(
    result
  )
  
  
  if (is.null(edge_df)) {
    
    return(NULL)
  }
  
  
  vertices <- result$ASVs %>%
    
    dplyr::transmute(
      
      name = as.character(
        ASV
      )
      
    ) %>%
    
    dplyr::distinct()
  
  
  graph_full <- igraph::graph_from_data_frame(
    
    edge_df %>%
      dplyr::select(
        from,
        to,
        distance_pct
      ),
    
    directed = FALSE,
    
    vertices = vertices
    
  )
  
  
  graph_mst <- igraph::mst(
    
    graph_full,
    
    weights =
      igraph::E(graph_full)$distance_pct
    
  )
  
  
  return(graph_mst)
}


# ============================================================
# 19.7. CALCULATE ALL SELECTED NETWORKS FIRST
# ============================================================

publication_networks <- list()

publication_all_distances <- numeric()


for (
  i in seq_len(
    nrow(selected_metrics_main)
  )
) {
  
  sp <- selected_metrics_main$species[i]
  
  mk <- selected_metrics_main$marker[i]
  
  
  message(
    "Preparing publication network: ",
    sp,
    " | ",
    mk
  )
  
  
  result_i <- get_publication_result(
    
    species_name = sp,
    
    marker_name = mk
    
  )
  
  
  mst_i <- build_publication_mst(
    result_i
  )
  
  
  key_i <- paste(
    sp,
    mk,
    sep = "___"
  )
  
  
  publication_networks[[key_i]] <- list(
    
    result = result_i,
    
    mst = mst_i
    
  )
  
  
  if (
    !is.null(mst_i) &&
    igraph::ecount(mst_i) > 0
  ) {
    
    publication_all_distances <- c(
      
      publication_all_distances,
      
      igraph::E(
        mst_i
      )$distance_pct
      
    )
    
  }
  
}


publication_all_distances <-
  publication_all_distances[
    is.finite(
      publication_all_distances
    )
  ]


if (
  length(
    publication_all_distances
  ) == 0
) {
  
  stop(
    "No genetic distances found."
  )
}


# ============================================================
# 19.8. GLOBAL GENETIC DISTANCE SCALE
# ============================================================

global_divergence_min <- 0

global_divergence_max <- max(
  publication_all_distances,
  na.rm = TRUE
)


# ------------------------------------------------------------
# Nice legend breaks
# ------------------------------------------------------------

divergence_breaks <- pretty(
  
  c(
    0,
    global_divergence_max
  ),
  
  n = 4
  
)


divergence_breaks <- divergence_breaks[
  
  divergence_breaks >= 0 &
    
    divergence_breaks <=
    global_divergence_max
  
]


# Make sure maximum is represented
if (
  length(divergence_breaks) < 2
) {
  
  divergence_breaks <- c(
    0,
    global_divergence_max
  )
  
}


cat(
  "\n========================================\n"
)

cat(
  "COMMON SEQUENCE DIVERGENCE SCALE\n"
)

cat(
  "Maximum divergence = ",
  round(
    global_divergence_max,
    3
  ),
  "%\n",
  sep = ""
)

cat(
  "Legend breaks = ",
  paste(
    divergence_breaks,
    collapse = ", "
  ),
  "%\n",
  sep = ""
)

cat(
  "========================================\n"
)


# ============================================================
# 19.9. COMMON EDGE-WIDTH FUNCTION
# ============================================================
#
# SAME transformation used:
#
#   1. in every network
#   2. in the custom distance legend
#
# ============================================================

distance_to_linewidth <- function(x) {
  
  if (
    global_divergence_max <= 0
  ) {
    
    return(
      rep(
        0.55,
        length(x)
      )
    )
    
  }
  
  
  scales::rescale(
    
    x,
    
    from = c(
      0,
      global_divergence_max
    ),
    
    to = c(
      0.35,
      2.20
    )
    
  )
  
}


# ============================================================
# 19.10. PREPARE PIE-CHART NODE DATA
# ============================================================

prepare_publication_nodes <- function(
    result,
    coordinates
) {
  
  site_abundance <- result$site_abundance %>%
    
    dplyr::filter(
      Site %in%
        names(site_colors)
    )
  
  
  site_pie <- site_abundance %>%
    
    tidyr::pivot_wider(
      
      names_from = Site,
      
      values_from = abundance,
      
      values_fill = 0
      
    )
  
  
  missing_sites <- setdiff(
    
    names(site_colors),
    
    names(site_pie)
    
  )
  
  
  for (ss in missing_sites) {
    
    site_pie[[ss]] <- 0
    
  }
  
  
  site_pie <- site_pie %>%
    
    dplyr::select(
      
      ASV,
      
      dplyr::all_of(
        names(site_colors)
      )
      
    ) %>%
    
    dplyr::mutate(
      
      total_abundance = rowSums(
        
        dplyr::across(
          
          dplyr::all_of(
            names(site_colors)
          )
          
        ),
        
        na.rm = TRUE
        
      )
      
    ) %>%
    
    dplyr::left_join(
      
      coordinates,
      
      by = c(
        "ASV" = "name"
      )
      
    )
  
  
  # ==========================================================
  # Node radius
  #
  # Adaptive maximum radius depending on network complexity.
  #
  # This is important because an 84-node network cannot use
  # the same node radius as a 3-node network.
  # ==========================================================
  
  n_nodes <- nrow(
    site_pie
  )
  
  
  radius_max <- dplyr::case_when(
    
    n_nodes <= 3  ~ 0.095,
    
    n_nodes <= 6  ~ 0.080,
    
    n_nodes <= 12 ~ 0.065,
    
    n_nodes <= 25 ~ 0.050,
    
    n_nodes <= 50 ~ 0.037,
    
    TRUE          ~ 0.028
    
  )
  
  
  radius_min <- radius_max * 0.50
  
  
  if (
    n_nodes == 1 ||
    length(
      unique(
        site_pie$total_abundance
      )
    ) == 1
  ) {
    
    site_pie$radius <-
      radius_max
    
  } else {
    
    site_pie$radius <-
      
      scales::rescale(
        
        sqrt(
          site_pie$total_abundance
        ),
        
        to = c(
          radius_min,
          radius_max
        )
        
      )
    
  }
  
  
  return(
    site_pie
  )
}


# ============================================================
# 19.11. NORMALISE NETWORK COORDINATES
# ============================================================
#
# THIS IS ONE OF THE MAIN VISUAL FIXES.
#
# Stress-layout coordinates vary enormously among networks.
#
# Here EVERY network is rescaled into:
#
#   x = 0.08 ... 0.92
#   y = 0.08 ... 0.92
#
# Therefore Peringia no longer becomes microscopic simply
# because its stress layout covers a larger numeric range.
#
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


get_publication_coordinates <- function(
    graph_mst
) {
  
  graph_tbl <- tidygraph::as_tbl_graph(
    graph_mst
  )
  
  
  layout_df <- ggraph::create_layout(
    
    graph_tbl,
    
    layout = "stress"
    
  ) %>%
    
    as.data.frame()
  
  
  coordinates <- tibble::tibble(
    
    name = as.character(
      layout_df$name
    ),
    
    x = normalise_coordinate(
      layout_df$x
    ),
    
    y = normalise_coordinate(
      layout_df$y
    )
    
  )
  
  
  return(
    coordinates
  )
}


# ============================================================
# 19.12. PREPARE EDGE COORDINATES
# ============================================================

prepare_publication_edges <- function(
    graph_mst,
    coordinates
) {
  
  edge_df <- igraph::as_data_frame(
    
    graph_mst,
    
    what = "edges"
    
  ) %>%
    
    dplyr::rename(
      
      distance_pct =
        distance_pct
      
    )
  
  
  edge_df <- edge_df %>%
    
    dplyr::left_join(
      
      coordinates %>%
        
        dplyr::rename(
          
          from = name,
          
          x = x,
          
          y = y
          
        ),
      
      by = "from"
      
    ) %>%
    
    dplyr::left_join(
      
      coordinates %>%
        
        dplyr::rename(
          
          to = name,
          
          xend = x,
          
          yend = y
          
        ),
      
      by = "to"
      
    ) %>%
    
    dplyr::mutate(
      
      edge_linewidth =
        distance_to_linewidth(
          distance_pct
        )
      
    )
  
  
  return(
    edge_df
  )
}


# ============================================================
# 19.13. CREATE ONE PUBLICATION PANEL
# ============================================================

make_publication_panel <- function(
    species_name,
    marker_name,
    status_value,
    metric_label
) {
  
  key_i <- paste(
    species_name,
    marker_name,
    sep = "___"
  )
  
  stored <- publication_networks[[key_i]]
  
  result <- stored$result
  graph_mst <- stored$mst
  
  n_nodes <- nrow(
    result$ASVs
  )
  
  
  # ==========================================================
  # SINGLE HAPLOTYPE
  # ==========================================================
  
  if (
    n_nodes <= 1 ||
    is.null(graph_mst)
  ) {
    
    coords <- tibble::tibble(
      
      name =
        as.character(
          result$ASVs$ASV[1]
        ),
      
      x = 0.5,
      
      y = 0.5
      
    )
    
    
    node_df <-
      prepare_publication_nodes(
        
        result = result,
        
        coordinates = coords
        
      )
    
    
    node_df$radius <- 0.10
    
    
    p <- ggplot2::ggplot() +
      
      scatterpie::geom_scatterpie(
        
        data = node_df,
        
        ggplot2::aes(
          x = x,
          y = y,
          r = radius
        ),
        
        cols =
          names(site_colors),
        
        colour = "black",
        
        linewidth = 0.45
        
      ) +
      
      ggplot2::scale_fill_manual(
        
        values = site_colors,
        
        limits =
          names(site_colors),
        
        drop = FALSE
        
      )
    
  } else {
    
    # ========================================================
    # MULTI-HAPLOTYPE NETWORK
    # ========================================================
    
    coords <-
      get_publication_coordinates(
        graph_mst
      )
    
    
    node_df <-
      prepare_publication_nodes(
        
        result = result,
        
        coordinates = coords
        
      )
    
    
    edge_df <-
      prepare_publication_edges(
        
        graph_mst = graph_mst,
        
        coordinates = coords
        
      )
    
    
    p <- ggplot2::ggplot() +
      
      # ------------------------------------------------------
    # Edges
    # ------------------------------------------------------
    
    ggplot2::geom_segment(
      
      data = edge_df,
      
      ggplot2::aes(
        
        x = x,
        y = y,
        
        xend = xend,
        yend = yend,
        
        linewidth =
          distance_pct
        
      ),
      
      colour = "grey35",
      
      alpha = 0.78,
      
      lineend = "round"
      
    ) +
      
      ggplot2::scale_linewidth_continuous(
        
        limits = c(
          0,
          global_divergence_max
        ),
        
        range = c(
          0.35,
          2.20
        ),
        
        guide = "none"
        
      ) +
      
      # ------------------------------------------------------
    # Haplotype nodes
    # ------------------------------------------------------
    
    scatterpie::geom_scatterpie(
      
      data = node_df,
      
      ggplot2::aes(
        
        x = x,
        
        y = y,
        
        r = radius
        
      ),
      
      cols =
        names(site_colors),
      
      colour = "black",
      
      linewidth = 0.40
      
    ) +
      
      ggplot2::scale_fill_manual(
        
        values = site_colors,
        
        limits =
          names(site_colors),
        
        drop = FALSE
        
      )
    
  }
  
  
  # ==========================================================
  # COMMON PANEL DESIGN
  # ==========================================================
  
  p <- p +
    
    ggplot2::coord_equal(
      
      xlim = c(
        -0.02,
        1.02
      ),
      
      ylim = c(
        -0.02,
        1.02
      ),
      
      expand = FALSE,
      
      clip = "off"
      
    ) +
    
    ggplot2::labs(
      
      title = species_name,
      
      subtitle = metric_label
      
    ) +
    
    # --------------------------------------------------------
  # Status in upper-left corner
  # --------------------------------------------------------
  
  ggplot2::annotate(
    
    "text",
    
    x = 0.01,
    
    y = 1.01,
    
    label = status_value,
    
    hjust = 0,
    
    vjust = 1,
    
    size = 3.0,
    
    fontface = "bold",
    
    colour = "grey25"
    
  ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      
      legend.position = "none",
      
      plot.title = ggplot2::element_text(
        
        face = "italic",
        
        size = 11.5,
        
        hjust = 0.5,
        
        margin = ggplot2::margin(
          b = 2
        )
        
      ),
      
      plot.subtitle =
        ggplot2::element_text(
          
          size = 7.5,
          
          hjust = 0.5,
          
          colour = "grey25",
          
          margin =
            ggplot2::margin(
              b = 5
            )
          
        ),
      
      plot.margin =
        ggplot2::margin(
          
          t = 5,
          
          r = 5,
          
          b = 5,
          
          l = 5
          
        )
      
    )
  
  
  return(p)
}


# ============================================================
# 19.14. BUILD ALL NETWORK PANELS
# ============================================================

publication_panels <- vector(
  
  mode = "list",
  
  length =
    nrow(
      selected_metrics_main
    )
  
)


for (
  i in seq_len(
    nrow(
      selected_metrics_main
    )
  )
) {
  
  sp <-
    selected_metrics_main$species[i]
  
  mk <-
    selected_metrics_main$marker[i]
  
  st <-
    selected_metrics_main$status[i]
  
  ml <-
    selected_metrics_main$metric_label[i]
  
  
  publication_panels[[i]] <-
    
    make_publication_panel(
      
      species_name = sp,
      
      marker_name = mk,
      
      status_value = st,
      
      metric_label = ml
      
    )
  
}


# ============================================================
# 19.15. CREATE ONE MANUAL SITE LEGEND
# ============================================================
#
# We deliberately DO NOT let patchwork collect ggplot legends.
#
# The legend is a real independent plot.
#
# Therefore it can appear only once.
#
# ============================================================

site_legend_df <- tibble::tibble(
  
  site = factor(
    
    names(site_colors),
    
    levels =
      names(site_colors)
    
  ),
  
  x =
    seq_along(
      site_colors
    ),
  
  y = 1
  
)


site_legend <- ggplot2::ggplot(
  
  site_legend_df,
  
  ggplot2::aes(
    x = x,
    y = y
  )
  
) +
  
  ggplot2::geom_point(
    
    ggplot2::aes(
      fill = site
    ),
    
    shape = 22,
    
    size = 5.2,
    
    colour = "black",
    
    stroke = 0.45
    
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = site
    ),
    
    y = 0.63,
    
    size = 3.0
    
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x = 0.35,
    
    y = 1.35,
    
    label = "Site",
    
    fontface = "bold",
    
    hjust = 0,
    
    size = 3.3
    
  ) +
  
  ggplot2::scale_fill_manual(
    
    values = site_colors,
    
    guide = "none"
    
  ) +
  
  ggplot2::scale_x_continuous(
    
    limits = c(
      0.2,
      length(site_colors) + 0.5
    )
    
  ) +
  
  ggplot2::coord_cartesian(
    
    ylim = c(
      0.40,
      1.45
    ),
    
    clip = "off"
    
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    
    plot.margin =
      ggplot2::margin(
        2,
        3,
        2,
        3
      )
    
  )


# ============================================================
# 19.16. CREATE ONE MANUAL GENETIC-DISTANCE LEGEND
# ============================================================

distance_legend_df <- tibble::tibble(
  
  distance_pct =
    divergence_breaks,
  
  row =
    seq_along(
      divergence_breaks
    )
  
) %>%
  
  dplyr::mutate(
    
    linewidth_value =
      distance_to_linewidth(
        distance_pct
      ),
    
    label = paste0(
      
      format(
        
        round(
          distance_pct,
          2
        ),
        
        trim = TRUE,
        
        scientific = FALSE
        
      ),
      
      "%"
      
    )
    
  )


distance_legend <- ggplot2::ggplot() +
  
  ggplot2::geom_segment(
    
    data =
      distance_legend_df,
    
    ggplot2::aes(
      
      x = 0,
      
      xend = 1,
      
      y = row,
      
      yend = row,
      
      linewidth =
        distance_pct
      
    ),
    
    colour = "grey30",
    
    lineend = "round"
    
  ) +
  
  ggplot2::scale_linewidth_continuous(
    
    limits = c(
      0,
      global_divergence_max
    ),
    
    range = c(
      0.35,
      2.20
    ),
    
    guide = "none"
    
  ) +
  
  ggplot2::geom_text(
    
    data =
      distance_legend_df,
    
    ggplot2::aes(
      
      x = 1.18,
      
      y = row,
      
      label = label
      
    ),
    
    hjust = 0,
    
    size = 3
    
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x = 0,
    
    y =
      max(
        distance_legend_df$row
      ) + 0.85,
    
    label =
      "Sequence divergence",
    
    hjust = 0,
    
    fontface = "bold",
    
    size = 3.3
    
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x = 0,
    
    y =
      max(
        distance_legend_df$row
      ) + 0.40,
    
    label =
      "(% differing aligned nucleotides)",
    
    hjust = 0,
    
    size = 2.6,
    
    colour = "grey35"
    
  ) +
  
  ggplot2::coord_cartesian(
    
    xlim = c(
      -0.05,
      1.75
    ),
    
    ylim = c(
      0.5,
      max(
        distance_legend_df$row
      ) + 1.05
    ),
    
    clip = "off"
    
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    
    plot.margin =
      ggplot2::margin(
        2,
        3,
        2,
        3
      )
    
  )


# ============================================================
# 19.17. COMBINE BOTH LEGENDS
# ============================================================

publication_legend <- (
  
  site_legend +
    
    distance_legend +
    
    patchwork::plot_layout(
      
      widths = c(
        2.6,
        1
      )
      
    )
  
)


# ============================================================
# 19.18. NETWORK GRID
# ============================================================
#
# 3 columns rather than 4.
#
# This gives substantially more space to every network.
#
# No main title here:
# for a manuscript, the figure caption provides that
# information and the graphical area should be prioritised.
#
# ============================================================

network_grid <- patchwork::wrap_plots(
  
  publication_panels,
  
  ncol = 3
  
) +
  
  patchwork::plot_annotation(
    
    tag_levels = "A",
    
    theme = ggplot2::theme(
      
      plot.tag =
        ggplot2::element_text(
          
          face = "bold",
          
          size = 12
          
        )
      
    )
    
  )


# ============================================================
# 19.19. FINAL PUBLICATION FIGURE
# ============================================================

main_haplotype_publication <- (
  
  network_grid /
    
    publication_legend
  
) +
  
  patchwork::plot_layout(
    
    heights = c(
      12,
      1.15
    )
    
  )


# ============================================================
# 19.20. SHOW
# ============================================================

print(
  main_haplotype_publication
)


# ============================================================
# 19.21. EXPORT
# ============================================================

publication_svg <- file.path(
  
  paper_haplo_dir,
  
  "MAIN_HAPLOTYPE_NETWORKS_PUBLICATION.svg"
  
)


publication_pdf <- file.path(
  
  paper_haplo_dir,
  
  "MAIN_HAPLOTYPE_NETWORKS_PUBLICATION.pdf"
  
)


publication_png <- file.path(
  
  paper_haplo_dir,
  
  "MAIN_HAPLOTYPE_NETWORKS_PUBLICATION.png"
  
)


# ------------------------------------------------------------
# SVG
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename =
    publication_svg,
  
  plot =
    main_haplotype_publication,
  
  width = 13.5,
  
  height = 19,
  
  units = "in",
  
  bg = "white"
  
)


# ------------------------------------------------------------
# PDF
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename =
    publication_pdf,
  
  plot =
    main_haplotype_publication,
  
  width = 13.5,
  
  height = 19,
  
  units = "in",
  
  device = cairo_pdf,
  
  bg = "white"
  
)


# ------------------------------------------------------------
# PNG
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename =
    publication_png,
  
  plot =
    main_haplotype_publication,
  
  width = 13.5,
  
  height = 19,
  
  units = "in",
  
  dpi = 500,
  
  bg = "white"
  
)


# ============================================================
# 19.22. EXPORT DISTANCE INFORMATION
# ============================================================

distance_information <- tibble::tibble(
  
  metric =
    c(
      "distance_definition",
      "global_maximum_percent"
    ),
  
  value =
    c(
      
      "ape::dist.dna model='raw'; pairwise deletion; multiplied by 100",
      
      as.character(
        global_divergence_max
      )
      
    )
  
)


readr::write_csv(
  
  distance_information,
  
  file.path(
    
    paper_haplo_dir,
    
    "MAIN_HAPLOTYPE_NETWORKS_distance_definition.csv"
    
  )
  
)


# ============================================================
# 19.23. FINAL MESSAGE
# ============================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "PUBLICATION HAPLOTYPE FIGURE FINISHED\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nSaved in:\n",
  normalizePath(
    paper_haplo_dir
  ),
  "\n",
  sep = ""
)

cat(
  "\nSVG:\n",
  publication_svg,
  "\n",
  sep = ""
)

cat(
  "\nPDF:\n",
  publication_pdf,
  "\n",
  sep = ""
)

cat(
  "\nPNG:\n",
  publication_png,
  "\n",
  sep = ""
)

cat(
  "\nSequence divergence range: 0 – ",
  round(
    global_divergence_max,
    3
  ),
  "%\n",
  sep = ""
)

cat(
  "============================================================\n"
)
