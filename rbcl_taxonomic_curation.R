# Load packages
library(dplyr) # used to manipulate data frames and tables
library(tidyr) # used to tidy data 
library(stringr) # string/text processing
library(purrr) # functional programming e.g. map
library(rgbif) # access the GBIF database from R
library(jsonlite)

# =============================================================================
# Semi-automated rbcL taxonomic curation
# =============================================================================
# Uses:
#   1. LULU x BASTA merged CSV
#   2. Eight-column BLAST results
#   3. GBIF occurrence and taxonomy information through rgbif
#   4. Preferred British common names from the NBN Atlas / UKSI
#
# The original merged-table columns are preserved. The requested curation
# columns are appended in the exact order defined near the end of this script.
# =============================================================================

# ------------------------------ 1. USER SETTINGS ------------------------------

# These paths allow the script to be opened in RStudio and run with Source.
lulu_basta_file <- "/home/ucbtns5/Scratch/rbcl/testing_tax_assign/rbcl_lulu98_bastai98l360_merge_test.csv"
blast_file <- "/home/ucbtns5/Scratch/rbcl/testing_tax_assign/rbcl_combined_results_test.txt"
basta_vout_file <- "/home/ucbtns5/Scratch/rbcl/testing_tax_assign/rbcl_bastai98l360m1p70_vout.txt"
output_file <- "/home/ucbtns5/Scratch/rbcl/testing_tax_assign/rbcl_automated_curation_v5.csv"
gbif_cache_file <- "/home/ucbtns5/Scratch/rbcl/testing_tax_assign/rbcl_gbif_cache_v5.csv"

# Optional command-line override:
# Rscript rbcl_taxonomic_curation_v5.R \
#   basta.vout blast.tsv lulu_x_basta.csv curated.csv gbif_cache.csv
# If no arguments are supplied, the USER SETTINGS paths above are used.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 5L) {
  basta_vout_file <- args[1]
  blast_file <- args[2]
  lulu_basta_file <- args[3]
  output_file <- args[4]
  gbif_cache_file <- args[5]
} else if (length(args) != 0L) {
  stop(
    "Supply either no command-line arguments (to use USER SETTINGS) or five: ",
    "<basta_vout.txt> <blast.txt> <lulu_x_basta.csv> ",
    "<curated_output.csv> <gbif_cache.csv>."
  )
}

# A richer private cache is retained beside the requested, compact cache. It
# avoids repeating GBIF taxonomy and occurrence calls on subsequent runs.
gbif_internal_cache_file <- sub(
  "([.]csv)?$", "_internal.csv", gbif_cache_file, ignore.case = TRUE
)
nbn_common_name_cache_file <- sub(
  "([.]csv)?$", "_nbn_common_names.csv", output_file, ignore.case = TRUE
)

country_code <- "GB"
minimum_identity <- 98
minimum_alignment_length <- 400
minimum_gbif_records <- 100

# Exact ties are used by default. Increase only if near-ties should be included.
best_pident_tolerance <- 0
best_alignment_tolerance <- 0

# Set FALSE to test the rest of the script without contacting GBIF.
query_gbif <- TRUE
gbif_pause_seconds <- 0.1

# Set FALSE to reuse an existing NBN cache without making new API requests.
query_nbn <- TRUE
nbn_pause_seconds <- 0.1

taxonomy_ranks <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)

required_packages <- c("dplyr", "stringr", "purrr", "rgbif", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these required packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1])) y else x
}

# ------------------------------ 2. HELPERS ------------------------------------

as_clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | tolower(x) %in% c("na", "nan", "none", "unknown", "unassigned")] <- NA_character_
  x
}

normalise_taxon <- function(x) {
  x <- as_clean_text(x)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

collapse_or_none <- function(x, none = "None") {
  x <- sort(unique(as_clean_text(x)))
  x <- x[!is.na(x)]
  if (!length(x)) none else paste(x, collapse = "; ")
}

# Convert genuinely empty character cells to NA before writing output files.
# This does not change meaningful labels such as "unknown" or "None".
blank_cells_to_na <- function(df) {
  character_columns <- vapply(df, is.character, logical(1))
  df[character_columns] <- lapply(df[character_columns], function(x) {
    x[!is.na(x) & trimws(x) == ""] <- NA_character_
    x
  })
  df
}

format_counts <- function(labels) {
  labels <- as_clean_text(labels)
  labels <- labels[!is.na(labels)]
  if (!length(labels)) return("None")
  counts <- sort(table(labels), decreasing = TRUE)
  paste0(names(counts), " (", as.integer(counts), ")", collapse = "; ")
}

rank_index <- function(x) match(tolower(x), tolower(taxonomy_ranks))

lowest_assigned_rank <- function(row) {
  values <- vapply(taxonomy_ranks, function(r) {
    normalise_taxon(row[[r]])
  }, character(1))
  assigned <- which(!is.na(values))
  if (!length(assigned)) "Unassigned" else taxonomy_ranks[max(assigned)]
}

# Extract the leading taxonomic label for display. This retains placeholders
# such as "Drimia sp." so they can be shown in top_BLAST_hits/removed_entries.
extract_display_label <- function(title) {
  title <- trimws(as.character(title))
  matched <- stringr::str_extract(
    title,
    "^[A-Z][A-Za-z.-]+[[:space:]]+(?:[a-z][A-Za-z0-9.-]*|sp\\.|cf\\.|aff\\.)"
  )
  ifelse(is.na(matched), title, matched)
}

# Accept a species only when the BLAST title begins with a conventional binomial.
# Anchoring at the beginning prevents voucher text such as "Bachtell et" from
# being mistaken for a scientific name.
extract_species <- function(title) {
  title <- trimws(as.character(title))
  species <- stringr::str_extract(
    title,
    "^[A-Z][a-z][A-Za-z.-]*(?:[[:space:]]+(?:x|×))?[[:space:]]+[a-z][A-Za-z.-]+"
  )
  normalise_taxon(species)
}

noninformative_pattern <- paste(
  c(
    "\\bsp\\.?\\b", "\\bcf\\.?\\b", "\\baff\\.?\\b",
    "uncultured", "unidentified", "environmental", "metagenome",
    "unknown", "unclassified", "clone"
  ),
  collapse = "|"
)

is_informative_species <- function(species, title) {
  !is.na(species) &
    grepl(
      "^[A-Z][a-z][A-Za-z.-]*(?:[[:space:]]+(?:x|×))?[[:space:]]+[a-z][A-Za-z.-]+$",
      species, perl = TRUE
    ) &
    !grepl(noninformative_pattern, title, ignore.case = TRUE, perl = TRUE)
}

genus_from_species <- function(species) {
  ifelse(is.na(species), NA_character_, sub("[[:space:]].*$", "", species))
}

# ------------------------------ 3. READ MERGED TABLE --------------------------

if (!file.exists(lulu_basta_file)) stop("Merged file not found: ", lulu_basta_file)

merged <- read.csv(
  lulu_basta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Give blank imported row-name columns a usable name while retaining their values.
bad_names <- is.na(names(merged)) | trimws(names(merged)) == ""
if (any(bad_names)) {
  names(merged)[bad_names] <- paste0("unnamed_column_", seq_len(sum(bad_names)))
}
names(merged) <- make.unique(names(merged))

required_merged_columns <- c("seqid", taxonomy_ranks)
missing_merged <- setdiff(required_merged_columns, names(merged))
if (length(missing_merged)) {
  stop("Merged table is missing: ", paste(missing_merged, collapse = ", "))
}
merged$seqid <- as.character(merged$seqid)
original_columns <- names(merged)

# ------------------------------ 4. READ BLAST ---------------------------------

if (!file.exists(blast_file)) stop("BLAST file not found: ", blast_file)

blast <- read.delim(
  blast_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fill = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(blast) != 8) {
  stop("Expected 8 BLAST columns but found ", ncol(blast), ".")
}

names(blast) <- c(
  "seqid", "qacc", "sseqid", "pident", "alignment_length",
  "mismatch", "evalue", "stitle"
)

blast <- blast %>%
  mutate(
    seqid = as.character(seqid),
    pident = suppressWarnings(as.numeric(pident)),
    alignment_length = suppressWarnings(as.numeric(alignment_length)),
    evalue = suppressWarnings(as.numeric(evalue)),
    display_label = extract_display_label(stitle),
    blast_species = extract_species(stitle),
    informative = is_informative_species(blast_species, stitle),
    blast_genus = if_else(informative, genus_from_species(blast_species), NA_character_)
  )

# Keep each subject accession once per ASV/score/name. This prevents a literal
# duplicate line from inflating number_best_hits while retaining distinct accessions.
blast_unique <- blast %>%
  distinct(
    seqid, sseqid, pident, alignment_length, display_label,
    blast_species, .keep_all = TRUE
  )

# Select the highest identity, then longest alignment among those highest hits.
top_hits <- blast_unique %>%
  group_by(seqid) %>%
  filter(pident >= max(pident, na.rm = TRUE) - best_pident_tolerance) %>%
  filter(
    alignment_length >= max(alignment_length, na.rm = TRUE) -
      best_alignment_tolerance
  ) %>%
  ungroup()

# ------------------------------ 5. READ BASTA VOUT ----------------------------

# The .vout is useful for tracing which lineages contributed to BASTA's merged
# assignment. It is supporting/audit evidence rather than an independent vote if
# it was produced from the same reference search as the BLAST table.
read_sectioned_basta <- function(path) {
  if (!file.exists(path)) stop("BASTA vout file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  current_qid <- NA_character_
  rows <- list()
  index <- 0L

  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next
    if (grepl("^###[[:space:]]*[^[:space:]]+", line)) {
      current_qid <- sub("^###[[:space:]]*", "", line)
      next
    }
    if (is.na(current_qid)) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 2) next
    index <- index + 1L
    rows[[index]] <- tibble(
      seqid = as.character(current_qid),
      support_count = suppressWarnings(as.integer(fields[1])),
      taxonomy_path = paste(fields[-1], collapse = "\t")
    )
  }
  if (!length(rows)) stop("No sectioned BASTA records were parsed from ", path)
  bind_rows(rows)
}

parse_basta_path <- function(path) {
  output <- setNames(rep(NA_character_, length(taxonomy_ranks)), taxonomy_ranks)
  if (is.na(path) || !nzchar(trimws(path))) return(output)
  parts <- trimws(unlist(strsplit(as.character(path), "[;|]")))
  parts <- gsub("^[dkpcofgs]__", "", parts, ignore.case = TRUE)
  parts[parts == "" | tolower(parts) %in% c("unknown", "unassigned", "na", "none")] <- NA_character_
  parts <- gsub("_", " ", parts, fixed = TRUE)
  if (length(parts) > length(output)) parts <- parts[seq_len(length(output))]
  output[seq_along(parts)] <- parts
  output
}

basta_vout <- read_sectioned_basta(basta_vout_file)
basta_vout_ranks <- t(vapply(
  basta_vout$taxonomy_path,
  parse_basta_path,
  FUN.VALUE = setNames(rep(NA_character_, length(taxonomy_ranks)), taxonomy_ranks)
))
basta_vout <- bind_cols(
  basta_vout,
  as.data.frame(basta_vout_ranks, stringsAsFactors = FALSE, check.names = FALSE)
)

vout_supports_taxon <- function(id, rank, taxon) {
  if (is.na(taxon) || !rank %in% taxonomy_ranks) return(FALSE)
  rows <- basta_vout[basta_vout$seqid == as.character(id), , drop = FALSE]
  if (!nrow(rows)) return(FALSE)
  values <- normalise_taxon(rows[[rank]])
  any(!is.na(values) & tolower(values) == tolower(normalise_taxon(taxon)))
}

# ------------------------------ 6. BLAST SUMMARY ------------------------------

summarise_one_blast <- function(df) {
  informative <- df %>% filter(informative)
  removed <- df %>% filter(!informative)

  species <- sort(unique(na.omit(informative$blast_species)))
  genera <- sort(unique(na.omit(informative$blast_genus)))

  tibble(
    top_pident = if (nrow(df)) max(df$pident, na.rm = TRUE) else NA_real_,
    top_alignment = if (nrow(df)) max(df$alignment_length, na.rm = TRUE) else NA_real_,
    top_BLAST_hits = format_counts(df$display_label),
    number_best_hits = nrow(df),
    number_unique_accessions = n_distinct(df$sseqid, na.rm = TRUE),
    number_unique_species = length(species),
    number_unique_genus = length(genera),
    best_species_list = collapse_or_none(species),
    best_genus_list = collapse_or_none(genera),
    removed_entries = format_counts(removed$display_label),
    number_removed_entries = nrow(removed)
  )
}

blast_summary <- top_hits %>%
  group_by(seqid) %>%
  group_modify(~ summarise_one_blast(.x)) %>%
  ungroup()

# Add explicit empty summaries for merged ASVs with no BLAST hits.
missing_blast_ids <- setdiff(merged$seqid, blast_summary$seqid)
if (length(missing_blast_ids)) {
  empty_blast <- tibble(
    seqid = missing_blast_ids,
    top_pident = NA_real_, top_alignment = NA_real_, top_BLAST_hits = "None",
    number_best_hits = 0L, number_unique_accessions = 0L,
    number_unique_species = 0L, number_unique_genus = 0L,
    best_species_list = "None", best_genus_list = "None",
    removed_entries = "None", number_removed_entries = 0L
  )
  blast_summary <- bind_rows(blast_summary, empty_blast)
}

# ------------------------------ 7. GBIF CACHE ---------------------------------

empty_gbif_cache <- tibble(
  queried_name = character(),
  gbif_name = character(),
  gbif_key = integer(),
  gbif_rank = character(),
  gbif_match_type = character(),
  gbif_status_taxonomy = character(),
  common_name = character(),
  kingdom = character(), phylum = character(), class = character(),
  order = character(), family = character(), genus = character(),
  species = character(),
  gbif_records = integer(),
  occurrences_GB_observations = integer(),
  occurrences_GB_specimens = integer(),
  occurrences_GB_all = integer(),
  gbif_record_status = character(),
  query_error = character(),
  query_date = character()
)

if (file.exists(gbif_internal_cache_file)) {
  gbif_cache <- read.csv(
    gbif_internal_cache_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  missing_cache_columns <- setdiff(names(empty_gbif_cache), names(gbif_cache))
  if (length(missing_cache_columns)) {
    # Migrate older caches by adding newly introduced columns. Rows missing the
    # new occurrence counts are automatically refreshed below.
    for (column in missing_cache_columns) {
      prototype <- empty_gbif_cache[[column]]
      gbif_cache[[column]] <- if (is.integer(prototype)) {
        NA_integer_
      } else {
        NA_character_
      }
    }
  }
  gbif_cache <- as_tibble(gbif_cache[, names(empty_gbif_cache), drop = FALSE])
  gbif_cache$gbif_records <- suppressWarnings(as.integer(gbif_cache$gbif_records))
  gbif_cache$occurrences_GB_observations <- suppressWarnings(
    as.integer(gbif_cache$occurrences_GB_observations)
  )
  gbif_cache$occurrences_GB_specimens <- suppressWarnings(
    as.integer(gbif_cache$occurrences_GB_specimens)
  )
  gbif_cache$occurrences_GB_all <- suppressWarnings(
    as.integer(gbif_cache$occurrences_GB_all)
  )
} else {
  gbif_cache <- empty_gbif_cache
}

get_common_name <- function(key) {
  tryCatch({
    result <- rgbif::name_vernaculars(key, limit = 100)$data
    if (is.null(result) || !nrow(result)) return(NA_character_)
    if ("language" %in% names(result)) {
      english <- result[tolower(result$language) %in% c("eng", "en", "english"), , drop = FALSE]
      if (nrow(english)) result <- english
    }
    name_column <- intersect(c("vernacularName", "vernacularname"), names(result))[1]
    if (is.na(name_column)) NA_character_ else as.character(result[[name_column]][1])
  }, error = function(e) NA_character_)
}

check_gbif <- function(species_name) {
  message("Checking GBIF: ", species_name)

  base <- tibble(
    queried_name = species_name,
    gbif_name = NA_character_, gbif_key = NA_integer_, gbif_rank = NA_character_,
    gbif_match_type = NA_character_, gbif_status_taxonomy = NA_character_,
    common_name = NA_character_, kingdom = NA_character_, phylum = NA_character_,
    class = NA_character_, order = NA_character_, family = NA_character_,
    genus = NA_character_, species = NA_character_, gbif_records = NA_integer_,
    occurrences_GB_observations = NA_integer_,
    occurrences_GB_specimens = NA_integer_,
    occurrences_GB_all = NA_integer_,
    gbif_record_status = "Query uncertain", query_error = NA_character_,
    query_date = as.character(Sys.Date())
  )

  if (!query_gbif) {
    base$query_error <- "GBIF querying disabled"
    return(base)
  }

  tryCatch({
    match_result <- rgbif::name_backbone(name = species_name, rank = "species")
    if ("data" %in% names(match_result)) match_result <- match_result$data
    if (is.data.frame(match_result)) {
      if (!nrow(match_result)) {
        base$query_error <- "No GBIF name match"
        return(base)
      }
      match_result <- as.list(match_result[1, , drop = FALSE])
      match_result <- lapply(match_result, function(z) z[[1]])
    }

    key <- match_result$usageKey %||% NA_integer_
    if (is.na(key)) {
      base$query_error <- "No GBIF usageKey"
      return(base)
    }

    base$gbif_key <- as.integer(key)
    base$gbif_name <- as.character(
      match_result$scientificName %||% match_result$canonicalName %||% species_name
    )
    base$gbif_rank <- as.character(match_result$rank %||% NA_character_)
    base$gbif_match_type <- as.character(match_result$matchType %||% NA_character_)
    base$gbif_status_taxonomy <- as.character(match_result$status %||% NA_character_)

    for (rank in c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
      base[[rank]] <- as.character(match_result[[rank]] %||% NA_character_)
    }
    base$common_name <- get_common_name(key)

    # Use the complete GBIF country count so occurrences_GB matches the count
    # displayed by a normal GBIF taxon/country query. The value is evidence of
    # records, not proof of native or established status.
    total <- as.numeric(rgbif::occ_count(
      taxonKey = as.integer(key),
      country = country_code
    ))
    observations <- as.numeric(rgbif::occ_count(
      taxonKey = as.integer(key),
      country = country_code,
      basisOfRecord = paste(
        c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OBSERVATION"),
        collapse = ";"
      ),
      occurrenceStatus = "PRESENT"
    ))
    specimens <- as.numeric(rgbif::occ_count(
      taxonKey = as.integer(key),
      country = country_code,
      basisOfRecord = paste(
        c(
          "PRESERVED_SPECIMEN", "LIVING_SPECIMEN", "FOSSIL_SPECIMEN",
          "MATERIAL_CITATION", "MATERIAL_SAMPLE"
        ),
        collapse = ";"
      ),
      occurrenceStatus = "PRESENT"
    ))
    base$gbif_records <- as.integer(total)
    base$occurrences_GB_observations <- as.integer(observations)
    base$occurrences_GB_specimens <- as.integer(specimens)
    base$occurrences_GB_all <- as.integer(total)
    base$gbif_record_status <- case_when(
      total >= minimum_gbif_records ~ "Above threshold",
      total > 0 ~ "Below threshold",
      total == 0 ~ "No qualifying records",
      TRUE ~ "Query uncertain"
    )
    base
  }, error = function(e) {
    base$query_error <- conditionMessage(e)
    base
  })
}

species_to_check <- sort(unique(top_hits$blast_species[top_hits$informative]))
species_to_check <- species_to_check[!is.na(species_to_check)]

# A cached row is reusable only when the occurrence query completed. Earlier
# versions used an incorrect rgbif argument, so rows with NA counts or a stored
# error must be queried again rather than silently treated as British absences.
valid_cached_names <- gbif_cache$queried_name[
  !is.na(gbif_cache$gbif_records) &
    !is.na(gbif_cache$occurrences_GB_observations) &
    !is.na(gbif_cache$occurrences_GB_specimens) &
    !is.na(gbif_cache$occurrences_GB_all) &
    (is.na(gbif_cache$query_error) | !nzchar(gbif_cache$query_error))
]
missing_from_cache <- setdiff(species_to_check, valid_cached_names)

if (length(missing_from_cache)) {
  gbif_cache <- gbif_cache %>%
    filter(!queried_name %in% missing_from_cache)

  new_gbif <- map_dfr(seq_along(missing_from_cache), function(i) {
    message("[", i, "/", length(missing_from_cache), "]")
    result <- check_gbif(missing_from_cache[i])
    if (gbif_pause_seconds > 0) Sys.sleep(gbif_pause_seconds)
    result
  })
  gbif_cache <- bind_rows(gbif_cache, new_gbif) %>%
    distinct(queried_name, .keep_all = TRUE)
}

gbif_cache <- blank_cells_to_na(gbif_cache)
write.csv(gbif_cache, gbif_internal_cache_file, row.names = FALSE, na = "NA")

# ------------------------------ 8. GBIF SUMMARY -------------------------------

top_informative_gbif <- top_hits %>%
  filter(informative) %>%
  distinct(seqid, blast_species) %>%
  left_join(gbif_cache, by = c("blast_species" = "queried_name"))

# Write the compact, ASV-specific GBIF evidence table requested for manual
# auditing. `support` is the BASTA vout support attached to that exact species;
# it is not treated as an independent BLAST vote. `present_GB` means at least
# one GBIF record of any type, while the main automated summary applies
# minimum_gbif_records to observation records only.
support_for_species <- function(id, taxon) {
  rows <- basta_vout[basta_vout$seqid == as.character(id), , drop = FALSE]
  if (!nrow(rows)) return(0L)
  matched <- !is.na(rows$Species) &
    tolower(normalise_taxon(rows$Species)) == tolower(normalise_taxon(taxon))
  if (!any(matched)) 0L else as.integer(sum(rows$support_count[matched], na.rm = TRUE))
}

gbif_cache_output <- top_informative_gbif %>%
  transmute(
    seqid = as.character(seqid),
    support = map2_int(seqid, blast_species, support_for_species),
    species = blast_species,
    gbif_name = gbif_name,
    usageKey = as.integer(gbif_key),
    occurrences_GB = as.integer(gbif_records),
    present_GB = !is.na(gbif_records) & gbif_records > 0
  ) %>%
  arrange(suppressWarnings(as.numeric(seqid)), seqid, species)

gbif_cache_output <- blank_cells_to_na(gbif_cache_output)
write.csv(gbif_cache_output, gbif_cache_file, row.names = FALSE, na = "NA")

gbif_summary <- top_informative_gbif %>%
  group_by(seqid) %>%
  summarise(
    number_species_checked = n_distinct(blast_species[!is.na(blast_species)]),
    number_species_present = n_distinct(
      blast_species[
        !is.na(occurrences_GB_observations) &
          occurrences_GB_observations >= minimum_gbif_records
      ]
    ),
    gbif_present_species = {
      keep <- !is.na(occurrences_GB_observations) &
        occurrences_GB_observations >= minimum_gbif_records
      if (!any(keep)) "None" else paste0(
        blast_species[keep], " (", occurrences_GB_observations[keep], ")",
        collapse = "; "
      )
    },
    .groups = "drop"
  )

missing_gbif_ids <- setdiff(merged$seqid, gbif_summary$seqid)
if (length(missing_gbif_ids)) {
  gbif_summary <- bind_rows(
    gbif_summary,
    tibble(
      seqid = missing_gbif_ids,
      number_species_checked = 0L,
      number_species_present = 0L,
      gbif_present_species = "None"
    )
  )
}

# ------------------------------ 9. EVIDENCE TABLE -----------------------------

evidence <- merged %>%
  left_join(blast_summary, by = "seqid") %>%
  left_join(gbif_summary, by = "seqid")

# Obtain GBIF taxonomy for the informative best species belonging to one ASV.
gbif_rows_for_asv <- function(id) {
  species <- top_hits$blast_species[top_hits$seqid == id & top_hits$informative]
  species <- unique(species[!is.na(species)])
  gbif_cache[match(species, gbif_cache$queried_name), , drop = FALSE]
}

taxonomy_consensus <- function(rows, rank) {
  if (is.null(rows) || !nrow(rows) || !rank %in% names(rows)) return(NA_character_)
  values <- unique(normalise_taxon(rows[[rank]]))
  values <- values[!is.na(values)]
  if (length(values) == 1 && all(!is.na(normalise_taxon(rows[[rank]])))) values else NA_character_
}

# ------------------------------ 10. CURATION RULES ----------------------------

curate_row <- function(i) {
  row <- evidence[i, , drop = FALSE]
  id <- as.character(row$seqid)
  basta_rank <- lowest_assigned_rank(row)
  basta_values <- setNames(
    vapply(taxonomy_ranks, function(r) normalise_taxon(row[[r]]), character(1)),
    taxonomy_ranks
  )

  gbif_tax <- gbif_rows_for_asv(id)
  final <- basta_values
  automated_rank <- basta_rank
  manual <- FALSE
  reason <- ""

  reliable <- !is.na(row$top_pident) && !is.na(row$top_alignment) &&
    row$top_pident >= minimum_identity &&
    row$top_alignment >= minimum_alignment_length

  n_species <- row$number_unique_species %||% 0L
  n_genera <- row$number_unique_genus %||% 0L
  informative_species <- unique(
    top_hits$blast_species[top_hits$seqid == id & top_hits$informative]
  )
  informative_species <- informative_species[!is.na(informative_species)]
  informative_genera <- unique(genus_from_species(informative_species))
  informative_genera <- informative_genera[!is.na(informative_genera)]

  # Recalculate consensus after applying the GB observation threshold. Use
  # GBIF-accepted species/genus names where available, so synonyms do not create
  # artificial cross-genus conflicts.
  asv_gbif_present <- top_informative_gbif[
    top_informative_gbif$seqid == id &
      !is.na(top_informative_gbif$occurrences_GB_observations) &
      top_informative_gbif$occurrences_GB_observations >=
        minimum_gbif_records,
    , drop = FALSE
  ]
  present_species <- if (nrow(asv_gbif_present)) {
    accepted <- normalise_taxon(asv_gbif_present$species)
    accepted[is.na(accepted)] <- normalise_taxon(
      asv_gbif_present$blast_species[is.na(accepted)]
    )
    unique(accepted[!is.na(accepted)])
  } else {
    character()
  }
  present_genera <- if (nrow(asv_gbif_present)) {
    accepted <- normalise_taxon(asv_gbif_present$genus)
    fallback <- genus_from_species(asv_gbif_present$blast_species)
    accepted[is.na(accepted)] <- fallback[is.na(accepted)]
    unique(accepted[!is.na(accepted)])
  } else {
    character()
  }
  decision_confidence <- NULL
  decision_evidence_basis <- NULL

  if (!reliable || n_species == 0) {
    automated_rank <- basta_rank
    manual <- TRUE
    reason <- paste0(
      "Retained the LULU x BASTA rank because no informative BLAST evidence met ",
      minimum_identity, "% identity and ", minimum_alignment_length,
      " bp alignment thresholds."
    )
  } else if (n_species == 1) {
    automated_rank <- "Species"
    final["Species"] <- informative_species[1]
    final["Genus"] <- genus_from_species(informative_species[1])
    for (r in c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
      value <- taxonomy_consensus(gbif_tax, r)
      if (!is.na(value)) final[taxonomy_ranks[match(r, tolower(taxonomy_ranks))]] <- value
    }

    species_agrees <- !is.na(basta_values["Species"]) &&
      tolower(basta_values["Species"]) == tolower(informative_species[1])
    genus_agrees <- !is.na(basta_values["Genus"]) &&
      tolower(basta_values["Genus"]) == tolower(genus_from_species(informative_species[1]))
    vout_species_support <- vout_supports_taxon(
      id, "Species", informative_species[1]
    )

    # Missing lower ranks are not conflicts. A clean species upgrade is allowed
    # when the vout contains that species and any already-assigned genus/species
    # agrees. Only an explicitly different assigned taxon is a conflict.
    assigned_species_conflict <- !is.na(basta_values["Species"]) &&
      tolower(basta_values["Species"]) != tolower(informative_species[1])
    assigned_genus_conflict <- !is.na(basta_values["Genus"]) &&
      tolower(basta_values["Genus"]) !=
        tolower(genus_from_species(informative_species[1]))

    # GBIF zero/low support is descriptive only and never triggers review.
    manual <- assigned_species_conflict || assigned_genus_conflict ||
      !vout_species_support
    reason <- if (species_agrees && vout_species_support) {
      "Retained species because all informative top BLAST hits support one species, the merged LULU x BASTA assignment agrees, and the species is represented in the BASTA vout. GBIF support is reported descriptively only."
    } else if (!assigned_species_conflict && !assigned_genus_conflict &&
               vout_species_support) {
      "Upgraded to species because the credible optimal BLAST hits support one species, that species occurs in the BASTA vout, and no assigned BASTA rank conflicts with its lineage. GBIF support is reported descriptively only."
    } else if (!vout_species_support) {
      "BLAST supports one species, but that species is not represented in the parsed BASTA vout evidence; retained the molecular proposal and flagged the discrepancy for review. GBIF support is descriptive only."
    } else {
      "BLAST supports one species but conflicts with the merged LULU x BASTA lineage; flagged for manual review. GBIF support is descriptive only."
    }
  } else if (length(present_species) == 1) {
    # Geographic tie-breaker: exactly one accepted species remains after the
    # GB observation filter.
    asv_gbif <- asv_gbif_present
    preferred_species <- present_species

    if (length(preferred_species) == 1) {
      preferred_species <- preferred_species[1]
      preferred_genus <- if (length(present_genera) == 1) {
        present_genera[1]
      } else {
        genus_from_species(preferred_species)
      }
      preferred_gbif <- asv_gbif_present

      automated_rank <- "Species"
      decision_confidence <- "Moderate"
      decision_evidence_basis <- "Geography-resolved species"
      final["Species"] <- preferred_species
      final["Genus"] <- preferred_genus
      for (r in c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
        value <- taxonomy_consensus(preferred_gbif, r)
        if (!is.na(value)) final[taxonomy_ranks[match(r, tolower(taxonomy_ranks))]] <- value
      }

      vout_species_support <- vout_supports_taxon(
        id, "Species", preferred_species
      )
      merged_species <- basta_values["Species"]
      merged_genus <- basta_values["Genus"]
      merged_lineage_support <- if (!is.na(merged_species)) {
        tolower(merged_species) == tolower(preferred_species)
      } else {
        !is.na(merged_genus) && tolower(merged_genus) == tolower(preferred_genus)
      }

      manual <- !vout_species_support || !merged_lineage_support
      reason <- if (!manual) {
        paste0(
          "Several species had equally best BLAST scores, but only ",
          preferred_species, " exceeded the GBIF threshold of ",
          minimum_gbif_records,
          " GB observation records. The merged BASTA lineage and BASTA vout ",
          "also support this taxon, so GBIF was used as a geographic tie-breaker."
        )
      } else {
        paste0(
          "Several species had equally best BLAST scores and only ",
          preferred_species, " exceeded the GB observation threshold. It was selected as ",
          "the geographic preference, but the merged BASTA lineage or BASTA vout ",
          "does not support it consistently, so manual review is required."
        )
      }
    }
  } else if (length(present_species) > 1 && length(present_genera) == 1) {
    # Several species remain geographically plausible, but all share one
    # accepted genus. This safely resolves genus, though not species.
    preferred_genus <- present_genera[1]
    automated_rank <- "Genus"
    final["Genus"] <- preferred_genus
    final["Species"] <- NA_character_
    for (r in c("kingdom", "phylum", "class", "order", "family", "genus")) {
      value <- taxonomy_consensus(asv_gbif_present, r)
      if (!is.na(value)) {
        final[taxonomy_ranks[match(r, tolower(taxonomy_ranks))]] <- value
      }
    }

    assigned_genus_conflict <- !is.na(basta_values["Genus"]) &&
      tolower(basta_values["Genus"]) != tolower(preferred_genus)
    manual <- assigned_genus_conflict
    decision_confidence <- "Moderate"
    decision_evidence_basis <- "Geography-resolved genus"
    reason <- if (!manual) {
      paste0(
        "Several species had equally best BLAST scores. After applying the GB ",
        "observation threshold and normalising GBIF taxonomy, all plausible ",
        "candidates belong to ", preferred_genus,
        "; assigned genus because the marker does not resolve species."
      )
    } else {
      paste0(
        "After applying the GB observation threshold, all plausible candidates ",
        "belong to ", preferred_genus,
        ", but this conflicts with the assigned BASTA genus; manual review is required."
      )
    }
  } else if (n_genera == 1) {
    automated_rank <- "Genus"
    final["Genus"] <- informative_genera[1]
    final["Species"] <- NA_character_
    for (r in c("kingdom", "phylum", "class", "order", "family", "genus")) {
      value <- taxonomy_consensus(gbif_tax, r)
      if (!is.na(value)) final[taxonomy_ranks[match(r, tolower(taxonomy_ranks))]] <- value
    }
    genus_agrees <- !is.na(basta_values["Genus"]) &&
      tolower(basta_values["Genus"]) == tolower(informative_genera[1])
    vout_genus_support <- vout_supports_taxon(id, "Genus", informative_genera[1])
    manual <- !genus_agrees || !vout_genus_support
    reason <- if (genus_agrees && vout_genus_support) {
      "Assigned genus because equally best informative BLAST hits contain multiple species from one genus, and that genus agrees with both the merged LULU x BASTA assignment and BASTA vout. GBIF support is descriptive only."
    } else {
      "Equally best BLAST hits support one genus, but the merged LULU x BASTA assignment or BASTA vout does not support the same genus; assigned genus and flagged for review. GBIF support is descriptive only."
    }
  } else {
    shared_family <- taxonomy_consensus(gbif_tax, "family")
    shared_order <- taxonomy_consensus(gbif_tax, "order")
    if (!is.na(shared_family)) {
      automated_rank <- "Family"
      final["Family"] <- shared_family
      final[c("Genus", "Species")] <- NA_character_
      manual <- TRUE
      reason <- "Equally best informative BLAST hits span multiple genera but share one GBIF family; assigned family and flagged for review."
    } else if (!is.na(shared_order)) {
      automated_rank <- "Order"
      final["Order"] <- shared_order
      final[c("Family", "Genus", "Species")] <- NA_character_
      manual <- TRUE
      reason <- "Equally best informative BLAST hits span multiple families but share one GBIF order; assigned order and flagged for review."
    } else {
      automated_rank <- basta_rank
      final <- basta_values
      manual <- TRUE
      reason <- "Equally best informative BLAST hits span conflicting higher taxa; retained the LULU x BASTA assignment and flagged for review."
    }
  }

  final_rank_index <- rank_index(automated_rank)
  if (!is.na(final_rank_index) && final_rank_index < length(taxonomy_ranks)) {
    final[(final_rank_index + 1):length(taxonomy_ranks)] <- NA_character_
  }
  basta_rank_index <- rank_index(basta_rank)

  curation_decision <- if (manual && (
    is.na(final_rank_index) ||
      (automated_rank == basta_rank && grepl("conflict|no informative", reason, ignore.case = TRUE))
  )) {
    "Manual review"
  } else if (is.na(basta_rank_index) || is.na(final_rank_index)) {
    "Manual review"
  } else if (final_rank_index == basta_rank_index) {
    paste("Retained", tolower(automated_rank))
  } else if (final_rank_index < basta_rank_index) {
    paste("Downgraded to", tolower(automated_rank))
  } else if (manual) {
    paste("Review possible upgrade to", tolower(automated_rank))
  } else {
    paste("Upgraded to", tolower(automated_rank))
  }

  final_taxon <- if (!is.na(final_rank_index)) final[final_rank_index] else NA_character_

  hit_consensus <- if (!reliable || n_species == 0) {
    "No informative BLAST consensus passed the configured quality thresholds."
  } else if (n_species == 1) {
    paste0(
      "All credible optimal BLAST hits support ", informative_species[1],
      " (", row$number_unique_accessions, " unique accession(s))."
    )
  } else if (n_genera == 1) {
    paste0(
      "Credible optimal hits contain ", n_species, " species but all belong to ",
      informative_genera[1], "."
    )
  } else {
    paste0(
      "Credible optimal hits contain ", n_species, " species across ",
      n_genera, " genera; the sequence does not resolve genus."
    )
  }

  gbif_status_text <- if ((row$number_species_checked %||% 0L) == 0L) {
    "No informative species was available for a GBIF check."
  } else if ((row$number_species_present %||% 0L) == 0L) {
    paste0(
      "No candidate reached the configured GB observation threshold of ",
      minimum_gbif_records, "; inspect the GBIF cache before geographic filtering."
    )
  } else {
    paste0(
      "GB candidates above the observation threshold: ",
      row$gbif_present_species, "."
    )
  }

  final_common_name <- NA_character_
  if (identical(automated_rank, "Species") && !is.na(final_taxon) && nrow(gbif_tax)) {
    common_match <- gbif_tax$common_name[
      !is.na(gbif_tax$queried_name) &
        tolower(gbif_tax$queried_name) == tolower(final_taxon)
    ]
    common_match <- common_match[!is.na(common_match) & nzchar(common_match)]
    if (length(common_match)) final_common_name <- common_match[1]
  }

  confidence <- if (!is.null(decision_confidence)) {
    decision_confidence
  } else if (!reliable || n_species == 0) {
    "Low"
  } else if (n_species == 1 && row$number_unique_accessions >= 2) {
    "Very high"
  } else if (n_species == 1 || n_genera == 1) {
    "High"
  } else {
    "Low"
  }

  evidence_basis <- if (!is.null(decision_evidence_basis)) {
    decision_evidence_basis
  } else if (!reliable || n_species == 0) {
    "Unresolved"
  } else if (n_species == 1 || n_genera == 1) {
    "Sequence-resolved"
  } else {
    "Manual geographic/ecological review required"
  }

  # Colour is based on structured evidence, not words in the reason text.
  # This prevents phrases such as "no assigned rank conflicts" from being
  # incorrectly classified as red merely because they contain "conflict".
  curation_colour <- if (isTRUE(manual)) {
    "Red"
  } else if (!is.na(basta_rank_index) && !is.na(final_rank_index) &&
             basta_rank_index == final_rank_index) {
    "Green"
  } else {
    "Yellow"
  }

  notes <- if (manual) {
    "Automated proposal requires manual confirmation; do not use BLAST row order or duplicate reference counts to break exact ties."
  } else {
    "Automated conservative assignment; GBIF occurrence is supporting context rather than molecular evidence."
  }

  tibble(
    BASTA_rank = basta_rank,
    final_kingdom = final["Kingdom"],
    final_phylum = final["Phylum"],
    final_class = final["Class"],
    final_order = final["Order"],
    final_family = final["Family"],
    final_genus = final["Genus"],
    final_species = final["Species"],
    final_taxon = final_taxon,
    final_rank = automated_rank,
    common_name = final_common_name,
    hit_consensus = hit_consensus,
    GBIF_UK_status = gbif_status_text,
    curation_decision = curation_decision,
    manual_review = manual,
    curation_reason = reason,
    confidence = confidence,
    evidence_basis = evidence_basis,
    curation_colour = curation_colour,
    notes = notes
  )
}

message("Curating ", nrow(evidence), " ASVs...")
decision_rows <- map_dfr(seq_len(nrow(evidence)), curate_row)

# ----------------------- 11. NBN BRITISH COMMON NAMES ------------------------

# The NBN Atlas species dictionary is based on the UK Species Inventory (UKSI).
# commonNameSingle is its preferred common name. These names are descriptive
# output only and never influence the taxonomic curation decision.
empty_nbn_cache <- tibble(
  queried_taxon = character(),
  queried_rank = character(),
  nbn_scientific_name = character(),
  nbn_rank = character(),
  nbn_guid = character(),
  british_common_name = character(),
  nbn_taxonomic_status = character(),
  nbn_match_status = character(),
  query_error = character(),
  query_date = character()
)

if (file.exists(nbn_common_name_cache_file)) {
  nbn_cache <- read.csv(
    nbn_common_name_cache_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "")
  )
  missing_nbn_columns <- setdiff(names(empty_nbn_cache), names(nbn_cache))
  for (column in missing_nbn_columns) nbn_cache[[column]] <- NA_character_
  nbn_cache <- as_tibble(nbn_cache[, names(empty_nbn_cache), drop = FALSE])
} else {
  nbn_cache <- empty_nbn_cache
}

check_nbn_common_name <- function(taxon, rank) {
  message("Checking NBN common name: ", taxon, " [", rank, "]")

  base <- tibble(
    queried_taxon = taxon,
    queried_rank = rank,
    nbn_scientific_name = NA_character_,
    nbn_rank = NA_character_,
    nbn_guid = NA_character_,
    british_common_name = NA_character_,
    nbn_taxonomic_status = NA_character_,
    nbn_match_status = "No exact accepted match",
    query_error = NA_character_,
    query_date = as.character(Sys.Date())
  )

  if (!query_nbn) {
    base$query_error <- "NBN querying disabled"
    return(base)
  }

  tryCatch({
    encoded_taxon <- utils::URLencode(taxon, reserved = TRUE)
    encoded_rank <- utils::URLencode(tolower(rank), reserved = TRUE)
    url <- paste0(
      "https://species-ws.nbnatlas.org/search?",
      "q=%2A%3A%2A",
      "&fq=idxtype%3ATAXON",
      "&fq=scientificName%3A%22", encoded_taxon, "%22",
      "&fq=rank%3A%22", encoded_rank, "%22",
      "&pageSize=20"
    )

    response <- jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
    matches <- response$searchResults$results
    if (is.null(matches) || !is.data.frame(matches) || !nrow(matches)) {
      return(base)
    }

    exact <- matches[
      !is.na(matches$scientificName) &
        tolower(trimws(matches$scientificName)) == tolower(trimws(taxon)) &
        !is.na(matches$rank) &
        tolower(matches$rank) == tolower(rank) &
        !is.na(matches$taxonomicStatus) &
        tolower(matches$taxonomicStatus) == "accepted",
      , drop = FALSE
    ]
    if (!nrow(exact)) return(base)

    match <- exact[1, , drop = FALSE]
    preferred_name <- as_clean_text(match$commonNameSingle %||% NA_character_)

    base$nbn_scientific_name <- as.character(match$scientificName[1] %||% NA_character_)
    base$nbn_rank <- as.character(match$rank[1] %||% NA_character_)
    base$nbn_guid <- as.character(match$guid[1] %||% NA_character_)
    base$british_common_name <- preferred_name[1] %||% NA_character_
    base$nbn_taxonomic_status <- as.character(
      match$taxonomicStatus[1] %||% NA_character_
    )
    base$nbn_match_status <- if (is.na(base$british_common_name)) {
      "Exact accepted taxon; no preferred common name"
    } else {
      "Exact accepted taxon with preferred common name"
    }
    base
  }, error = function(e) {
    base$query_error <- conditionMessage(e)
    base$nbn_match_status <- "Query uncertain"
    base
  })
}

taxa_to_check <- decision_rows %>%
  transmute(
    queried_taxon = normalise_taxon(final_taxon),
    queried_rank = as_clean_text(final_rank)
  ) %>%
  filter(!is.na(queried_taxon), !is.na(queried_rank)) %>%
  distinct() %>%
  arrange(queried_taxon, queried_rank)

nbn_key <- function(taxon, rank) {
  paste(tolower(trimws(taxon)), tolower(trimws(rank)), sep = "||")
}

valid_nbn_cache <- is.na(nbn_cache$query_error) | !nzchar(trimws(nbn_cache$query_error))
cached_nbn_keys <- nbn_key(
  nbn_cache$queried_taxon[valid_nbn_cache],
  nbn_cache$queried_rank[valid_nbn_cache]
)
required_nbn_keys <- nbn_key(taxa_to_check$queried_taxon, taxa_to_check$queried_rank)
new_nbn_taxa <- taxa_to_check[!required_nbn_keys %in% cached_nbn_keys, , drop = FALSE]

if (nrow(new_nbn_taxa)) {
  retry_nbn_keys <- nbn_key(new_nbn_taxa$queried_taxon, new_nbn_taxa$queried_rank)
  existing_nbn_keys <- nbn_key(nbn_cache$queried_taxon, nbn_cache$queried_rank)
  nbn_cache <- nbn_cache[!existing_nbn_keys %in% retry_nbn_keys, , drop = FALSE]
  new_nbn_rows <- map2_dfr(
    new_nbn_taxa$queried_taxon,
    new_nbn_taxa$queried_rank,
    function(taxon, rank) {
      answer <- check_nbn_common_name(taxon, rank)
      if (nbn_pause_seconds > 0) Sys.sleep(nbn_pause_seconds)
      answer
    }
  )
  nbn_cache <- bind_rows(nbn_cache, new_nbn_rows)
}

nbn_cache <- nbn_cache %>%
  distinct(queried_taxon, queried_rank, .keep_all = TRUE) %>%
  arrange(queried_taxon, queried_rank)
write.csv(
  blank_cells_to_na(nbn_cache),
  nbn_common_name_cache_file,
  row.names = FALSE,
  na = "NA"
)

decision_nbn_keys <- nbn_key(decision_rows$final_taxon, decision_rows$final_rank)
cache_nbn_keys <- nbn_key(nbn_cache$queried_taxon, nbn_cache$queried_rank)
nbn_common_names <- nbn_cache$british_common_name[
  match(decision_nbn_keys, cache_nbn_keys)
]
use_nbn_name <- !is.na(nbn_common_names) & nzchar(trimws(nbn_common_names))
decision_rows$common_name[use_nbn_name] <- nbn_common_names[use_nbn_name]

# ------------------------------ 12. ORDER AND SAVE ----------------------------

requested_columns <- c(
  "top_pident",
  "top_alignment",
  "top_BLAST_hits",
  "number_best_hits",
  "number_unique_accessions",
  "number_unique_species",
  "number_unique_genus",
  "best_species_list",
  "best_genus_list",
  "removed_entries",
  "number_removed_entries",
  "number_species_checked",
  "number_species_present",
  "gbif_present_species"
)

result <- bind_cols(
  evidence[, original_columns, drop = FALSE],
  evidence[, requested_columns, drop = FALSE],
  decision_rows
)

result <- blank_cells_to_na(result)
write.csv(result, output_file, row.names = FALSE, na = "NA")

message("Curation complete.")
message("Output: ", normalizePath(output_file, mustWork = FALSE))
message("GBIF cache: ", normalizePath(gbif_cache_file, mustWork = FALSE))
message("Internal reusable GBIF cache: ", normalizePath(gbif_internal_cache_file, mustWork = FALSE))
message("NBN common-name cache: ", normalizePath(nbn_common_name_cache_file, mustWork = FALSE))
message(
  "Manual-review ASVs: ",
  sum(result$manual_review %in% TRUE, na.rm = TRUE),
  " of ", nrow(result)
)
