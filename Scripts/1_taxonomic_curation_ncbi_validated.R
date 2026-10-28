# Load packages
library(dplyr) # used to manipulate data frames and tables
library(tidyr) # used to tidy data
library(stringr) # string/text processing
library(purrr) # functional programming e.g. map
library(rgbif) #access the GBIF database from R
library(jsonlite)

# =============================================================================
# Semi-automated LgV2 taxonomic curation
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

lulu_basta_file <- "path/to/merged_LULU_BASTA_ASV_table.csv"
blast_file <- "path/to/combined_BLAST_results.txt"
basta_vout_file <- "path/to/BASTA_lineage_output.txt"
output_file <- "path/to/curated_ASV_table.csv"
gbif_cache_file <- "path/to/GBIF_evidence_cache.csv"

# Optional command-line override:
# Rscript taxonomic_curation_ncbi_validated.R \
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
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(gbif_cache_file), recursive = TRUE, showWarnings = FALSE)

gbif_internal_cache_file <- sub(
  "([.]csv)?$", "_internal.csv", gbif_cache_file, ignore.case = TRUE
)
nbn_common_name_cache_file <- sub(
  "([.]csv)?$", "_nbn_common_names.csv", output_file, ignore.case = TRUE
)

# One execution mode for BLAST diagnostics and GBIF threshold exploration.
analysis_mode <- "curate"  # "explore" or "curate"

# Taxonomic scope: Orange means outside this scope, not necessarily incorrect.
target_phyla <- c("Arthropoda")

# Marker-specific parameters (provisional examples; calibrate per marker).
expected_amplicon_length <- 400
species_min_identity <- 99
genus_min_identity <- 95
family_min_identity <- 90
order_min_identity <- 85
species_min_alignment <- 300
genus_min_alignment <- 260
family_min_alignment <- 190
order_min_alignment <- 140
species_pident_margin <- 0.2
genus_pident_margin <- 0.3
family_pident_margin <- 0.4
order_pident_margin <- 0.6
species_alignment_margin <- 10
genus_alignment_margin <- 20
family_alignment_margin <- 30
order_alignment_margin <- 40


# Only distinct subject accessions count as reference support. Support is
# secondary; a majority never overrides a materially superior alignment.
min_accession_advantage <- 2L

country_code <- "GB"
minimum_gbif_records <- 20
gbif_thresholds_to_test <- c(5, 10, 15, 20, 25, 50, 100)

if (!analysis_mode %in% c("explore", "curate")) {
  stop("analysis_mode must be 'explore' or 'curate'.")
}
for (rank in c("species", "genus", "family", "order")) {
  ident <- get(paste0(rank, "_min_identity"))
  length_min <- get(paste0(rank, "_min_alignment"))
  pi_margin <- get(paste0(rank, "_pident_margin"))
  len_margin <- get(paste0(rank, "_alignment_margin"))
  if (!is.finite(ident) || ident < 0 || ident > 100 ||
      !is.finite(length_min) || length_min <= 0 ||
      !is.finite(pi_margin) || pi_margin < 0 ||
      !is.finite(len_margin) || len_margin < 0) {
    stop("Invalid taxonomic thresholds for ", rank)
  }
}
if (!(species_min_identity >= genus_min_identity &&
      genus_min_identity >= family_min_identity &&
      family_min_identity >= order_min_identity &&
      species_min_alignment >= genus_min_alignment &&
      genus_min_alignment >= family_min_alignment &&
      family_min_alignment >= order_min_alignment)) {
  stop("Minimum identity and alignment thresholds must decrease from species to order.")
}

# Control columns are detected from their sample names. Edit these expressions
# if another primer run uses different naming conventions. Matching is case
# insensitive. A count at or above control_detection_minimum is reported in
# `notes`, but control presence does not automatically change the taxonomy.
# A full stop is accepted as a boundary because R may rename duplicate sample
# headings as POS.1, POS.2, and so on.
positive_control_pattern <- "(^|[_.-])(pos|positive|posctrl|positivectrl)([_.-]|$)"
negative_control_pattern <- "(pcrneg|negative|negctrl|negativectrl|ntc|blank)"
control_detection_minimum <- 1

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

format_control_note <- function(row, columns, label) {
  if (!length(columns)) return(NA_character_)
  counts <- suppressWarnings(as.numeric(unlist(row[1, columns, drop = FALSE])))
  names(counts) <- columns
  detected <- !is.na(counts) & counts >= control_detection_minimum
  if (!any(detected)) return(NA_character_)
  paste0(
    label, ": ",
    paste0(names(counts)[detected], "=", counts[detected], collapse = "; ")
  )
}

control_notes_for_row <- function(row) {
  notes <- c(
    format_control_note(row, positive_control_columns, "Positive-control detection"),
    format_control_note(row, negative_control_columns, "Negative-control detection")
  )
  notes <- notes[!is.na(notes)]
  if (!length(notes)) NA_character_ else paste(notes, collapse = " | ")
}

rank_index <- function(x) match(tolower(x), tolower(taxonomy_ranks))

lowest_assigned_rank <- function(row) {
  values <- vapply(taxonomy_ranks, function(r) {
    normalise_taxon(row[[r]])
  }, character(1))
  assigned <- which(!is.na(values))
  if (!length(assigned)) "Unassigned" else taxonomy_ranks[max(assigned)]
}

# Remove database-description prefixes that appear before an otherwise usable
# scientific name. These records are treated exactly like ordinary named
# GenBank records and are not separately flagged.
remove_informative_record_prefix <- function(title) {
  title <- trimws(as.character(title))
  sub(
    "^(TPA_asm:|PREDICTED:)[[:space:]]*",
    "",
    title,
    ignore.case = TRUE
  )
}

# Keep the leading scientific-name label for reporting; a qualifier does not
# erase a useful genus. For example: Entomobrya cf. lanuginosa 14FR7_1.
extract_display_label <- function(title) {
  title <- remove_informative_record_prefix(title)
  matched <- stringr::str_extract(
    title,
    "^[A-Z][A-Za-z.-]+[[:space:]]+(?:(?:cf|aff|nr|sp|spp)\\.[[:space:]]+)?[a-z][A-Za-z0-9._-]*"
  )
  ifelse(is.na(matched), title, matched)
}

# Rank-aware parsing: uncertain species labels (cf., aff., nr., sp., spp.)
# are never used as confirmed species, but their genus can still be used.
# Avoid deriving genera from 'uncultured', 'unidentified', etc.
# Family-level labels in GenBank titles (e.g. "Aphididae sp. voucher")
# must never be interpreted as genera or confirmed species.
extract_family_label <- function(title) {
  title <- remove_informative_record_prefix(title)
  first <- stringr::str_extract(title, "^[A-Z][A-Za-z.-]*")
  first <- normalise_taxon(first)
  is_family <- !is.na(first) & grepl("idae$", first, ignore.case = TRUE)
  first[!is_family] <- NA_character_
  first
}

extract_genus <- function(title) {
  title <- remove_informative_record_prefix(title)
  genus <- stringr::str_extract(title, "^[A-Z][a-z][A-Za-z.-]*")
  genus <- normalise_taxon(genus)
  invalid <- c("Uncultured", "Unidentified", "Unknown", "Environmental",
               "Metagenome", "Unclassified", "Clone", "Bacterium",
               "Arthropod", "Insect", "Mitochondrion")
  is_family <- !is.na(genus) & grepl("idae$", genus, ignore.case = TRUE)
  genus[!is.na(genus) & (tolower(genus) %in% tolower(invalid) | is_family)] <- NA_character_
  genus
}

extract_species <- function(title) {
  title <- remove_informative_record_prefix(title)
  # Only conventional, unqualified binomials count as species-level evidence.
  # Do not silently promote 'cf. lanuginosa' to 'lanuginosa'.
  qualified <- grepl(
    "^[A-Z][a-z][A-Za-z.-]*[[:space:]]+(cf|aff|nr|sp|spp)\\.",
    title, ignore.case = TRUE, perl = TRUE
  )
  species <- stringr::str_extract(
    title,
    "^[A-Z][a-z][A-Za-z.-]*(?:[[:space:]]+(?:x|×))?[[:space:]]+[a-z][A-Za-z.-]+"
  )
  species <- normalise_taxon(species)
  species[qualified | !is.na(extract_family_label(title))] <- NA_character_
  species
}

noninformative_pattern <- paste(
  c("\\bsp\\.?\\b", "\\bcf\\.?\\b", "\\baff\\.?\\b",
    "\\bnr\\.?\\b", "uncultured", "unidentified", "environmental",
    "metagenome", "unknown", "unclassified", "clone"), collapse = "|"
)
is_informative_species <- function(species) {
  !is.na(species) &
    grepl("^[A-Z][a-z][A-Za-z.-]*(?:[[:space:]]+(?:x|×))?[[:space:]]+[a-z][A-Za-z.-]+$",
          species, perl = TRUE) &
    !grepl(noninformative_pattern, species, ignore.case = TRUE, perl = TRUE)
}
genus_from_species <- function(species) {
  ifelse(is.na(species), NA_character_, sub("[[:space:]].*$", "", species))
}

# ------------------------------ 3. READ MERGED TABLE --------------------------

if (!file.exists(lulu_basta_file)) stop("Merged file not found: ", lulu_basta_file)

merged <- read.csv(
  lulu_basta_file,
  sep = ",",
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

positive_control_columns <- names(merged)[grepl(
  positive_control_pattern, names(merged), ignore.case = TRUE, perl = TRUE
)]
negative_control_columns <- names(merged)[grepl(
  negative_control_pattern, names(merged), ignore.case = TRUE, perl = TRUE
)]

message(
  "Positive-control columns detected: ",
  if (length(positive_control_columns)) {
    paste(positive_control_columns, collapse = ", ")
  } else {
    "none"
  }
)
message(
  "Negative-control columns detected: ",
  if (length(negative_control_columns)) {
    paste(negative_control_columns, collapse = ", ")
  } else {
    "none"
  }
)

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
    informative = is_informative_species(blast_species),
    blast_genus = extract_genus(stitle),
    blast_family_label = extract_family_label(stitle),
    genus_informative = !is.na(blast_genus)
  )

# --------------------- NCBI ACCESSION TAXONOMY VALIDATION ---------------------
# The validation_status column summarizes the consistency of the
# final taxonomic assignment with NCBI Taxonomy references.
#
# Possible values:
#
# lineage_consistent: At least one supporting NCBI lineage is compatible with all
#   populated final taxonomic ranks.
#
# lineage_inconsistent: No supporting NCBI lineage is compatible with all populated
#   final taxonomic ranks. Manual review is recommended.
#
# lineage_disagreement: The final taxon disagrees with the NCBI assignment at the
#   corresponding taxonomic rank. Manual review is recommended.
#
# rank_absent_in_ncbi: The required taxonomic rank is unavailable in the retrieved
#   NCBI references. Validation is incomplete.
#
# no_verified_accession: No NCBI-verified BLAST accession is available for the ASV.
#   Validation is incomplete.
#
# unassigned: The ASV has no valid final taxonomic assignment to evaluate.
#
# IMPORTANT:
#   Validation checks taxonomic consistency, not the accuracy of
#   the molecular identification itself.
#   Validation results do not modify final taxonomic assignments.
#
# ================================================================

if (!requireNamespace("rentrez", quietly = TRUE) ||
    !requireNamespace("xml2", quietly = TRUE)) {
  stop("Install 'rentrez' and 'xml2' for accession-based taxonomy validation.")
}
ncbi_cache_file <- sub("[.]csv$", "_ncbi_accession_cache.rds", output_file)
ncbi_audit_file <- sub("[.]csv$", "_ncbi_accession_audit.csv", output_file)
ncbi_final_audit_file <- sub("[.]csv$", "_ncbi_final_audit.csv", output_file)
ncbi_cache <- if (file.exists(ncbi_cache_file)) readRDS(ncbi_cache_file) else list()
ncbi_session_cache <- new.env(parent = emptyenv())
ncbi_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
empty_lineage <- function() setNames(rep(NA_character_, length(ncbi_ranks)), ncbi_ranks)
xml_value <- function(node, path) {
  x <- xml2::xml_find_first(node, path)
  if (inherits(x, "xml_missing")) NA_character_ else xml2::xml_text(x)
}
get_ncbi_accession <- function(accession) {
  if (is.na(accession) || !nzchar(accession)) {
    return(list(status = "invalid_accession", taxid = NA_character_,
                lineage = empty_lineage(), detail = "Unrecognized accession format"))
  }
  if (exists(accession, envir = ncbi_session_cache, inherits = FALSE)) {
    return(get(accession, envir = ncbi_session_cache, inherits = FALSE))
  }
  if (!is.null(ncbi_cache[[accession]])) return(ncbi_cache[[accession]])
  result <- list(status = "query_error", taxid = NA_character_,
                 lineage = empty_lineage(), detail = NA_character_)
  result <- tryCatch({
    es <- rentrez::entrez_summary(db = "nuccore", id = accession)
    taxid <- as.character(es$taxid %||% NA_character_)
    if (is.na(taxid) || !nzchar(taxid)) {
      result$status <- "no_taxid"
      result
    } else {
      xml <- rentrez::entrez_fetch(db = "taxonomy", id = taxid,
                                   rettype = "xml", retmode = "xml")
      doc <- xml2::read_xml(xml)
      taxon <- xml2::xml_find_first(doc, ".//Taxon[TaxId]")
      if (inherits(taxon, "xml_missing")) stop("No taxonomy XML taxon")
      nodes <- c(list(taxon), as.list(xml2::xml_find_all(taxon, "./LineageEx/Taxon")))
      lin <- empty_lineage()
      for (node in nodes) {
        r <- tolower(xml_value(node, "./Rank"))
        nm <- xml_value(node, "./ScientificName")
        if (!is.na(r) && r %in% ncbi_ranks && !is.na(nm)) lin[r] <- nm
      }
      list(status = "matched", taxid = taxid, lineage = lin, detail = NA_character_)
    }
  }, error = function(e) {
    result$detail <- conditionMessage(e)
    result
  })
  # Cache every result during this run to avoid repeated failed requests.
  # Persist successful lookups only; transient errors remain retryable next run.
  assign(accession, result, envir = ncbi_session_cache)
  # Cache successful lookups only: transient errors must be retryable.
  if (result$status %in% c("matched", "no_taxid")) {
    ncbi_cache[[accession]] <<- result
    saveRDS(ncbi_cache, ncbi_cache_file)
  }
  result
}
get_accession <- function(x) {
  x <- trimws(as.character(x))
  # Accept BLAST identifiers such as gi|...|tpg|BK012650.1|,
  # gi|...|gb|PQ885887.1|, ref|NC_029727.1| and bare accessions.
  tokens <- strsplit(x, "|", fixed = TRUE)[[1]]
  prefixes <- c("gb", "ref", "emb", "dbj", "tpg", "tpe", "tpd", "pdb")
  positions <- which(tolower(tokens) %in% prefixes)
  if (length(positions)) {
    idx <- positions[length(positions)] + 1L
    if (idx <= length(tokens) && nzchar(tokens[idx])) return(tokens[idx])
  }
  # An unrecognized pipe-delimited identifier is not a valid nuccore UID.
  if (grepl("|", x, fixed = TRUE)) return(NA_character_)
  if (nzchar(x)) x else NA_character_
}
accessions <- unique(na.omit(vapply(blast$sseqid, get_accession, character(1))))
message("Validating ", length(accessions), " distinct BLAST accessions with NCBI...")
for (j in seq_along(accessions)) {
  get_ncbi_accession(accessions[j])
  if (j %% 25L == 0L) message("NCBI accessions checked: ", j, "/", length(accessions))
}
blast$ncbi_accession <- vapply(blast$sseqid, get_accession, character(1))
blast$ncbi_taxid <- vapply(blast$ncbi_accession, function(a)
  get_ncbi_accession(a)$taxid, character(1))
blast$ncbi_status <- vapply(blast$ncbi_accession, function(a)
  get_ncbi_accession(a)$status, character(1))
for (rank in ncbi_ranks) {
  blast[[paste0("ncbi_", rank)]] <- vapply(blast$ncbi_accession,
    function(a) get_ncbi_accession(a)$lineage[[rank]], character(1))
}
name_equal <- function(a, b) !is.na(a) & !is.na(b) &
  tolower(trimws(a)) == tolower(trimws(b))
# Species requires exact title/NCBI binomial agreement; a synonym is flagged,
# never silently promoted. Genus requires exact agreement with NCBI genus.
blast$ncbi_species_valid <- name_equal(blast$blast_species, blast$ncbi_species)
blast$ncbi_genus_valid <- name_equal(blast$blast_genus, blast$ncbi_genus)
blast$ncbi_family_label_valid <- name_equal(blast$blast_family_label, blast$ncbi_family)
blast$blast_species[!blast$ncbi_species_valid] <- NA_character_
blast$blast_genus[!blast$ncbi_genus_valid] <- NA_character_
blast$blast_family_label[!blast$ncbi_family_label_valid] <- NA_character_
blast$informative <- is_informative_species(blast$blast_species)
blast$genus_informative <- !is.na(blast$blast_genus)
write.csv(blast[, c("seqid", "sseqid", "stitle", "ncbi_accession", "ncbi_taxid",
  "ncbi_status", "ncbi_species_valid", "ncbi_genus_valid", "ncbi_family_label_valid",
  paste0("ncbi_", ncbi_ranks))], ncbi_audit_file, row.names = FALSE, na = "NA")
message("NCBI BLAST audit: ", ncbi_audit_file)

# Keep each subject accession once per ASV/score/name. This prevents a literal
# duplicate line from inflating number_best_hits while retaining distinct accessions.
blast_unique <- blast %>%
  distinct(
    seqid, sseqid, pident, alignment_length, display_label,
    blast_species, .keep_all = TRUE
  )

# Restrict downstream BLAST/GBIF evidence to ASVs in the final LULU x BASTA
# dataset. Preserve the original `blast` and `blast_unique` objects unchanged.
# This keeps explore and curate on the same ASV population.
curation_asv_ids <- unique(merged$seqid[!is.na(merged$seqid) & nzchar(merged$seqid)])
blast_asv_ids <- unique(blast_unique$seqid[!is.na(blast_unique$seqid) &
                                             nzchar(blast_unique$seqid)])
message("BLAST ASVs (all): ", length(blast_asv_ids),
        " | LULU x BASTA ASVs: ", length(curation_asv_ids),
        " | Shared: ", length(intersect(blast_asv_ids, curation_asv_ids)),
        " | BLAST-only (excluded from curation): ",
        length(setdiff(blast_asv_ids, curation_asv_ids)))

blast_curation <- blast_unique %>%
  filter(seqid %in% curation_asv_ids)

# Audit BLAST records identified only at family level, across all ASVs.
family_label_audit <- blast_curation %>%
  filter(!is.na(blast_family_label)) %>%
  select(seqid, sseqid, pident, alignment_length, stitle,
         blast_family_label, blast_genus, blast_species) %>%
  arrange(seqid, desc(pident), desc(alignment_length))
family_label_audit_file <- sub("([.]csv)?$", "_family_label_audit.csv",
                               output_file, ignore.case = TRUE)
write.csv(family_label_audit, family_label_audit_file,
          row.names = FALSE, na = "NA")
message("Family-only BLAST title audit: ", nrow(family_label_audit),
        " hit(s) across ", dplyr::n_distinct(family_label_audit$seqid), " ASV(s).")


# Select the highest identity, then longest alignment among those highest hits.
top_hits <- blast_curation %>%
  group_by(seqid) %>%
  filter(!is.na(pident), !is.na(alignment_length)) %>%
  filter(pident >= max(pident) - best_pident_tolerance) %>%
  filter(
    alignment_length >= max(alignment_length, na.rm = TRUE) -
      best_alignment_tolerance
  ) %>%
  ungroup()

# One representative best hit per retained ASV for dataset diagnostics.
# A deterministic tie-break selects longer alignments, then accession ID.
# IMPORTANT: use the same LULU x BASTA ASVs as the final curation output.
blast_best_per_asv <- blast_curation %>%
  filter(!is.na(pident), !is.na(alignment_length)) %>%
  arrange(seqid, desc(pident), desc(alignment_length), sseqid) %>%
  group_by(seqid) %>% slice_head(n = 1) %>% ungroup()

blast_diagnostic_file <- sub("([.]csv)?$", "_blast_summary.csv",
                             output_file, ignore.case = TRUE)
blast_stat <- function(x, fun) if (!length(x)) NA_real_ else fun(x)
blast_identity <- blast_best_per_asv$pident
blast_lengths <- blast_best_per_asv$alignment_length
blast_dataset_summary <- tibble(
  Metric = c("Total ASVs", "ASVs with BLAST hits",
             "Mean alignment length (bp)", "Median alignment length (bp)",
             "Min alignment length (bp)", "Max alignment length (bp)",
             "Mean identity (%)", "Median identity (%)",
             "Min identity (%)", "Max identity (%)"),
  Value = c(length(curation_asv_ids),
            sum(curation_asv_ids %in% blast_best_per_asv$seqid),
            blast_stat(blast_lengths, mean), blast_stat(blast_lengths, median),
            blast_stat(blast_lengths, min), blast_stat(blast_lengths, max),
            blast_stat(blast_identity, mean), blast_stat(blast_identity, median),
            blast_stat(blast_identity, min), blast_stat(blast_identity, max))
)
if (analysis_mode == "explore") {
  message("\nBLAST dataset summary (one best hit per retained LULU x BASTA ASV):")
  print(blast_dataset_summary, n = Inf)
  write.csv(blast_dataset_summary, blast_diagnostic_file,
            row.names = FALSE, na = "NA")
  message("BLAST summary saved: ", blast_diagnostic_file)
}

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

# Return a rank only when all informative BASTA-vout lineages that contain that
# rank agree. This provides a GenBank-derived fallback without importing the
# alternative rank hierarchy used by GBIF.
vout_taxonomy_consensus <- function(id, rank) {
  if (!rank %in% taxonomy_ranks) return(NA_character_)
  rows <- basta_vout[basta_vout$seqid == as.character(id), , drop = FALSE]
  if (!nrow(rows)) return(NA_character_)
  values <- normalise_taxon(rows[[rank]])
  values <- unique(values[!is.na(values)])
  if (length(values) == 1) values[1] else NA_character_
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
  occurrences_GB_human_observation = integer(),
  occurrences_GB_machine_observation = integer(),
  occurrences_GB_observation = integer(),
  occurrences_GB_preserved_specimen = integer(),
  occurrences_GB_living_specimen = integer(),
  occurrences_GB_fossil_specimen = integer(),
  occurrences_GB_material_citation = integer(),
  occurrences_GB_material_sample = integer(),
  occurrences_GB_occurrence = integer(),
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
  gbif_integer_columns <- names(empty_gbif_cache)[
    vapply(empty_gbif_cache, is.integer, logical(1))
  ]
  for (column in gbif_integer_columns) {
    gbif_cache[[column]] <- suppressWarnings(as.integer(gbif_cache[[column]]))
  }
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
    occurrences_GB_human_observation = NA_integer_,
    occurrences_GB_machine_observation = NA_integer_,
    occurrences_GB_observation = NA_integer_,
    occurrences_GB_preserved_specimen = NA_integer_,
    occurrences_GB_living_specimen = NA_integer_,
    occurrences_GB_fossil_specimen = NA_integer_,
    occurrences_GB_material_citation = NA_integer_,
    occurrences_GB_material_sample = NA_integer_,
    occurrences_GB_occurrence = NA_integer_,
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
    # Retrieve every basisOfRecord count in one faceted query. This provides the
    # detailed internal audit columns without sending one request per category.
    basis_levels <- c(
      "HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OBSERVATION",
      "PRESERVED_SPECIMEN", "LIVING_SPECIMEN", "FOSSIL_SPECIMEN",
      "MATERIAL_CITATION", "MATERIAL_SAMPLE", "OCCURRENCE"
    )
    basis_counts <- setNames(rep(0L, length(basis_levels)), basis_levels)
    basis_result <- rgbif::occ_search(
      taxonKey = as.integer(key),
      country = country_code,
      occurrenceStatus = "PRESENT",
      facet = "basisOfRecord",
      facetLimit = 100,
      limit = 0
    )
    basis_facet <- basis_result$facets$basisOfRecord
    if (!is.null(basis_facet) && nrow(basis_facet)) {
      facet_names <- toupper(as.character(basis_facet$name))
      facet_counts <- suppressWarnings(as.integer(basis_facet$count))
      matched <- match(facet_names, names(basis_counts))
      usable <- !is.na(matched) & !is.na(facet_counts)
      basis_counts[matched[usable]] <- facet_counts[usable]
    }

    base$occurrences_GB_human_observation <- basis_counts["HUMAN_OBSERVATION"]
    base$occurrences_GB_machine_observation <- basis_counts["MACHINE_OBSERVATION"]
    base$occurrences_GB_observation <- basis_counts["OBSERVATION"]
    base$occurrences_GB_preserved_specimen <- basis_counts["PRESERVED_SPECIMEN"]
    base$occurrences_GB_living_specimen <- basis_counts["LIVING_SPECIMEN"]
    base$occurrences_GB_fossil_specimen <- basis_counts["FOSSIL_SPECIMEN"]
    base$occurrences_GB_material_citation <- basis_counts["MATERIAL_CITATION"]
    base$occurrences_GB_material_sample <- basis_counts["MATERIAL_SAMPLE"]
    base$occurrences_GB_occurrence <- basis_counts["OCCURRENCE"]

    observations <- sum(basis_counts[c(
      "HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OBSERVATION"
    )])
    specimens <- sum(basis_counts[c(
      "PRESERVED_SPECIMEN", "LIVING_SPECIMEN", "FOSSIL_SPECIMEN",
      "MATERIAL_CITATION", "MATERIAL_SAMPLE"
    )])
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
gbif_basis_count_columns <- c(
  "occurrences_GB_human_observation",
  "occurrences_GB_machine_observation",
  "occurrences_GB_observation",
  "occurrences_GB_preserved_specimen",
  "occurrences_GB_living_specimen",
  "occurrences_GB_fossil_specimen",
  "occurrences_GB_material_citation",
  "occurrences_GB_material_sample",
  "occurrences_GB_occurrence"
)
has_complete_basis_counts <- rowSums(is.na(
  gbif_cache[, gbif_basis_count_columns, drop = FALSE]
)) == 0

valid_cached_names <- gbif_cache$queried_name[
  !is.na(gbif_cache$gbif_records) &
    !is.na(gbif_cache$occurrences_GB_observations) &
    !is.na(gbif_cache$occurrences_GB_specimens) &
    !is.na(gbif_cache$occurrences_GB_all) &
    has_complete_basis_counts &
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

# ---------------------- OPTIONAL GBIF THRESHOLD EXPLORATION -------------------

if (analysis_mode == "explore") {

  # Keep only ASV-species candidates for which a GBIF observation count
  # was successfully obtained. Missing GBIF results (NA) are excluded.
  gbif_evaluable <- top_informative_gbif %>%
    filter(!is.na(occurrences_GB_observations))

  # Summarise GBIF observation support at the ASV level.
  # For ASVs with multiple informative BLAST species candidates, the highest
  # GBIF observation count among those candidates is retained.
  gbif_asv_summary <- gbif_evaluable %>%
    group_by(seqid) %>%
    summarise(
      max_gbif_observations = max(occurrences_GB_observations),
      .groups = "drop"
    )

  # Count ASVs with no GBIF observation support in the selected country.
  n_zero_observations <- gbif_asv_summary %>%
    filter(max_gbif_observations == 0) %>%
    nrow()

  # Count ASVs with at least one GBIF observation in the selected country.
  # This group is used as the baseline for threshold sensitivity analysis.
  n_with_observations <- gbif_asv_summary %>%
    filter(max_gbif_observations >= 1) %>%
    nrow()

  # Evaluate increasingly stringent GBIF observation thresholds.
  # Threshold comparison starts at >= 5, using >= 1 observation as the baseline.
  # Therefore, the first change represents the number of ASVs lost when moving
  # from >= 1 to >= 5 observations.
  gbif_threshold_exploration <- map_dfr(
    gbif_thresholds_to_test,
    function(threshold) {

      n_passing <- gbif_asv_summary %>%
        filter(max_gbif_observations >= threshold) %>%
        nrow()

      tibble(
        GBIF_threshold = paste0(">= ", threshold),
        ASVs_passing = n_passing
      )
    }
  ) %>%
    mutate(
      change_from_previous = ASVs_passing -
        lag(ASVs_passing, default = n_with_observations)
    )

  # Display the threshold sensitivity table.
  print(gbif_threshold_exploration)

  # Save the exploratory table in the same folder as the GBIF cache.
  gbif_threshold_exploration_file <- sub(
    "([.]csv)?$",
    "_threshold_exploration.csv",
    gbif_cache_file,
    ignore.case = TRUE
  )

  write.csv(
    gbif_threshold_exploration,
    gbif_threshold_exploration_file,
    row.names = FALSE,
    na = "NA"
  )

  # Add GBIF coverage information below the threshold table.
  # These values are provided as context and are not additional thresholds.
  cat(
    "\n",
    "GBIF coverage summary:\n",
    paste0("ASVs with no observations (= 0): ", n_zero_observations, "\n"),
    paste0("ASVs with observations (>= 1): ", n_with_observations, "\n"),
    "Note: Threshold comparison starts at >= 5, using >= 1 as the baseline.\n",
    file = gbif_threshold_exploration_file,
    append = TRUE,
    sep = ""
  )

  # Interpretation:
  # "change_from_previous" shows the number of ASVs gained or lost relative
  # to the preceding threshold. For the first threshold (>= 5), the comparison
  # is made against the >= 1 baseline.
  # This analysis is intended only to help select a reasonable GBIF threshold.

  message("")
  message("GBIF threshold exploration saved to:")
  message(normalizePath(gbif_threshold_exploration_file, mustWork = FALSE))
}

if (analysis_mode == "explore") {

  message("")
  message("Exploratory GBIF threshold analysis complete.")
  message("Choose a threshold, set analysis_mode = \"curate\",")
  message("set minimum_gbif_records to the selected value, and run the script again.")
  message("")

} else {
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

  # Obtain GBIF matches for occurrence filtering and accepted-name comparison.
  # These rows must not be used to populate the final BASTA/GenBank lineage.
  gbif_rows_for_asv <- function(id) {
    species <- top_hits$blast_species[top_hits$seqid == id & top_hits$informative]
    species <- unique(species[!is.na(species)])
    gbif_cache[match(species, gbif_cache$queried_name), , drop = FALSE]
  }

  # ------------------------ BLAST TAXONOMIC RESOLUTION --------------------------
  # Compare the strongest accession per taxon, then use other distinct accessions
  # only as secondary support. Identity and length must point to the same winner
  # when differences are within the configured identity margin.
  # Family is obtained from matching BASTA-vout genus paths or GBIF taxonomy;
  # it is NEVER guessed from a genus string.
  # Resolve higher ranks only from a lineage linked to the BLAST genus.
  # An unrelated BASTA consensus must not be transferred to a BLAST hit.
  higher_for_genus <- function(id, genus, rank) {
    if (is.na(genus) || !rank %in% c("Family", "Order")) return(NA_character_)
    rows <- basta_vout[basta_vout$seqid == as.character(id), , drop = FALSE]
    matches <- !is.na(rows$Genus) &
      tolower(normalise_taxon(rows$Genus)) == tolower(normalise_taxon(genus))
    values <- unique(normalise_taxon(rows[[rank]][matches]))
    values <- values[!is.na(values)]
    if (length(values) > 1L) return(NA_character_)
    if (length(values) == 1L) return(values[1])
    gb <- gbif_cache[!is.na(gbif_cache$genus) &
                       tolower(normalise_taxon(gbif_cache$genus)) ==
                       tolower(normalise_taxon(genus)), , drop = FALSE]
    values <- unique(normalise_taxon(gb[[tolower(rank)]]))
    values <- values[!is.na(values)]
    if (length(values) == 1L) values[1] else NA_character_
  }
  family_for_genus <- function(id, genus) {
    h <- blast_curation[blast_curation$seqid == id &
      !is.na(blast_curation$blast_genus) &
      tolower(blast_curation$blast_genus) == tolower(genus), , drop = FALSE]
    vals <- unique(na.omit(h$ncbi_family))
    if (length(vals) == 1L) vals[1] else NA_character_
  }
  order_for_genus <- function(id, genus) {
    h <- blast_curation[blast_curation$seqid == id &
      !is.na(blast_curation$blast_genus) &
      tolower(blast_curation$blast_genus) == tolower(genus), , drop = FALSE]
    vals <- unique(na.omit(h$ncbi_order))
    if (length(vals) == 1L) vals[1] else NA_character_
  }

  resolve_blast_rank <- function(id, rank) {
    cfg <- tolower(rank)
    identity_floor <- get(paste0(cfg, "_min_identity"))
    alignment_floor <- get(paste0(cfg, "_min_alignment"))
    pi_margin <- get(paste0(cfg, "_pident_margin"))
    len_margin <- get(paste0(cfg, "_alignment_margin"))
    hits <- blast_curation %>%
      filter(seqid == id, !is.na(pident),
             !is.na(alignment_length), pident >= identity_floor,
             alignment_length >= alignment_floor)
    if (!nrow(hits)) return(list(taxon = NA_character_, reason = paste0(
      rank, ": no BLAST hit passed ", identity_floor,
      "% identity and ", alignment_floor, " bp."),
      support = 0L, best = NULL, ambiguity = FALSE))
    hits$taxon <- if (rank == "Species") {
      ifelse(hits$informative, hits$blast_species, NA_character_)
    } else if (rank == "Genus") {
      hits$blast_genus
    } else if (rank == "Family") {
      ifelse(!is.na(hits$blast_genus) | !is.na(hits$blast_family_label),
             hits$ncbi_family, NA_character_)
    } else {
      ifelse(!is.na(hits$blast_genus) | !is.na(hits$blast_family_label),
             hits$ncbi_order, NA_character_)
    }
    hits <- hits[!is.na(hits$taxon), , drop = FALSE]
    if (!nrow(hits)) return(list(taxon = NA_character_, reason = paste0(
      rank, ": taxonomy unavailable for quality-filtered BLAST hits."),
      support = 0L, best = NULL, ambiguity = FALSE))
    # One best alignment per accession/taxon; no duplicated accession votes.
    hits <- hits %>% arrange(desc(pident), desc(alignment_length), sseqid) %>%
      distinct(taxon, sseqid, .keep_all = TRUE)
    groups <- split(hits, tolower(hits$taxon))
    candidates <- lapply(groups, function(g) {
      g <- g[order(-g$pident, -g$alignment_length), , drop = FALSE]
      list(taxon = g$taxon[1], pident = g$pident[1],
           alignment = g$alignment_length[1], accessions = n_distinct(g$sseqid))
    })
    candidates <- candidates[order(
      -vapply(candidates, `[[`, numeric(1), "pident"),
      -vapply(candidates, `[[`, numeric(1), "alignment"))]
    leader <- candidates[[1]]
    competing <- if (length(candidates) > 1L) candidates[-1] else list()
    ambiguity <- FALSE
    comparison <- character()
    for (challenger in competing) {
      dpi <- leader$pident - challenger$pident
      dlen <- leader$alignment - challenger$alignment
      # A challenger with a substantially longer alignment cannot be ignored
      # even if its identity is marginally lower. Conflicting signals => ambiguity.
      close_pi <- dpi <= pi_margin
      if (close_pi) {
        if (dlen >= len_margin && len_margin > 0) {
          comparison <- c(comparison, paste0(challenger$taxon,
                                             " (identity gap ", round(dpi, 3), " pp; length advantage +",
                                             round(dlen), " bp)"))
        } else if (abs(dlen) < len_margin &&
                   leader$accessions >= challenger$accessions + min_accession_advantage &&
                   leader$accessions >= 2L) {
          comparison <- c(comparison, paste0(challenger$taxon,
                                             " (near tie; accession support ", leader$accessions,
                                             " vs ", challenger$accessions, ")"))
        } else {
          ambiguity <- TRUE
          comparison <- c(comparison, paste0(challenger$taxon,
                                             " (unresolved: identity gap ", round(dpi, 3),
                                             " pp; length difference ", round(dlen), " bp)"))
        }
      } else if (dlen < -len_margin && len_margin > 0) {
        ambiguity <- TRUE
        comparison <- c(comparison, paste0(challenger$taxon,
                                           " (conflicting identity and alignment-length evidence)"))
      }
    }
    desc <- paste0(rank, " BLAST candidate ", leader$taxon, " (",
                   round(leader$pident, 3), "% identity; ", round(leader$alignment),
                   " bp; ", leader$accessions, " distinct accession(s))",
                   if (length(comparison)) paste0(". Compared with ",
                                                  paste(comparison, collapse = "; ")) else ". No other qualifying taxon",
                   if (ambiguity) "; taxonomic resolution insufficient." else "; resolved.")
    list(taxon = if (ambiguity) NA_character_ else leader$taxon,
         reason = desc, support = leader$accessions, best = leader,
         ambiguity = ambiguity)
  }

  # ------------------------------ 10. CURATION RULES ----------------------------
  curate_row <- function(i) {
    row <- evidence[i, , drop = FALSE]
    id <- as.character(row$seqid)
    basta_rank <- lowest_assigned_rank(row)
    basta_values <- setNames(vapply(taxonomy_ranks, function(r)
      normalise_taxon(row[[r]]), character(1)), taxonomy_ranks)
    final <- basta_values
    for (r in taxonomy_ranks) if (is.na(final[r])) {
      value <- vout_taxonomy_consensus(id, r)
      if (!is.na(value)) final[r] <- value
    }
    proposed <- lapply(c("Species", "Genus", "Family", "Order"), function(r)
      resolve_blast_rank(id, r))
    names(proposed) <- c("Species", "Genus", "Family", "Order")
    accepted <- names(proposed)[vapply(proposed, function(x)
      !is.na(x$taxon), logical(1))]
    manual <- FALSE
    if (length(accepted)) {
      automated_rank <- accepted[1]
      selected <- proposed[[automated_rank]]
      final[automated_rank] <- selected$taxon
      if (automated_rank == "Species") {
        final["Genus"] <- genus_from_species(selected$taxon)
        f <- family_for_genus(id, final["Genus"])
        if (!is.na(f)) {
          if (!is.na(basta_values["Family"]) &&
              tolower(basta_values["Family"]) != tolower(f)) {
            final[c("Order", "Class", "Phylum", "Kingdom")] <- NA_character_
          }
          final["Family"] <- f
        }
      }
      if (automated_rank == "Genus") {
        f <- family_for_genus(id, selected$taxon)
        if (!is.na(f)) {
          if (!is.na(basta_values["Family"]) &&
              tolower(basta_values["Family"]) != tolower(f)) {
            final[c("Order", "Class", "Phylum", "Kingdom")] <- NA_character_
          }
          final["Family"] <- f
        }
      }
      # For an assigned family/order, clear conflicting higher-rank BASTA
      # labels unless the linked BASTA lineage supports them. The output must
      # never present an unsupported combination of ranks as one lineage.
      if (automated_rank == "Family") {
        matched <- basta_vout[basta_vout$seqid == id &
                                !is.na(basta_vout$Family) &
                                tolower(normalise_taxon(basta_vout$Family)) ==
                                tolower(normalise_taxon(selected$taxon)), , drop = FALSE]
        for (r in c("Order", "Class", "Phylum", "Kingdom")) {
          vals <- unique(normalise_taxon(matched[[r]]))
          vals <- vals[!is.na(vals)]
          final[r] <- if (length(vals) == 1L) vals[1] else NA_character_
        }
      }
      if (automated_rank == "Order") {
        matched <- basta_vout[basta_vout$seqid == id &
                                !is.na(basta_vout$Order) &
                                tolower(normalise_taxon(basta_vout$Order)) ==
                                tolower(normalise_taxon(selected$taxon)), , drop = FALSE]
        for (r in c("Class", "Phylum", "Kingdom")) {
          vals <- unique(normalise_taxon(matched[[r]]))
          vals <- vals[!is.na(vals)]
          final[r] <- if (length(vals) == 1L) vals[1] else NA_character_
        }
      }
      # Never leave a contradictory BASTA lineage above the BLAST rank.
      if (automated_rank %in% c("Species", "Genus")) {
        f <- family_for_genus(id, final["Genus"])
        if (is.na(f) && !is.na(basta_values["Genus"]) &&
            tolower(basta_values["Genus"]) != tolower(final["Genus"])) {
          final[c("Family", "Order", "Class", "Phylum", "Kingdom")] <- NA_character_
        }
      }
      reason <- paste(vapply(proposed[seq_len(match(automated_rank,
                                                    names(proposed)))], `[[`, character(1), "reason"), collapse = " | ")
      reason <- paste0("Assigned ", tolower(automated_rank), ". ", reason)
      assigned_taxon <- normalise_taxon(basta_values[automated_rank])
      vout_taxon <- vout_taxonomy_consensus(id, automated_rank)
      if (!is.na(assigned_taxon) &&
          tolower(assigned_taxon) != tolower(selected$taxon)) {
        reason <- paste0(reason, " BASTA merged ", automated_rank,
                         " disagrees (", assigned_taxon,"); BLAST-supported assignment retained.")
      } else if (!is.na(assigned_taxon)) {
        reason <- paste0(reason, " BASTA merged assignment agrees.")
      } else if (!is.na(vout_taxon) &&
                 tolower(vout_taxon) == tolower(selected$taxon)) {
        reason <- paste0(reason, " BASTA vout supports the selected taxon.")
      } else {
        reason <- paste0(reason, " BASTA provides no concordant evidence at this rank.")
      }
      confidence <- if (selected$support >= 2L && automated_rank == "Species")
        "High" else "Moderate"
      evidence_basis <- "BLAST molecular resolution; BASTA audited"
    } else {
      # Retain only a rank supported by a coherent BASTA/vout lineage.
      consensus <- vapply(taxonomy_ranks, function(r)
        vout_taxonomy_consensus(id, r), character(1))
      names(consensus) <- taxonomy_ranks
      usable <- which(!is.na(consensus))
      if (length(usable)) {
        automated_rank <- taxonomy_ranks[max(usable)]
        final <- consensus
        manual <- TRUE
        reason <- paste0("BLAST could not resolve species, genus, family or order with the configured criteria. ",
                         paste(vapply(proposed, `[[`, character(1), "reason"), collapse = " | "),
                         " Retained unanimous BASTA-vout ", tolower(automated_rank),
                         " as provisional; manual review required.")
      } else {
        automated_rank <- basta_rank
        final <- basta_values
        manual <- TRUE
        reason <- paste0("BLAST did not resolve species, genus, family or order; BASTA-vout has no unanimous lineage. ",
                         paste(vapply(proposed, `[[`, character(1), "reason"), collapse = " | "),
                         " Retained original merged assignment provisionally; manual review required.")
      }
      confidence <- "Low"
      evidence_basis <- "Provisional BASTA; manual review required"
    }
    final_rank_index <- rank_index(automated_rank)
    if (!is.na(final_rank_index) && final_rank_index < length(taxonomy_ranks)) {
      final[(final_rank_index + 1):length(taxonomy_ranks)] <- NA_character_
    }
    # Rank-specific taxon names cannot be inherited from unrelated BASTA lineages.
    if (automated_rank == "Family" && !is.na(final["Family"]) &&
        !is.na(basta_values["Family"]) &&
        tolower(basta_values["Family"]) != tolower(final["Family"])) {
      final[c("Order", "Class", "Phylum", "Kingdom")] <- NA_character_
    }
    basta_rank_index <- rank_index(basta_rank)
    final_taxon <- if (!is.na(final_rank_index)) final[final_rank_index] else NA_character_
    curation_decision <- if (manual || is.na(final_rank_index)) {
      "Manual review"
    } else if (is.na(basta_rank_index)) {
      paste("Assigned", tolower(automated_rank))
    } else if (final_rank_index == basta_rank_index) {
      paste("Retained", tolower(automated_rank))
    } else if (final_rank_index < basta_rank_index) {
      paste("Downgraded to", tolower(automated_rank))
    } else paste("Upgraded to", tolower(automated_rank))

    # Orange marks scope, while manual_review independently records uncertainty.
    final_phylum <- normalise_taxon(final["Phylum"])
    outside_target <- !is.na(final_phylum) &&
      !tolower(final_phylum) %in% tolower(target_phyla)
    if (outside_target) reason <- paste0(reason, " Outside target phyla: ",
                                         final_phylum, " (targets: ", paste(target_phyla, collapse = ", "), ").")
    same_rank_agrees <- !is.na(basta_rank_index) &&
      basta_rank_index == final_rank_index &&
      !is.na(basta_values[automated_rank]) &&
      !is.na(final_taxon) &&
      tolower(basta_values[automated_rank]) == tolower(final_taxon)
    curation_colour <- if (outside_target) "Orange" else if (manual)
      "Red" else if (same_rank_agrees) "Green" else "Yellow"

    gbif_tax <- gbif_rows_for_asv(id)
    final_common_name <- NA_character_
    if (identical(automated_rank, "Species") && !is.na(final_taxon) && nrow(gbif_tax)) {
      common_match <- gbif_tax$common_name[
        !is.na(gbif_tax$queried_name) &
          tolower(gbif_tax$queried_name) == tolower(final_taxon)]
      common_match <- common_match[!is.na(common_match) & nzchar(common_match)]
      if (length(common_match)) final_common_name <- common_match[1]
    }
    gbif_status_text <- if ((row$number_species_checked %||% 0L) == 0L) {
      "No informative species available for GBIF check."
    } else if ((row$number_species_present %||% 0L) == 0L) {
      paste0("No candidate reached ", minimum_gbif_records,
             " GBIF observation records; geographic support is descriptive only.")
    } else paste0("GBIF candidates above observation threshold: ",
                  row$gbif_present_species, ".")
    hit_consensus <- if (length(accepted)) proposed[[accepted[1]]]$reason else
      "No molecular assignment resolved at species, genus, family or order."
    tibble(
      BASTA_rank = basta_rank,
      final_kingdom = final["Kingdom"], final_phylum = final["Phylum"],
      final_class = final["Class"], final_order = final["Order"],
      final_family = final["Family"], final_genus = final["Genus"],
      final_species = final["Species"], final_taxon = final_taxon,
      final_rank = automated_rank, common_name = final_common_name,
      hit_consensus = hit_consensus, GBIF_UK_status = gbif_status_text,
      curation_decision = curation_decision, manual_review = manual,
      curation_reason = reason, confidence = confidence,
      evidence_basis = evidence_basis, curation_colour = curation_colour,
      outside_target_phylum = outside_target,
      notes = control_notes_for_row(row)
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
  # Post-curation lineage audit: check all seven final_* ranks against
  # complete NCBI lineages supporting the final-rank taxon. Never edit results.
  # Other BLAST hits may belong to different taxa and must not be used to
  # declare a contradiction with the selected final assignment.
  normalize_audit_name <- function(x) {
    x <- as.character(x)
    x <- trimws(gsub("[[:space:]]+", " ", x))
    x[is.na(x) | !nzchar(x) | tolower(x) %in% c("na", "unassigned")] <- NA_character_
    tolower(x)
  }
  final_audit <- lapply(seq_len(nrow(result)), function(i) {
    id <- as.character(result$seqid[i])
    rank <- normalize_audit_name(result$final_rank[i])
    taxon <- normalize_audit_name(result$final_taxon[i])
    hits <- blast[as.character(blast$seqid) == id &
                    !is.na(blast$ncbi_status) & blast$ncbi_status == "matched", , drop = FALSE]
    final_values <- setNames(lapply(ncbi_ranks, function(r)
      normalize_audit_name(result[[paste0("final_", r)]][i])), ncbi_ranks)
    observed <- vapply(final_values, function(x) !is.na(x), logical(1))
    rank_status <- setNames(rep("not_assigned", length(ncbi_ranks)), ncbi_ranks)
    if (is.na(rank) || !(rank %in% ncbi_ranks) || is.na(taxon)) {
      status <- "unassigned"
    } else if (!nrow(hits)) {
      status <- "no_verified_accession"
      rank_status[observed] <- "no_verified_accession"
    } else {
      # Restrict to NCBI reference lineages that actually support the
      # selected final taxon at its declared rank.
      anchor <- normalize_audit_name(hits[[paste0("ncbi_", rank)]])
      supporting <- hits[!is.na(anchor) & anchor == taxon, , drop = FALSE]
      if (!nrow(supporting)) {
        candidate_names <- normalize_audit_name(hits[[paste0("ncbi_", rank)]])
        status <- if (all(is.na(candidate_names))) "rank_absent_in_ncbi" else "lineage_disagreement"
        rank_status[observed] <- "no_supporting_lineage"
      } else {
        # Compare each populated rank against the SAME supporting lineage.
        # A lineage can have absent NCBI ranks; absence is not disagreement.
        compatible <- rep(TRUE, nrow(supporting))
        for (r in ncbi_ranks) {
          if (!observed[[r]]) next
          ref <- normalize_audit_name(supporting[[paste0("ncbi_", r)]])
          compatible <- compatible & (is.na(ref) | ref == final_values[[r]])
          rank_status[[r]] <- if (any(!is.na(ref) & ref == final_values[[r]])) {
            "matched"
          } else if (all(is.na(ref))) {
            "rank_absent_in_ncbi"
          } else {
            "lineage_disagreement"
          }
        }
        # Even if each rank is individually supported by different accessions,
        # a full lineage must be consistent with at least one accession.
        status <- if (any(compatible)) "lineage_consistent" else "lineage_inconsistent"
      }
    }
    out <- data.frame(seqid = id, final_rank = rank, final_taxon = taxon,
                      validation_status = status,
                      supporting_accessions = if (exists("supporting", inherits = FALSE))
                        length(unique(supporting$ncbi_accession)) else 0L,
                      stringsAsFactors = FALSE)
    for (r in ncbi_ranks) {
      out[[paste0("final_", r)]] <- result[[paste0("final_", r)]][i]
      out[[paste0("audit_", r)]] <- rank_status[[r]]
    }
    out
  })
  final_audit <- dplyr::bind_rows(final_audit)
  write.csv(final_audit, ncbi_final_audit_file, row.names = FALSE, na = "NA")
  message("Full-lineage audit: ", ncbi_final_audit_file)
  print(table(final_audit$validation_status, useNA = "ifany"))
  message("Audit is diagnostic only; final assignments are unchanged.")

  # Add only the NCBI validation status immediately after final_taxon.
  audit_idx <- match(as.character(result$seqid), as.character(final_audit$seqid))
  if (anyNA(audit_idx) || anyDuplicated(final_audit$seqid)) {
    stop("Cannot safely join the NCBI audit to curation results: missing or duplicate seqid.")
  }
  result$validation_status <- final_audit$validation_status[audit_idx]
  result$validation_flag <- NULL
  result <- dplyr::relocate(result, validation_status, .after = final_taxon)
  message("NCBI validation statuses in main curation CSV:")
  print(table(result$validation_status, useNA = "ifany"))
  write.csv(result, output_file, row.names = FALSE, na = "NA")

  message("Curation complete.")
  message("Output: ", normalizePath(output_file, mustWork = FALSE))
  message("GBIF cache: ", normalizePath(gbif_cache_file, mustWork = FALSE))
  message("Internal reusable GBIF cache: ", normalizePath(gbif_internal_cache_file, mustWork = FALSE))
  message("NBN common-name cache: ", normalizePath(nbn_common_name_cache_file, mustWork = FALSE))
  message("Resolved ranks: ", paste(names(table(result$final_rank)),
                                    as.integer(table(result$final_rank)), collapse = "; "))
  message(
    "Manual-review ASVs: ",
    sum(result$manual_review %in% TRUE, na.rm = TRUE),
    " of ", nrow(result)
  )

}  # End analysis_mode == "curate"
