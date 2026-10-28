# ============================================================
# COMPLETE MISSING FINAL TAXONOMIC RANKS USING NCBI TAXONOMY
# Independent post-curation step; never overwrites the input.
# Preserves validation_status and all existing final_* values.
# NCBI may not provide a formal kingdom rank for many animals;
# missing kingdom values are not guessed or replaced with superkingdom.
# Run from the RStudio project root.
# ============================================================

# Input and output file paths
input_file  <- "path/to/curated_ASV_table.csv"
output_file <- "path/to/completed_ASV_table.csv"
audit_file  <- "path/to/NCBI_completion_audit.csv"
cache_file  <- "path/to/NCBI_completion_cache.rds"

# NCBI request settings
# Delay (in seconds) between NCBI requests.
request_pause_seconds <- 0.40

# Target taxonomic group(s)
# Assignments outside these phyla are flagged as Orange, not necessarily incorrect.
# Use the same target phyla configured in the taxonomic curation script.
target_phyla <- c("Arthropoda")

required_packages <- c("rentrez", "xml2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install packages: install.packages(c(", paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
if (normalizePath(input_file, mustWork = TRUE) == normalizePath(output_file, mustWork = FALSE)) stop("Output must differ from input.")

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
final_cols <- paste0("final_", ranks)
clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in% c("", "na", "n/a", "none", "unknown", "unassigned")] <- NA_character_
  x
}
canonical <- function(x) tolower(gsub("[[:space:]]+", " ", clean(x)))

# Keep original column names and all original values unchanged.
df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE,
               fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
if (!"seqid" %in% names(df)) stop("Missing seqid column.")
if (!"final_rank" %in% names(df)) stop("Missing final_rank column.")
for (col in final_cols) if (!col %in% names(df)) df[[col]] <- NA_character_
if (!"final_taxon" %in% names(df)) stop("Missing final_taxon column.")
if (anyDuplicated(as.character(df$seqid))) stop("Duplicate seqid values detected.")

# Reuse successful and unsuccessful NCBI lookups on subsequent runs.
cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
if (!is.list(cache)) stop("Invalid NCBI cache file.")

xml_text <- function(node, path) {
  value <- xml2::xml_find_first(node, path)
  if (inherits(value, "xml_missing")) return(NA_character_)
  clean(xml2::xml_text(value))
}

query_ncbi <- function(name, expected_rank) {
  key <- paste(expected_rank, canonical(name), sep = "::")
  if (!is.null(cache[[key]])) return(cache[[key]])
  answer <- list(status = "not_found", taxid = NA_character_, lineage = setNames(rep(NA_character_, length(ranks)), ranks), detail = "")
  tryCatch({
    # Scientific Name field prevents broad, fuzzy matches.
    term <- paste0('"', gsub('"', '', name, fixed = TRUE), '"[Scientific Name]')
    search <- rentrez::entrez_search(db = "taxonomy", term = term, retmax = 25)
    ids <- unique(as.character(search$ids))
    if (!length(ids)) {
      answer$detail <- "No scientific-name match in NCBI Taxonomy"
    } else if (as.numeric(search$count) > 25) {
      answer$status <- "ambiguous"
      answer$detail <- "Too many NCBI matches; not resolved automatically"
    } else {
      fetched <- rentrez::entrez_fetch(db = "taxonomy", id = ids, rettype = "xml", retmode = "xml")
      doc <- xml2::read_xml(fetched)
      taxa <- xml2::xml_find_all(doc, "/TaxaSet/Taxon")
      exact <- taxa[vapply(taxa, function(t) {
        sci <- xml_text(t, "./ScientificName")
        rank <- xml_text(t, "./Rank")
        !is.na(sci) && !is.na(rank) && identical(canonical(sci), canonical(name)) &&
          identical(tolower(rank), expected_rank)
      }, logical(1))]
      if (length(exact) != 1L) {
        answer$status <- if (length(exact) > 1L) "ambiguous" else "not_found"
        answer$detail <- paste("Exact scientific-name/rank matches:", length(exact))
      } else {
        taxon <- exact[[1]]
        lineage_nodes <- xml2::xml_find_all(taxon, "./LineageEx/Taxon")
        nodes <- c(as.list(lineage_nodes), list(taxon))
        lineage <- answer$lineage
        for (node in nodes) {
          rank <- tolower(xml_text(node, "./Rank"))
          sci <- xml_text(node, "./ScientificName")
          if (!is.na(rank) && rank %in% ranks && !is.na(sci)) lineage[[rank]] <- sci
        }
        answer$status <- "matched"
        answer$taxid <- xml_text(taxon, "./TaxId")
        answer$lineage <- lineage
        answer$detail <- "Exact NCBI scientific-name and rank match"
      }
    }
  }, error = function(e) {
    answer$status <<- "query_error"
    answer$detail <<- conditionMessage(e)
  })
  # Do not cache temporary network errors, so they can be retried.
  if (answer$status != "query_error") {
    cache[[key]] <<- answer
    saveRDS(cache, cache_file)
  }
  if (request_pause_seconds > 0) Sys.sleep(request_pause_seconds)
  answer
}

# Only use the taxon explicitly selected by the curation pipeline.
# A lower-rank name must never be inferred from a higher-rank name.
choose_anchor <- function(row) {
  assigned_rank <- tolower(clean(row[["final_rank"]]))
  if (is.na(assigned_rank) || !assigned_rank %in% ranks) return(NULL)
  anchor_name <- clean(row[[paste0("final_", assigned_rank)]])
  chosen_name <- clean(row[["final_taxon"]])
  if (is.na(anchor_name) || is.na(chosen_name) ||
      canonical(anchor_name) != canonical(chosen_name)) return(NULL)
  # Exclude provisional species names and non-binomial species labels.
  if (assigned_rank == "species" &&
      (!grepl("^[A-Z][a-zA-Z.-]+ [a-z][a-zA-Z.-]+$", anchor_name) ||
       grepl("\\b(cf|aff|nr|sp|spp)\\.", anchor_name, ignore.case = TRUE))) return(NULL)
  list(rank = assigned_rank, name = anchor_name)
}

missing_before <- vapply(final_cols, function(col) sum(is.na(clean(df[[col]]))), integer(1))
audit <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  row <- df[i, , drop = FALSE]
  anchor <- choose_anchor(row)
  old <- vapply(final_cols, function(col) clean(row[[col]]), character(1))
  names(old) <- ranks
  status <- "skipped"
  detail <- "No valid final-rank/final-taxon anchor"
  taxid <- NA_character_
  filled <- character()
  if (!is.null(anchor)) {
    anchor_index <- match(anchor$rank, ranks)
    higher <- if (anchor_index > 1L) ranks[seq_len(anchor_index - 1L)] else character()
    missing_higher <- higher[is.na(old[higher])]
    if (!length(missing_higher)) {
      status <- "complete"
      detail <- "No missing higher ranks"
    } else {
      response <- query_ncbi(anchor$name, anchor$rank)
      status <- response$status
      detail <- response$detail
      taxid <- response$taxid
      if (identical(status, "matched")) {
        ncbi <- response$lineage
        # Require all already populated ranks (including the anchor) to agree.
        # If any conflict exists, leave the entire ASV unchanged.
        known <- ranks[!is.na(old)]
        conflicts <- known[!is.na(ncbi[known]) & canonical(old[known]) != canonical(ncbi[known])]
        if (length(conflicts)) {
          status <- "conflict"
          detail <- paste0("Existing vs NCBI disagreement at: ", paste(conflicts, collapse = ", "))
        } else {
          filled <- missing_higher[!is.na(ncbi[missing_higher])]
          for (rank in filled) df[[paste0("final_", rank)]][i] <- ncbi[[rank]]
          status <- if (length(filled)) "filled" else "no_ncbi_higher_ranks"
          detail <- if (length(filled)) paste("Filled:", paste(filled, collapse = ", ")) else "NCBI has no requested higher ranks"
        }
      }
    }
  }
  audit[[i]] <- data.frame(seqid = as.character(row$seqid),
                           final_rank = as.character(row$final_rank),
                           final_taxon = as.character(row$final_taxon),
                           ncbi_taxid = taxid, status = status,
                           fields_filled = paste(filled, collapse = "; "),
                           detail = detail, stringsAsFactors = FALSE)
  if (i %% 25L == 0L || i == nrow(df)) message("[PROGRESS] ", i, "/", nrow(df), " ASVs")
}

# Reconcile curation colours AFTER NCBI rank completion.
# Orange has absolute priority when the final phylum is known to be out of scope.
# Preserve the original colour for traceability; never promote Red/Yellow to Green.
if (!"curation_colour" %in% names(df)) stop("Missing curation_colour column.")
if (!"outside_target_phylum" %in% names(df)) df$outside_target_phylum <- NA
if (!"curation_colour_before_ncbi" %in% names(df)) {
  df$curation_colour_before_ncbi <- df$curation_colour
}
phylum <- canonical(df$final_phylum)
target <- canonical(target_phyla)
known_phylum <- !is.na(phylum)
outside <- known_phylum & !(phylum %in% target)
df$outside_target_phylum <- ifelse(known_phylum, outside, NA)
old_colours <- as.character(df$curation_colour)
df$curation_colour[outside] <- "Orange"

# A previously Orange ASV is not automatically promoted if the final phylum
# is missing or has changed; such cases require explicit review.
colour_changed <- which(!is.na(old_colours) & old_colours != df$curation_colour)
colour_audit <- data.frame(
  seqid = as.character(df$seqid),
  final_phylum = df$final_phylum,
  final_class = df$final_class,
  colour_before = old_colours,
  colour_after = df$curation_colour,
  outside_target_phylum = df$outside_target_phylum,
  changed = seq_len(nrow(df)) %in% colour_changed,
  stringsAsFactors = FALSE
)
colour_audit_file <- sub("[.]csv$", "_colour_audit.csv", output_file)
write.csv(colour_audit, colour_audit_file, row.names = FALSE, na = "NA")
if (any(outside & df$curation_colour != "Orange")) {
  stop("Colour consistency check failed: out-of-scope ASV is not Orange.")
}
message("[COLOUR] Out-of-scope ASVs enforced as Orange: ", sum(outside))
message("[COLOUR] Colours changed: ", length(colour_changed))
message("[COLOUR AUDIT] ", colour_audit_file)

audit_df <- do.call(rbind, audit)
missing_after <- vapply(final_cols, function(col) sum(is.na(clean(df[[col]]))), integer(1))
summary_df <- data.frame(rank = ranks, missing_before = as.integer(missing_before),
                         missing_after = as.integer(missing_after),
                         values_filled = as.integer(missing_before - missing_after))
summary_file <- "LgV2/LgV2_ncbi_completion_summary.csv"
write.csv(summary_df, summary_file, row.names = FALSE)
write.csv(df, output_file, row.names = FALSE, na = "NA")
write.csv(audit_df, audit_file, row.names = FALSE, na = "NA")
message("\n[DONE] Original curation file preserved: ", input_file)
message("[OUTPUT] ", output_file)
message("[AUDIT] ", audit_file)
message("[CACHE] ", cache_file)
message("[SUMMARY] ", summary_file)
print(table(audit_df$status))
print(summary_df, row.names = FALSE)
