# ============================================================
# MARKER COMPARISON - BASIC DESCRIPTIVE STATISTICS
# Minimum abundance: 10 reads per ASV across all samples
# Run from the RStudio project root:
# "Testing primers on subset samples"
# ============================================================

options(stringsAsFactors = FALSE)
min_total_reads <- 10L
marker_folders <- c("16s", "LgV1reamp", "LgV2", "LgV2reamp", "ShV1reamp")
output_dir <- "Marker_comparison"
figures_dir <- file.path(output_dir, "Figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Check that the working directory is the project root.
missing_folders <- marker_folders[!dir.exists(marker_folders)]
if (length(missing_folders) > 0L) {
  stop("Marker folders not found: ", paste(missing_folders, collapse = ", "),
       "\nOpen the RStudio project and run this script from its root directory.")
}

# Read a CSV, preserving original column names and detecting common delimiters.
read_curation_csv <- function(path) {
  header <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  delimiters <- c(",", ";", "\t")
  counts <- vapply(delimiters, function(sep) {
    length(strsplit(header, sep, fixed = TRUE)[[1L]])
  }, integer(1L))
  sep <- delimiters[which.max(counts)]
  df <- read.csv(path, sep = sep, check.names = FALSE,
                 stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  required <- c("seqid", "Kingdom", "final_rank", "final_taxon",
                "final_kingdom", "final_phylum", "final_class", "final_order",
                "final_family", "final_genus", "final_species",
                "curation_colour", "validation_status", "manual_review",
                "outside_target_phylum")
  if (!all(required %in% names(df))) {
    stop("Missing required columns in ", path, ": ",
         paste(setdiff(required, names(df)), collapse = ", "))
  }
  df
}

# Sample columns in the curation output occur between seqid and Kingdom.
# This excludes the row index, seqid, and all taxonomy/evidence columns.
get_sample_columns <- function(df, path) {
  first <- match("seqid", names(df)) + 1L
  last <- match("Kingdom", names(df)) - 1L
  if (first > last) stop("No sample columns found in: ", path)
  cols <- names(df)[seq.int(first, last)]
  if (length(cols) == 0L) stop("No sample columns found in: ", path)
  cols
}

# Convert numeric sample counts safely; missing counts are treated as zero.
get_read_totals <- function(df, sample_cols, path) {
  sample_data <- lapply(df[sample_cols], function(x) {
    y <- suppressWarnings(as.numeric(as.character(x)))
    invalid <- is.na(y) & !is.na(x) & nzchar(trimws(as.character(x)))
    if (any(invalid)) {
      stop("Non-numeric sample read counts in ", path)
    }
    if (any(y < 0, na.rm = TRUE)) stop("Negative read counts in ", path)
    y[is.na(y)] <- 0
    y
  })
  rowSums(as.data.frame(sample_data), na.rm = TRUE)
}

clean_taxon <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in% c("", "na", "n/a", "none", "unknown", "unassigned")] <- NA_character_
  x
}

rank_levels <- c("Species", "Genus", "Family", "Order", "Class", "Unassigned", "Other")
all_summaries <- list()
all_ranks <- list()
all_filtered <- list()
all_orders <- list()
all_exclusions <- list()
all_categories <- list()
all_validation <- list()
all_richness <- list()
all_class_composition <- list()
all_insect_orders <- list()

for (marker in marker_folders) {
  message("[READ] ", marker)
  folder <- marker
  # Select the exact final curation CSV for this marker.
  # Do not select auxiliary summaries or previous curation versions.
  expected_file <- paste0(marker, "_automated_curation_v5g_ncbi_validated_completed.csv")
  selected_file <- file.path(folder, expected_file)
  if (!file.exists(selected_file)) {
    stop("Expected final curation file not found: ", selected_file,
         "\nCheck the file name and its marker folder. ",
         "The script requires the exact *_automated_curation_v5g_ncbi_validated_completed.csv file.")
  }
  message("  File: ", selected_file)
  df <- read_curation_csv(selected_file)
  if (anyDuplicated(as.character(df$seqid))) {
    stop("Duplicate seqid values in: ", selected_file)
  }
  sample_cols <- get_sample_columns(df, selected_file)
  df$total_reads <- get_read_totals(df, sample_cols, selected_file)
  original_n <- nrow(df)
  df <- df[df$total_reads >= min_total_reads, , drop = FALSE]
  retained_n <- nrow(df)

  # Record curation categories BEFORE excluding Red and Orange ASVs.
  # NCBI validation is a separate diagnostic, not an additional acceptance rule.
  colour <- tolower(trimws(as.character(df$curation_colour)))
  colour[is.na(colour) | !nzchar(colour)] <- "missing"
  allowed_colours <- c("green", "yellow", "red", "orange")
  if (any(!colour %in% allowed_colours)) {
    stop("Unexpected or missing curation_colour in ", selected_file, ": ",
         paste(unique(colour[!colour %in% allowed_colours]), collapse = ", "))
  }
  category_levels <- c("Green", "Yellow", "Red - manual review", "Orange - outside target")
  categories <- c(green = category_levels[1], yellow = category_levels[2],
                  red = category_levels[3], orange = category_levels[4])
  category <- unname(categories[colour])
  category_counts <- as.integer(table(factor(category, levels = category_levels)))
  all_categories[[marker]] <- data.frame(
    Marker = marker, Category = category_levels, ASVs = category_counts,
    Percentage = if (retained_n > 0L) round(100 * category_counts / retained_n, 1) else NA_real_,
    stringsAsFactors = FALSE
  )
  validation <- clean_taxon(df$validation_status)
  validation[is.na(validation)] <- "Missing"
  validation_levels <- sort(unique(c("lineage_consistent", "lineage_inconsistent",
                                     "lineage_disagreement", "rank_absent_in_ncbi",
                                     "no_verified_accession", "unassigned",
                                     "Missing", validation)))
  validation_table <- table(factor(validation, levels = validation_levels),
                            factor(category, levels = category_levels))
  all_validation[[marker]] <- data.frame(
    Marker = marker, validation_status = rep(validation_levels, times = length(category_levels)),
    Curation_category = rep(category_levels, each = length(validation_levels)),
    ASVs = as.integer(validation_table), stringsAsFactors = FALSE
  )
  all_exclusions[[marker]] <- data.frame(
    Marker = marker, seqid = df$seqid, total_reads = df$total_reads,
    curation_colour = df$curation_colour,
    inclusion_status = ifelse(colour %in% c("green", "yellow"),
                              "Included", ifelse(colour == "red", "Excluded - manual review",
                                                  "Excluded - outside target")),
    final_rank = df$final_rank, final_taxon = df$final_taxon,
    validation_status = df$validation_status,
    manual_review = df$manual_review,
    outside_target_phylum = df$outside_target_phylum,
    stringsAsFactors = FALSE
  )
  green_n <- sum(colour == "green")
  yellow_n <- sum(colour == "yellow")
  red_n <- sum(colour == "red")
  orange_n <- sum(colour == "orange")
  # The taxonomic denominator includes only Green and Yellow ASVs.
  # Never use validation_status alone to override a curation colour.
  df <- df[colour %in% c("green", "yellow"), , drop = FALSE]
  eligible_n <- nrow(df)

  ranks <- clean_taxon(df$final_rank)
  ranks <- ifelse(is.na(ranks), "Unassigned", ranks)
  ranks[!ranks %in% rank_levels] <- "Other"
  phylum <- clean_taxon(df$final_phylum)
  class <- if ("final_class" %in% names(df)) clean_taxon(df$final_class) else rep(NA_character_, nrow(df))
  if (!("final_class" %in% names(df))) {
    warning("final_class is absent in ", selected_file, "; class values will be missing.")
  }
  arthropoda <- !is.na(phylum) & tolower(phylum) == "arthropoda"
  missing_phylum <- is.na(phylum)
  n_species <- sum(ranks == "Species")
  n_genus <- sum(ranks == "Genus")
  n_family <- sum(ranks == "Family")
  n_order <- sum(ranks == "Order")
  n_class <- sum(ranks == "Class")
  n_unassigned <- sum(ranks == "Unassigned")
  pct <- function(n, total) if (total == 0L) NA_real_ else round(100 * n / total, 1)

  all_summaries[[marker]] <- data.frame(
    Marker = marker,
    ASVs_before_filter = original_n,
    ASVs_after_filter = retained_n,
    ASVs_removed = original_n - retained_n,
    Retained_pct = pct(retained_n, original_n),
    Eligible_ASVs_n = eligible_n,
    Eligible_pct_of_read_filtered = pct(eligible_n, retained_n),
    Green_n = green_n,
    Yellow_n = yellow_n,
    Red_excluded_n = red_n,
    Orange_excluded_n = orange_n,
    Excluded_total_n = red_n + orange_n,
    Excluded_pct_of_read_filtered = pct(red_n + orange_n, retained_n),
    Species_n = n_species,
    Species_pct = pct(n_species, eligible_n),
    Genus_n = n_genus,
    Genus_pct = pct(n_genus, eligible_n),
    Family_n = n_family,
    Family_pct = pct(n_family, eligible_n),
    Order_n = n_order,
    Order_pct = pct(n_order, eligible_n),
    Class_n = n_class,
    Class_pct = pct(n_class, eligible_n),
    Classes_distinct_n = length(unique(class[!is.na(class)])),
    Missing_final_class_n = sum(is.na(class)),
    Missing_final_class_pct = pct(sum(is.na(class)), eligible_n),
    Unassigned_n = n_unassigned,
    Unassigned_pct = pct(n_unassigned, eligible_n),
    Species_or_genus_n = n_species + n_genus,
    Species_or_genus_pct = pct(n_species + n_genus, eligible_n),
    Arthropoda_explicit_n = sum(arthropoda),
    Arthropoda_explicit_pct = pct(sum(arthropoda), eligible_n),
    Missing_final_phylum_n = sum(missing_phylum),
    Missing_final_phylum_pct = pct(sum(missing_phylum), eligible_n),
    stringsAsFactors = FALSE
  )

  all_ranks[[marker]] <- data.frame(
    Marker = marker,
    Rank = rank_levels,
    ASVs = as.integer(table(factor(ranks, levels = rank_levels))),
    Percentage = vapply(as.integer(table(factor(ranks, levels = rank_levels))),
                        pct, numeric(1L), total = eligible_n),
    stringsAsFactors = FALSE
  )
  all_filtered[[marker]] <- data.frame(
    Marker = marker, seqid = df$seqid, total_reads = df$total_reads,
    curation_colour = df$curation_colour,
    final_rank = df$final_rank, final_taxon = df$final_taxon,
    validation_status = df$validation_status,
    manual_review = df$manual_review,
    final_kingdom = df$final_kingdom, final_phylum = df$final_phylum,
    final_class = class,
    final_order = if ("final_order" %in% names(df)) df$final_order else NA_character_,
    final_family = if ("final_family" %in% names(df)) df$final_family else NA_character_,
    final_genus = if ("final_genus" %in% names(df)) df$final_genus else NA_character_,
    final_species = if ("final_species" %in% names(df)) df$final_species else NA_character_,
    stringsAsFactors = FALSE
  )
  # Order composition uses final_order, not final_rank: ASVs identified at
  # species, genus, or family still contribute to their resolved order.
  order_values <- clean_taxon(df$final_order)
  order_values[is.na(order_values)] <- "Unresolved"
  order_levels <- c(sort(unique(order_values[order_values != "Unresolved"])), "Unresolved")
  order_counts <- as.integer(table(factor(order_values, levels = order_levels)))
  all_orders[[marker]] <- data.frame(
    Marker = marker,
    Order = order_levels,
    ASVs = order_counts,
    Percentage = vapply(order_counts, pct, numeric(1L), total = eligible_n),
    stringsAsFactors = FALSE
  )

  # Count distinct taxon names at each rank, using all populated final_* fields.
  # A species-level ASV also contributes its verified higher-rank lineage.
  # Do not sum richness values across ranks: these are separate taxonomic units.
  richness_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (scope in c("Arthropoda", "Insecta")) {
    scope_mask <- if (scope == "Insecta") {
      !is.na(clean_taxon(df$final_class)) &
        tolower(clean_taxon(df$final_class)) == "insecta"
    } else {
      !is.na(clean_taxon(df$final_phylum)) &
        tolower(clean_taxon(df$final_phylum)) == "arthropoda"
    }
    for (rank in richness_ranks) {
      taxa <- clean_taxon(df[[paste0("final_", tolower(rank))]][scope_mask])
      all_richness[[length(all_richness) + 1L]] <- data.frame(
        Marker = marker, Scope = scope, Rank = rank,
        Distinct_taxa = length(unique(tolower(taxa[!is.na(taxa)]))),
        Eligible_ASVs_in_scope = sum(scope_mask), stringsAsFactors = FALSE
      )
    }
  }

  # Class composition: proportions are based on ASV counts, not read counts.
  class_values <- clean_taxon(df$final_class)
  class_values[is.na(class_values)] <- "Unresolved"
  class_counts <- table(class_values)
  all_class_composition[[marker]] <- data.frame(
    Marker = marker, Class = names(class_counts), ASVs = as.integer(class_counts),
    Percentage = if (eligible_n > 0L) 100 * as.integer(class_counts) / eligible_n else NA_real_,
    stringsAsFactors = FALSE
  )

  # Insect orders only; unresolved insect orders are counted separately.
  insect_mask <- !is.na(clean_taxon(df$final_class)) &
    tolower(clean_taxon(df$final_class)) == "insecta"
  insect_order <- clean_taxon(df$final_order[insect_mask])
  insect_order[is.na(insect_order)] <- "Unresolved"
  insect_counts <- table(insect_order)
  all_insect_orders[[marker]] <- data.frame(
    Marker = marker, Order = names(insect_counts), ASVs = as.integer(insect_counts),
    stringsAsFactors = FALSE
  )

  message("  Samples: ", length(sample_cols),
          " | ASVs: ", original_n, " -> ", retained_n,
           " | Eligible: ", eligible_n, " | Red: ", red_n, " | Orange: ", orange_n,
          " | Species: ", n_species,
          " | Explicit Arthropoda: ", sum(arthropoda))
}

summary_df <- do.call(rbind, all_summaries)
rank_df <- do.call(rbind, all_ranks)
filtered_df <- do.call(rbind, all_filtered)
order_df <- do.call(rbind, all_orders)
category_df <- do.call(rbind, all_categories)
exclusions_df <- do.call(rbind, all_exclusions)
validation_df <- do.call(rbind, all_validation)
rownames(summary_df) <- NULL
rownames(rank_df) <- NULL
rownames(filtered_df) <- NULL
rownames(order_df) <- NULL
rownames(category_df) <- NULL
rownames(exclusions_df) <- NULL
rownames(validation_df) <- NULL
richness_df <- do.call(rbind, all_richness)
class_composition_df <- do.call(rbind, all_class_composition)
insect_orders_df <- do.call(rbind, all_insect_orders)
rownames(richness_df) <- NULL
rownames(class_composition_df) <- NULL
rownames(insect_orders_df) <- NULL

write.csv(summary_df, file.path(output_dir, "01_marker_summary.csv"), row.names = FALSE, na = "")
write.csv(rank_df, file.path(output_dir, "02_taxonomic_resolution.csv"), row.names = FALSE, na = "")
write.csv(filtered_df, file.path(output_dir, "03_filtered_ASV_taxonomy.csv"), row.names = FALSE, na = "")
write.csv(order_df, file.path(output_dir, "04_order_distribution.csv"), row.names = FALSE, na = "")
write.csv(category_df, file.path(output_dir, "05_curation_categories.csv"), row.names = FALSE, na = "")
write.csv(exclusions_df, file.path(output_dir, "06_ASV_inclusion_audit.csv"), row.names = FALSE, na = "")
write.csv(validation_df, file.path(output_dir, "07_NCBI_validation_by_category.csv"), row.names = FALSE, na = "")
write.csv(richness_df, file.path(output_dir, "08_taxonomic_richness_by_rank.csv"), row.names = FALSE, na = "")
write.csv(class_composition_df, file.path(output_dir, "09_class_composition.csv"), row.names = FALSE, na = "")
write.csv(insect_orders_df, file.path(output_dir, "10_insect_order_coverage.csv"), row.names = FALSE, na = "")

# Figure 1: ASVs before and after the minimum-read filter.
png(file.path(figures_dir, "01_ASV_retention.png"), width = 1300, height = 850, res = 150)
par(mar = c(8, 5, 4, 1))
values <- rbind(summary_df$ASVs_before_filter, summary_df$ASVs_after_filter,
                summary_df$Eligible_ASVs_n)
bp <- barplot(values, beside = TRUE, names.arg = summary_df$Marker,
              col = c("grey75", "steelblue3", "#287C71"), las = 2,
              ylab = "Number of ASVs", main = "ASV retention and taxonomic eligibility",
              ylim = c(0, max(1, max(values)) * 1.18))
text(bp, values, labels = values, pos = 3, cex = 0.8)
legend("topright", legend = c("Before filtering", "At least 10 reads", "Green + Yellow"),
       fill = c("grey75", "steelblue3", "#287C71"), bty = "n")
dev.off()

# Figure 2: Taxonomic resolution with percentage labels inside each segment.
rank_matrix <- sapply(marker_folders, function(m) {
  z <- rank_df[rank_df$Marker == m, ]
  z$ASVs[match(rank_levels, z$Rank)]
})
rownames(rank_matrix) <- rank_levels
rank_colors <- c("#287C71", "#76B5A5", "#E7B85A", "#C78B59", "#8E91C3", "#B5B5B5", "#7C80A0")
rank_pct_matrix <- sweep(rank_matrix, 2, colSums(rank_matrix), FUN = "/") * 100
rank_pct_matrix[!is.finite(rank_pct_matrix)] <- 0

png(file.path(figures_dir, "02_taxonomic_resolution.png"),
    width = 1280, height = 900, res = 150)
par(mar = c(8, 5, 4, 10), xpd = FALSE)
bar_centers <- barplot(
  rank_pct_matrix,
  names.arg = marker_folders,
  col = rank_colors,
  las = 2,
  ylab = "Percentage of eligible ASVs",
  main = "Taxonomic resolution by marker (Green + Yellow)",
  ylim = c(0, 105),
  border = "white"
)

# Label segments large enough for readable text.
# Smaller segments remain represented in the legend and output CSV.
min_label_pct <- 3
for (j in seq_len(ncol(rank_pct_matrix))) {
  cumulative <- 0
  for (i in seq_len(nrow(rank_pct_matrix))) {
    segment_pct <- rank_pct_matrix[i, j]
    if (segment_pct >= min_label_pct) {
      label_y <- cumulative + segment_pct / 2
      label_color <- if (i == 1L) "white" else "#222222"
      text(bar_centers[j], label_y,
           labels = sprintf("%.1f%%", segment_pct),
           col = label_color, cex = 0.83, font = 2)
    }
    cumulative <- cumulative + segment_pct
  }
}
par(xpd = NA)
legend("topright", inset = c(-0.36, 0),
       legend = rank_levels, fill = rank_colors,
       bty = "n", cex = 0.9)
par(xpd = FALSE)
dev.off()


# Figure 3: Show excluded categories separately without contaminating resolution.
category_matrix <- sapply(marker_folders, function(m) {
  z <- category_df[category_df$Marker == m, , drop = FALSE]
  z$ASVs[match(category_levels, z$Category)]
})
rownames(category_matrix) <- category_levels
category_percent <- sweep(category_matrix, 2, colSums(category_matrix), FUN = "/") * 100
category_percent[!is.finite(category_percent)] <- 0
category_colors <- c("#287C71", "#E7B85A", "#B95D65", "#D48A45")
png(file.path(figures_dir, "03_curation_categories.png"),
    width = 1400, height = 900, res = 150)
par(mar = c(8, 5, 4, 12), xpd = FALSE)
category_centers <- barplot(category_percent, names.arg = marker_folders,
                            col = category_colors, las = 2, border = "white",
                            ylim = c(0, 105),
                            ylab = "Percentage of ASVs with >= 10 reads",
                            main = "Curation categories by marker")
for (j in seq_len(ncol(category_percent))) {
  cumulative <- 0
  for (i in seq_len(nrow(category_percent))) {
    value <- category_percent[i, j]
    if (value >= 3) {
      text(category_centers[j], cumulative + value / 2,
           labels = sprintf("%.1f%%", value), cex = 0.83, font = 2,
           col = if (i == 1L) "white" else "#222222")
    }
    cumulative <- cumulative + value
  }
}
par(xpd = NA)
legend("topright", inset = c(-0.48, 0), legend = category_levels,
       fill = category_colors, bty = "n", cex = 0.85)
par(xpd = FALSE)
dev.off()

# Figure 4: Distinct-taxon richness across ranks for Arthropoda and Insecta.
# Each rank counts unique taxa independently, not a cumulative sum across ranks.
marker_colors <- c("#287C71", "#C78B59", "#6A8DC1", "#9A76B5", "#D1A441")
richness_scopes <- list(
  Arthropoda = c("Class", "Order", "Family", "Genus", "Species"),
  Insecta = c("Order", "Family", "Genus", "Species")
)
png(file.path(figures_dir, "04_taxonomic_richness_by_rank.png"),
    width = 1900, height = 1050, res = 150)
# Reserve a dedicated strip for the shared legend; no legend overlaps the data.
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
       heights = c(5.5, 1.15))
for (scope in names(richness_scopes)) {
  ranks_to_plot <- richness_scopes[[scope]]
  z <- richness_df[richness_df$Scope == scope, , drop = FALSE]
  mat <- sapply(marker_folders, function(m) {
    v <- z[z$Marker == m, , drop = FALSE]
    v$Distinct_taxa[match(ranks_to_plot, v$Rank)]
  })
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(marker_folders))
  mat[is.na(mat)] <- 0
  par(mar = c(6.5, 5.5, 3.5, 1.2), mgp = c(3.2, 0.8, 0), las = 1)
  plot(seq_along(ranks_to_plot), mat[, 1], type = "n",
       xlim = c(0.85, length(ranks_to_plot) + 0.15),
       ylim = c(0, max(1, mat) * 1.12), xaxt = "n",
       xlab = "", ylab = "Distinct taxa", main = paste(scope, "richness"),
       cex.main = 1.25)
  axis(1, at = seq_along(ranks_to_plot), labels = ranks_to_plot,
       las = 1, cex.axis = 0.88, padj = 1.1)
  for (j in seq_along(marker_folders)) {
    lines(seq_along(ranks_to_plot), mat[, j], type = "b", pch = 16,
          col = marker_colors[j], lwd = 2.4, cex = 1.05)
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = marker_folders, col = marker_colors,
       lty = 1, lwd = 2.4, pch = 16, horiz = TRUE,
       cex = 1.0, bty = "n", x.intersp = 0.8)
dev.off()

# Figure 5: Class composition by marker, including unresolved ASVs.
# Classes are not inferred from orders when final_class is missing.
class_names <- sort(unique(class_composition_df$Class))
class_names <- c(setdiff(class_names, "Unresolved"), intersect("Unresolved", class_names))
class_matrix <- sapply(marker_folders, function(m) {
  z <- class_composition_df[class_composition_df$Marker == m, , drop = FALSE]
  v <- z$ASVs[match(class_names, z$Class)]
  v[is.na(v)] <- 0L
  v
})
if (is.null(dim(class_matrix))) class_matrix <- matrix(class_matrix, ncol = length(marker_folders))
rownames(class_matrix) <- class_names
class_pct <- sweep(class_matrix, 2, pmax(1, colSums(class_matrix)), "/") * 100
class_palette <- setNames(grDevices::hcl.colors(length(class_names), "Dark 3"), class_names)
if ("Insecta" %in% class_names) class_palette["Insecta"] <- "#287C71"
if ("Unresolved" %in% class_names) class_palette["Unresolved"] <- "#B5B5B5"
# Give the legend its own full-width panel to prevent cropping.
legend_rows <- max(1, ceiling(length(class_names) / 4))
png(file.path(figures_dir, "05_arthropod_class_composition.png"),
    width = 1750, height = 1100 + 90 * legend_rows, res = 150)
layout(matrix(c(1, 2), ncol = 1), heights = c(5, 1.1 + 0.45 * legend_rows))
par(mar = c(6.2, 5.8, 3.5, 1.5), mgp = c(3.5, 0.85, 0))
barplot(class_pct, names.arg = marker_folders, col = class_palette[class_names],
        las = 1, border = "white", ylim = c(0, 100),
        ylab = "Percentage of eligible ASVs", main = "ASV composition by class",
        cex.names = 0.86, cex.main = 1.3)
par(mar = c(0, 1, 0, 1))
plot.new()
legend("center", legend = class_names, fill = class_palette[class_names],
       ncol = min(4L, length(class_names)), bty = "n", cex = 0.92,
       x.intersp = 0.8, y.intersp = 1.2)
dev.off()

# Figure 6: Insect order heatmap. Each cell is the number of eligible ASVs.
# Only orders actually detected are shown; unresolved orders are separate in CSV.
observed_orders <- sort(setdiff(unique(insect_orders_df$Order), "Unresolved"))
if (length(observed_orders) > 0L) {
  order_matrix <- sapply(marker_folders, function(m) {
    z <- insect_orders_df[insect_orders_df$Marker == m, , drop = FALSE]
    v <- z$ASVs[match(observed_orders, z$Order)]
    v[is.na(v)] <- 0L
    v
  })
  if (is.null(dim(order_matrix))) {
    order_matrix <- matrix(order_matrix, nrow = length(observed_orders),
                           ncol = length(marker_folders))
  }
  rownames(order_matrix) <- observed_orders
  colnames(order_matrix) <- marker_folders
  # Prioritize orders found by more markers, then by total ASV counts.
  order_idx <- order(-rowSums(order_matrix > 0), -rowSums(order_matrix), rownames(order_matrix))
  order_matrix <- order_matrix[order_idx, , drop = FALSE]
  n_orders <- nrow(order_matrix)
  # Publication-ready heatmap: horizontal column headers and no axis lines.
  png(file.path(figures_dir, "06_insect_order_heatmap.png"),
      width = 1600, height = max(850, 470 + 75 * n_orders), res = 180)
  par(mar = c(4.2, 11.5, 7.0, 2.0), xpd = NA, family = "sans")
  palette <- grDevices::colorRampPalette(c("#F2F5F4", "#9CCFC5", "#287C71"))(101)
  max_count <- max(1, order_matrix)
  plot(NA, xlim = c(0.5, ncol(order_matrix) + 0.5),
       ylim = c(n_orders + 0.5, 0.5), axes = FALSE,
       xlab = "", ylab = "", main = "", type = "n", bty = "n")
  # Position labels in the top margin, leaving room for the title above.
  axis(3, at = seq_along(marker_folders), labels = marker_folders,
       las = 1, tick = FALSE, line = -0.4, cex.axis = 0.92, font = 2)
  axis(2, at = seq_len(n_orders), labels = rownames(order_matrix),
       las = 1, tick = FALSE, line = -0.4, cex.axis = 0.95)
  for (i in seq_len(n_orders)) for (j in seq_along(marker_folders)) {
    value <- order_matrix[i, j]
    color_idx <- 1L + round(100 * value / max_count)
    rect(j - 0.47, i - 0.45, j + 0.47, i + 0.45,
         col = palette[color_idx], border = NA)
    text(j, i, labels = if (value == 0L) "–" else as.character(value),
         col = if (value / max_count > 0.65) "white" else "#172A28",
         cex = 0.95, font = 2)
  }
  mtext("Insect order coverage by marker", side = 3, line = 4.4,
        cex = 1.2, font = 2)
  mtext("Cell values: number of ASVs (not read counts)", side = 1,
        line = 2.2, cex = 0.85, col = "#555555")
  dev.off()
} else {
  message("[INFO] No resolved insect orders: insect heatmap not generated.")
}

message("\n[DONE] Outputs saved in: ", normalizePath(output_dir))
message("  01_marker_summary.csv")
message("  02_taxonomic_resolution.csv")
message("  03_filtered_ASV_taxonomy.csv")
message("  04_order_distribution.csv")
message("  05_curation_categories.csv")
message("  06_ASV_inclusion_audit.csv")
message("  07_NCBI_validation_by_category.csv")
message("  08_taxonomic_richness_by_rank.csv")
message("  09_class_composition.csv")
message("  10_insect_order_coverage.csv")
message("  Figures/01_ASV_retention.png")
message("  Figures/02_taxonomic_resolution.png")
message("  Figures/03_curation_categories.png")
message("  Figures/04_taxonomic_richness_by_rank.png")
message("  Figures/05_arthropod_class_composition.png")
message("  Figures/06_insect_order_heatmap.png (if insect orders detected)\n")
print(summary_df[, c("Marker", "ASVs_before_filter", "ASVs_after_filter",
                     "Eligible_ASVs_n", "Red_excluded_n", "Orange_excluded_n",
                      "Species_n", "Genus_n", "Family_n", "Order_n", "Class_n",
                     "Classes_distinct_n", "Missing_final_class_n",
                     "Unassigned_n", "Arthropoda_explicit_n")], row.names = FALSE)
