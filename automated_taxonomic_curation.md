# Automated Taxonomic Curation and NCBI Lineage Completion

This document describes a two-stage R workflow for automated taxonomic curation and NCBI lineage completion of amplicon sequence variants (ASVs). It covers input requirements, configuration, taxonomic decision rules, external database validation, output files, and execution.

## 1. Overview

### 1.1. Purpose and scope

The workflow reduces the time and manual effort required to curate taxonomic assignments in metabarcoding datasets, particularly those containing thousands of ASVs. It integrates molecular evidence with external taxonomic resources, applies standardized decision rules, and identifies assignments requiring manual review.

The workflow comprises two sequential stages:

1. **Automated taxonomic curation and validation:** integrates LULU–BASTA, BLAST, GBIF, and NCBI information to generate curated taxonomic assignments.
2. **NCBI lineage completion:** fills supported missing higher taxonomic ranks and reconciles curation categories.

The workflow begins after ASV inference, LULU–BASTA processing, and BLAST searches have been completed. It does not generate the initial molecular assignments.

### 1.2. Scripts and execution order

| Script | Function |
|---|---|
| `taxonomic_curation_ncbi_validated.R` | Performs exploratory diagnostics or automated taxonomic curation, depending on the selected execution mode. |
| `complete_final_taxonomy_ncbi_color_revision.R` | Completes missing higher taxonomic ranks using NCBI Taxonomy and reconciles curation categories. |

The recommended execution sequence is:

1. Run the first script in `explore` mode.
2. Examine the BLAST diagnostics and GBIF threshold exploration results.
3. Select the most appropriate curation parameters for the dataset and molecular marker.
4. Run the first script in `curate` mode using the selected parameters.
5. Run the second script to complete missing taxonomic ranks and reconcile curation categories.

## 2. Requirements and Dependencies

### 2.1. R packages

The workflow uses the following R packages:

`dplyr`, `stringr`, `purrr`, `rgbif`, `jsonlite`, `rentrez`, and `xml2`.

They can be installed using:

```r
install.packages(c(
  "dplyr", "stringr", "purrr", "rgbif",
  "jsonlite", "rentrez", "xml2"
))
```

An internet connection is required when retrieving information not already available in local caches.

### 2.2. External databases and services

The workflow interacts with several external resources:

- **NCBI Nucleotide (`nuccore`):** retrieves taxonomic identifiers associated with BLAST reference accessions.
- **NCBI Taxonomy (`taxonomy`):** provides scientific names, taxonomic ranks, and lineages for validation and completion.
- **GBIF:** provides taxonomic matching and geographic occurrence evidence.
- **NBN Atlas / UK Species Inventory:** provides common names when available.

These resources contribute different forms of evidence. Geographic occurrence information is supplementary and does not independently confirm a molecular identification.

## 3. Script 1 — Automated Taxonomic Curation

**Script:** `taxonomic_curation_ncbi_validated.R`

The first script supports two execution modes:

- `explore`: generates BLAST diagnostics and evaluates GBIF observation thresholds.
- `curate`: performs automated taxonomic curation and NCBI validation.

The selected mode is controlled by `analysis_mode`.

### 3.1. Input files

The script requires three input files generated during the preceding LULU–BASTA and BLAST processing steps.

In the filenames below, `<marker_name>` represents the molecular marker or dataset identifier (e.g., `LgV2` or `LgV1reamp`).

| Setting | Expected filename | Description |
|---|---|---|
| `lulu_basta_file` | `<marker_name>_asvs_taxonomy_gbif.csv` | Merged ASV abundance and taxonomic assignment table produced by the preceding LULU–BASTA workflow. |
| `blast_file` | `<marker_name>_combined_results.txt` | Combined BLAST results in tab-separated format, without a header. |
| `basta_vout_file` | `<marker_name>_basta_vout.txt` | BASTA `.vout` output containing supporting taxonomic lineage information for each ASV. |

**Input format requirements**

- **LULU–BASTA table:** CSV containing `seqid`, sample abundance columns, and taxonomic assignment information.
- **BLAST results:** Tab-separated text file without a header, containing eight fields in this order: `seqid`, `qacc`, `sseqid`, `pident`, `alignment_length`, `mismatch`, `evalue`, and `stitle`.
- **BASTA lineage file:** Text file containing ASV-specific sections identified by `### ASV_IDENTIFIER`, followed by supporting taxonomic lineage information.

ASV identifiers must be consistent across the input files.

The filename patterns reflect the original dataset organization. The script uses the configured file paths and does not require these exact filenames.

### 3.2. Configuration

#### File paths and execution settings

| Setting | Description |
|---|---|
| `lulu_basta_file` | Path to the merged LULU–BASTA ASV table. |
| `blast_file` | Path to the eight-column BLAST results. |
| `basta_vout_file` | Path to the BASTA `.vout` file. |
| `output_file` | Path to the main curated CSV output. |
| `gbif_cache_file` | Path to the compact GBIF evidence cache. |
| `analysis_mode` | Selects `explore` for diagnostics or `curate` for automated taxonomic curation. |
| `target_phyla` | Defines the phyla included in the target taxonomic scope. |

#### BLAST identity and alignment parameters

| Parameter | Description |
|---|---|
| `expected_amplicon_length` | Expected amplicon length, used as a marker-specific reference. |
| `species_min_identity` | Minimum percentage identity for species-level evaluation. |
| `genus_min_identity` | Minimum percentage identity for genus-level evaluation. |
| `family_min_identity` | Minimum percentage identity for family-level evaluation. |
| `order_min_identity` | Minimum percentage identity for order-level evaluation. |
| `species_min_alignment` | Minimum alignment length for species-level evaluation. |
| `genus_min_alignment` | Minimum alignment length for genus-level evaluation. |
| `family_min_alignment` | Minimum alignment length for family-level evaluation. |
| `order_min_alignment` | Minimum alignment length for order-level evaluation. |

Minimum identity and alignment requirements must decrease consistently from species to order.

The script validates this ordering and stops when the configured thresholds violate it.

These parameters are marker-dependent and should be evaluated before transferring the workflow to another amplicon.

#### Candidate comparison and support

| Parameter | Description |
|---|---|
| `species_pident_margin` | Percentage-identity margin used to distinguish competing species-level candidates. |
| `genus_pident_margin` | Percentage-identity margin used to distinguish competing genus-level candidates. |
| `family_pident_margin` | Percentage-identity margin used to distinguish competing family-level candidates. |
| `order_pident_margin` | Percentage-identity margin used to distinguish competing order-level candidates. |
| `species_alignment_margin` | Alignment-length margin used to distinguish competing species-level candidates. |
| `genus_alignment_margin` | Alignment-length margin used to distinguish competing genus-level candidates. |
| `family_alignment_margin` | Alignment-length margin used to distinguish competing family-level candidates. |
| `order_alignment_margin` | Alignment-length margin used to distinguish competing order-level candidates. |
| `min_accession_advantage` | Minimum advantage in distinct reference-accession support required to resolve certain competitive candidates. |
| `best_pident_tolerance` | Maximum allowed difference in percentage identity from the highest BLAST identity when selecting best hits for each ASV. |
| `best_alignment_tolerance` | Maximum allowed difference in alignment length (bp) from the longest alignment among the identity-filtered hits. |

Only distinct subject accessions count toward reference support.

Reference support is secondary: an accession majority does not override materially superior alignment evidence.

The two `best_*_tolerance` settings define the initial best-hit set. They are distinct from the rank-specific margins used to resolve competing taxonomic candidates.

By default, both tolerances are set to zero, requiring exact ties in the corresponding BLAST statistics.

#### GBIF and external query parameters

| Parameter | Description |
|---|---|
| `country_code` | Country code used for geographic GBIF occurrence queries. |
| `minimum_gbif_records` | Minimum number of qualifying GBIF observation records used in geographic evidence evaluation. |
| `gbif_thresholds_to_test` | Candidate GBIF observation thresholds evaluated during exploration. |
| `query_gbif` | Enables or disables new GBIF queries. |
| `gbif_pause_seconds` | Delay between GBIF requests. |
| `query_nbn` | Enables or disables new NBN Atlas queries. |
| `nbn_pause_seconds` | Delay between NBN Atlas requests. |

#### Control detection parameters

| Parameter | Description |
|---|---|
| `positive_control_pattern` | Positive-control sample columns. |
| `negative_control_pattern` | Negative-control sample columns. |
| `control_detection_minimum` | Minimum count required to report control detection in the ASV notes. |

Control detection is recorded as supplementary information and does not automatically modify the taxonomic assignment or curation category.

### 3.3. Explore mode

Set:

```r
analysis_mode <- "explore"
```

This mode characterizes the molecular and geographic evidence available in the dataset to support informed parameter selection before automated curation.

It includes two complementary exploratory components.

#### BLAST diagnostics

The script identifies a representative best BLAST hit for each ASV represented in the retained input dataset.

Representative hits are selected deterministically using percentage identity, alignment length, and subject identifier where needed.

The diagnostic summary reports:

- Total number of ASVs.
- Number of ASVs with BLAST hits.
- Mean, median, minimum, and maximum alignment lengths.
- Mean, median, minimum, and maximum percentage identities.

These statistics help characterize the molecular evidence and evaluate whether the configured identity and alignment-length requirements are appropriate for the marker and dataset. 

#### GBIF threshold exploration

The script evaluates candidate thresholds specified in:

```r
gbif_thresholds_to_test
```

For each threshold, it summarizes the number of ASVs with sufficient geographic occurrence evidence.

It also distinguishes ASVs without observation records from those with at least one record.

Comparing thresholds helps assess how restrictive different GBIF observation cutoffs would be for the dataset.

#### Interpreting exploratory results

The BLAST and GBIF results should be considered together when selecting curation parameters.

BLAST statistics support evaluation of marker-specific identity and alignment thresholds, while GBIF exploration informs the selection of a geographic observation cutoff.

The exploratory outputs guide parameter selection; they do not automatically optimize the configuration.

After reviewing the results, update the relevant settings and run the script in `curate` mode.

### 3.4. Curate mode

Set:

```r
analysis_mode <- "curate"
```

This mode evaluates the available taxonomic evidence for each ASV, applies automated decision rules, and produces a curated classification with supporting diagnostic and validation information.


#### 3.4.1. Taxonomic evidence integration

The script integrates five complementary sources of taxonomic evidence:

| Evidence source | Role in the curation workflow |
|---|---|
| LULU–BASTA merged table | Provides the original taxonomic assignments and ASV abundance information. |
| BLAST results | Provide molecular evidence through percentage identity, alignment length, candidate taxon names, and distinct reference-accession support. |
| BASTA `.vout` | Provides supporting taxonomic lineages for concordance assessment and provisional assignments. |
| GBIF | Provides taxonomic name matching and geographic occurrence evidence based on the configured country and observation threshold. |
| NCBI | Provides accession-associated taxonomic information and reference lineages for validation. |

GBIF is used to match candidate scientific names against its taxonomic backbone and to retrieve geographic occurrence records for the configured country. The occurrence evidence is evaluated using `minimum_gbif_records` to distinguish taxa that meet the selected observation threshold from those with insufficient or no retrieved records.

Geographic occurrence evidence is supplementary and does not independently confirm or invalidate a molecular taxonomic assignment.


#### 3.4.2. Automated curation decision rules

The script applies a hierarchical decision process to select the most specific taxonomic assignment supported by the molecular evidence.

BLAST-based resolution is prioritized. BASTA assignments contribute to concordance assessment, lineage consistency checks, and provisional fallback decisions.

**A. Rank-specific candidate resolution**

At each rank, candidate hits must satisfy the corresponding minimum identity and alignment-length thresholds.

Candidates are grouped by taxon, retaining the strongest alignment for each distinct reference accession.

The leading candidate is evaluated against competing taxa using the rank-specific identity and alignment margins.

The decision follows three general cases:

1. **Clear molecular advantage:** a candidate is accepted when its alignment evidence sufficiently distinguishes it from competitors.
2. **Near ties:** when competing candidates have similar identity and alignment statistics, an advantage in distinct reference-accession support may resolve the comparison.
3. **Unresolved competition:** when the evidence does not sufficiently distinguish competing taxa, the rank is treated as ambiguous.

Accession support cannot override materially superior alignment evidence.

**B. Selection of the final taxonomic rank**

The script evaluates species, genus, family, and order, selecting the most specific successfully resolved rank.

| BLAST evidence | Automated decision |
|---|---|
| Species-level criteria are satisfied and competing species are resolved. | Assign species. |
| Species-level evidence is ambiguous, but genus-level evidence is resolved. | Assign genus. |
| Species and genus cannot be resolved, but family-level evidence is resolved. | Assign family. |
| Only order-level evidence can be resolved. | Assign order. |
| No evaluated rank can be resolved. | Apply the BASTA fallback procedure. |

This prevents the workflow from assigning a taxonomic rank more specific than the available molecular evidence supports.

**C. Comparison with BASTA**

When BLAST successfully resolves a taxon, the selected BLAST assignment is retained even if it disagrees with the original merged BASTA assignment.

The script examines whether:

- The original BASTA assignment agrees with the selected taxon.
- The original BASTA assignment disagrees with the selected taxon.
- The BASTA `.vout` consensus supports the selected taxon.
- No concordant BASTA evidence is available at the selected rank.

The outcome is documented in the curation explanation.

BASTA agreement contributes to the subsequent curation category but does not override a successfully resolved BLAST assignment.

**D. Taxonomic lineage consistency**

After selecting a taxon, the script constructs the final lineage while preventing unsupported combinations of taxonomic ranks.

The implemented safeguards depend on the selected rank:

| Selected rank | Lineage handling |
|---|---|
| Species | The genus is derived from the selected species name. Additional higher-rank information may be recovered from matching, NCBI-validated BLAST evidence. |
| Genus | Higher-rank information may be recovered from NCBI-validated BLAST evidence associated with the selected genus. |
| Family | Higher ranks are retained only when supported by a matching BASTA `.vout` family lineage. |
| Order | Higher ranks are retained only when supported by a matching BASTA `.vout` order lineage. |

Higher-rank assignments inherited from BASTA are cleared when they conflict with the selected BLAST-supported lineage.

Taxonomic information from an unrelated lineage is not automatically transferred to the selected taxon.

Ranks below the final assigned rank are cleared to avoid retaining unsupported taxonomic resolution.

**E. Fallback when BLAST resolution fails**

If none of the four evaluated ranks can be resolved, the script applies a conservative fallback.

**Case 1 — Unanimous BASTA `.vout` lineage available**

The script retains the most specific rank supported by the unanimous BASTA `.vout` consensus.

**Case 2 — No unanimous BASTA `.vout` lineage available**

The original merged LULU–BASTA assignment is retained provisionally.

Both cases require manual review.

The fallback procedure records low confidence, identifies the provisional evidence basis, and preserves an explanation of why BLAST resolution failed.

A provisional assignment should not be interpreted as a confidently resolved molecular identification.

**F. Assignment outcome and confidence**

After selecting the final assignment, the script records its outcome in `curation_decision` by comparing the final taxonomic rank with the original BASTA rank. The possible outcomes are `Assigned`, `Retained`, `Upgraded`, `Downgraded`, or `Manual review`.

The `confidence` field is assigned as `High` for resolved species-level assignments supported by at least two distinct reference accessions, `Moderate` for other BLAST-resolved assignments, and `Low` for provisional BASTA assignments requiring manual review.

These fields describe the outcome and supporting evidence of the automated decision; they do not introduce additional taxonomic selection rules.

#### 3.4.3. Curation categories and manual review

Each ASV receives a category in `curation_colour`.

| Category | Meaning |
|---|---|
| `Green` | The final assignment agrees with BASTA at the same taxonomic rank, without an overriding review or out-of-target condition. |
| `Yellow` | An automated assignment was obtained, but the conditions for Green were not met. |
| `Red` | The assignment requires manual review because the available evidence remains unresolved. |
| `Orange` | The final assignment belongs to a known phylum outside the configured target taxonomic scope. |

The implemented priority order is:

```text
Orange > Red > Green > Yellow
```

#### 3.4.4. NCBI validation and lineage audit

The first script performs NCBI-based checks using BLAST reference accessions and their associated taxonomic information.

**Accession-associated validation**

BLAST reference accessions are linked to NCBI taxonomic identifiers using nucleotide-record information.

These identifiers allow the script to retrieve NCBI Taxonomy records and examine the taxonomic information associated with the reference sequences.

**Post-curation lineage audit**

After selecting the final assignment, the script evaluates its consistency with supporting NCBI lineages.

The audit considers the seven principal taxonomic ranks:

```text
kingdom, phylum, class, order, family, genus, species
```

The outcome is recorded in `validation_status`.

| Validation status | Meaning |
|---|---|
| `lineage_consistent` | At least one supporting NCBI lineage is compatible with all populated final taxonomic ranks. |
| `lineage_inconsistent` | No supporting NCBI lineage is compatible with all populated final taxonomic ranks. Manual review is recommended. |
| `lineage_disagreement` | The final taxon disagrees with the NCBI assignment at the corresponding rank. Manual review is recommended. |
| `rank_absent_in_ncbi` | The required rank is unavailable in the retrieved NCBI references, preventing complete validation. |
| `no_verified_accession` | No NCBI-verified BLAST accession is available for the ASV. |
| `unassigned` | No valid final taxonomic assignment is available for evaluation. |

Validation evaluates taxonomic consistency, not the biological accuracy of the molecular identification itself.

The procedure does not modify the selected final assignment.

`validation_status` is independent of `curation_colour` and should not be treated as an additional curation category or automatic exclusion criterion.

### 3.5. Output files

The script generates a main curated ASV table and supporting diagnostic, audit, and cache files. Some outputs depend on the selected execution mode.

In the filenames below, `<marker_name>` represents the molecular marker or dataset identifier.

| Output filename | Description |
|---|---|
| `<marker_name>_automated_curation_v5g_ncbi_validated.csv` | **Main output:** curated ASV table and input for Script 2. |
| `<marker_name>_automated_curation_v5g_ncbi_validated_family_label_audit.csv` | Family-name extraction audit. |
| `<marker_name>_automated_curation_v5g_ncbi_validated_nbn_common_names.csv` | NBN common-name cache. |
| `<marker_name>_automated_curation_v5g_ncbi_validated_ncbi_accession_audit.csv` | NCBI accession validation audit. |
| `<marker_name>_automated_curation_v5g_ncbi_validated_ncbi_accession_cache.rds` | NCBI accession lookup cache. |
| `<marker_name>_automated_curation_v5g_ncbi_validated_ncbi_final_audit.csv` | Final NCBI lineage consistency audit. |
| `<marker_name>_gbif_cache_v5g_rankfix.csv` | Compact GBIF evidence cache. |
| `<marker_name>_gbif_cache_v5g_rankfix_internal.csv` | Internal GBIF query cache. |
| `<marker_name>_gbif_cache_v5g_threshold_exploration.csv` | GBIF observation-threshold exploration results. |

The main curated CSV is the only file passed directly to Script 2.

The filename patterns reflect the original configuration. Actual filenames depend on the paths configured in the script.

## 4. Script 2 — NCBI Lineage Completion

**Script:** `complete_final_taxonomy_ncbi_color_revision.R`

This script performs a post-curation procedure that completes missing higher taxonomic ranks using NCBI Taxonomy.

It preserves existing populated assignments, prevents inconsistent lineage completion, and reconciles curation categories.

### 4.1. Input file

The second script requires the main curated ASV table generated by Script 1 in `curate` mode.

| Setting | Expected filename | Description |
|---|---|---|
| `input_file` | `<marker_name>_automated_curation_v5g_ncbi_validated.csv` | Main output of Script 1, containing curated taxonomic assignments, supporting evidence, curation categories, and NCBI validation results. |

The input must be a CSV containing unique `seqid` identifiers, `final_taxon`, `final_rank`, the seven `final_*` taxonomic columns, and the curation fields required for category reconciliation.

Audit files and caches generated by Script 1 are not required as direct inputs to Script 2.

### 4.2. Configuration

#### File paths

| Setting | Description |
|---|---|
| `input_file` | Path to the main curated CSV generated by Script 1. |
| `output_file` | Path to the completed taxonomic CSV. Must differ from `input_file`. |
| `audit_file` | Path to the ASV-level NCBI lineage completion audit. |
| `cache_file` | Path to the NCBI Taxonomy lookup cache. |

#### Additional parameters

| Parameter | Description |
|---|---|
| `request_pause_seconds` | Delay between NCBI taxonomy requests. |
| `target_phyla` | Target phyla used when reconciling final curation categories. |

The script evaluates the seven principal taxonomic ranks:

```text
kingdom, phylum, class, order, family, genus, species
```

The completion summary filename is defined separately through `summary_file` and should be adapted when changing datasets.

### 4.3. NCBI taxonomy retrieval and caching

The script queries NCBI Taxonomy using the selected scientific name and expected taxonomic rank.

It requires an exact scientific-name and rank match.

This prevents broad or ambiguous search results from being accepted automatically.

For a successfully matched taxon, the script retrieves its NCBI taxonomic identifier and available lineage information.

Successful and unsuccessful non-error lookups are cached to avoid unnecessary repeated queries.

Temporary query errors are not cached, allowing them to be retried during subsequent executions.

### 4.4. Lineage completion criteria

The procedure uses the taxon selected during the original curation as its **taxonomic anchor**.

A valid anchor requires:

1. A supported taxonomic rank in `final_rank`.
2. A usable scientific name in `final_taxon`.
3. Agreement between `final_taxon` and the corresponding `final_*` field.
4. Compliance with the implemented name-format checks for species-level anchors.

Provisional or non-binomial species labels are not accepted as valid species anchors.

**Higher-rank completion only**

The script fills only missing taxonomic ranks above the selected anchor.

For example, a genus-level assignment may support completion of family, order, class, phylum, and kingdom.

It does not infer species from a genus-level assignment.

**Preservation of existing assignments**

Populated `final_*` fields are never overwritten.

The script does not replace existing assignments based solely on the retrieved NCBI lineage.

**Lineage consistency**

Before completing missing ranks, the script compares the existing populated lineage against the retrieved NCBI lineage.

If a conflict is detected, no ranks are filled for that ASV.

This prevents the creation of inconsistent hybrid taxonomic lineages.

**Missing formal ranks**

NCBI does not always provide a formal kingdom rank for animal taxa.

The script does not replace missing kingdom information with a superkingdom or an inferred value.

### 4.5. Conflict detection and unresolved assignments

Each ASV receives an audit status describing the completion outcome.

| Status | Meaning |
|---|---|
| `skipped` | No valid taxonomic anchor was available. |
| `complete` | No higher taxonomic ranks were missing. |
| `filled` | One or more missing higher ranks were completed. |
| `no_ncbi_higher_ranks` | NCBI did not provide the missing higher ranks requested. |
| `conflict` | An existing assignment disagreed with the retrieved NCBI lineage. |
| `not_found` | No suitable exact NCBI taxonomic match was found. |
| `ambiguous` | The NCBI search could not be resolved to a single matching taxon. |
| `query_error` | The NCBI request failed. |

An exact NCBI match is an intermediate lookup outcome; the final ASV-level status reflects whether the requested lineage completion was possible.

Unresolved and conflicting assignments remain unchanged.

### 4.6. Curation category reconciliation

After lineage completion, the script reevaluates the final phylum against `target_phyla`.

If `final_phylum` is known and falls outside the target group, the ASV is assigned:

```r
curation_colour <- "Orange"
```

Orange has priority over the previous curation category.

The script does not automatically promote Red or Yellow assignments to Green.

A previously Orange assignment is not automatically promoted if the final phylum is missing or changes.

The previous curation category is preserved in `curation_colour_before_ncbi`, allowing changes introduced during lineage completion to be traced.

### 4.7. Output files

The script produces the final completed ASV table and supporting NCBI lineage-completion audits.

| Output filename | Description |
|---|---|
| `<marker_name>_automated_curation_v5g_ncbi_validated_completed.csv` | **Main output:** final ASV table after NCBI lineage completion and curation category reconciliation. |
| `<marker_name>_automated_curation_v5g_ncbi_validated_completed_colour_audit.csv` | Audit of curation category changes. |
| `<marker_name>_ncbi_completion_audit.csv` | ASV-level lineage completion audit. |
| `<marker_name>_ncbi_completion_cache.rds` | NCBI Taxonomy lookup cache. |
| `<marker_name>_ncbi_completion_summary.csv` | Summary of missing and completed taxonomic ranks. |

The completed CSV is the final output of the two-stage workflow and is intended for downstream analyses after quality-control review.

The original curated CSV remains unchanged.

## 5. Running the Workflow

### 5.1. Explore mode

Configure the required input paths in the first script and set:

```r
analysis_mode <- "explore"
```

Review the marker-specific thresholds and candidate GBIF thresholds before execution.

Run:

```bash
Rscript taxonomic_curation_ncbi_validated.R
```

Alternatively, supply the five supported path arguments:

```bash
Rscript taxonomic_curation_ncbi_validated.R \
  basta.vout \
  blast_results.txt \
  lulu_basta.csv \
  curated_output.csv \
  gbif_cache.csv
```

The arguments must be supplied in exactly this order.

The first script accepts either zero or five positional arguments.

### 5.2. Selecting parameters

Use the exploratory results to review the configuration before automated curation.

Consider:

- BLAST identity and alignment-length distributions.
- Expected amplicon length and marker-specific taxonomic resolution.
- Rank-specific identity and alignment requirements.
- GBIF occurrence coverage and the effects of different observation thresholds.
- Target taxonomic group and geographic scope.

The selected parameters should reflect the characteristics of the dataset and molecular marker.

The values defined in the original script are provisional configurations rather than universal recommendations.

### 5.3. Curate mode

Set:

```r
analysis_mode <- "curate"
```

Confirm that the selected thresholds and taxonomic settings have been configured.

Run:

```bash
Rscript taxonomic_curation_ncbi_validated.R
```

The script generates the curated ASV table and supporting diagnostic and validation files.

After execution, inspect the curation categories, unresolved assignments, and NCBI validation results.

### 5.4. NCBI lineage completion

Configure the second script using the curated CSV from Script 1.

Update the input, output, audit, and cache paths, as well as the target taxonomic scope.

Run:

```bash
Rscript complete_final_taxonomy_ncbi_color_revision.R
```

The script completes supported higher ranks, checks lineage consistency, and reconciles curation colours.

Inspect the completion audit, colour audit, and rank-level summary before using the final dataset.

### 5.5. Adapting the workflow to another marker

To apply the workflow to another molecular marker:

1. Generate the corresponding LULU–BASTA and BLAST inputs.
2. Verify consistent ASV identifiers and input formats.
3. Update the input and output paths.
4. Review marker-specific identity, alignment, and candidate comparison thresholds.
5. Configure the target taxonomic group and geographic scope.
6. Run `explore` mode and examine its diagnostic outputs.
7. Select appropriate curation parameters.
8. Run `curate` mode.
9. Configure and run NCBI lineage completion.
10. Inspect the final outputs and audit files.

Dataset-specific configurations and empirical results should be documented separately.

## 6. Final Output Data Dictionary

The final output of the workflow is the completed ASV table:

`<marker_name>_automated_curation_v5g_ncbi_validated_completed.csv`

This file integrates the curated taxonomic assignments generated by the first script with the NCBI lineage completion and curation category reconciliation performed by the second script.

The table preserves the original ASV identifiers, sample abundance data, and relevant input information while incorporating the final taxonomic assignments and supporting evidence.

The following sections describe the principal columns of the completed dataset, organized by their analytical purpose rather than by the script that generated them.

### 6.1. BLAST and GBIF evidence

| Column | Meaning |
|---|---|
| `top_pident` | Percentage identity of the top BLAST evidence. |
| `top_alignment` | Alignment length associated with the top BLAST evidence. |
| `top_BLAST_hits` | Summary of top BLAST hits considered for the ASV. |
| `number_best_hits` | Number of hits retained in the best-hit set. |
| `number_unique_accessions` | Number of distinct reference accessions represented in the selected evidence. |
| `number_unique_species` | Number of distinct species represented in the selected evidence. |
| `number_unique_genus` | Number of distinct genera represented in the selected evidence. |
| `best_species_list` | Species names represented among the best candidate hits. |
| `best_genus_list` | Genus names represented among the best candidate hits. |
| `removed_entries` | Entries removed during evidence filtering. |
| `number_removed_entries` | Number of removed entries. |
| `number_species_checked` | Number of species candidates evaluated using GBIF information. |
| `number_species_present` | Number of evaluated species meeting the geographic observation criterion. |
| `gbif_present_species` | Species names meeting the applicable GBIF observation threshold. |
| `hit_consensus` | Summary of molecular evidence or consensus supporting the assignment. |
| `GBIF_UK_status` | Textual summary of geographic GBIF evidence. |

### 6.2. Final taxonomic assignments

| Column | Meaning |
|---|---|
| `BASTA_rank` | Taxonomic rank of the original BASTA assignment. |
| `final_kingdom` | Final kingdom assignment. |
| `final_phylum` | Final phylum assignment. |
| `final_class` | Final class assignment. |
| `final_order` | Final order assignment. |
| `final_family` | Final family assignment. |
| `final_genus` | Final genus assignment. |
| `final_species` | Final species assignment. |
| `final_taxon` | Selected taxonomic name at the final assigned rank. |
| `final_rank` | Rank of the selected taxon. |
| `common_name` | Common name retrieved from available reference information, when present. |

The `final_*` fields reflect the curated classification after both scripts have been executed.

The second script may fill missing higher-rank fields but does not overwrite populated assignments.

### 6.3. Curation decisions and confidence

| Column | Meaning |
|---|---|
| `curation_decision` | Description of assignment retention, upgrade, downgrade, or manual review. |
| `manual_review` | Indicates whether further manual assessment is required. |
| `curation_reason` | Explanation of the automated decision. |
| `confidence` | Confidence category assigned by the curation procedure. |
| `evidence_basis` | Description of the evidence supporting the selected assignment. |
| `curation_colour` | Final Green, Yellow, Red, or Orange category. |
| `outside_target_phylum` | Indicates whether the known final phylum falls outside the target group. |
| `notes` | Additional observations, including relevant control-detection information. |
| `curation_colour_before_ncbi` | Curation category assigned before NCBI lineage completion and colour reconciliation. |

### 6.4. NCBI validation and lineage completion

| Column | Meaning |
|---|---|
| `validation_status` | Outcome of the NCBI taxonomic consistency audit performed during automated curation. |

NCBI lineage completion does not introduce a separate completion-status column into the final ASV table. Detailed completion outcomes are recorded in the accompanying NCBI completion audit.

The original `validation_status` is preserved without modification by the second script.

## 7. Reproducibility Recommendations

For reproducibility, retain:

- Original input files.
- Exact scripts used.
- R and package versions.
- Configured parameters.
- Exploratory diagnostic outputs.
- Curated and completed taxonomic tables.
- GBIF and NCBI caches.
- Validation, completion, and colour audit files.
- Execution date.

External databases may change over time. Preserving caches and audit files helps document the evidence available during the analysis.
