# Semi-automated taxonomic curation - created by NEHLIN SAYED nehlin.sayed.24@ucl.ac.uk

This repository contains an R workflow for semi-automated curation of metabarcoding ASVs. The current script was developed for an rbcL dataset (but can work for any primer set if rules are followed carefully) and combines:

- a LULU × BASTA merged ASV table `.csv` file
- BASTA `.vout` .txt file
- BLAST tabular results `.txt` file
- GBIF occurrence information for geographic plausibility
- preferred British common names from the NBN Atlas/UK Species Inventory (UKSI)

The workflow produces a proposed final taxonomic assignment for every ASV and clearly identifies assignments that still require manual review. It is deliberately semi-automated: GBIF occurrence counts provide supporting geographic context and are not proof that a species is present or absent

## Main script

`semi_automated_taxonomic_curation.R`

If the downloaded script is currently named `rbcl_taxonomic_curation_v5.R`, either use that filename in the job script or rename it to match the name above.

## Requirements

The script requires R and the following packages:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "rgbif",
  "jsonlite"
))
```

The libraries are loaded at the beginning of the script.

## Input files

### 1. LULU × BASTA merged table

A CSV containing one row per ASV. It must contain:

- `seqid`, which identifies the ASV
- sample abundance columns containing read counts
- the taxonomy columns `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus` and `Species`.

The original columns and read abundances are retained in the curated output.

### 2. BLAST results

A tab-delimited file without a header. Exactly eight columns are expected, in this order:

```text
seqid  qacc  sseqid  pident  alignment_length  mismatch  evalue  stitle
```

The script selects the greatest percentage identity and then the longest alignment among those hits. Literal duplicate lines for the same subject accession are not allowed to inflate the number of best hits.

Named clone records are retained. Named records beginning with `TPA_asm:` or `PREDICTED:` are also retained and processed as ordinary taxonomic matches.

Entries without a reliable species name, such as `uncultured eukaryote`, environmental sequences, `cf.` identifications and names ending in `sp.`, are excluded from the species-level consensus and reported in `removed_entries`.

### 3. BASTA .vout

A sectioned text file containing the alternative BASTA-supported lineages for each ASV. Sections are expected in this form:

```text
###1
10	unknown;Streptophyta;
8	unknown;Streptophyta;Magnoliopsida;...
```

The `.vout` is used to confirm whether the BLAST-supported taxon is represented in the BASTA evidence and to fill missing ranks only when the BASTA lineages agree.

## User settings

Edit the `USER SETTINGS` section before running the script:

```r
lulu_basta_file <- "/path/to/lulu_x_basta.csv"
blast_file <- "/path/to/blast_results.txt"
basta_vout_file <- "/path/to/basta.vout.txt"
output_file <- "/path/to/curated_output.csv"
gbif_cache_file <- "/path/to/gbif_cache.csv"

country_code <- "GB"
minimum_identity <- 98
minimum_alignment_length <- 400
minimum_gbif_records <- 100
```

Required changes:

- Change all five file paths.
- Change `country_code` when analysing a country other than Great Britain. See https://rs.gbif.org/areas/ for more codes. 
- `minimum_identity`: minimum BLAST percentage identity used as reliable evidence. Change depending on what primer you are using.
- `minimum_alignment_length`: minimum alignment length in base pairs. Change depending on what primer you are using.
- `minimum_gbif_records`: minimum number of accepted GBIF observation records used for geographic filtering. This threshold may need to be changed separately for different taxonomic groups. A low count or zero records must not automatically be interpreted as biological absence.

Optional changes:

- `best_pident_tolerance` and `best_alignment_tolerance`: include near-ties instead of exact best hits.
- `query_gbif` and `query_nbn`: disable new online queries and reuse existing caches.
- `gbif_pause_seconds` and `nbn_pause_seconds`: pauses between online queries.
- control-name patterns and `control_detection_minimum`.


## Control detection

Positive- and negative-control columns are detected from their sample names using case-insensitive patterns.

Positive controls can contain names such as:

```text
POS
POS_pilot
positive_control
posctrl
```

Negative controls can contain names such as:

```text
PCRNEG
negative_control
negctrl
NTC
blank
```

The `notes` column reports the control column and read count when an ASV reaches `control_detection_minimum`. If the ASV is not detected in either type of control, `notes` is `NA`. Control detection is reported for laboratory quality control but does not automatically change the taxonomic assignment.

## Curation rules

For each ASV, the script applies the following rules:

1. Retain BLAST hits with the greatest identity and then the longest alignment, subject to the configured tolerances.
2. Require the best evidence to meet the identity and alignment-length thresholds.
3. Remove non-informative species labels from the species-level consensus.
4. If all informative best hits support one species, propose that species and compare it with the merged BASTA assignment and `.vout` evidence.
5. If equally best hits contain several species from one genus, assign the genus because the marker does not resolve species.
6. When several equally good species are present, GBIF observation counts may resolve the result only if the geographic filtering leaves one species or one genus.
7. If equally best hits span multiple genera and cannot be resolved geographically, retain the lowest shared BASTA/GenBank rank, normally family or order, and request manual review.
8. Preserve the merged BASTA/GenBank hierarchy. Missing ranks are filled only from unanimous BASTA `.vout` evidence, not from the GBIF hierarchy.
9. Obtain an NBN/UKSI preferred common name for the final taxon when an exact accepted name and rank are available. Common names never affect the taxonomic decision.
10. Record `NA` in empty output cells.

## GBIF records

The main geographic threshold uses only the combined count of:

- `HUMAN_OBSERVATION`;
- `OBSERVATION`; and
- `MACHINE_OBSERVATION`.

The internal cache additionally records excluded records:

- preserved specimens;
- living specimens;
- fossil specimens;
- material citations;
- material samples;
- generic occurrences; and
- all GBIF records.

The compact GBIF cache retains the original audit format and reports the total GBIF country count. Consequently, `occurrences_GB` in the compact cache and the observation count used by the automated threshold are not necessarily the same number.

## Output files

### Curated CSV

The main output contains all original ASV and abundance columns followed by BLAST, GBIF and curation columns.

Important columns include:

| Column | Meaning |
|---|---|
| `top_pident` | Highest percentage identity among the selected BLAST hits. |
| `top_alignment` | Longest alignment among the highest-identity hits. |
| `top_BLAST_hits` | Taxonomic labels and counts among the selected optimal hits. |
| `number_best_hits` | Number of selected best BLAST rows after literal duplicate removal. |
| `number_unique_accessions` | Number of distinct subject accessions among the best hits. | This column can be removed as it is not working.
| `number_unique_species` | Number of informative species among the best hits. |
| `number_unique_genus` | Number of genera represented by those informative species. |
| `best_species_list` | Unique informative species among the best hits. |
| `best_genus_list` | Unique genera represented by the best species. |
| `removed_entries` | Best-hit labels excluded from species-level consensus. |
| `number_species_checked` | Number of informative candidate species checked against GBIF. |
| `number_species_present` | Candidates reaching the configured GBIF observation threshold. |
| `gbif_present_species` | Candidate names and observation counts above that threshold. |
| `BASTA_rank` | Lowest rank assigned in the merged LULU × BASTA table. |
| `final_kingdom`–`final_species` | Proposed final taxonomic hierarchy. |
| `final_taxon` | Name at the proposed final rank. |
| `final_rank` | Rank of `final_taxon`. |
| `common_name` | Preferred British common name returned by NBN/UKSI, when available. |
| `hit_consensus` | Plain-language summary of agreement among informative best hits. |
| `GBIF_UK_status` | Summary of candidates reaching the GB observation threshold. |
| `curation_decision` | Whether the assignment was retained, upgraded, downgraded or requires review. | 
| `manual_review` | `TRUE` when human examination is required. |
| `curation_reason` | Evidence-based explanation of the proposed decision. |
| `confidence` | `Very high`, `High`, `Moderate` or `Low`, based on the decision pathway. | This column can be removed as it is not working.
| `evidence_basis` | Whether the proposal was sequence-resolved, geography-resolved or unresolved. |
| `curation_colour` | Quick review category described below. | This column can be removed. 
| `notes` | Positive- or negative-control detection and read counts; otherwise `NA`. |

### Colour categories

- Green: no manual review is required and the final rank matches the original BASTA rank.
- Yellow: no manual review is required, but the automated result changes the assigned rank.
- Red: `manual_review` is `TRUE`.

Colour is a review aid rather than an independent measure of taxonomic certainty. This column can be removed. 

### Compact GBIF cache

The filename is specified by `gbif_cache_file`. It contains one row per ASV/candidate-species combination:

```text
seqid, support, species, gbif_name, usageKey,
occurrences_GB, present_GB
```

### Internal GBIF cache

This is created automatically by adding `_internal` to the compact-cache filename. It stores name matches, taxonomy, detailed counts by basis of record, query errors and query dates. It is reused on later runs to reduce repeated GBIF requests.

### NBN common-name cache

This is created automatically by adding `_nbn_common_names` to the curated-output filename. It records exact NBN/UKSI matches and preferred British common names and is reused on later runs.

## Running interactively

After changing the paths in `USER SETTINGS`, open the script in R or RStudio and run the complete file.

The script also accepts five optional command-line arguments in this order:

```text
basta.vout blast.txt lulu_x_basta.csv curated_output.csv gbif_cache.csv
```

If no arguments are supplied, the paths in `USER SETTINGS` are used.

## Running on UCL Myriad

Example job file:

```bash
#!/bin/bash -l

#$ -l h_rt=24:00:00
#$ -l mem=8G
#$ -pe smp 4
#$ -N rbcl_semiauto_curation
#$ -wd /home/ucbtns5/Scratch/rbcl/testing_tax_assign

module -f unload compilers mpi gcc-libs
module load r/4.5.1-openblas/gnu-10.2.0

R --no-save < /home/ucbtns5/Scratch/scripts/semi_automated_taxonomic_curation.R \
  > rbcl_semiauto_curation.out 2>&1
```

Submit it using:

```bash
qsub semi_automated_taxonomic_curation.sh
```

The first run can take considerably longer because GBIF and NBN must be queried. Later runs reuse the cache files. If GBIF reports too many requests, increase `gbif_pause_seconds` and rerun; successfully cached species will not need to be queried again.

## Reproducibility and manual review

For reproducible analyses:

- retain the original input files
- retain the curated output and both online-service caches
- record the script version or Git commit used
- record all threshold settings
- inspect every ASV for which `manual_review` is `TRUE`
- report that GBIF was used as geographic supporting evidence rather than proof of absence.

Database content and accepted names can change over time. Retaining the caches and query dates provides an audit trail of the information used for a particular analysis.

## Scope

The defaults were developed for an rbcL plant metabarcoding dataset. The workflow can be adapted to other primer sets, including Leray COI, but identity thresholds, alignment-length thresholds, taxonomic groups, reference-title parsing, geographic occurrence thresholds, (and potentially reference databases in the future) must be validated for the marker and target organisms before use.


