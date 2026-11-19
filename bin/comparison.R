#!/usr/bin/env Rscript

# =====================================================================
# Tiara vs FCS-GX comparison and visualization script version 7
# =====================================================================

# written by @NiallG1 on 27/03/2026

#install package
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)

################################################################################
# Input Parameters from Nextflow
################################################################################

# ------------------------------
# Get command-line arguments
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)

input_tsv <- args[1]  # FCS-GX TSV file
input_txt <- args[2]  # TIARA TXT file (if you need it)
output_prefix <- args[3]  # Sample ID

# ------------------------------
# Read input files
# ------------------------------
cat("Reading FCS-GX results from:", input_tsv, "\n")

fcs <- read.table(
  input_tsv,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

cat("Loaded", nrow(fcs), "rows from FCS-GX\n")

# ------------------------------
# Read TIARA input file
# ------------------------------
cat("Reading TIARA results from:", input_txt, "\n")

tiara <- read.table(
  input_txt,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

cat("Loaded", nrow(tiara), "rows from TIARA\n")

# ------------------------------
# Step 3: Process FCS-GX results
# ------------------------------
df_fcs <- fcs %>%
  select(
    seq_id = seq.id,
    seq_len = seq.len,
    div1 = div,
    species_fcs = tax.name.1,
    taxid_fcs = tax.id.1
  ) %>%
  mutate(
    seq_id = as.character(seq_id),
    seq_len = as.numeric(seq_len),
    div1 = ifelse(
      is.na(div1) | div1 == "" | tolower(div1) %in% c("unassigned", "unknown"),
      "unknown",
      div1
    ),
    species_fcs = ifelse(
      is.na(species_fcs) | species_fcs == "" |
      tolower(species_fcs) == "unknown",
      "unknown",
      species_fcs
    ),
    div1 = str_trim(div1),
    species_fcs = str_trim(species_fcs),
    domain_fcs = case_when(
      grepl("^prok", div1, ignore.case = TRUE) |
      grepl("^bact", div1, ignore.case = TRUE) ~ "bacteria",
      grepl("^arch", div1, ignore.case = TRUE) ~ "archaea",
      grepl("^fung", div1, ignore.case = TRUE) |
      grepl("^plnt", div1, ignore.case = TRUE) |
      grepl("^anml", div1, ignore.case = TRUE) ~ "eukarya",
      tolower(div1) == "unknown" ~ "unknown",
      TRUE ~ "unknown"
    )
  )

# --------------------------------------------
# Step 3B: Process FCS-GX results for chimeras
# --------------------------------------------
df_fcs_collapsed <- df_fcs %>%
  mutate(contig_base = str_remove(seq_id, "~.*")) %>%
  group_by(contig_base) %>%
  summarise(
    seq_len = sum(seq_len),

    n_domains = n_distinct(domain_fcs[domain_fcs != "unknown"]),
    n_species = n_distinct(species_fcs[species_fcs != "unknown"]),

    species_list = paste(
      unique(species_fcs[species_fcs != "unknown"]),
      collapse = "; "
    ),

    species_fcs = case_when(
      n_domains > 1 ~ "possible cross-domain chimera",
      n_species == 1 ~ first(species_fcs[species_fcs != "unknown"]),
      n_species > 1 ~ "possible chimera",
      TRUE ~ "unknown"
    ),

    chimera_species = case_when(
      species_fcs %in% c("possible chimera", "possible cross-domain chimera") ~ species_list,
      TRUE ~ NA_character_
    ),

    taxid = case_when(
      species_fcs %in% c("possible chimera","possible cross-domain chimera","unknown") ~ "32644",
      TRUE ~ as.character(first(taxid_fcs[!is.na(taxid_fcs)]))
    ),

    domain_fcs = case_when(
      n_domains == 1 ~ first(domain_fcs[domain_fcs != "unknown"]),
      n_domains > 1 ~ "mixed",
      TRUE ~ "unknown"
    ),

    div1 = first(div1),
    .groups = "drop"
  ) %>%
  rename(seq_id = contig_base)


# ------------------------------
# Step 4: Process Tiara results
# ------------------------------
df_tiara <- tiara %>%
  mutate(
    seq_id = as.character(sequence_id),
    domain_tiara = ifelse(
      class_fst_stage == "" | tolower(class_fst_stage) == "unknown",
      "unknown",
      class_fst_stage
    ),
    domain_tiara = str_trim(domain_tiara)
  ) %>%
  select(seq_id, domain_tiara)

# ------------------------------
# Step 5: Merge and compare new
# ------------------------------

df_compare <- df_fcs_collapsed %>%
  left_join(
    df_tiara %>% select(seq_id, domain_tiara),
    by = "seq_id"
  ) %>%
  replace_na(list(domain_tiara = "unknown")) %>%
  mutate(
    taxid = case_when(
       domain_tiara == "unknown" ~ "32644",
      TRUE ~ as.character(taxid)
    ),
    match = case_when(
      species_fcs %in% c("possible chimera", "possible cross-domain chimera") ~ "chimera",
      domain_tiara == domain_fcs ~ "match",
      TRUE ~ "mismatch"
    )
  )

# ------------------------
# step 6: create  blobtags
# ------------------------
df_compare <- df_compare %>%
  mutate(
    # canonicalize species unknowns (case-insensitive) and trim
    species_fcs = ifelse(is.na(species_fcs) | tolower(species_fcs) == "unknown", "unknown", species_fcs),
    species_fcs = stringr::str_trim(species_fcs),
    species_fcs = gsub(" +", "_", species_fcs),        # replace spaces with underscores

    # canonicalize domains (case-insensitive) and sanitize for tags
    domain_fcs  = ifelse(is.na(domain_fcs)  | tolower(domain_fcs)  == "unknown", "unknown", domain_fcs),
    domain_tiara = ifelse(is.na(domain_tiara) | tolower(domain_tiara) == "unknown", "unknown", domain_tiara),
    domain_fcs = tolower(gsub(" +", "_", stringr::str_trim(domain_fcs))),
    domain_tiara = tolower(gsub(" +", "_", stringr::str_trim(domain_tiara))),

    # Create blob_tag depending on match status (use lowercase "unknown" consistently)
    blob_tag = case_when(
      seq_len < 1000 ~ "unknown",                # override small contigs
      match == "match" ~ species_fcs,            # agreement → species
      match == "chimera" ~ "possible_chimera",   # chimera case
      TRUE ~ paste0(domain_tiara, "_", domain_fcs)  # mismatch (both normalized)
    )
  )

# ----------------------------------------
# step 7: Save the merged comparison table
# ----------------------------------------
out_file <- paste0(output_prefix, "_tiara_vs_fcs_compare.tsv")
write_tsv(df_compare, out_file)
cat("Saved comparison table to:", out_file, "\n")

# -----------------------------
# Step 8: Export comparison CSV
# -----------------------------
#now all contigs >1kbp are labelled as "unknown"
#this is as tiara will not test contigs below 1kbp and FCS-GX is not accurate below 1kbp.

blob_taxonomy <- df_compare %>%
  select(seq_id, blob_tag) %>%
  rename(taxonomy = blob_tag)

# Write it to the output directory
blob_file <- paste0(output_prefix, "_blobtools_taxonomy.tsv")
write_tsv(blob_taxonomy, blob_file)
cat("Saved blobtools taxonomy to:", blob_file, "\n")


# ------------------------------
# Step 9: create a fake blast tsv with taxid
# ------------------------------
# 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=qstart,11=qend,14=evalue

blob_hits <- df_compare %>%
  select(seq_id, taxid) %>%
  mutate(bitscore = 1, four = "NA", sseqid = "NA", six = 1, seven = 1, eight = 1, nine = 1, ten = 1, eleven = 1, twele = 1, thirteen = 1, foureen = 1)

# Write it to the output directory
hits_file <- paste0(output_prefix, "_blobtools_hits.tsv")
write_tsv(blob_hits, hits_file, col_names = FALSE)
cat("Saved blobtools hits to:", hits_file, "\n")
