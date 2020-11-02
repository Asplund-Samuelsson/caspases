options(width=130)

# Load libraries
library(tidyverse)

# Define infiles
tax_file = "results/uniprot.p20.filtered.taxonomy.tab"
dyad_file = "results/uniprot.p20.CH_dyad_classification.tab"
p10_file = "results/uniprot.p20.filtered.p10_CD_classification.tab"
para_file = "results/uniprot.p20.paracaspase.txt"
len_file = "results/uniprot.p20_only.complete_dyad.seq_lengths.tab"
meta_file = "results/uniprot.p20.metacaspase.tab"
acid_file = "results/uniprot.p20.DD_dyad_classification.tab"
pfam_file = "results/uniprot.p20.pfam_architecture.tab"
bp20_file = "results/uniprot.p20.basic_pocket.tab"
bp10_file = "results/uniprot.p10.basic_pocket.tab"

# Load data
tax = read_tsv(tax_file, quote="", comment="")
dyd = read_tsv(dyad_file, col_names=c("UniProt_ID", "HC_dyad"))
p10 = read_tsv(p10_file, col_names=c("UniProt_ID", "p10_CD"))
para = scan(para_file, character())
len = read_tsv(len_file, col_names=c("UniProt_ID", "p20_length"))
meta = read_tsv(meta_file)
acid = read_tsv(acid_file, col_names=c("UniProt_ID", "Acidic_pocket"))
pfam = read_tsv(pfam_file)
bp20 = read_tsv(bp10_file, col_names=c("UniProt_ID", "bp20"))
bp10 = read_tsv(bp20_file, col_names=c("UniProt_ID", "bp10"))

tax = select(tax, -full_lineage)
tax = rename(
  tax,
  UniProt_ID = identifier, Tax_ID = taxid,
  Group = group, Superkingdom = superkingdom, Kingdom = kingdom,
  Phylum = phylum, Class = class, Order = order, Family = family,
  Genus = genus, Species = species, Lowest_rank = lowest_rank_name
)

# NA in dyd is literally "NA"
dyd[is.na(dyd$HC_dyad),]$HC_dyad = "NA"

# Combine data
dyd = mutate(
  dyd,
  Predicted_activity = ifelse(HC_dyad == "HC", "Active", "Inactive"),
  p10_domain = ifelse(UniProt_ID %in% p10$UniProt_ID, "Yes", "No")
)

len_tax_dyd = inner_join(inner_join(tax, len), dyd)

# Add paracaspase information
len_tax_dyd = mutate(
  len_tax_dyd, Paracaspase = ifelse(UniProt_ID %in% para, "Yes", "No")
)

# Rename Lowest_rank to Organism
len_tax_dyd = rename(len_tax_dyd, Organism = Lowest_rank)

# Add p10 domain Cys Asp (CD) information instead of Yes/No
len_tax_dyd = rbind(
  mutate(filter(len_tax_dyd, p10_domain == "No"), p10_domain = ""),
  inner_join(
    filter(len_tax_dyd, p10_domain == "Yes"), p10
  ) %>% select(-p10_domain) %>% rename(p10_domain = p10_CD)
)

# Add metacaspase and architecture information
tb = inner_join(len_tax_dyd, meta)

# Change Metacaspase NA to No
tb = mutate(tb, Metacaspase = ifelse(is.na(Metacaspase), "No", Metacaspase))

# Add acidic pocket classification
tb = tb %>%
  left_join(acid) %>%
  mutate(Acidic_pocket = replace_na(Acidic_pocket, "")) %>%
  # Classify as caspase if Type I metacaspase in Metazoa and not acidic pocket
  mutate(
    Caspase = ifelse(
      Metacaspase != "No" & Acidic_pocket != "DD" & Kingdom == "Metazoa",
      "Yes", "No"
    )
  )

# Add Pfam architecture
tb = tb %>%
  left_join(pfam) %>%
  mutate(Pfam_Architecture = replace_na(Pfam_Architecture, ""))

# Add basic pocket annotation
tb = tb %>%
  left_join(bp20) %>%
  left_join(bp10) %>%
  mutate(
    Basic_pocket = paste(replace_na(bp20, ""), replace_na(bp10, ""), sep="")
  )

# Re-order columns
tb = select(
  tb,
  UniProt_ID, HC_dyad, Predicted_activity, p20_length,
  Acidic_pocket, Basic_pocket,
  p10_domain, Interdomain_dist, n_p20, n_p10,
  Architecture, Length, Metacaspase, Paracaspase, Caspase,
  Organism, Group, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus,
  Species, Tax_ID, Pfam_Architecture
)

# Write to outfile
outfile = "results/uniprot.p20.classification_and_taxonomy.tab"
write_tsv(tb, outfile)
