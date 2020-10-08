# Load libraries
library(tidyverse)

# Define infile and outfiles
infile = "results/uniprot.p20.filtered.domains.tab"
out_para = "results/uniprot.p20.filtered.paracaspase.txt"
out_p10 = "results/uniprot.p20.filtered.p10.txt"

# Load data
dom = read_tsv(infile, col_names=F)
colnames(dom) = c("seqid", "domain", "start", "end")

# Identify sequences with a p10 domain
write(unique(filter(dom, startsWith(domain, "p10"))$seqid), out_p10)

# Identify paracaspases (at least one ig, Ig_2, or Ig_3 domain prior to p20)
para = dom %>% filter(!startsWith(domain, "p10")) %>% select(-end) %>%
  mutate(
    domain = ifelse(domain == "p20", "p20", "ig"),
    start = ifelse(domain == "p20", start, -start)
  ) %>% group_by(seqid, domain) %>% summarise(start = abs(max(start))) %>%
  spread(domain, start) %>% filter(!is.na(ig) & p20 > ig)

write(unique(para$seqid), out_para)
