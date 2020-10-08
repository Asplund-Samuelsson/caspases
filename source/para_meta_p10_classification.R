options(width=150)

# Load libraries
library(tidyverse)
library(foreach)
library(doMC)
registerDoMC(16)

# Define infile and outfiles
infile = "results/uniprot.p20.domains.2.tab"
in_len = "results/uniprot.p20.filtered.seq_lengths.tab"
out_para = "results/uniprot.p20.paracaspase.2.txt"
out_meta = "results/uniprot.p20.metacaspase.2.tab"

# Load data
dom = read_tsv(infile, col_names=c("seqid", "domain", "start", "end"))
len = read_tsv(in_len, col_names=c("seqid", "length"))

# Identify paracaspases (at least one ig, Ig_2, or Ig_3 domain prior to p20)
para = dom %>% filter(!startsWith("p10", domain)) %>% select(-end) %>%
  mutate(
    domain = ifelse(domain == "p20", "p20", "ig"),
    start = ifelse(domain == "p20", start, -start)
  ) %>% group_by(seqid, domain) %>% summarise(start = abs(max(start))) %>%
  spread(domain, start) %>% filter(!is.na(ig) & p20 > ig)

write(unique(para$seqid), out_para)

# Classify metacaspases
# Type I: min(p10_start - p20_end) >= 67 for all p10_start > p20_start
# Type II: min(p10_start - p20_end) <= 66 for all p10_start > p20_start
# Type III: min(p10_start) < min(p20_start)

# Interpret any p10 as p10
dom = mutate(dom, domain = ifelse(startsWith(domain, "p10"), "p10", domain))

# If the p10 domains start and end within 10 aa of one another, merge them
domp = foreach(sid=unique(filter(dom, domain == "p10")$seqid)) %dopar% {
    domp = filter(dom, seqid == sid, domain == "p10") %>% apply(1, function(x){
        # Display all the positions inside the proteins
        tibble(seqid = x[1], pos = x[3]:x[4])
      }) %>% bind_rows() %>% distinct() %>% arrange(pos)
    # Determine start and end of domains
    x = c(
      1,
      rbind(
        which(diff(domp$pos) != 1),
        which(diff(domp$pos) != 1)+1
      ),
      nrow(domp)
    )
    # Consolidate p10 identification
    tibble(
      seqid = sid,
      domain = "p10",
      start = domp$pos[x[c(T,F)]],
      end = domp$pos[x[c(F,T)]]
    )
  } %>% bind_rows()

# Replace p10 in domain table
dom = bind_rows(filter(dom, domain != "p10"), domp)

# Determine architecture and add full protein length
arch = filter(dom, domain %in% c("p20", "p10")) %>% arrange(seqid, start) %>%
  mutate(
    dom_str = paste(domain,"[", paste(start, end, sep=":"), "]", sep="")
  ) %>% group_by(seqid) %>% summarise(
    architecture = paste(dom_str, collapse="-"),
    n_p20 = sum(domain == "p20"),
    n_p10 = sum(domain == "p10")
  ) %>% inner_join(len)

meta = full_join(
  rename(
    select(filter(dom, domain == "p20"), -domain),
    p20_start = start, p20_end = end
  ),
  rename(
    select(filter(dom, domain == "p10"), -domain),
    p10_start = start, p10_end = end
  )
) %>% mutate(
  p10_p20_dist = p20_start - p10_end - 1,
  p20_p10_dist = p10_start - p20_end - 1
)

meta = meta %>% group_by(seqid) %>% summarise(
    p10_p20_dist = min(p10_p20_dist),
    p20_p10_dist = min(p20_p10_dist)
  ) %>% inner_join(arch) %>% mutate(
    Type = ifelse(p20_p10_dist >= 67, "II", "I")
  ) %>% mutate(
    Type = ifelse(p10_p20_dist >= 0, "III", Type)
  ) %>% mutate(
    Type = ifelse(is.na(Type), "", Type)
  ) %>% mutate(
    Type = ifelse(Type != "" & (n_p20 > 1 | n_p10 > 1), "Ambiguous", Type)
  )

# Create one column for the interdomain distance
meta = meta %>% mutate(
    Interdomain_dist = ifelse(
      Type %in% c("I","II"), p20_p10_dist, p10_p20_dist
    )
  ) %>% mutate(
    Interdomain_dist = ifelse(
      Type != "" & (n_p20 > 1 | n_p10 > 1), NA, Interdomain_dist
    )
  )

# Format and write metacaspase classification to outfile
meta %>% select(
    seqid, n_p20, n_p10, architecture, length, Type, Interdomain_dist
  ) %>% rename(
    UniProt_ID = seqid, Architecture = architecture,
    Length = length, Metacaspase = Type
  ) %>% write_tsv(out_meta)
