options(width=130)

# Load libraries
library(tidyverse)
library(foreach)
library(doMC)
registerDoMC(16)

# Define infiles
pfam_file = "results/uniprot.p20.pfam.domains.tab"

# Load data
pfam = read_tsv(
  pfam_file, col_names=c("UniProt_ID", "Domain", "iEval", "Start", "End")
)

# pfam %>% group_by(UniProt_ID) %>% summarise(Count = length(UniProt_ID)) %>% arrange(-Count)

arch = foreach(x=unique(pfam$UniProt_ID)) %dopar% {
  psub = pfam %>%
    filter(UniProt_ID == x) %>%
    apply(1, function(x){
      # Display all the positions inside the proteins
      tibble(Domain = x[2], iEval = as.numeric(x[3]), Position = x[4]:x[5])
    }) %>%
    bind_rows() %>%
    # Select lowest iEval domain for each position, without ties
    group_by(Position) %>%
    slice_min(order_by=iEval, n=1, with_ties=F) %>%
    ungroup() %>%
    arrange(Position) %>%
    select(-iEval)

  # Determine ends of domains
  ends = c(which(
    diff(psub$Position) != 1 |
    psub$Domain[2:length(psub$Domain)] != psub$Domain[1:(length(psub$Domain)-1)]
  ), length(psub$Domain))

  # Determine beginnings of domains
  begs = c(1, ends[1:(length(ends)-1)]+1)

  # Reduce to beginning and end positions
  psub[c(rbind(begs, ends)),] %>%
    # Add domain order and state whether position is start or end
    mutate(
      Order = rep(1:length(begs), each=2),
      Location = rep(c("Start", "End"), length(begs))
    ) %>%
    # Spread our Start and End
    spread(Location, Position) %>%
    # Remove domains shorter than 7 amino acids
    filter(End - Start + 1 > 6) %>%
    # Arrange by start position
    arrange(Start) %>%
    # Create architecture
    mutate(
      dom_str = paste(Domain,"[", paste(Start, End, sep=":"), "]", sep="")
    ) %>%
    summarise(Pfam_Architecture = paste(dom_str, collapse="-")) %>%
    # Add protein ID
    mutate(UniProt_ID = x)

} %>% bind_rows()

# filter(arch, !grepl("Peptidase_C14", Pfam_Architecture))

# Save architectures
write_tsv(arch, "results/uniprot.p20.pfam_architecture.tab")
