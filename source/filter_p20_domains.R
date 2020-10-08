library(tidyverse)

# Load data
infile = "results/uniprot.p20.filtered.p20_scores.tab"
tb = read_tsv(infile, col_names=c("UniProt_ID", "Score", "Start", "End"))

# Store number of p20
n_p20 = length(unique(tb$UniProt_ID))

# Find min-max Score
max_score = tb %>% group_by(UniProt_ID) %>% summarise(max_score = max(Score))
min_max = min(max_score$max_score)

# Filter out domains below min_max score
tb = filter(tb, Score >= min_max)

# Ensure that each sequence has at least one domain
n_p20 <= nrow(tb)

# Write the filtered data to outfile
outfile = str_replace(infile, "p20_scores", "minmax_domains")
write_tsv(tb, outfile)
