# 1) HMMsearch for p20 and p10:
hmmsearch --cpu 16 --tblout results/uniprot.p20.hmmsearch.out \
data/hmms/Pfam_p20.hmm data/uniprot/uniprot.fasta

hmmsearch --cpu 16 --tblout results/uniprot.p10.hmmsearch.out \
data/hmms/p10.hmm data/uniprot/uniprot.fasta


# 2) Extract unique significant hits (filters by DOMAIN E-VALUE < 0.0001):
source/filter_hmmer_output.R results/uniprot.p20.hmmsearch.out \
results/uniprot.p20.hmmsearch.significant_IDs.txt

source/filter_hmmer_output.R results/uniprot.p10.hmmsearch.out \
results/uniprot.p10.hmmsearch.significant_IDs.txt


# 3) Extraction of p20-containing sequences:
source/filter_fasta_by_id.py data/uniprot/uniprot.fasta \
results/uniprot.p20.hmmsearch.significant_IDs.txt \
results/uniprot.p20.fasta


# 4) Alignment to the p20 HMM:
hmmalign --trim data/hmms/Pfam_p20.hmm results/uniprot.p20.fasta | \
seqmagick convert --input-format stockholm - results/uniprot.p20_only.ali.fasta


# 5) Check position of catalytic dyad (H and C) in alignment:
seqmagick convert --output-format fasta --sample 1000 \
results/uniprot.p20_only.ali.fasta - | belvu -


# 6) Cut down p20 alignment so that it only includes catalytic dyad:
seqmagick convert --cut 793:1207 --drop 2:414 \
results/uniprot.p20_only.ali.fasta \
results/uniprot.p20_only.CH_dyad_only.ali.fasta


# 7) Create a table with dyad classification for all sequences:
cut -f 1 -d \  results/uniprot.p20_only.CH_dyad_only.ali.fasta | \
tr "\n" "\t" | tr ">" "\n" | sed -e 's/$/\n/' | sed -e 's/\t$//' | \
grep -v "^$" > results/uniprot.p20.CH_dyad_classification.tab


# 8) Select only p20 sequences without catalytic dyad gaps:
cat results/uniprot.p20.CH_dyad_classification.tab | \
grep -vP "\t(\-.)|(.\-)$" | cut -f 1 > results/uniprot.p20.complete_dyad.txt


# 9) Measure length of each p20 sequence:
seqmagick convert --ungap --include-from-file \
results/uniprot.p20.complete_dyad.txt --output-format fasta \
results/uniprot.p20_only.ali.fasta - | source/lengths_of_sequences.py > \
results/uniprot.p20_only.complete_dyad.seq_lengths.tab


# 10) Filter out sequences whose length is not within two standard deviations of
# the mean, rounding to nearest whole amino acid inclusively:
source/filter_p20_lengths.R \
results/uniprot.p20_only.complete_dyad.seq_lengths.tab \
results/uniprot.p20.filtered.txt


# 11) Extract the filtered sequences
source/filter_fasta_by_id.py results/uniprot.p20.fasta \
results/uniprot.p20.filtered.txt results/uniprot.p20.filtered.fasta


# 12) Identification of Ig and p10 domains.

hmmsearch --domT 21.8 --domtblout results/uniprot.p20.filtered.ig.domtblout \
data/hmms/ig.hmm results/uniprot.p20.filtered.fasta

hmmsearch --domT 27.0 --domtblout results/uniprot.p20.filtered.Ig_2.domtblout \
data/hmms/Ig_2.hmm results/uniprot.p20.filtered.fasta

hmmsearch --domT 30.0 --domtblout results/uniprot.p20.filtered.Ig_3.domtblout \
data/hmms/Ig_3.hmm results/uniprot.p20.filtered.fasta

hmmsearch --domT 15.0 --domtblout results/uniprot.p20.filtered.p10.domtblout \
data/hmms/Pfam_p10.hmm results/uniprot.p20.filtered.fasta

hmmsearch --domtblout results/uniprot.p20.filtered.p20.domtblout \
data/hmms/Pfam_p20.hmm results/uniprot.p20.filtered.fasta


# 13) Obtain taxonomy IDs and taxonomy information:
(grep ">" results/uniprot.p20.filtered.fasta | tr " " "\n" | \
grep -P "^>|^OX=" | sed -e 's/^OX=//' | tr ">" "&" | tr "\n" "\t" | \
tr "&" "\n" | sed -e 's/\t$//' | grep -v "^$") > \
results/uniprot.p20.filtered.taxids_from_fasta.tab

# Getting full taxonomy information for the taxonomy IDs:
source/taxid-to-taxonomy.py \
-i results/uniprot.p20.filtered.taxids_from_fasta.tab \
-n data/ncbi/taxonomy/names.dmp \
-d data/ncbi/taxonomy/nodes.dmp \
-o results/uniprot.p20.filtered.taxonomy.tab


# 14) Getting the length of each sequence:
cat results/uniprot.p20.filtered.fasta | \
source/lengths_of_sequences.py > results/uniprot.p20.filtered.seq_lengths.tab


# 15) Filter p20 domains with minmax

# Creating a table with the scores for each sequence:
cat results/uniprot.p20.filtered.p20.domtblout | grep -v "^#" | \
sed -e 's/ \+/\t/g' | cut -f 1,14,18,19 > \
results/uniprot.p20.filtered.p20_scores.tab

# Filter
Rscript source/filter_p20_domains.R

# Create new domains list
(cat results/uniprot.p20.filtered.*.domtblout | grep -v "^#" | \
sed -e 's/ \+/\t/g' | cut -f 1,4,18,19 | grep -P "\t[Ii]g"; \
cat results/uniprot.p20.filtered.p10.domtblout | grep -v "^#" | \
sed -e 's/ \+/\t/g' | cut -f 1,4,18,19 | sed -e 's/PF00656_seed.p20.ali/p20/' \
-e 's/PF00656_seed.p10.ali/p10/'; \
cat results/uniprot.p20.filtered.minmax_domains.tab | tail -n +2 | \
sed -e 's/\t/\tp20\t/' | cut -f 1,2,4,5) > results/uniprot.p20.domains.tab


# 16) Perform domain architecture analysis
Rscript source/para_meta_p10_classification.R


# 17) Align to the p10 HMM:
seqmagick convert --output-format fasta --include-from-file \
results/uniprot.p20.filtered.p10.txt results/uniprot.p20.filtered.fasta - | \
hmmalign data/hmms/Pfam_p10.hmm - | \
seqmagick convert --input-format stockholm - \
results/uniprot.p20.filtered.p10_ali.fasta

# Checking the positions of the Cys and Asp residues for a
# random subset of 500 sequences:
seqmagick convert --output-format fasta --sample 500 \
results/uniprot.p20.filtered.p10_ali.fasta - | belvu -

# Extracting those two positions:
seqmagick convert --cut 3409:3438 --drop 2:29 \
results/uniprot.p20.filtered.p10_ali.fasta \
results/uniprot.p20.filtered.p10_CD_only.p10_ali.fasta

# Counting the number of conserved sites and the number of replacements:
grep -v ">" results/uniprot.p20.filtered.p10_CD_only.p10_ali.fasta | sort | \
uniq -c | sort -rn > results/uniprot.p20.p10_CD_count.txt

# Creating a table with a classification for all sequences:
cut -f 1 -d \  results/uniprot.p20.filtered.p10_CD_only.p10_ali.fasta | \
tr "\n" "\t" | tr ">" "\n" | sed -e 's/$/\n/' | sed -e 's/\t$//' | \
grep -v "^$" > results/uniprot.p20.filtered.p10_CD_classification.tab


# 18) Check position of acidic pocket dyad (D and D) in alignment:
seqmagick convert --output-format fasta --sample 1000 \
results/uniprot.p20_only.ali.fasta - | belvu -

# Cut down p20 alignment so that it only includes acidic dyad:
seqmagick convert --cut 127:1200 --drop 2:1073 \
results/uniprot.p20_only.ali.fasta \
results/uniprot.p20_only.DD_dyad_only.ali.fasta

# Create a table with DD dyad classification for all sequences:
cut -f 1 -d \  results/uniprot.p20_only.DD_dyad_only.ali.fasta | \
tr "\n" "\t" | tr ">" "\n" | sed -e 's/$/\n/' | sed -e 's/\t$//' | \
grep -v "^$" > results/uniprot.p20.DD_dyad_classification.tab


# 19) Identify Pfam domain architecture

# Run hmmsearch versus Pfam
hmmsearch --cut_tc --cpu 16 --noali --domtblout results/uniprot.p20.pfam.domtblout \
data/Pfam-A.hmm.gz results/uniprot.p20.filtered.fasta

# Get sequence IDs, i-Eval, and the positions of different domains:
cat results/uniprot.p20.pfam.domtblout | grep -v "^#" | \
sed -e 's/ \+/\t/g' | cut -f 1,4,13,18,19 > \
results/uniprot.p20.pfam.domains.tab

# Compute domain architectures
Rscript source/domain_architectures.R


# 20) Align to caspase-1 and p20 and p10 HMMs
hmmalign --trim data/hmms/Pfam_p20.hmm \
<(cat data/caspase-1.fasta results/uniprot.p20.fasta) | \
seqmagick convert --input-format stockholm - results/uniprot.p20_cas1.ali.fasta

hmmalign --trim data/hmms/Pfam_p10.hmm \
<(cat data/caspase-1.fasta results/uniprot.p20.fasta) | \
seqmagick convert --input-format stockholm - results/uniprot.p10_cas1.ali.fasta

# Check position of basic pocket triad (R, Q, and R) in alignment:
seqmagick convert --output-format fasta --head 1000 \
results/uniprot.p20_cas1.ali.fasta - | belvu -

seqmagick convert --output-format fasta --head 1000 \
results/uniprot.p10_cas1.ali.fasta - | belvu -

# Cut down p20 and p10 alignments so that they only includes basic pocket triad:
seqmagick convert --cut 87:1200 --drop 2:1113 \
results/uniprot.p20_cas1.ali.fasta \
results/uniprot.p20_cas1.basic_pocket.ali.fasta

seqmagick convert --cut 43:43 \
results/uniprot.p10_cas1.ali.fasta \
results/uniprot.p10_cas1.basic_pocket.ali.fasta

# Create tables with triad classification for all sequences:
cut -f 1 -d \  results/uniprot.p20_cas1.basic_pocket.ali.fasta | \
tr "\n" "\t" | tr ">" "\n" | sed -e 's/$/\n/' | sed -e 's/\t$//' | \
grep -v "^$" > results/uniprot.p20.basic_pocket.tab

cut -f 1 -d \  results/uniprot.p10_cas1.basic_pocket.ali.fasta | \
tr "\n" "\t" | tr ">" "\n" | sed -e 's/$/\n/' | sed -e 's/\t$//' | \
grep -v "^$" > results/uniprot.p10.basic_pocket.tab


# FINAL STEP) Create the big table:
Rscript source/combine_info.R

# Gzip all results files that are not used in the final step
(
  cat source/combine_info.R | grep "results/" | cut -f 2 -d \"
  ls results/*
) | sort | uniq -u | parallel --no-notice --jobs 16 'gzip {}'
