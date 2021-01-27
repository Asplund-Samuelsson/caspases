## Caspase homolog analysis

### Description
Analysis of caspase homologs from UniProt as outlined in `source/Identification.sh`.

The final output is this file (generated by `source/combine_info.R`):
```
results/uniprot.p20.classification_and_taxonomy.tab
```

The final output file contains the following information:

| Column | Description |
| ------ | ----------- |
| UniProt_ID | UniProt identifier |
| HC_dyad | Histidine-cysteine dyad residues |
| Predicted_activity | "Active" if HC_dyad == HC, otherwise "Inactive" |
| p20_length | Length of p20 domain |
| Acidic_pocket | Acidic pocket residues |
| Basic_pocket | Basic pocket residues |
| p10_domain | Critical residues in p10 domain (if present) |
| Interdomain_dist | Distance between p20 and p10 domains |
| n_p20 | Number of p20 domains detected |
| n_p10 | Number of p10 domains detected |
| Architecture | Architecture of p20 and p10 domains |
| Length | Total length of protein |
| Metacaspase | Metacaspase classification |
| Paracaspase | Paracaspase classification |
| Caspase | Caspase classification |
| Organism | Lowest NCBI taxonomy rank name |
| Group | NCBI Phylum (Class for Proteobacteria) |
| Superkingdom | NCBI taxonomy Superkingdom name |
| Kingdom | NCBI taxonomy Kingdom name |
| Phylum | NCBI taxonomy Phylum name |
| Class | NCBI taxonomy Class name |
| Order | NCBI taxonomy Order name |
| Family | NCBI taxonomy Family name |
| Genus | NCBI taxonomy Genus name |
| Species | NCBI taxonomy Species name |
| Tax_ID | NCB taxonomy identifier |
| Pfam_Architecture | Pfam domain architecture |

### Author
Johannes Asplund-Samuelsson, KTH (johannes.asplund.samuelsson@scilifelab.se)
