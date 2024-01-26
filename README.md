# Summary
This script takes as input the Phylogenetic Hierarchical Orthogroups (HOGs) from orthofinder and calcuates the core/shell/cloud of the investiaged genomes. In addition a list of all protein IDs included in the orthofinder analysis must be provided. 

# requirements
pandas

# Definitions
core genes: genes present in all genomes in any ratio (N:N)
singleCopy core genes: genes present in all genomes as one copy (1:1)
shell genes: genes present all but one genome
cloud genes: genes present in only one genome, single copy or multiple copies

