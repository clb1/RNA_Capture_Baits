## Computational Design of RNA Baits for Capturing Full-length mRNA Isoforms
The latest generation of DNA sequencing technology (e.g. PacBio, Oxford Nanopore) are capable of sequencing the complete structure mRNA isoforms. Given the low throughput of these technologies and the diversity and copy number skew of mRNA isoforms that exist in a cell, it is not yet feasible to use these technologies for targeted analysis of particular subpopulations of isoforms.

Herein is included an algorithm and software for designing capture baits for any subpopulation of genes and isoforms. The algorithm biases the baits for capturing full length molecules (minimizing capture of molecular fragments), optimizes the capture of on-target vs off-target molecules, and incorporates exhaustive RNA thermodynamic calculations to maximize the hybridization efficiency of the baits selected. Importantly, the approach embodied in herein is for direct capture of mRNA molecules and not cDNA.