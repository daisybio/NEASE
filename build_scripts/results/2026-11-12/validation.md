# Validation of /nfs/roberts/scratch/pi_mjh24/raa79/NEASE/build_scripts/work_nease/output/nease_data

## Human

- Joint graph, domain level: 88732 shipped, 271577 new, 60643 shared (68.3% of shipped kept, Jaccard 0.20)
- Joint graph, gene level: 52467 shipped, 124012 new, 35396 shared (67.5% of shipped kept, Jaccard 0.25)
- Protein interaction network: 416722 shipped, 996631 new, 358289 shared (86.0% of shipped kept, Jaccard 0.34)
- Evidence of the shipped domain edges that are recovered: {'ddi_3did2020': 49083, 'ddi_ppidm_gold': 44363, 'ddi_pdb': 39579, 'i3d_model_pair': 4692, 'pdb_pair': 4097, 'humanppi_pair': 4026, 'pdb_homolog_pair': 3643, 'i3d_structure_pair': 3441, 'af_pair': 2241, 'boltz_pair': 226}
- Evidence of all new domain edges: {'ddi_ppidm_gold': 195821, 'ddi_3did2020': 128062, 'ddi_pdb': 121736, 'pdb_homolog_pair': 16238, 'pdb_pair': 13534, 'humanppi_pair': 13146, 'i3d_structure_pair': 10688, 'i3d_model_pair': 10237, 'af_pair': 9238, 'boltz_pair': 751}
- Motif exons (exon, ELM class): 2919 shipped, 2358 new, 1662 shared (56.9% of shipped kept, Jaccard 0.46)
- Motif-domain interactions: 4846 shipped, 11667 new, 4092 shared (84.4% of shipped kept, Jaccard 0.33)
- Residue-level interfaces (exon, partner): 49835 shipped rows, 276084 new rows; 41308 shipped rows (82.9%) have a new row with the same gene and partner and overlapping coordinates (how NEASE matches them; exon IDs differ between Ensembl releases)
- Network sizes (edges) shipped: {'Degree in the PPI/DDI': 52467, 'Degree in the structural PPI': 60235}, new: {'Degree in the PPI/DDI': 124012, 'Degree in the structural PPI': 135097}
- Network sizes NEASE uses for this network (network_size_human.json): {'Degree in the PPI/DDI': 124012, 'Degree in the structural PPI': 135097}
- Pathway degrees "Degree in the PPI/DDI": Spearman correlation shipped vs new 0.973
- Pathway degrees "Degree in the structural PPI": Spearman correlation shipped vs new 0.981

## Mouse

- Joint graph, domain level: 26558 shipped, 259804 new, 20770 shared (78.2% of shipped kept, Jaccard 0.08)
- Joint graph, gene level: 13348 shipped, 115072 new, 10111 shared (75.7% of shipped kept, Jaccard 0.09)
- Protein interaction network: 62791 shipped, 883687 new, 58136 shared (92.6% of shipped kept, Jaccard 0.07)
- Evidence of the shipped domain edges that are recovered: {'ddi_3did2020': 16555, 'ddi_ppidm_gold': 15731, 'ddi_pdb': 12166, 'pdb_pair_ortholog': 1164, 'pdb_homolog_pair': 1163, 'humanppi_pair_ortholog': 1150, 'i3d_structure_pair_ortholog': 1007, 'i3d_model_pair_ortholog': 850, 'af_pair': 618, 'i3d_model_pair': 538, 'af_pair_ortholog': 481, 'pdb_pair': 376, 'i3d_structure_pair': 297, 'boltz_pair_ortholog': 69}
- Evidence of all new domain edges: {'ddi_ppidm_gold': 190637, 'ddi_3did2020': 122550, 'ddi_pdb': 113963, 'pdb_homolog_pair': 13180, 'pdb_pair_ortholog': 10495, 'humanppi_pair_ortholog': 10242, 'i3d_structure_pair_ortholog': 8231, 'af_pair': 7832, 'i3d_model_pair_ortholog': 7716, 'af_pair_ortholog': 5833, 'pdb_pair': 2088, 'i3d_structure_pair': 1830, 'i3d_model_pair': 1636, 'boltz_pair_ortholog': 524}
- Motif exons (exon, ELM class): 366 shipped, 2326 new, 234 shared (63.9% of shipped kept, Jaccard 0.10)
- Motif-domain interactions: 257 shipped, 11888 new, 222 shared (86.4% of shipped kept, Jaccard 0.02)
- Residue-level interfaces (exon, partner): 68432 shipped rows, 230333 new rows; 12834 shipped rows (18.8%) have a new row with the same gene and partner and overlapping coordinates (how NEASE matches them; exon IDs differ between Ensembl releases)
- Network sizes (edges) shipped: {'Degree in the PPI/DDI': 13348, 'Degree in the structural PPI': 31144}, new: {'Degree in the PPI/DDI': 115072, 'Degree in the structural PPI': 124543}
- Network sizes NEASE uses for this network (network_size_mouse.json): {'Degree in the PPI/DDI': 115072, 'Degree in the structural PPI': 124543}
- Pathway degrees "Degree in the PPI/DDI": Spearman correlation shipped vs new 0.960
- Pathway degrees "Degree in the structural PPI": Spearman correlation shipped vs new 0.963
