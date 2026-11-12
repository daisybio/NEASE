# Shipped residue-level interactions missing from /nfs/roberts/scratch/pi_mjh24/raa79/NEASE/build_scripts/work_nease/output/nease_data

See 16_diagnose_losses.py for the reasons. Per-row tables: losses_<species>.tsv.gz

## Human

- Shipped residue rows (exon, partner): 49835; kept 41308 (82.9%), lost 8527
- Of these, rows whose gene pair is in the shipped PPI network: 47785; kept 41180 (86.2%)
- New residue rows: 276084

| reason | lost rows | lost rows in shipped PPI | lost gene pairs | examples (gene - partner) |
|---|---|---|---|---|
| pair_kept_other_exons | 4484 | 3907 | 2196 | F8-F8, GNAI1-RGS14, ACY1-ABHD14A-ACY1, PLEC-ITGB4 |
| no_pdb_entry_with_both | 1663 | 1663 | 813 | TNPO1-HNRNPA1L2, NPAS3-ARNT, RAP1GDS1-RHOC, ARNT-NPAS3 |
| no_uniprot_partner | 1350 | 30 | 304 | CHMP3-RNF103-CHMP3 |
| below_min_contacts | 490 | 490 | 320 | PTCH1-PTCH1, H3-3A-TONSL, GNG2-GNAI1, GNAI1-GNG2 |
| same_entry_no_contact | 248 | 248 | 181 | KPNA1-KPNA1, SSBP4-SSBP4, CAPN8-CAPN8, LATS1-LATS1 |
| not_coding_in_ensembl | 189 | 187 | 182 | PAXX-PAXX, SOD2-SOD2, FLT1-FLT1, H2AC18-CENPA |
| entry_excluded:resolution | 49 | 49 | 30 | SHMT2-BABAM2, BABAM2-SHMT2, DYNC2I2-DYNC2H1, H2BC11-CENPN |
| no_uniprot_gene | 36 | 13 | 21 | CBS-CBS, EIF1AX-RAN, H2BS1-CENPA, H2BS1-CENPN |
| entry_excluded:single chain | 18 | 18 | 5 | RALGPS1-RALGPS1, DIXDC1-DIXDC1, MAPK6-MAPK6, TLN2-TLN2 |

## Mouse

- Shipped residue rows (exon, partner): 68432; kept 12834 (18.8%), lost 55598
- Of these, rows whose gene pair is in the shipped PPI network: 2880; kept 1587 (55.1%)
- New residue rows: 230333

| reason | lost rows | lost rows in shipped PPI | lost gene pairs | examples (gene - partner) |
|---|---|---|---|---|
| no_pdb_entry_with_both | 45058 | 560 | 16759 | Myh2-Myh2, Myh13-Myh13, Polr3a-Polr3a, Sspo-Sspo |
| pair_kept_other_exons | 6516 | 630 | 1801 | Myo5a-Myo5a, Cdh23-Cdh23, Ank1-Ank1, Camk2b-Camk2b |
| same_entry_no_contact | 1142 | 46 | 224 | Myo7a-Myo7a, Anks1b-Anks1b, Myo3b-Myo3b, Kif13b-Kif13b |
| no_uniprot_gene | 1064 | 0 | 396 | Atad2b-Atad2b, Stard9-Kifc5b, Mylk4-Mylk4, Gm20547-Mst1 |
| entry_excluded:single chain | 629 | 19 | 135 | Pik3cd-Pik3cd, Scn8a-Scn8a, Cfh-Cfh, Gnal-Gnal |
| not_coding_in_ensembl | 469 | 25 | 341 | Nnt-Nnt, Septin2-Septin11, Gm42715-Ruvbl2, Gm21985-Slc12a4 |
| human_interface_not_transferred | 301 | 7 | 42 | Muc5ac-Muc5ac, Itgam-Itgam, Ly75-Ly75, Mst1-Mst1 |
| entry_excluded:resolution | 218 | 0 | 44 | Atp6v0a1-Atp6v0a1, Cct4-Cct4, Actr3-Actr3, Arpc1b-Arpc1b |
| no_uniprot_partner | 146 | 0 | 74 | Kif21a-Kif28, Grk6-Mylk4, Sptan1-Sptbn5, Musk-Hmcn2 |
| interface_no_identical_transcript | 37 | 6 | 4 | Gadl1-Gadl1, Lgals9-Lgals9, Ahr-Ahr, Sycp3-Sycp3 |
| below_min_contacts | 18 | 0 | 5 | Ankrd31-Ankrd31, Endov-Endov, Rpl13a-Rpl9, Mrpl28-Mrpl15 |
