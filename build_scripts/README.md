# Building the NEASE networks from current data

These scripts rebuild NEASE's interaction networks (human and mouse) from current public resources. The
network has the same structure and file formats as the one shipped with NEASE, so the output can replace
`nease/data/` directly. Each edge also records the evidence it comes from.

## What is built

A NEASE edge says that an exon is likely to affect a protein interaction. The build creates four layers:

| Layer | Question it answers | Sources |
|---|---|---|
| Physical interactions (PPI) | Do proteins A and B bind? | BioGRID (physical), IntAct (IMEx), HuRI, BioPlex 3.0, optionally STRING (physical, experimental) |
| Domain level (joint graph) | Which domains of A and B make the contact? | Structures/models of the pair; domain-domain catalogues (computed from the PDB, 3did 2020 snapshot, PPIDM) |
| Residue level ("pdb" table) | Which exons contain interface residues? | PDB structures and confident AlphaFold DB complex models of the pair; optionally Interactome3D structures and homology models, and PDB structures of closely related proteins (other species, paralogs) |
| Motifs (ELM) | Which exons contain a binding motif, and who binds it? | ELM instances, curated ELM interactions, ELM motif-binding domains |

Evidence tiers of the domain-level network (`network.evidence` in `config.yaml`):

| Tier | Edge between domain dA of A and dB of B when |
|---|---|
| `pdb_pair` | A and B are in one PDB structure and dA, dB are in contact |
| `af_pair` | an AlphaFold DB model of A–B is confident (ipSAE ≥ 0.6, pDockQ2 ≥ 0.23) and dA, dB are in contact |
| `i3d_structure_pair` / `i3d_model_pair` | dA, dB are in contact in an Interactome3D experimental structure (also of closely related proteins) / homology model of A–B (step 15) |
| `humanppi_pair` | dA, dB are in contact in a confident RoseTTAFold2-PPI/AlphaFold2 model of the predicted interactions of Zhang et al., Science 2025 (step 18) |
| `boltz_pair` | dA, dB are in contact in a Boltz-2 model of a measured direct interaction (Ille et al., bioRxiv 2025; step 19) |
| `pdb_homolog_pair` | dA, dB are in contact in a PDB structure of proteins ≥ 95% identical to A and B (another species or a paralog; step 17) |
| `ddi_pdb` | A and B interact physically and the Pfam pair (dA, dB) is in contact in any PDB structure (any species; a self-made, current replacement of 3did) |
| `ddi_3did2020` | the same with 3did's 2020 catalogue (contained in the PPIDM data) |
| `ddi_ppidm_gold` / `_silver` / `_bronze` | the same with the PPIDM predicted domain-domain interactions |
| `*_ortholog` (mouse) | the tier for the human orthologs, transferred to mouse |

Mouse: few mouse interactions are measured directly. Human interactions are transferred through one-to-one orthologs
(Ensembl Compara), and residue positions through alignments of the orthologous proteins. These edges are flagged
as `*_ortholog` (domain level and residue level) or listed under `human_sources` (PPI).

Domains are those of NEASE's exon → domain table (Ensembl 116, `nease/data/database/Human|Mouse`), so
every domain node can be hit by an exon. Pfam families that were killed since are mapped to their successors
(Pfam-A.dead).

## Steps

```
1 download ─┬─ 2 ensembl ─┬─ 10 orthologs ─ 11 ppi ─ 12 elm ─┐
            │             └──────────────────┘                │
            ├─ 3 pdb_select ─ 4 pdb_fetch ─ 5 pdb_contacts* ─ 6 pdb_merge ─┤
            ├─ 7 af_select ─ 8 af_fetch* ─ 9 af_contacts* ─────────────────┤
            ├─ (2, 5) ─ 17 pdb_homologs (optional) ───────────────────────────┤
            └─ (2) ─ 15 interactome3d (optional) ─────────────────────────────┴─ 13 assemble ─┬─ 14 validate
                                                                                                 └─ 16 losses
                                                     (* array jobs)
```

| Step | Script | What it does | Internet | Resources (full run, rough) |
|---|---|---|---|---|
| 1 | `01_download.py` | tables: SIFTS, Pfam, PPIDM, BioGRID, IntAct, HuRI, BioPlex, STRING, ELM, UniProt (ID mapping, reference proteomes, Swiss-Prot sequences), Ensembl, AFDB metadata (~6 GB) | yes | 1–3 h |
| 2 | `02_ensembl.py` | exon coding coordinates; UniProt ↔ identical Ensembl transcripts; UniProt entries without an identical translation aligned to their transcripts (≥ 95% of the UniProt sequence identical) | | 10–20 min, 8 cores |
| 3 | `03_pdb_select.py` | PDB entries with Pfam domains or human/mouse proteins (~180,000), shards | | 2 min |
| 4 | `04_pdb_fetch.py` | download the entries (~60 GB), or use a local mirror (`pdb.mirror`) | yes | hours (or rsync) |
| 5 | `05_pdb_contacts.py` | per shard: biological assembly, residue contacts (heavy atoms ≤ 5 Å), domain contacts | | ~2 CPU-hours in total |
| 6 | `06_pdb_merge.py` | domain-domain catalogue and residue interfaces of human/mouse protein pairs | | 10 min, 32 GB |
| 7 | `07_af_select.py` | confident AFDB heterodimers (optional homodimers) of human/mouse proteins (5,791 human, 4,924 mouse) | | 1 min |
| 8 | `08_af_fetch.py` | download the selected models: single files from the AlphaFold DB website (~10 GB), or the FTP tar files that contain them (`alphafold.download: tars`; 2,066 files, ~3 TB, days) | yes | download-bound |
| 9 | `09_af_contacts.py` | interface residues of the models (pLDDT ≥ 70) | | ~3 ms per model |
| 10 | `10_orthologs.py` | human–mouse one-to-one orthologs (Ensembl Compara) and alignments of the orthologous proteins | | 1 min, 16 cores |
| 11 | `11_ppi.py` | physical interactions per species (NCBI gene IDs) | | 2 min, 4 GB |
| 12 | `12_elm.py` | motif exons and motif-domain interactions | | 1 min |
| 13 | `13_assemble.py` | NEASE files, pathway degrees and network sizes | | 1 h, 64–128 GB |
| 14 | `14_validate.py` | comparison with the network shipped with NEASE, optional benchmark | | 10 min |
| 16 | `16_diagnose_losses.py` | why residue-level rows of the shipped network are missing from the new one (per row, e.g. no identical transcript, below the contact threshold, structure excluded, no structure with both proteins) | | 10 min, 32 GB |
| 18 | `18_humanppi.py` | optional, runs before 13: interfaces of the confident models of the predicted human interactions of Zhang et al. (Science 2025; 90% or 80% precision set), read from the 67 GB model archive in one pass | | ~30 min, 16 cores |
| 19 | `19_boltz.py` | optional, runs before 13: interfaces of the Boltz-2 models of 1,394 measured direct human interactions (GitHub repository zip, ~0.2 GB, downloaded by step 1) | yes (download) | minutes |
| 17 | `17_pdb_homologs.py` | optional, runs before 13: PDB chains of other species and paralogs mapped to human/mouse proteins at ≥ 95% identity (word search + alignment), their interfaces transferred | | ~1 h, 16 cores, 96 GB |
| 15 | `15_interactome3d.py` | optional, runs before 13: Interactome3D structures and homology models of interactions (representative set: human 1.4 GB, mouse 0.3 GB), residues mapped to UniProt by alignment, contacts as for the PDB | yes (download) | minutes, 8 cores |

## Running it

1. Create the environment: `conda env create -f environment.yml`
2. Copy `config.yaml` and set `work_dir` (a large shared file system, ~150 GB) and `nease_repo` (this repository).
   The rest of the configuration has working defaults; see the comments in the file.
3. Slurm, from the `build_scripts` folder: `sbatch slurm/submit_all.sh my_config.yaml` (or `bash slurm/submit_all.sh
   my_config.yaml` on a login node). It only submits one job per step, with the dependencies, and ends.
   - Set `PARTITION`, `DOWNLOAD_PARTITION` and `NEASE_BUILD_ACTIVATE` (how to activate the environment)
     at the top of `slurm/submit_all.sh` or as environment variables.
   - Steps 1, 4, 8 and the download of step 15 need internet access. If your compute nodes have none, run them on a
     transfer node (`python 01_download.py --config my_config.yaml`, ..., `python 15_interactome3d.py --config
     my_config.yaml --download-only`) and submit the other steps by number:
     `bash slurm/submit_all.sh my_config.yaml 2 3 5 6 7 9 10 11 12 15 17 18 19 13 14 16`.
   - Array tasks can be rerun individually: `sbatch --array=17 ...`, or `python 05_pdb_contacts.py --config my_config.yaml --task 17`.
   - Logs: `<work_dir>/logs/`; the log of the submission job is `nease_submit_<job ID>.out`.
4. Without Slurm: `JOBS=16 bash run_local.sh my_config.yaml` (any subset of steps can be given).

Every step reads its inputs from `work_dir` and can be rerun on its own, e.g. to try other settings:
`python 13_assemble.py --config my_config.yaml --name structures_only --evidence pdb_pair af_pair ddi_pdb`.

For a quick test run, `03_pdb_select.py --max-entries 400` and `07_af_select.py --max-tars 2` restrict the
structures to a sample, and fewer shards/tasks (`pdb.shards`, `alphafold.tasks`) keep the runs small.

## Output

`<work_dir>/output/<name>/` (default name `nease_data`):

- `nease_data/` in the layout of `nease/data/`:
  `network/graph_*.pkl`, `network/PPI_*.pkl`, `database/elm*`, `database/ELM_interactions*`, `database/pdb*`,
  `pathways/pathways_*` (the NEASE pathway tables with recomputed degrees)
- `nease_data/network/network_size_<species>.json`: the network sizes (edges of the PPI/DDI and structural
  networks) for NEASE's enrichment test, computed exactly as NEASE builds these networks
- `network_info_<species>.json`: network sizes, evidence used, counts per evidence
- `details/<species>/`: all tables with their evidence (`joint_edges`, `ppi`, `pdb`, `elm`, `elm_interactions`)
- `validation.md` (step 14), `losses.md` and `losses_<species>.tsv.gz` (step 16), `downloads_manifest.tsv`
  (sources and download dates)

To use a network in NEASE (this repository's version, which reads `network_size_<species>.json`), point NEASE to it
before importing it:

```bash
export NEASE_DATA_DIR=<work_dir>/output/nease_data/nease_data
python my_analysis.py        # prints "NEASE: using the network data in ..." on import
```

Files that are not in that folder (the exon tables, gene sets) come from NEASE's own data. Alternatively copy
`nease_data/*` over `nease/data/`; always copy the network files and `network_size_<species>.json` together. NEASE
uses the size file only if it belongs to the loaded joint graph (same number of edges) and otherwise warns and
falls back to the sizes of the shipped network, which are hard-coded in `builtin_network_size()`
(`nease/functions.py`). `enrich_msigdb()` computes the sizes from the loaded network in any case. Enrichment results
change with the network, because the background of the test changes, so keep the network version with your results.

## Validation and benchmark

`14_validate.py` compares the new network with the shipped one: how many shipped edges are recovered and by
which evidence, layer sizes, network sizes and pathway degree correlation.

With `--benchmark file.tsv` it also tests both networks against measured effects of exon skipping on interactions,
for example from isoform-resolved interaction screens (Yang et al., Cell 2016). The table needs the columns
`species, gene (Ensembl gene ID), start, end (exon), partner (NCBI gene ID), lost (1 lost / 0 kept)`, and it has to be
prepared from such a study's supplementary tables. A useful network predicts lost interactions more often than
kept ones.

## Sources and licences

| Source | Licence / terms |
|---|---|
| PDB, SIFTS (EBI) | CC0 / free |
| Pfam (InterPro) | CC0 |
| AlphaFold DB complexes | CC BY 4.0 |
| Interactome3D | no licence stated on the site; cite Mosca et al., Nat Methods 2013, and ask the authors before redistributing |
| BioGRID | MIT |
| IntAct | CC BY 4.0 |
| HuRI, BioPlex | free for research; cite the papers |
| STRING | CC BY 4.0 |
| PPIDM (Zenodo) | CC BY 4.0 |
| humanPPI, Zhang et al. 2025 (Dryad doi:10.5061/dryad.15dv41p84) | CC0; cite the paper |
| Boltz-2 models, Ille et al. 2025 (github.com/structural-interactome/human-interactome) | MIT; cite the preprint |
| ELM | free for academic use; check the terms before redistributing it (e.g. in the pip package) |
| UniProt, Ensembl | CC BY 4.0 / free |

## Status and test notes

These scripts are a first version. Each step was run on real data on one machine (October 2026): all
downloads, the Ensembl, ortholog, PPI, motif, assembly and validation steps at full size, and the structure
steps on samples (400 random PDB entries, 2 AlphaFold tar files). NEASE loads the resulting files unchanged,
and `run`, `enrich`, `path_analysis` and `Vis_path` work with them. A full-size run of the structure steps
on a cluster has not been done yet.

Numbers from the test run (structure layers incomplete because of the samples):

| | Human | Mouse |
|---|---|---|
| Physical interactions | 983,408 (BioGRID, IntAct, BioPlex) | 80,842 measured + 789,888 from human orthologs |
| Joint graph, gene level (shipped) | 96,929 edges (52,467) | 90,095 edges (13,348) |
| Shipped gene-level edges recovered | 65% | 73% |
| Motif-domain interactions (shipped) | 11,171 (4,846; 84% recovered) | 9,700 (257) |
| Pathway degree correlation with shipped (Spearman) | 0.98 | 0.96 |

Things to decide or watch:

- **Mouse orthologs**: with `ppi.mouse_from_human: true` (default) most mouse interactions come from human,
  which makes the mouse network about as large as the human one. Set it to `false` for a mouse network of
  measured mouse interactions only. Either way the edges are flagged.
- **PPIDM gold** contributes the most domain-level edges, and they are predicted. For a stricter network, leave it
  out of `network.evidence` (`13_assemble.py --name strict --evidence pdb_pair af_pair ddi_pdb ddi_3did2020`).
- **AlphaFold download**: step 8 downloads the confident human/mouse models one by one from the AlphaFold DB website
  (`alphafold.download: models`, ~10 GB). If the website refuses the downloads, the step stops and says so; the
  alternative, the FTP tar files (`alphafold.download: tars`), is ~3 TB because the models are spread over ~2,000 tar
  files of 0.1–11 GB (two of them, shards 1037 and 1050, are missing on the server). Set `alphafold.enabled: false`
  to skip AlphaFold.
- **HuRI**: interactome-atlas.org currently serves an invalid HTTPS certificate, so the download fails. The step
  is optional, because HuRI is also part of IntAct.
- **UniProt entries without an identical Ensembl translation** (e.g. one residue differs, or UniProt's canonical
  isoform is not an Ensembl translation; KCNN2's Q9H2S1 is an example) are aligned to their transcripts in step 2 and
  used when ≥ `uniprot.min_identity` % (95) of the UniProt sequence is identical. This affects every residue-level
  source (PDB, AlphaFold, Interactome3D, ELM) and the ortholog transfer.
- **UniProt isoform entries** (Q674X7-2): Ensembl links them to transcripts without identity values; step 2 compares
  the sequences (Swiss-Prot isoforms, `uniprot_sprot_varsplic.fasta.gz`) and uses the identical ones like other
  entries, so isoform-specific models and motifs are placed on the isoform's own exons.
- **Boltz-2 models** (`boltz*`, step 19): 1,394 measured direct interactions from IntAct (pull-down with purified
  proteins), about half with a confident interface. The model files repeat a metric table that strict mmCIF parsers
  reject; step 19 drops the second copy.
- **Structures of related proteins** (`pdb_homolog*`, step 17): a complex solved for the rat or cow proteins, or for
  a near-identical paralog, is used for the human/mouse proteins at ≥ 95% identity. These are physical structures,
  but not of the proteins themselves (several human paralogs can receive the same interface, e.g. actins, histones):
  leave the tiers out of `network.evidence` / `interface_evidence` if you want structures of the proteins themselves only.
- **humanPPI predictions** (`humanppi*`, step 18): the two files come from Dryad (doi:10.5061/dryad.15dv41p84), which
  does not allow scripted downloads without an API token; set their paths in `humanppi:`. About 80% of the predicted
  pairs were screened because they were already in BioGRID/STRING/UniProt, so the tier mostly adds predicted
  interfaces to known interactions; `humanppi.require_known: true` keeps only those. The step needs no internet and
  only the UniProt sequences of step 1, so it can run where the archive is and its output
  (`intermediate/humanppi/interfaces.parquet`, a few MB) can be copied into the build.
- **Which structures count as interactions** (`network.structure_implies_ppi`): a list of tiers whose
  structures/models make a protein pair a physical interaction. The default leaves out `humanppi`: its pairs then
  give domain edges (`humanppi_pair`) and residue interfaces, but enter the PPI network and the catalogue tiers only
  when another source has them.
- **Network sizes**: NEASE reads them from `network/network_size_<species>.json` next to the network (see Output);
  without that file it uses the hard-coded sizes of the shipped network.
