[![PyPI version](https://badge.fury.io/py/nease.svg)](https://badge.fury.io/py/nease)



# NEASE
NEASE  (Network Enrichment method for Alternative Splicing Events) a tool for functional enrichment of alternative splicing exons/events. 


## General info
The python package NEASE (Network-based Enrichment method for Alternative Splicing Events) first detects protein features affected by AS such as domains, motifs and residues. Next, NEASE uses  a protein-protein interactions integrated with domain-domain interactions, residue-level and domain-motif interactions to identify interaction partners likely affected by AS

Next, NEASE performs gene set overrepresentation analysis and identifies enriched pathways based on affected edges. Furthermore, since the statistical approach is network-based, it also prioritizes (differentially) spliced genes and finds new disease biomarkers candidates in case of aberrant splicing.


![alt text](https://i.imgur.com/5QPhUw3.png)



## Installation

NEASE requires Python 3.10 or newer (tested up to Python 3.15).

To install the package from PyPI please run:

`pip install nease` 

To install the package from git:

`git clone https://github.com/louadi/NEASE.git  && cd NEASE`

`pip install .`


Enjoy your instance of NEASE




## Data input

The standard input of the package (also the recommended) is a DataFrame object with Ensembl IDs of the genes and the exon coordinates from the human genome build GRCh38 (hg38) or the mouse genome build GRCm39 (mm39). NEASE's exon-domain annotation comes from Ensembl release 116 for both organisms.

Coordinates from any Ensembl release (or GENCODE version) on these genome builds can be used: the builds have not changed, only the annotation of transcripts and exons. Exons that only exist in a much older or newer annotation may not be matched. Mouse data on the older GRCm38 (mm10) build, used by Ensembl up to release 102, has to be lifted over to GRCm39 first.

- First column  - genes IDs (Only Ensembl gene IDs can be used).
- Second column - start of the exon coordinate.
- Third column  - end of the exon coordinate.
- Fourth column - dPSI (optionally)



| Gene              |   Start   |   End     |dpsi  | 
|-------------------|-----------|-----------|------|
| ENSG00000154263   | 69314431  | 69315425  |-0.10 | 
| ENSG00000154265   | 87411893  | 87412033  | 0.13 | 


In the case of differential splicing analysis, please make sure to filter events by significance beforehand, for example using a p-value cutoff. 

The package also supports the output of multiple AS differential detection tools such as rMATs, Whippet and MAJIQ. For the newer version of MAJIQ, pleasee check [the issue](https://github.com/louadi/NEASE/issues/3). 

If you need help with your data or need to add support for another tools, open an issue or contact us.



## Main functions and examples

Please note, that all functions are annotated with dockstrings for more details.

Import NEASE package and pandas:

```python
import nease
import pandas as pd
```



### Run NEASE 

table: Data input as DataFrame object as explained in "Data input".

organism: Either 'Human' or 'Mouse'.

input_type: Either 'Standard', 'Whippet', 'rmats'or "MAJIQ".

remove_non_in_frame: Remove exons that are predicted to disturb the ORF ([Prediction source](https://vastdb.crg.eu/wiki/FAQ#:~:text=positions%20were%20analyzed.-,How%20is%20the%20impact%20on%20the%20ORF%20predicted%3F,-The%20pipeline%20to)). 

only_divisible_by_3: Remove exons not divisible by 3.

only_DDIs: Only use DDI annotations (No PDB and ELM).

 p_value_cutoff: The enrichment p-value cutoff.
 

```python
events=nease.run(table, organism='Human',input_type='Standard',
                 remove_non_in_frame=True, only_divisible_by_3=False)
```


###  General functions
Get statistics of your data.

```python
events.get_stats()
```


Get a list of all affected domains.
```python
events.get_domains()
```



Get a list of all affected linear motifs.

```python
events.get_elm()
```



Get a list of all affected residues and their interactions from the PDB.

```python
events.get_pdb()
```


List of affected interactions from domains and linear motif binding.

```python
events.get_edges()
```

###  NEASE enrichment 

The main function of NEASE

database: a list of pathway databases to run enrichment on it. 


```python
# Supported databases:
database=  ['PharmGKB', 'HumanCyc', 'Wikipathways', 'Reactome','KEGG', 'SMPDB',
            'Signalink','NetPath', 'EHMN', 'INOH','BioCarta','PID']

# Run enrichment on Reactome only
events.enrich(database=['Reactome'])
```


###  NEASE enrichment with MSigDB gene sets and KEGG

`enrich_msigdb()` works like `enrich()`, with the gene set collections of [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) (release 2026.2) and the current KEGG pathways from the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html).
The gene sets are downloaded the first time a database is used and kept in a local cache (`~/.cache/nease`, or the folder set in the `NEASE_CACHE` environment variable). NEASE computes the degree of every gene set in its networks when the gene sets are loaded.

database: list of gene set databases (all databases of the organism if not given). `events.supported_gene_sets` lists them:

| Database | Human (MSigDB collection) | Mouse (MSigDB collection) |
|---|---|---|
| Hallmark | H | MH |
| Positional | C1 | M1 |
| CGP | C2:CGP | M2:CGP |
| BioCarta | C2:CP:BIOCARTA | M2:CP:BIOCARTA |
| KEGG_MEDICUS | C2:CP:KEGG_MEDICUS | – |
| PID | C2:CP:PID | – |
| Reactome | C2:CP:REACTOME | M2:CP:REACTOME |
| WikiPathways | C2:CP:WIKIPATHWAYS | M2:CP:WIKIPATHWAYS |
| miRDB | C3:MIR:MIRDB | M3:MIRDB |
| GTRD | C3:TFT:GTRD | M3:GTRD |
| 3CA | C4:3CA | – |
| CGN | C4:CGN | – |
| CM | C4:CM | – |
| GO:BP | C5:GO:BP | M5:GO:BP |
| GO:CC | C5:GO:CC | M5:GO:CC |
| GO:MF | C5:GO:MF | M5:GO:MF |
| HPO | C5:HPO | – |
| MP:General | – | M5:MPO:MP |
| MPT:Tumor | – | M5:MPO:MPT |
| C6:Onco | C6 | – |
| C7:ImmuneSigDB | C7:IMMUNESIGDB | – |
| C7:VAX | C7:VAX | – |
| M7:Immuno | – | M7 |
| C8:Cell_type | C8 | – |
| M8:Cell_type | – | M8 |
| C9:Comp_Perturb | C9 | – |
| KEGG | KEGG REST API (hsa) | KEGG REST API (mmu) |

cutoff: the p value cutoff.

```python
# Hallmark, Reactome and KEGG
events.enrich_msigdb(database=['Hallmark', 'Reactome', 'KEGG'])

# Gene Ontology
events.enrich_msigdb(database=['GO:BP', 'GO:CC', 'GO:MF'])

# all databases
events.enrich_msigdb()
```

The results have the same format as `enrich()`; 'Source' is the gene set database. The pathway IDs are the MSigDB systematic names (e.g. M5890) or the KEGG pathway IDs (e.g. hsa04010, mmu04010) and can be used with `path_analysis()` and `Vis_path()`.

Gene symbols of MSigDB are mapped to NCBI gene IDs with the mapping of MSigDB 2026.2 bundled with NEASE. Without internet access, or to use another MSigDB release, you can pass MSigDB JSON bundles downloaded from the [MSigDB download page](https://www.gsea-msigdb.org/gsea/downloads.jsp) (human (Hs) bundles for human data, mouse (Mm) bundles for mouse data) and, for another release, the matching "NCBI_Gene_ID" chip file:

```python
events.enrich_msigdb(database=['Reactome', 'GO:BP'], gene_sets='msigdb.v2026.2.Hs.json')
events.enrich_msigdb(database=['Reactome'], gene_sets='c2.cp.reactome.v2025.1.Hs.json', chip='Human_NCBI_Gene_ID_MSigDB.v2025.1.Hs.chip')
```

The KEGG pathways are downloaded once and then read from the cache; delete the `kegg` folder of the cache to get the current KEGG pathways. NEASE prints the date of the KEGG download with the results.

MSigDB is available under the [Creative Commons Attribution 4.0 International License](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp); some collections (KEGG MEDICUS, BioCarta) have additional terms. Please cite MSigDB when using its gene sets. The KEGG REST API is provided for academic use ([KEGG terms](https://www.kegg.jp/kegg/legal.html)); please cite KEGG when using KEGG pathways.


###  Pathway specific analysis


#### Get list of genes affecting pathways and their statistics
path_id: a string representing the Pathway ID. You can find pathways id in the enrichment table results.


```python
events.path_analysis('R-HSA-112314')
```


#### Visualize a pathway in the PPI:

Generate an HTML file for visualizing the network [example](https://tender-elion-977996.netlify.app/).

path_id: a string representing the Pathway ID.

file: a string representing a local file path for the HTML file.

save_pdf: If True, also save the figure as PDF.

k: a Float for the algorithm  Fruchterman-Reingold force-directed for nodes positions, to be tuned by the user. You might need to run the following function multiple times for an optimal visualizations.
        [more details](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html).
       

```python
events.Vis_path("R-HSA-5674135",file='AS data/enrichment/',save_pdf=False,k=0.8)
```



###  Classic gene set enrichment (Gene level)

gseapy_databases: gseapy pathways databases

non_symmetrical_only: Run classical gene set enrichment for only non-symmetrical alternative exons, that are likely to cause NMD

```python
# Run on KEGG gene set
events.classic_enrich(gseapy_databases=['KEGG_2019_Human'],non_symmetrical_only=True)
```


## Annotation data

The tables that map exons to Pfam domains (`nease/data/database/Human` and `Mouse`) are built from Ensembl release 116 (human GRCh38.p14, mouse GRCm39) through Ensembl BioMart with (needs the `requests` package):

```
python scripts/build_exon_domain_table.py Human nease/data/database/Human
python scripts/build_exon_domain_table.py Mouse nease/data/database/Mouse
```

Each row is an exon of a protein-coding transcript whose coding part overlaps a Pfam domain of that transcript's protein. Earlier NEASE versions used older annotation (mouse: Ensembl 109). Run with `--mart https://feb2023.archive.ensembl.org`, the script reproduces that Ensembl 109 mouse table exactly. Because the new annotation contains more transcripts, more exons can be matched to domains, and results can differ slightly from earlier versions.

The linear motif (ELM) and residue-level (PDB) interface annotations, the list of exons that disturb the ORF, and the interaction networks are unchanged.


## Tutorials


A step-by-step guide to use NEASE is available [here](https://github.com/louadi/NEASE-tutorials).


A simple example for running NEASE on a standard input:
([Notebook](https://github.com/louadi/NEASE-tutorials/blob/main/DCM_analysis.ipynb)/[Google Colab](https://colab.research.google.com/github/louadi/NEASE-tutorials/blob/main/DCM_analysis.ipynb))



## Development


First released by: Zakaria Louadi

Maintainer: Elias Albrecht 


## Cite

If you use NEASE, please cite:

Louadi, Z., Elkjaer, M.L., Klug, M. et al. Functional enrichment of alternative splicing events with NEASE reveals insights into tissue identity and diseases. Genome Biol 22, 327 (2021). https://doi.org/10.1186/s13059-021-02538-1



