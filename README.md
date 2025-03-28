[![PyPI version](https://badge.fury.io/py/nease.svg)](https://badge.fury.io/py/nease)



# NEASE
NEASE  (Network Enrichment method for Alternative Splicing Events) a tool for functional enrichment of alternative splicing exons/events. 


## General info
The python package NEASE (Network-based Enrichment method for Alternative Splicing Events) first detects protein features affected by AS such as domains, motifs and residues. Next, NEASE uses  a protein-protein interactions integrated with domain-domain interactions, residue-level and domain-motif interactions to identify interaction partners likely affected by AS

Next, NEASE performs gene set overrepresentation analysis and identifies enriched pathways based on affected edges. Furthermore, since the statistical approach is network-based, it also prioritizes (differentially) spliced genes and finds new disease biomarkers candidates in case of aberrant splicing.


![alt text](https://i.imgur.com/5QPhUw3.png)



## Installation

To install the package from PyPI please run:

`pip install nease` 

To install the package from git:

`git clone https://github.com/louadi/NEASE.git  && cd NEASE`

`pip install .`


Enjoy your instance of NEASE




## Data input

The standard input of the package (also the recommended) is a DataFrame object with Ensembl IDs of the genes and the exon coordinates from human genome build hg38 (GRCh38). For the mouse build, we used Ensembl version 109. 

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



