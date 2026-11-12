"""Shared helpers of the NEASE network build: configuration, paths, downloads and species tables."""
import gzip
import os
from functools import lru_cache
import pickle
import sys
import time
import urllib.parse
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests
import yaml

SPECIES = {
    'Human': dict(taxid=9606, ensembl='homo_sapiens', ensembl_file='Homo_sapiens.GRCh38', uniprot='HUMAN_9606',
                  elm_taxon='Homo sapiens', proteome='UP000005640',
                  nease=dict(graph='network/graph_human.pkl', ppi='network/PPI_Human.pkl', elm='database/elm',
                             elm_interactions='database/ELM_interactions', pdb='database/pdb',
                             pathways='pathways/pathways_human', exon_domains='database/Human')),
    'Mouse': dict(taxid=10090, ensembl='mus_musculus', ensembl_file='Mus_musculus.GRCm39', uniprot='MOUSE_10090',
                  elm_taxon='Mus musculus', proteome='UP000000589',
                  nease=dict(graph='network/graph_mouse.pkl', ppi='network/PPI_mouse.pkl', elm='database/elm_mouse',
                             elm_interactions='database/ELM_interactions_mouse', pdb='database/pdb_mouse',
                             pathways='pathways/pathways_mouse', exon_domains='database/Mouse')),
}

# Local file of every downloaded source (relative to <work_dir>/downloads); {fields} are filled per species
FILES = {
    'sifts_pfam': 'sifts/pdb_pfam_mapping.tsv.gz',
    'sifts_uniprot': 'sifts/pdb_chain_uniprot.tsv.gz',
    'pfam_dead': 'pfam/Pfam-A.dead.gz',
    'ppidm': 'ppidm/PPIDMdata_2021.zip',
    'biogrid': 'biogrid/BIOGRID-ALL-LATEST.tab3.zip',
    'intact': 'intact/intact.zip',
    'huri': 'huri/HuRI.tsv',
    'bioplex_293t': 'bioplex/BioPlex_293T_Network_10K_Dec_2019.tsv',
    'bioplex_hct116': 'bioplex/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv',
    'string_physical': 'string/{taxid}.protein.physical.links.detailed.txt.gz',
    'string_aliases': 'string/{taxid}.protein.aliases.txt.gz',
    'elm_instances': 'elm/instances_{species}.tsv',
    'elm_interactions': 'elm/interactions.tsv',
    'elm_domains': 'elm/interaction_domains.tsv',
    'uniprot_idmapping': 'uniprot/{uniprot}_idmapping_selected.tab.gz',
    'uniprot_proteome': 'uniprot/{uniprot}_reference_proteome.fasta.gz',
    'uniprot_sprot': 'uniprot/uniprot_sprot.fasta.gz',
    'uniprot_varsplic': 'uniprot/uniprot_sprot_varsplic.fasta.gz',
    'boltz': 'boltz/human-interactome-main.zip',
    'ensembl_gtf': 'ensembl/{ensembl_file}.{release}.gtf.gz',
    'ensembl_pep': 'ensembl/{ensembl_file}.pep.all.fa.gz',
    'ensembl_uniprot': 'ensembl/{ensembl_file}.{release}.uniprot.tsv.gz',
    'ensembl_entrez': 'ensembl/{ensembl_file}.{release}.entrez.tsv.gz',
    'compara_homologies': 'ensembl/Compara.{release}.protein_default.homologies.mus_musculus.tsv.gz',
    'afdb_heterodimers': 'alphafold/heterodimer_metadata.csv',
    'afdb_homodimers': 'alphafold/homodimer_metadata.csv',
}

# Sources that are only used for some species (the human cell line screens are transferred to mouse via orthologs)
HUMAN_ONLY = {'huri', 'bioplex_293t', 'bioplex_hct116'}


def log(*message):
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), *message, flush=True)


def load_config(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg['work_dir'] = Path(os.path.expandvars(os.path.expanduser(cfg['work_dir'])))
    cfg['nease_repo'] = Path(os.path.expandvars(os.path.expanduser(cfg['nease_repo'])))
    if not (cfg['nease_repo'] / 'nease' / 'data').is_dir() and (cfg['nease_repo'] / 'data').is_dir():
        cfg['nease_repo'] = cfg['nease_repo'].parent  # the nease package folder was given
    if not (cfg['nease_repo'] / 'nease' / 'data').is_dir():
        sys.exit(f'nease_repo: no NEASE repository at {cfg["nease_repo"]} (nease/data not found)')
    for species in cfg['species']:
        if species not in SPECIES:
            sys.exit(f'Unknown species {species}; supported: {list(SPECIES)}')
    if 'Mouse' in cfg['species'] and 'Human' not in cfg['species']:
        sys.exit('The mouse network uses the human data (orthologs): add Human to species')
    return cfg


def paths(cfg):
    work = cfg['work_dir']
    p = dict(downloads=work / 'downloads', intermediate=work / 'intermediate', output=work / 'output',
             logs=work / 'logs', pdb=Path(cfg['pdb']['mirror']) if cfg['pdb'].get('mirror') else work / 'pdb_mmcif',
             alphafold=work / 'alphafold')
    for path in p.values():
        path.mkdir(parents=True, exist_ok=True)
    return p


def species_dir(cfg, species):
    path = paths(cfg)['intermediate'] / species
    path.mkdir(parents=True, exist_ok=True)
    return path


def pdb_file(cfg, pdb):
    '''mmCIF file of a PDB entry (divided layout of the wwPDB archive: ab/1abc.cif.gz).'''
    return paths(cfg)['pdb'] / pdb[1:3] / f'{pdb}.cif.gz'


def fields(cfg, species=None):
    out = {'release': cfg['ensembl_release']}
    if species is not None:
        out.update(SPECIES[species], species=species, elm_taxon=urllib.parse.quote(SPECIES[species]['elm_taxon']))
    return out


def per_species(name):
    return any(f'{{{key}}}' in FILES[name] for key in ('taxid', 'species', 'uniprot', 'ensembl_file'))


def input_path(cfg, name, species=None):
    '''Local path of a downloaded source.'''
    return paths(cfg)['downloads'] / FILES[name].format(**fields(cfg, species))


def source_url(cfg, name, species=None):
    return cfg['sources'][name].format(**fields(cfg, species))


def download(url, dest, retries=5, timeout=1800, force=False):
    '''Download url to dest (once; written to a temporary file first so that interrupted downloads are not used).'''
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(f'{dest.name}.{os.getpid()}.part')
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=timeout, headers={'User-Agent': 'nease-network-build'}) as r:
                r.raise_for_status()
                with open(part, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            if part.stat().st_size == 0:
                raise IOError('empty download')
            os.replace(part, dest)
            return dest
        except (requests.RequestException, IOError) as e:
            log(f'  download failed ({e}), attempt {attempt + 1}/{retries}: {url}')
            part.unlink(missing_ok=True)
            time.sleep(min(30 * (attempt + 1), 300))
    raise RuntimeError(f'download failed: {url}')


def task_id(args):
    '''Array task index: --task, otherwise $SLURM_ARRAY_TASK_ID.'''
    if args.task is not None:
        return args.task
    if 'SLURM_ARRAY_TASK_ID' in os.environ:
        return int(os.environ['SLURM_ARRAY_TASK_ID'])
    sys.exit('Give --task or run as a Slurm array job')


def object_columns(table):
    '''Text columns (and the column index) as object dtype, so that pickles written with pandas 3 load with pandas 2.'''
    table = table.astype({c: object for c in table.columns if pd.api.types.is_string_dtype(table[c])})
    table.columns = table.columns.astype(object)
    return table


def save_pickle(obj, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(obj, pd.DataFrame):
        object_columns(obj).to_pickle(path)
    else:
        with open(path, 'wb') as f:
            pickle.dump(obj, f)


def nease_file(cfg, species, key):
    return cfg['nease_repo'] / 'nease' / 'data' / SPECIES[species]['nease'][key]


def gene_domains(cfg, species, forward=None):
    '''NCBI gene ID -> Pfam domains of NEASE's exon -> domain table (the domains NEASE can find affected).'''
    table = pd.read_pickle(nease_file(cfg, species, 'exon_domains'), compression='gzip')
    table = pd.DataFrame({'entrez': table['NCBI gene ID'].astype(str), 'pfam': table['Pfam ID'].astype(str)}).drop_duplicates()
    if forward:
        table = forward_pfam(table, ['pfam'], forward).drop_duplicates()
    return table.reset_index(drop=True)


def read_uniprot_entrez(cfg, species):
    '''UniProt accession -> NCBI gene ID (one row per pair) from the UniProt ID mapping of the species.'''
    table = pd.read_csv(input_path(cfg, 'uniprot_idmapping', species), sep='\t', header=None, usecols=[0, 2],
                        names=['uniprot', 'entrez'], dtype=str).dropna()
    table['entrez'] = table['entrez'].str.split('; ')
    return table.explode('entrez').drop_duplicates().reset_index(drop=True)


def threads():
    '''CPUs to use: the Slurm allocation, otherwise the machine.'''
    return int(os.environ.get('SLURM_CPUS_PER_TASK') or os.cpu_count() or 1)


def read_uniprot_fasta(path, accessions=None):
    '''UniProt FASTA -> {accession: sequence} (only the given accessions, if any).'''
    sequences, accession, parts = {}, None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if accession and (accessions is None or accession in accessions):
                    sequences[accession] = ''.join(parts)
                accession, parts = line[1:].split('|')[1], []
            else:
                parts.append(line.strip())
    if accession and (accessions is None or accession in accessions):
        sequences[accession] = ''.join(parts)
    return sequences


def uniprot_sequences(cfg, accessions=None):
    '''
    Sequences of UniProt entries: the reference proteomes of the species of the build, Swiss-Prot of all species
    and Swiss-Prot isoforms. Returns ({accession: sequence}, set of Swiss-Prot accessions).
    '''
    sequences = {}
    for species in cfg['species']:
        sequences.update(read_uniprot_fasta(input_path(cfg, 'uniprot_proteome', species), accessions))
    swissprot = read_uniprot_fasta(input_path(cfg, 'uniprot_sprot'), accessions)
    sequences.update(swissprot)
    # Swiss-Prot isoforms (Q674X7-2); only when the accessions asked for include isoforms
    if accessions is None or any('-' in a for a in accessions):
        path = input_path(cfg, 'uniprot_varsplic')
        if path.exists():
            isoforms = read_uniprot_fasta(path, accessions)
            sequences.update(isoforms)
            swissprot.update(isoforms)
    return sequences, set(swissprot)


@lru_cache(maxsize=1)
def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'))
    aligner.open_gap_score, aligner.extend_gap_score = -11, -1
    return aligner


def matrix_letters(sequence):
    '''Letters outside BLOSUM62 (U selenocysteine, O pyrrolysine, J, ...) as X; positions are unchanged.'''
    alphabet = _aligner().substitution_matrix.alphabet
    return ''.join(c if c in alphabet else 'X' for c in sequence.upper())


def align_positions(seq_a, seq_b):
    '''Local alignment (BLOSUM62): aligned positions (1-based) of a and b and the number of identical residues.'''
    seq_a, seq_b = matrix_letters(seq_a), matrix_letters(seq_b)
    alignment = _aligner().align(seq_a, seq_b)[0]
    pos_a, pos_b, identical = [], [], 0
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        for a, b in zip(range(a0, a1), range(b0, b1)):
            pos_a.append(a + 1)
            pos_b.append(b + 1)
            identical += seq_a[a] == seq_b[b]
    return pos_a, pos_b, identical


def chain_positions(chain, sequence, min_identity):
    '''
    Residues of a model chain (gemmi) on a UniProt sequence: {index along the chain (1-based): UniProt position}.
    By the residue numbers if they fit the sequence (>= 95% of the residues), otherwise by local alignment
    (>= min_identity % identity of the aligned residues); None if neither works.
    '''
    import gemmi
    letters = ''.join((gemmi.find_tabulated_residue(r.name).one_letter_code.upper() if gemmi.find_tabulated_residue(r.name) else 'X')
                      for r in chain)
    numbers = [r.seqid.num for r in chain]
    by_number = sum(1 for n, a in zip(numbers, letters) if 0 < n <= len(sequence) and sequence[n - 1] == a)
    if by_number >= 0.95 * len(letters):
        return {i: n for i, n in enumerate(numbers, 1) if 0 < n <= len(sequence)}
    chain_pos, uniprot_pos, identical = align_positions(letters, sequence)
    if not chain_pos or 100 * identical / len(chain_pos) < min_identity:
        return None
    return dict(zip(chain_pos, uniprot_pos))


def placeable_proteins(cfg, species):
    '''
    UniProt entries whose positions can be placed on exons of the species: entries identical to an Ensembl
    translation (proteins.parquet) and entries aligned to one (uniprot_maps.parquet, step 2). One row per entry and
    gene: uniprot, gene_id, canonical, swissprot, length, sequence (UniProt numbering), identical.
    '''
    folder = species_dir(cfg, species)
    proteins = pd.read_parquet(folder / 'proteins.parquet', columns=['uniprot', 'gene_id', 'canonical', 'swissprot', 'length', 'sequence'])
    proteins = proteins.sort_values('canonical', ascending=False, kind='stable').drop_duplicates(['uniprot', 'gene_id']).assign(identical=True)
    path = folder / 'uniprot_maps.parquet'
    if path.exists():
        maps = pd.read_parquet(path, columns=['uniprot', 'gene_id', 'canonical', 'swissprot', 'sequence'])
        maps = maps[~maps['uniprot'].isin(set(proteins['uniprot']))]
        maps = maps.sort_values('canonical', ascending=False, kind='stable').drop_duplicates(['uniprot', 'gene_id'])
        proteins = pd.concat([proteins, maps.assign(length=maps['sequence'].str.len(), identical=False)], ignore_index=True)
    return proteins.reset_index(drop=True)


def pfam_forward(cfg):
    '''Killed Pfam families -> the families they were merged into or split into (Pfam-A.dead).'''
    forward, accession = {}, None
    with gzip.open(input_path(cfg, 'pfam_dead'), 'rt') as f:
        for line in f:
            if line.startswith('#=GF AC'):
                accession = line.split()[2]
            elif line.startswith('#=GF FW') and len(line.split()) > 2:
                forward.setdefault(accession, []).extend(x.strip(';') for x in line.split()[2:])
    return forward


def overlapping(rows, other, keys, start='Genomic coding start', end='Genomic coding end',
                other_start=None, other_end=None):
    '''
    For every row of rows: is there a row of other with the same keys and overlapping coordinates? This is how NEASE
    matches exons to its motif and residue tables (Ensembl gene ID and coordinate overlap; exon IDs are not used,
    so tables made with other Ensembl releases can be compared).
    '''
    other_start, other_end = other_start or start, other_end or end
    left = rows[keys + [start, end]].reset_index(drop=True).reset_index().rename(columns={start: '_s', end: '_e'})
    right = other[keys + [other_start, other_end]].drop_duplicates().rename(columns={other_start: '_os', other_end: '_oe'})
    pairs = left.merge(right, on=keys)
    hit = pairs.loc[(pairs[['_s', '_os']].max(axis=1) <= pairs[['_e', '_oe']].min(axis=1)), 'index'].unique()
    out = pd.Series(False, index=range(len(rows)))
    out[hit] = True
    return out.set_axis(rows.index)


def forward_pfam(table, columns, forward):
    '''Replace killed Pfam families in the given columns by the families they were merged or split into.'''
    for column in columns:
        table = table.assign(**{column: table[column].map(lambda x: forward.get(x, [x]))}).explode(column)
    return table.reset_index(drop=True)
