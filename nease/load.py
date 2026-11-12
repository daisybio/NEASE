import pandas as pd
import pickle
import json
import os
import warnings


# Silence SettingWithCopyWarning on pandas < 3. The option is deprecated (and has
# no effect) since pandas 3.0, where Copy-on-Write removed that warning.
if int(pd.__version__.split('.')[0]) < 3:
    pd.options.mode.chained_assignment = None  # default='warn'


here=os.path.dirname(__file__)

# A rebuilt network (e.g. the nease_data folder written by build_scripts/13_assemble.py) can be used instead of
# the built-in one by setting the environment variable NEASE_DATA_DIR to its folder before importing nease.
# The folder has the layout of nease/data; files missing there are taken from the built-in data.
data_dir=os.environ.get('NEASE_DATA_DIR')
if data_dir:
    data_dir=os.path.abspath(os.path.expanduser(data_dir))
    if not os.path.isdir(data_dir):
        raise FileNotFoundError('NEASE_DATA_DIR is not a folder: '+data_dir)
    print('NEASE: using the network data in '+data_dir)


def data_path(relative):
    '''Path of a data file: in NEASE_DATA_DIR if it is there, otherwise in the built-in data.'''
    if data_dir and os.path.exists(os.path.join(data_dir, relative)):
        return os.path.join(data_dir, relative)
    return os.path.join(here, 'data', relative)


# Helper function

def load_obj(data_folder):
    data_folder=os.path.join(os.path.dirname(__file__), data_folder)
    with open(data_folder + '.pkl', 'rb') as f:
        return pickle.load(f)


def load_pickle(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


# Save data as dictionary
# Databases
database_mapping={}
# Gene-level ID table (Ensembl gene ID, gene name, NCBI gene ID) used for ID conversions
gene_mapping={}
Pathways={}
# Join graph
network={}
# The PPI
PPI={}

elm={}
elm_interactions={}

pdb={}


# Human
# exon -> Pfam domain table from Ensembl 116 (gzip-compressed, see scripts/build_exon_domain_table.py)
database_mapping['Human']= pd.read_pickle(data_path("database/Human"), compression='gzip')
gene_mapping['Human']= database_mapping['Human'][['Gene stable ID','Gene name','NCBI gene ID']].drop_duplicates()
Pathways['Human']= pd.read_pickle(data_path("pathways/pathways_human"))
network['Human']=load_pickle(data_path('network/graph_human.pkl'))
PPI['Human']=load_pickle(data_path('network/PPI_Human.pkl'))
elm['Human']= pd.read_pickle(data_path("database/elm"))
elm_interactions['Human']= pd.read_pickle(data_path("database/ELM_interactions"))
pdb['Human']= pd.read_pickle(data_path("database/pdb"))
non_coding=load_pickle(data_path('database/non_coding.pkl'))

# exon -> Pfam domain table from Ensembl 116 (gzip-compressed, see scripts/build_exon_domain_table.py)
database_mapping['Mouse']= pd.read_pickle(data_path("database/Mouse"), compression='gzip')
gene_mapping['Mouse']= database_mapping['Mouse'][['Gene stable ID','Gene name','NCBI gene ID']].drop_duplicates()
Pathways['Mouse']= pd.read_pickle(data_path("pathways/pathways_mouse"))
network['Mouse']=load_pickle(data_path('network/graph_mouse.pkl'))
PPI['Mouse']=load_pickle(data_path('network/PPI_mouse.pkl'))
elm['Mouse']= pd.read_pickle(data_path("database/elm_mouse"))
elm_interactions['Mouse']= pd.read_pickle(data_path("database/ELM_interactions_mouse"))
pdb['Mouse']= pd.read_pickle(data_path("database/pdb_mouse"))


# Network sizes (number of edges) for the enrichment test of the built-in pathway tables. A rebuilt network comes
# with the file network/network_size_<organism>.json (13_assemble.py); without it, the sizes of the built-in
# network are used (functions.builtin_network_size). The file is only used if it belongs to the loaded network.
network_sizes={}
for organism in ['Human', 'Mouse']:
    path=data_path('network/network_size_'+organism.lower()+'.json')
    if not os.path.exists(path):
        continue
    with open(path) as f:
        info=json.load(f)
    if info.get('joint_graph_edges')==network[organism].number_of_edges():
        network_sizes[organism]=info['network_edges']
    else:
        warnings.warn(path+' does not belong to the loaded '+organism+' network (different number of edges); '
                      'the network sizes of the built-in network are used. Copy the network files and this file together.')
