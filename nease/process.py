from .load import *


import numpy as np 
import networkx as nx
import csv
import re
import collections
import mygene



pd.set_option('display.max_colwidth', 1000)



# DIGGER Exon-level link: 
"""Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein interactions and domain-domain interactions
into a joint graph and maps interacting residues to exons. DIGGER allows the users to query exons or isoforms individually or as a set 
to visually explore their interactions. The following modes of DIGGER can be used interchangeably:
"""
DIGGER='https://exbio.wzw.tum.de/digger/ID/exon/'






def splitDataFrameList(df,target_column):

        ''' 
        Efficiently split Pandas Dataframe cells containing lists into multiple rows,
        duplicating the other column's values.
        Original code from https://gist.github.com/jlln/338b4b0b55bd6984f883

        df = dataframe to split,
        target_column = the column containing the values to split
        separator = the symbol used to perform the split
        returns: a dataframe with each entry for the target column separated, with each element moved into a new row. 
        The values in the other columns are duplicated across the newly divided rows.
        '''

        def splitListToRows(row,row_accumulator,target_column):
            split_row = row[target_column]
            for s in split_row:
                new_row = row.to_dict()
                new_row[target_column] = s
                row_accumulator.append(new_row)
        new_rows = []
        df.apply(splitListToRows,axis=1,args = (new_rows,target_column))
        new_df = pd.DataFrame(new_rows)
        return new_df

    
    
    
# Standard data processing 
# Gene ID     Start       End        Delta  (optionally) 


def process_standard (data,
                      mapping,
                      min_delta,
                      only_DDIs,
                      nease_data,
                      remove_non_in_frame,
                      only_divisible_by_3):
    
        # verify the input data format
        columns=data.columns

        # make sure you have at least two columns
        if len(columns)<3:
                raise ValueError('Make sure your table have at least 3 columns:    Gene ensembl ID    EXON START    EXON END    dPSI (optional)')

        else:
                genes=list(data[columns[0]])

                # Verify the gene IDs
                if not  all(x.startswith('ENS') for x in genes):
                     raise ValueError(' Could not recognize Ensembl gene ID. Please make sure that the first column corresponds to gene IDs.')

                # Verify the start and end
                try:
                        
                        # Remove non numerical 
                        if  not data[columns[1]].dtype=='int64':
                            data=data[data[columns[1]].apply(lambda x: x.isnumeric())]
                            # convert to int
                            data[columns[1]]=data[columns[1]].astype(int)
                            
                            
                        if not data[columns[2]].dtype=='int64':
                            data=data[data[columns[2]].apply(lambda x: x.isnumeric())]
                            data[columns[2]]=data[columns[2]].astype(int)
                except:
                    raise ValueError('Could not find exons coordinates. Please make sure that the second column corresponds to the exon start and the third to the exon end (hg38).')
        
        
        data['symetric']=data.apply(lambda x: abs(x[columns[1]]- x[columns[2]])%3==0, axis=1)
        
                
        # get all coding genes affected by splicing
        spliced_genes=list(data[columns[0]].unique())
        
        

        data['symetric']=data.apply(lambda x: (abs(x[columns[1]]- x[columns[2]])+1)%3==0, axis=1)
        
        # only non-symetric genes                                            
        non_symetric_genes=list(data[~data['symetric']][columns[0]].unique())
        
                
        #################
        # filtering steps
        initial=len(data)
        
        if only_divisible_by_3:
            data=data[data.symetric]
            print("You selected to remove all exons that are not divisible by 3 [(abs(start-end)+1) / 3].\nIf you would like to change this option, and include all exons, please change the parameter only_divisible_by_3 to False.")
                
        if remove_non_in_frame:
            # parse gene coord as one set 
            data['coord']=data.apply(lambda x:{x.iloc[0],x.iloc[1],x.iloc[2]}  ,axis = 1)
            # only keep coding exons
            coding=[ (e not in nease_data.non_coding)  for e in data['coord']]
            data=data[coding]
            
            
            data=data.drop(columns=[ 'coord'])
            
            print("\nYou selected to remove all exons that are predicted to disturb the ORF or known to result in a non-coding gene.\nIf you would like to change this option and include all exons, please change the parameter remove_non_in_frame to False. ")
            
        
        filtred=initial-len(data)            
        if filtred!=0:
                print('\n'+str(filtred)+' Exons out of '+str(initial)+' were removed.')
        
        

        # map to domains by calculating the overlap of exon coordinate and domain
        mapping_tb=pd.merge(data, mapping,   left_on=columns[0], right_on='Gene stable ID').drop_duplicates()   
        
        elm_affected,pdb_affected=[],[]
        # map to pdb and elm    
        if not only_DDIs:
                
                elm_affected,pdb_affected=get_interfaces(data,nease_data,min_delta,overlap_cal=True)
        


        if len(mapping_tb)==0:
            return []
        
        try:
                #try to get the delta PSI from the user input
                mapping_tb['max_change']=mapping_tb[columns[3]].astype(float)
                
                # significance filter (delta PSI)
                mapping_tb=mapping_tb[mapping_tb['max_change'].abs()>=min_delta]
        except:
                # No delta PSI given
                
                print('Delta PSI column was not found. Proceeding with all events (no filtering)')
                mapping_tb['max_change']='-'
                


        
        
        # any overlap is considered
        mapping_tb['overl']=mapping_tb[["Genomic coding start", columns[1]]].max(axis=1) <= mapping_tb[["Genomic coding end", columns[2]]].min(axis=1)
        mapping_tb=mapping_tb[mapping_tb['overl']].drop_duplicates()
        
        
        
        mapping_tb=mapping_tb[['Gene name','NCBI gene ID','Gene stable ID','Exon stable ID','Pfam ID','max_change']]
        #mapping_tb=mapping_tb.sort_values(['max_change'].abs(), ascending=False)
        try:
            mapping_tb=mapping_tb.reindex(mapping_tb['max_change'].abs().sort_values(ascending=False, kind='stable').index)
        except:
            pass
        mapping_tb=mapping_tb[mapping_tb['NCBI gene ID'].notnull()]
        mapping_tb['NCBI gene ID']=mapping_tb['NCBI gene ID'].astype('int').astype('str')
        #mapping_tb=mapping_tb.drop_duplicates(['Gene name','NCBI gene ID','Gene stable ID','Pfam ID'],keep= 'first')
        
        
        
        #mapping_tb=mapping_tb.groupby(['Gene name','NCBI gene ID','Gene stable ID','Exon stable ID','Pfam ID']).max()['max_change']

        
        return mapping_tb,spliced_genes,elm_affected,pdb_affected,non_symetric_genes



def process_spycone (spycone,
                      mapping,
                      ):
    
        '''
        # This function is used for output of spycone only
        # example of input
            domains
            64225/PF02263
            64225/PF02841
            64864/PF18326
            6188/PF07650
        
        
       '''
        
        mapping['id']=mapping['NCBI gene ID'].astype('str')+'/'+mapping['Pfam ID'].astype('str')

        # map to domains by calculating the overlap of exon coordinate and domain
        mapping_tb=pd.merge(mapping,spycone,left_on='id',right_on=spycone.columns[0]).drop_duplicates()
        
        
        # get all coding genes affected by splicing
        # Only genes with Pfam domain will be considred here
        # NCBI id used here for the network visualization
        spliced_genes=list(mapping_tb['NCBI gene ID'].unique())

        if len(mapping_tb)==0:
            return []
        

        mapping_tb=mapping_tb[['Gene name','NCBI gene ID','Gene stable ID','Exon stable ID','Pfam ID']]
         
        mapping_tb=mapping_tb[mapping_tb['NCBI gene ID'].notnull()]
        mapping_tb['NCBI gene ID']=mapping_tb['NCBI gene ID'].astype('int').astype('str')
    

        
        
      
        return mapping_tb,spliced_genes,elm_affected,pdb_affected




# Majiq output
def process_MAJIQ(data,
                  mapping,
                  Majiq_confidence, 
                  min_delta,
                  only_DDIs,
                  nease_data):
        try:     
            
            # extract exon skipping events:
            if data['ES'].dtype=='bool':
                data=data[ data['ES']==True]
            else:
                data=data[ data['ES']=='True']
                
                
            # helper functions:
            print('Processing MAJIQ format...')
            f = lambda x: [abs(float(y)) for y in x.split(';')]
            l=lambda x:  any( y >= min_delta for y in x )

            # get Delta PSI values for each junction
            data["delta"] = data['E(dPSI) per LSV junction'].apply(f)
            data["P(|dPSI|>=0.20)"] = data['P(|dPSI|>=0.20) per LSV junction'].apply(f)
            data["delta_sif"] = data['delta'].apply(l)
            data=data[ data['delta_sif']]
        except :
            raise Exception  ("Could not recognize MAJIQ format. Please make sure your table matches MAJIQ output or use the standard format.")

        #filter for significant of diff. AS events
        # only keep diff used junction with confidence higher than 'Majiq_confidence' (for instance: 0.95)
        
        if data.empty:
              raise Exception  ("No significant events with the chosen cutoff.")


        l=lambda x:  any( y >=Majiq_confidence for y in x )
        data["delta_sif"] = data['P(|dPSI|>=0.20)'].apply(l)
        data=data[ data['delta_sif']]
        
        if data.empty:
              raise Exception  ("No significant events with the chosen cutoff.")
        
        z= lambda x: x.split(':')[1]
        data['Gene ID'] = data['Gene ID'].apply(z)
        
        
            
        #############################
        # get all spliced genes here
        # get all coding genes affected by splicing
        # Only genes with Pfam domain will be considred here
        # NCBI id used here for the network visualization
        spliced_genes=list(data['Gene ID'].unique())
        #spliced_genes=[Ensemb_to_entrez(x,mapping) for x in spliced_genes ]

        # Define a custom function to handle NaN values during join
        def custom_join(row):
            if pd.notna(row["Junctions coords"]) and pd.notna(row["IR coords"]):
                return f"{row['Junctions coords']};{row['IR coords']}"
            else:
                return row['Junctions coords']


        junctions=list(data['Junctions coords'])
        ir=list(data['IR coords'])
        data['test'] = data.apply(custom_join, axis=1)
        junctions = list(data['test'])
        confidence=list(data['P(|dPSI|>=0.20)'])
        junc_confid=dict(zip(junctions, confidence))

        complexx=0
        jun_to_target={}
        jun_to_source={}
        for j in junc_confid.keys():
            #find the source of the junction
            x=re.findall(r"[\w']+", j)
            source=list(set([l for l in x if x.count(l) == len(x)/2]))

            #complex events : skip for now
            # TO DO later
            if len(source)==0:      complexx+=1

            #source found
            #search for all possible targets
            else: 
                targets=[ y for y in x if y not in source ]

                test = len(targets)
                test2 = len(junc_confid[j])
                #check if we have correct confidence for every target
                if len(targets)==len(junc_confid[j]):

                    #filter low confident diff used event
                    targets=[ x  for x,y in zip(targets, junc_confid[j]) if y>= Majiq_confidence]

                    #save
                    jun_to_target[j]=targets
                    jun_to_source[j]=int(source[0])
        t=lambda x:  jun_to_target[x] if x in jun_to_target.keys() else False
        data["targets"] = data['test'].apply(t)
        
        
        g=lambda x:  jun_to_source[x] if x in jun_to_source.keys() else False
        data["source"] = data['test'].apply(g)
        
        
        data=data[data['targets']!=False]
        
        data=splitDataFrameList(data,'targets')
        
        
        mapping_tb=[]
        
        # check mapped exons
        if len(data)==0:
                    print('None of MAJIQ junctions maps to annotated exons (Ensembl Exons). Try to use the standard input instead.')

        else:
                data=data  [ [ 'Gene ID','E(dPSI) per LSV junction' , 'Junctions coords','P(|dPSI|>=0.20) per LSV junction','delta','targets','source' ]]

                data['targets']=data['targets'].astype(int)


                #Map exons to domain
                mapping_tb=pd.merge(mapping, data,  left_on='Genomic coding start', right_on='targets')


                #check if the mapping based on coordinate match the gene ID provided from Majiq
                mapping_tb=mapping_tb[mapping_tb['Gene ID']==mapping_tb['Gene stable ID']]
                
                
                
                
                
                
                elm_affected,pdb_affected=[],[]                
                if not only_DDIs:
                        # map to pdb and elm
                        elm_affected,pdb_affected=get_interfaces(data,nease_data,overlap_cal=False)
                        
                        
                        
                        
                    
                

                m= lambda x: max(list(x))
                mapping_tb['max_change']=mapping_tb['delta'].apply(m)

                mapping_tb=mapping_tb[['Exon stable ID','Gene name',
                               'NCBI gene ID','Gene stable ID','Genomic coding start','Genomic coding end','max_change','Pfam ID','source','targets']]


                mapping_tb=mapping_tb.sort_values(['max_change'], ascending=False, kind='stable')
                mapping_tb=mapping_tb[mapping_tb['NCBI gene ID'].notnull()]
                mapping_tb['NCBI gene ID']=mapping_tb['NCBI gene ID'].astype('int').astype('str')
                mapping_tb=mapping_tb.drop_duplicates()

                # an extra filtering step
                # make sure that the source belong to an annotated exon 
                #mapping_tb=pd.merge(mapping[['Genomic coding end']], mapping_tb,  left_on='Genomic coding end', right_on='source')



    
                mapping_tb=mapping_tb[['Gene name','NCBI gene ID','Gene stable ID','Exon stable ID','Pfam ID','max_change']].drop_duplicates()


                print('MAJIQ output converted successfully to NEASE format.')


        return mapping_tb,spliced_genes,elm_affected,pdb_affected
    
    
    
    
    
    
    
def get_interfaces(data,nease_data,min_delta=.05,overlap_cal=True):
    
    # data: list of exons to map
    if overlap_cal:
        
        
            columns=data.columns
            
            #merge with elm table
            elm_mapped=pd.merge(data,nease_data.elm,   left_on=columns[0], right_on='Gene stable ID').drop_duplicates()
            
            
            #merge with pdb table
            pdb_mapped=pd.merge(data,nease_data.pdb,   left_on=columns[0], right_on='Gene stable ID').drop_duplicates()
            
        
            # any overlap is considered
            elm_mapped['overl']=elm_mapped[["Genomic coding start", columns[1]]].max(axis=1) <= elm_mapped[["Genomic coding end",columns[2] ]].min(axis=1)
            pdb_mapped['overl']=pdb_mapped[["Genomic coding start", columns[1]]].max(axis=1) <= pdb_mapped[["Genomic coding end",columns[2] ]].min(axis=1)



            elm_mapped=elm_mapped[elm_mapped['overl']].drop_duplicates()
            pdb_mapped=pdb_mapped[pdb_mapped['overl']].drop_duplicates()

            try:
                        #try to get the delta PSI from the user input
                elm_mapped['dPSI']=elm_mapped[columns[3]].astype(float)
                pdb_mapped['dPSI']=pdb_mapped[columns[3]].astype(float)


                        # significance filter (delta PSI)
                elm_mapped=elm_mapped[elm_mapped['dPSI'].abs()>=min_delta]
                pdb_mapped=pdb_mapped[pdb_mapped['dPSI'].abs()>=min_delta]

            except:
                    # No delta PSI given

                    elm_mapped['dPSI']='-'
                    pdb_mapped['dPSI']='-'




            elm_mapped=elm_mapped[['Gene stable ID','entrezgene','ELMIdentifier','dPSI']].drop_duplicates()
                
                

                  
    
    else:
            #no overal: direct mapping of junctions (Majiq output)
            # already filtered by signifcant
            elm_mapped=pd.merge(nease_data.elm, data,  left_on='Genomic coding start', right_on='targets')
            pdb_mapped=pd.merge(nease_data.pdb, data,  left_on='Genomic coding start', right_on='targets')
                
                

            #check if the mapping based on coordinate match the gene ID provided from Majiq
            elm_mapped=elm_mapped[elm_mapped['Gene ID']==elm_mapped['Gene stable ID']]
            pdb_mapped=pdb_mapped[pdb_mapped['Gene ID']==pdb_mapped['Gene stable ID']]


            m= lambda x: max(list(x))
            elm_mapped['dPSI']=elm_mapped['delta'].apply(m)
            pdb_mapped['dPSI']=pdb_mapped['delta'].apply(m)

                
                
    #convert 
    elm_mapped['Gene name']=elm_mapped['entrezgene'].apply(lambda x: Entrez_to_name(str(x),nease_data.gene_mapping,nease_data.organism))
    elm_mapped['ID']=elm_mapped['entrezgene'].astype("str")+"/"+elm_mapped['ELMIdentifier'].astype("str")
    elm_mapped=elm_mapped[['Gene name','entrezgene','Gene stable ID','ELMIdentifier','dPSI','ID']].drop_duplicates()
        
    
    
    pdb_mapped=(pdb_mapped.dropna().groupby(["symbol",'Gene stable ID'])["entrezgene"].agg(list)).sort_values(ascending=False).reset_index(name='entrezgene')
    pdb_mapped['entrezgene']=pdb_mapped['entrezgene'].apply(lambda x: list(set(x)))
    pdb_mapped['NCBI gene ID']=pdb_mapped['symbol'].apply(lambda x: name_to_entrez(str(x),nease_data.gene_mapping,nease_data.organism))
    pdb_mapped['NCBI gene ID']=pdb_mapped['NCBI gene ID'].astype('str')
    

    
    return elm_mapped,pdb_mapped
    
    

    
# Functions to convert IDs
# organism: 'Human' or 'Mouse', used as the species ('human'/'mouse') of the online (mygene) lookup

# Answers of the online lookups, so that a gene missing from the local table is only queried once per session
online_lookups={}

mg = mygene.MyGeneInfo()
# the client sleeps (1 s) after every batch to space out large batch queries; these queries are a single batch (one gene)
mg.delay = 0

def lookup_online(gene,scope,field,organism):
    key=(str(gene),scope,field,organism.lower())
    if key not in online_lookups:
        out = mg.querymany(str(gene), scopes=scope, fields=field, species=organism.lower(), verbose=False, as_dataframe = True,
              df_index = True)
        # not found: no column for the field
        online_lookups[key]=out[field].to_numpy() if field in out else None
    if online_lookups[key] is None:
        raise KeyError(gene)
    return online_lookups[key]


def Entrez_to_name(gene,mapping,organism='Human'):
    try:
        name=mapping[mapping['NCBI gene ID']==int(gene)]['Gene name'].unique()

        if len(name)==0:

            # try to convert it online
            name=lookup_online(gene,"entrezgene",'symbol',organism)
        
        return name[0]
    
    except :
        return gene
    
    
def Ensemb_to_entrez(gene,mapping,organism='Human'):
    try:
        entrez=mapping[mapping['Gene stable ID']==gene]['NCBI gene ID'].unique()
        if len(entrez)==0:
            # try to convert it online
            entrez=lookup_online(gene,"ensembl.gene",'entrezgene',organism)
            

        return int((entrez[0]))
    
    except :
        return gene
    
def Ensemb_to_name(gene,mapping,organism='Human'):
    try:
        name=mapping[mapping['Gene stable ID']==gene]['Gene name'].unique()

        if len(name)==0:
            # try to convert it online
            name=lookup_online(gene,"ensembl.gene",'symbol',organism)
            
        return name[0]
    
    except :
        return gene
    
    
# Functions to convert IDs 
def name_to_entrez(gene,mapping,organism='Human'):
    try:
        name=mapping[mapping['Gene name']==str(gene)]['NCBI gene ID'].unique()


        if len(name)==0:
            # try to convert it online
            name=lookup_online(gene,"symbol",'entrezgene',organism)
        
        return name[0]
    
    except :
        return str(gene)