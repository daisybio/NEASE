
from .functions import *
from .gene_sets import load_gene_sets, GENE_SET_DATABASES, DEGREE_DDI, DEGREE_STRUCTURAL
import statsmodels.api as sm
import gseapy as gp


# Main Code 

class run(object):

    
    
    
    def __init__(self, 
                 data ,
                 organism='Human',
                 input_type='Standard',
                 p_value_cutoff=0.05,
                 remove_non_in_frame=True,
                 only_divisible_by_3=False,
                 min_delta=0.05,
                 Majiq_confidence=0.95,
                 only_DDIs=False
                
                ):
        
        
        
        """
            Create an instance of NEASE object.
            
            
            Parameters
            ---------- 
            data: dataframe 
                Splicing (or diff. splicing) exons  or events. 
               
                Standard input: Ensemble gene ID    Start of exon  End of exon     dPsi (optional)
                for MAJOQ output, Please change the  input_type to "MAJIQ"
                
            organism : str, optional
                NEASE 1.0 support only Human.
            
            input_type: str, optional
                
                Either "Standard","Spycone",'Whippet','rmats','DEXSeq'or "MAJIQ", If you need support of more types of outputs. Please contact: louadi@wzw.tum.de
                
            p_value_cutoff: float, optional
               The p value cutoff used to compute NEASE scores. (default is 0.05)
               
            remove_non_in_frame: boolean, default: True
                   If True remove all exons that are predicted to disturb the ORF or known to result in a non-coding gene. If you would like to change this option and include all exons, please change the parameter remove_non_in_frame to False.
           
           only_divisible_by_3: boolean, default: False
                    If True remove exons not divisible by 3
             
            min_delta: float, optional
                min delta to consider in case your input has dPsi column. The value also corresponds  to logfold change in case of Dexeq.   (default is 0.05)
            
            Majiq_confidence: float, optional
                In case of input_type='MAJIQ'. the parameter P(dPSI > 20%) is needed. Check MAJIQ paper for details about this  (default is 0.95)  
                
            only_DDIs: boolean, optional, default: False
               if True only use DDI and exclude PDB and linear motif interactions
            
                
         """
        self.data=[]
        self.organism=organism
        self.p_value_cutoff=p_value_cutoff
        # gene sets loaded with enrich_msigdb()
        self.gene_sets=None

        if organism!='Human' and organism!='Mouse':
            print('Error: Please choose one of the  supported  organism: "Human" and "Mouse".')
        
         
        if input_type!='MAJIQ' and input_type!='Standard' and input_type!='Spycone' and input_type!='Whippet' and input_type!='rmats' and input_type!='DEXSeq':
            raise ValueError('Input type not supported')
            
        else:
        
            
            #Open the Join graph and databases of the selected organism:
            Join=network[organism]
            self.mapping=database_mapping[organism]
            self.gene_mapping=gene_mapping[organism]
            self.path=Pathways[organism]
            self.ppi=PPI[organism]
            self.only_DDIs=only_DDIs
            self.elm_interactions=elm_interactions[organism]
            self.input_type=input_type

            if not only_DDIs:
                self.elm=elm[organism]
                self.elm_interactions=elm_interactions[organism]
                self.pdb=pdb[organism]
            self.non_coding=non_coding
            self.data=[]
            self.spliced_genes=[]
            
            
            
            if input_type=='MAJIQ':

                # Processing Majiq output
                        self.data,self.spliced_genes,self.elm_affected,self.pdb_affected=process_MAJIQ(data,self.mapping, Majiq_confidence, min_delta ,self.only_DDIs,self)
                        if len(self.data)==0:
                            print('Found no overlap with protein domains.')

                
            
            
            elif input_type=='Standard':
                    
                try:
                    self.data,self.spliced_genes,self.elm_affected,self.pdb_affected,self.symetric_genes=process_standard(data,self.mapping,min_delta ,self.only_DDIs,self,remove_non_in_frame,only_divisible_by_3)
                    if len(self.data)==0:
                        print('Found no overlap with protein domains.')
                        print('Make sure that the genomic coordinates of the exons correspond to the human genome build hg38 (GRCh38).')

                    
                except:
                        print('Could not recognize the standard format. Please make sure your table matches the standard format.')
                        print('Gene ensembl ID          EXON START        EXON END          dPSI (optional)')
                        print('Make sure that the genomic coordinates of the exons correspond to the human genome build hg38 (GRCh38).')
                    
            elif input_type=='Whippet':
                
                try:
                    data = data[data['DeltaPsi'] >= 0.90]
                    data=data.rename_axis('Gene ID').reset_index()
                    data['tmp']= data['Node'].apply(lambda x: x.split(':')[1])

                    data['start']=data['tmp'].apply(lambda x: x.split('-')[0])
                    data['end']=data['tmp'].apply(lambda x: x.split('-')[1])
                    data=data[['Gene ID' , 'start','end','DeltaPsi']]

                    self.data,self.spliced_genes,self.elm_affected,self.pdb_affected,self.symetric_genes=process_standard(data,self.mapping,min_delta ,self.only_DDIs,self,remove_non_in_frame,only_divisible_by_3)
                    
                except:
                        print('process failed....Try to use the Standard input')
                
                
            elif input_type=='rmats':  
                
                try:
                    
                            data=data[data['FDR']<=p_value_cutoff]
                            data=data[['GeneID','exonStart_0base','exonEnd','IncLevelDifference']]
                            data['GeneID']=data['GeneID'].apply(lambda x: x.split('.')[0])

                            self.data,self.spliced_genes,self.elm_affected,self.pdb_affected,self.symetric_genes=process_standard(data,self.mapping,min_delta ,self.only_DDIs,self,remove_non_in_frame,only_divisible_by_3)

                except:
                        print('process failed....Try to use the Standard input')
                        
                        
            elif input_type=='DEXSeq':  
                
                try:
                            data=data[data['padj']<=p_value_cutoff]
                            data=data[[data.columns[1],'genomicData.start','genomicData.end','log2fold_control_case']]
                            print('proceding with log2fold threshold: '+str(min_delta))
                            self.data,self.spliced_genes,self.elm_affected,self.pdb_affected,self.symetric_genes=process_standard(data,self.mapping,min_delta ,self.only_DDIs,self,remove_non_in_frame,only_divisible_by_3)

                except:
                        print('process failed....Try to use the Standard input')
            
            elif input_type=='Spycone':
                    
                    self.data,self.spliced_genes=process_spycone(data,self.mapping)
                    
                    
                    
                    #spycone only uses DDI for now
                    self.only_DDIs=True
                    
                    if len(self.data)==0:
                        print('Found no overlap with protein domains.')
                        
  
                    
            

            

            if len(self.data)==0:#
                print('process canceled...')
                
                
                


            else : 
                self.data=self.data.drop_duplicates(['Gene name','NCBI gene ID','Gene stable ID','Pfam ID'],keep= 'first')
                
                
                
                # check interaction of the domains
                 # prepare elm to doamin interactions file
                self.elm_interactions['interactor 1']=self.elm_interactions['Interator gene 1'].astype('str')+"/"+self.elm_interactions['Elm id of gene 1']
                self.elm_interactions['interactor 2']=self.elm_interactions['Interator gene 2'].astype('str')+"/"+self.elm_interactions['Domain of gene 2']
                
                self.data=exons_to_edges(self.data,Join,self.elm_interactions,self.organism)
                
                
                
                print('\n\t\tData Summary')
                print('**************************************************')

                print(str(len(self.data['Domain ID'].unique()))+' protein domains are affected by AS.\n')
                if not self.only_DDIs:
                    print(str(len(self.elm_affected['ELMIdentifier'].unique()))+" linear motifs are affected by AS.\n"
                          +str(len(self.pdb_affected))+ ' interacting residue are affected by AS.\n')
                      
                print(str(len(self.data[self.data['Interacting domain']]['Domain ID'].unique()))+' of the affected domains/motifs have known interactions.') 
                
                
                


               
                
                
              
                
                
                # Identify binding of affected domains = Edges in the PPI
                # get DMI adn DDI for spliced domains

                self.interacting_domains=affected_edges(self,Join,self.only_DDIs)
            
                             


                #get all edges of a gene from DDIs
                # self.interacting_domains,: DMI and DDI
                # self.pdb_affected: co-resolved interactions
                self.g2edges=gene_to_edges(self.interacting_domains,self.pdb_affected,self.only_DDIs)
                
                print(str(len([item for sublist in self.g2edges.values() for item in sublist]))+' protein interactions/binding affected.')
                
                
                #get affected edges from ELM and PDB
                
                
                
                
                

                
                # Runing Enrichment analysis
                
                print('\n**************************************************')
                print('Running enrichment analysis...')
                
                self.supported_database=  list(self.path['source'].unique())
                self.supported_gene_sets=GENE_SET_DATABASES[organism]
                self.enrichment=pathway_enrichment(self.g2edges,self.path, self.gene_mapping,organism,p_value_cutoff,self.only_DDIs).reset_index(drop=True)
                print('NEASE enrichment done.')
                
                
                
    def get_stats(self, file_path=''):
        
        """
            Display the stats of affected domains by splicing.
            
            Parameters
            ---------- 
            file_path : str
                Path for saving the statistics figure.
            
        """
                
        if len(self.data)==0 :
            print('Processing failed')
        
        else:
            
            # initial number of genes with known affected features
            
            affecting_number=0
            self.spliced_genes=list(set(self.spliced_genes))    

            
            # number of genes with affected domains/number of all events (genes)
            genes_with_domain=self.data['Gene stable ID'].unique()
            domain_number=len(genes_with_domain)
            
            
            affecting_number=domain_number
            number_of_features=domain_number

            # no linear motif / residue data when running with only_DDIs=True
            elm_number=0
            pdb_number=0

            if not self.only_DDIs:
              
                    # genes with affected elm
                    genes_with_elm=list(self.elm_affected['Gene stable ID'].unique())
                    
                    # number of genes with elm affected
                    elm_number=len(genes_with_elm)
                    
                    
                    
                    
                    # number of genes with affected pdb
                    genes_with_pdb=list(self.pdb_affected['Gene stable ID'].unique())
                    pdb_number=len(genes_with_pdb)
                    
                    
                    
                    affecting_number=len(list(set(list(genes_with_elm)+list(genes_with_pdb)+list(genes_with_domain))))
                    
                    number_of_features=number_of_features+elm_number+pdb_number
            
      


                    

            affecting_percentage= round(affecting_number/len(self.spliced_genes), 2)



            
            stats_domains(affecting_percentage,number_of_features,domain_number,elm_number,pdb_number,file_path)
        return
        
        
        
        
    def get_domains(self):
        
        """
            Display the list of AS events affectting domains in NEASE format.
            
        Returns
        -------
        pd.DataFrame Object
        
        """
                
        if len(self.data)==0 :
            print('Processing failed')


        # elif self.organism=="Mouse":
        #         #no visualization available for mouse in DIGGER
        #         return self.data.drop(columns=['Domain ID','Visualization link'])
        # else:
            
        #DIGGER visualization available for Human and Mouse
        return self.data.drop(columns=[ 'Domain ID','DDI','elm']).reset_index(drop=True)
    
    
    
    
    def get_elm(self):
        
        """
            Display list of linear motifs affected by AS.
            
        Returns
        -------
        pd.DataFrame Object
        
        """
        
        if self.only_DDIs:
            print('You ran NEASE with the option:  only_DDIs=True. No ELM data available.')
            
        elif self.elm_affected.empty:
            
            print('No known linear motif are affected by AS')
            
        else:
            self.elm_affected['ELM link']='http://elm.eu.org/elms/'+self.elm_affected.ELMIdentifier.astype('str')
            
            return self.elm_affected.drop(columns=[ 'ID','Affected binding (NCBI)']).reset_index(drop=True)
           
        
    def get_pdb(self):
        
        """
            Display list of linear motifs affected by AS.
            
        Returns
        -------
        pd.DataFrame Object
        
        """
        if self.only_DDIs:
            print('You ran NEASE with the option: only_DDIs=True. No pdb data available.')
     
        elif self.pdb_affected.empty:
            
            print('No residue from the PDB database motif are affected by AS')
            
        else:
            pdb_affected=self.pdb_affected.rename(columns={"symbol": "Gene name",'entrezgene':'Co-resolved interactions'}).copy()
            
            #Convert IDs to names
            c=lambda x: [ Entrez_to_name(gene,self.gene_mapping,self.organism) for gene in list(set(x))]
            pdb_affected['Co-resolved interactions symbol']=pdb_affected['Co-resolved interactions'].apply(c)
            
            
            
            a=lambda x: ",".join(x)
            pdb_affected['Co-resolved interactions']=pdb_affected['Co-resolved interactions'].apply(a)
            pdb_affected['Co-resolved interactions symbol']=pdb_affected['Co-resolved interactions symbol'].apply(a)
            
            return pdb_affected[['Gene name','NCBI gene ID','Gene stable ID','Co-resolved interactions symbol','Co-resolved interactions']].reset_index(drop=True)
        
        
        
        
            
    def get_edges(self):
        
        """   
            Display affected interactions from AS. 
            
        Returns
        -------
        pd.DataFrame Object
        
        """
        if len(self.data)==0:
            print('Processing failed')
        elif len(self.interacting_domains)==0:
            print('No affected edges identified.')
        else:
            
            
            edges=self.interacting_domains.copy()
            #Convert IDs to names
            c=lambda x: [ Entrez_to_name(gene,self.gene_mapping,self.organism) for gene in list(set(x))]
            edges['Affected binding']=edges['Affected binding (NCBI)'].apply(c)

            # count number of affected PPI for every domain
            count=lambda x: len(x)
            edges['Number of affected interactions']=edges['Affected binding'].apply(count)
            
            
            a=lambda x: ",".join(x)
            edges['Affected binding']=edges['Affected binding'].apply(a)
            edges['Affected binding (NCBI)']=edges['Affected binding (NCBI)'].apply(a)
            edges=edges.drop_duplicates()
            edges=edges.sort_values('Number of affected interactions', ascending=False)

            return edges[['Gene name','NCBI gene ID','Identifier','dPSI','Number of affected interactions','Affected binding','Affected binding (NCBI)']].reset_index(drop=True)
        
        
    
          
    
    
    
    def classic_enrich(self,
                       gseapy_databases,
                       outdir='Enrichr_gseapy',
                       non_symmetrical_only=False):
        
        '''
        Classic gene level enrichement using the python library gseapy.
        
        Parameters
        ---------- 
        gseapy_databases: str, list, tuple of Enrichr Library name(s). 
                  or custom defined gene_sets (dict, or gmt file). For more details, please check: https://pypi.org/project/gseapy/ 
                  
        non_symmetrical_only: Bool 
                    Run classic gene set enrichment only on non_symmetrical exons.

         Returns
        -------
        pd.DataFrame Object    
        
        Example:
            gseapy_databases=['KEGG_2019_Human', 'Reactome_2016','WikiPathways_2019_Human']
        
        '''

        if not self.spliced_genes:
            print('No genes found on your input.')
            
        else:
            # gseapy >= 1.3 only accepts lower case organism names in enrichr ('human', 'mouse')
            organism=self.organism.lower()

            # get gene sets supported in gseapy
            gseapy_library = gp.get_library_name(organism=organism)

            if isinstance(gseapy_databases,str):
                gseapy_databases=[gseapy_databases]

            # compare with user input with gseapy_library
            gene_set_database=list(set(gseapy_library).intersection(gseapy_databases))
            
            if not gene_set_database:
                print('none of the gene set databases provided is supported in gseapy library.')
                print('please check https://pypi.org/project/gseapy/')
                return
                
            else:
            
                # run gene set enrichment
                
                if non_symmetrical_only:
                    if self.input_type=="MAJIQ":
                        print('Non-symmetrical exons enrichment is not available for Majiq output. Please use standard input. ')
                        return

                    else:
                        # run only on non-symteric
                        gene_list=[ Ensemb_to_name(x,self.gene_mapping,self.organism) for x in self.symetric_genes]
                else:     
                     # run on all genes
                     gene_list=[ Ensemb_to_name(x,self.gene_mapping,self.organism) for x in self.spliced_genes]
                        
                        
                try:
                    enr=gp.enrichr(gene_list=gene_list,organism=organism,gene_sets=gene_set_database,outdir=outdir)
                except ValueError as e:
                    # gseapy fails while drawing its bar plot when no term passes its cutoff (0.05)
                    if 'No enrich terms' not in str(e):
                        raise
                    print('No gene set passes the gseapy cutoff (adjusted p value 0.05); no plot was saved.')
                    enr=gp.enrichr(gene_list=gene_list,organism=organism,gene_sets=gene_set_database,outdir=outdir,no_plot=True)

                return enr.results.sort_values('Adjusted P-value')
            
            
            
            
            
        
    def enrich(self, database=  ['PharmGKB',
                                 'HumanCyc',
                                 'Wikipathways',
                                 'Reactome',
                                 'KEGG',
                                 'SMPDB',
                                 'Signalink',
                                 'NetPath',
                                 'EHMN',
                                 'INOH',
                                 'BioCarta',
                                 'PID'], cutoff=0.05 ):
        
        """ 
        Run enrichement analysis
        
        
        Parameters
        ---------- 
        database: List 
            gene set sources for enrichment.
        cutoff: float
            The p value cutoff.
        
        
            
        Returns
        -------
        pd.DataFrame Object
            
        Example: 
                events=nease.run(table, organism='Human')
                events.enrich(database=['Reactome'])
        
        """
        
        if len(self.data)==0:
            print('Processing failed')
        elif len(self.interacting_domains)==0:
            print('No affected edges identified.') 
        
        else:
            # Check if user input matches the available databases
            database=[ x for x in database if x in self.supported_database]

            if len(database)==0: 
                print('Please select a supported pathway database as an argument. ')
                print('supported databases for '+self.organism+' are :',[ x for x in  self.supported_database],".")
                print('\n')
            else:


                enrich_results=self.enrichment[self.enrichment['Source'].isin(database)]

                print('NEASE enrichment for the pathway databases:\n',[ x for x in database])
                return self._multiple_testing(enrich_results,cutoff)



    def _multiple_testing(self, enrich_results, cutoff):

        # Correct for multiple testing
        # fdr_bh : Benjamini/Hochberg (non-negative)
        enrich_results['adj p_value']=sm.stats.multipletests(list(enrich_results['p_value']),method='fdr_bh',alpha=0.05)[1]
        #enrich_results['adj p_value']=sm.stats.fdrcorrection(list(enrich_results['p_value']),alpha=0.01)[1]

        # shift column 'score to the last position
        scores = enrich_results.pop('Nease score')
        enrich_results.insert(len(enrich_results.columns), 'Nease score', scores)

        num=len(enrich_results[enrich_results['adj p_value']<=cutoff])
        if num==0:
            print('No enrichment found with the cutoff '+str(cutoff)+'.')
        else:
            print("Found "+str(num)+" enriched pathways after multiple testing correction.\n")

        return enrich_results.sort_values(['Nease score','p_value'],ascending=[False,True]).reset_index(drop=True)



    def enrich_msigdb(self, database=None, cutoff=0.05, gene_sets=None, chip=None):

        """
        Run NEASE enrichment with gene sets from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb)
        and KEGG pathways from the KEGG REST API (https://www.kegg.jp/kegg/rest/keggapi.html).

        The gene sets are downloaded the first time a database is used (MSigDB release 2026.2, current KEGG pathways)
        and read from a local cache afterwards (~/.cache/nease, or the folder in the NEASE_CACHE environment variable).
        The degree of every gene set in the NEASE networks is computed when the gene sets are loaded.


        Parameters
        ----------
        database: List
            Gene set databases for enrichment (all databases if not given). The databases of the organism are listed in
            events.supported_gene_sets:
            Human: 'Hallmark', 'Positional', 'CGP', 'BioCarta', 'KEGG_MEDICUS', 'PID', 'Reactome', 'WikiPathways', 'miRDB',
                   'GTRD', '3CA', 'CGN', 'CM', 'GO:BP', 'GO:CC', 'GO:MF', 'HPO', 'C6:Onco', 'C7:ImmuneSigDB', 'C7:VAX',
                   'C8:Cell_type', 'C9:Comp_Perturb', 'KEGG'
            Mouse: 'Hallmark', 'Positional', 'CGP', 'BioCarta', 'Reactome', 'WikiPathways', 'miRDB', 'GTRD', 'GO:BP',
                   'GO:CC', 'GO:MF', 'MP:General', 'MPT:Tumor', 'M7:Immuno', 'M8:Cell_type', 'KEGG'
        cutoff: float
            The p value cutoff.
        gene_sets: str or list of str, optional
            MSigDB JSON bundle(s) to use instead of downloading the MSigDB databases, e.g. the complete msigdb.v2026.2.Hs.json
            (human (Hs) bundles for human data, mouse (Mm) bundles for mouse data).
        chip: str, optional
            MSigDB "NCBI_Gene_ID" chip file matching the MSigDB release of gene_sets, e.g. Human_NCBI_Gene_ID_MSigDB.v2026.2.Hs.chip,
            used to map gene symbols to NCBI gene IDs. By default the mapping of MSigDB 2026.2 bundled with NEASE is used.


        Returns
        -------
        pd.DataFrame Object
            Same columns as enrich(). 'Pathway ID' is the MSigDB systematic name (e.g. M5890) or the KEGG pathway ID
            (e.g. mmu04010), 'Source' the gene set database. The IDs can be used with path_analysis() and Vis_path().

        Example:
                events=nease.run(table, organism='Mouse')
                events.enrich_msigdb(database=['Hallmark','Reactome','KEGG'])

        """

        if len(self.data)==0:
            print('Processing failed')
        elif len(self.interacting_domains)==0:
            print('No affected edges identified.')

        else:
            if database is None:
                database=self.supported_gene_sets
            elif isinstance(database,str):
                database=[database]

            # Check if user input matches the available databases
            not_supported=[ x for x in database if x not in self.supported_gene_sets]
            database=[ x for x in self.supported_gene_sets if x in database]

            if len(database)==0:
                print('Please select a supported gene set database as an argument. ')
                print('supported gene set databases for '+self.organism+' are :',[ x for x in self.supported_gene_sets],".")
                print('\n')
                return
            if not_supported:
                print('Not supported for '+self.organism+' (ignored):',not_supported)

            table,n,summary=load_gene_sets(database,self.organism,gene_sets,chip)

            if len(table)==0:
                print('No gene sets found for the selected databases.')
                return

            print(str(summary['gene sets'])+' gene sets loaded.')
            if 'gene symbols' in summary:
                print(str(summary['gene symbols in the network'])+' of '+str(summary['gene symbols'])
                      +' MSigDB gene symbols are part of the NEASE network.')
            if 'KEGG download date' in summary:
                print('KEGG pathways downloaded on '+summary['KEGG download date']+'.')

            ppi_type=DEGREE_DDI if self.only_DDIs else DEGREE_STRUCTURAL
            enrichment=pathway_enrichment(self.g2edges,table,self.gene_mapping,self.organism,self.p_value_cutoff,self.only_DDIs,
                                          n=n[ppi_type],ppi_type=ppi_type).reset_index(drop=True)

            # keep the gene sets for path_analysis() and Vis_path()
            self.gene_sets={'table':table,'enrichment':enrichment,'n':n[ppi_type],'degree column':ppi_type}

            print('NEASE enrichment for the gene set databases:\n',summary['databases'])
            return self._multiple_testing(enrichment.copy(),cutoff)



    def _pathway_source(self, path_id):

        # Pathway table, enrichment, number of network edges and degree column for a pathway ID:
        # built-in pathways first, then the gene sets loaded with enrich_msigdb()
        if self.gene_sets is not None and not (self.enrichment['Pathway ID']==path_id).any():
            return self.gene_sets['table'],self.gene_sets['enrichment'],self.gene_sets['n'],self.gene_sets['degree column']
        return self.path,self.enrichment,None,None




    def path_analysis(self,path_id):
        
        '''
        Run enrichment analysis on a specific pathway with details of impact of AS.
        
        Parameters
        ---------- 
        path_id: str 
            The ID of the pathway f interest.
            Run enrich() or enrich_msigdb() to find enriched pathways and their IDs.
            
            
        Returns
        -------
        pd.DataFrame Object
        
        exanple:
                # Run general enrichment
                events=nease.run(table, organism='Human')
                events.enrich(database=['Reactome']
                
                # Pathway specific analysis
                events.path_analysis('R-HSA-388396')

        
        '''
            
        if self.data.empty:
            print('Processing failed')
        elif  self.interacting_domains.empty:
            print('No affected edges identified.') 
            
            
        paths,enrichment,n,ppi_type=self._pathway_source(path_id)
        path_info=enrichment[enrichment['Pathway ID']==path_id]
        
        if len(path_info)==0:
            print('No pathway with the given ID found.')

        elif len(path_info)>1:
            # e.g. INOH pathways have no ID ('None')
            print('The given ID is shared by several pathways. Pathway specific analysis is not possible for this ID.')

        else:
            path_name=list(path_info['Pathway name'])[0]
            print('Enrichment of the pathway: '+path_name+'.\n')
            print('Overall p_value: ',list(path_info['p_value'])[0])
            print('\n')
            # run enrichment
            enrich,_=single_path_enrich(path_id,paths,self.g2edges,self.gene_mapping,self.organism,self.only_DDIs,n,ppi_type)
            
            if len(enrich)==0:
                print('No enrichment or genes found for the selected pathway.')
                return
            else:
                return enrich.sort_values(['p_value']).reset_index(drop=True)
            
            


    def Vis_path(self,
                 path_id,
                 file='', 
                 k=0.8,
                 save_pdf=False,
                 auto_open=True):

            '''
               Visualize the network module of a specific pathway.
               
               
            Parameters
            ---------- 
            path_id: str 
                The ID of the pathway f interest.
                Run enrich() or enrich_msigdb() to find enriched pathways and their IDs.
                
            file: str
                 A string representing a local file path for the html file.

            k: float 
                 -  Float (default=None))is a parameter to be tuned by the user:
                            Position nodes using Fruchterman-Reingold force-directed algorithm.
                            Optimal distance between nodes. If None the distance is set to 1/sqrt(n) where n is the number of nodes. 
                            Increase this value to move nodes farther apart.
                    Link: networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html.
                
            auto_open: Boolean 
            
            
            Returns
            -------
            pd.DataFrame Object

            exanple:
                    # Run general enrichment
                    events=nease.run(table, organism='Human')
                    events.enrich(database=['Reactome']

                    # Pathway specific analysis
                    events.path_analysis('R-HSA-388396')

               
            path_id:    -  str: An unique pathway id. Run  enrich() to get all pathways and their ids.
            K:          -  Float (default=None))is a parameter to be tuned by the user:
                            Position nodes using Fruchterman-Reingold force-directed algorithm.
                            Optimal distance between nodes. If None the distance is set to 1/sqrt(n) where n is the number of nodes. 
                            Increase this value to move nodes farther apart.
                        Link: networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html.
            file         - A string representing a local file path.
            '''
            if self.data.empty:
                print('Processing failed')
            elif self.interacting_domains.empty:
                print('No affected edges identified.') 


            paths,enrichment,n,ppi_type=self._pathway_source(path_id)
            path_info=enrichment[enrichment['Pathway ID']==path_id]

            if len(path_info)==0:
                print('No pathway with the given ID found.')

            elif len(path_info)>1:
                # e.g. INOH pathways have no ID ('None')
                print('The given ID is shared by several pathways. Pathway specific analysis is not possible for this ID.')

            else:
                path_name=list(path_info['Pathway name'])[0]
                print('Enrichment of the pathway: '+path_name+'.\n')
                print('Overall p_value: ',list(path_info['p_value'])[0])
                print('\n')
                # run enrichment

                enrich,affected_graph=single_path_enrich(path_id,paths,self.g2edges,self.gene_mapping,self.organism,self.only_DDIs,n,ppi_type)

                if len(enrich)==0:
                    print('No enrichment or genes found for the selected pathway.')

                    return

                else:

                    # Get genes of the pathway (Entrez IDs)
                    path_genes=list(paths[ paths['external_id']==path_id ]['entrez_gene_ids'])[0]

                    significant= list(enrich [enrich ['p_value']<=0.05] ['NCBI gene ID'].unique())

                    graph_data= extract_subnetwork(path_genes,
                                                   self.ppi,
                                                   list(enrich['NCBI gene ID'].unique()),
                                                   self.spliced_genes,
                                                   k,
                                                   self.gene_mapping,
                                                   affected_graph,
                                                   significant,
                                                   self.organism)




                    path_info=enrichment[enrichment['Pathway ID']==path_id]
                    path_name=list(path_info['Pathway name'])[0]

                    fig = go.Figure(data=graph_data,
                                     layout=go.Layout(
                                        title=dict(text='<br>'+path_name, font=dict(size=16)),
                                        showlegend=False,
                                        hovermode='closest',
                                        margin=dict(b=20,l=5,r=5,t=40),
                                        annotations=[ dict(
                                            text="<br> The large nodes have p_value<=0.05 (affecting the pathway).<br> 🔴 Spliced gene and known to be part of the patwhay.<br> 🟠 Spliced gene but not known to be in the pathway.",
                                            showarrow=False,
                                            font=dict(size=20),
                                            xref="paper", yref="paper",
                                            x=0.005, y=-0.002 ) ],
                                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
                    

                    # pathway names can contain '/' (e.g. 'PI3K/AKT Signaling in Cancer')
                    file_name=path_name.replace('/','_')
                    file_path=os.path.join(os.path.dirname(file),file_name+'.html')
                    fig.write_html(file_path, auto_open=auto_open)
                    
                    if save_pdf:
                        fig.write_image(os.path.join(os.path.dirname(file),file_name+'.pdf'),format='pdf',width=2544,height=1313)
                    

                    print('Visualization of the pathway generated in: '+file_path)
                    
                    return
