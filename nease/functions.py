from .process import *
import scipy.stats as stats
import plotly.graph_objects as go

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np




# main functions for nease
      
def exons_to_edges(mapped,G,elm_interactions,organism):
    
        # check if domains have known interactions/binding:
        mapped['domain']=mapped['NCBI gene ID']+'/'+mapped['Pfam ID']
        
        # domain in ddi
        mapped['DDI']=mapped['domain'].apply(lambda x: G.has_node(x))
        
        
        # domaimn in elm
        mapped['elm']=mapped['domain'].apply(lambda x: x in  list(elm_interactions['interactor 2'].unique()))

        
        mapped['Interacting domain']=mapped.apply(lambda x: bool ( x['elm'] + x['DDI']) , axis=1)
        
        def interactiontypye(ddi,elm):
            if (ddi & elm):
                return "DDI and DMI"
            elif ddi:
                return "DDI"
            elif elm:
                 return "DMI"
            
                                                  
        mapped['Interaction type']=mapped.apply(lambda x: interactiontypye(x['DDI'], x['elm']), axis=1)             
                                                 
        #mapped['Interacting domain']=mapped['domain'].apply(lambda x: G.has_node(x))
        
        mapped=mapped.rename(columns={"max_change": "dPSI",
                                      "domain": "Domain ID"}).reset_index(drop=True)
        mapped['Visualization link'] = ''
        if organism == 'Human':
            mapped.loc[mapped['Interacting domain'], ['Visualization link']] = DIGGER + 'human/' + mapped['Exon stable ID']
        elif organism == 'Mouse':
            mapped.loc[mapped['Interacting domain'], ['Visualization link']] = DIGGER + 'mouse/' + mapped['Exon stable ID']
        return mapped
    
    

def affected_edges(nease_data,Join,only_DDIs):
    
            data=nease_data.data
            mapping=nease_data.mapping
            elm_interactions=nease_data.elm_interactions

            
            
            
            # get domains with DDIs
            interacting_domains=data[data['Interacting domain']]



            # Identify binding of affected domains = Edges in the PPI
            
            
            # helper function to search in ddi
            #retun nothing in case there is no node in joint graph
            def get_neighbors(graph,node):
                try:
                    return graph.neighbors(node)
                except :
                    return []
            
            
            def get_elm(node):
                
                interactors=list(elm_interactions[elm_interactions['interactor 2']==node]['Interator gene 1'].unique())
                return [str(x) for x in interactors]
            
                    
            # search in DDI and ELM        
            t=lambda node:  list(set(    [x.split('/')[0] for x in [n for n in get_neighbors(Join,node)]]     
                                                      +  get_elm(node)    
                                             ))
            
            interacting_domains['Affected binding (NCBI)']=interacting_domains['Domain ID'].apply(t)
            interacting_domains=interacting_domains.rename(columns={"Pfam ID": "Identifier"}).reset_index(drop=True) 
            
            
            # get edges from elm nad pdb
            if not only_DDIs:
                
                elm_affected=nease_data.elm_affected
                pdb_affected=nease_data.pdb

                
                # get elm edges
                t=lambda x: [ str(x) for x in list(elm_interactions[elm_interactions['interactor 1']==x ]['Interator gene 2'].unique())]
                elm_affected['Affected binding (NCBI)']=elm_affected['ID'].apply(t)
                
                # only elm with interactions
                elm_affected=elm_affected[elm_affected['Affected binding (NCBI)'].map(lambda d: len(d)) > 0]
                
                
                # get elm edges

                
                
                
                elm_affected=elm_affected.rename(columns={"ELMIdentifier": "Identifier",
                                          "entrezgene": "NCBI gene ID"}).drop(columns=[ 'Gene stable ID','ID']).reset_index(drop=True)





                interacting_domains=pd.concat([interacting_domains, elm_affected], ignore_index=True)



            return interacting_domains
    
    
    
    
    
    
def gene_to_edges(data,pdb,only_DDIs):
    
    #For every gene get all edges
        gene_edges={}
        for gene in data['NCBI gene ID'].unique():
            edges=data[data['NCBI gene ID']==gene]['Affected binding (NCBI)']
            edges=[item for sublist in edges for item in sublist]
            
            gene_edges[gene]=list(set(edges))
            
        # get pdb interactions    
        if not only_DDIs:
            
            for gene in pdb['NCBI gene ID'].unique():
                    edges=pdb[pdb['NCBI gene ID']==gene]['entrezgene']
                    edges=[item for sublist in edges for item in sublist]
                    
                    if gene in gene_edges:
                        
                        gene_edges[gene]=list(set(gene_edges[gene]+edges))
                    else:
                        gene_edges[gene]=list(set(edges))

            
        return gene_edges
        
    
    
    
    
def builtin_network_size(organism,only_DDIs):
    # Totat degree of structural network for human (pre-computer)
    # For statistical test: edge enrichment
    # Used with the built-in pathway tables; returns (number of edges, pathway degree column)

    # rebuilt network (NEASE_DATA_DIR or replaced data files): sizes from its network_size_<organism>.json
    if organism in network_sizes:
        ppi_type='Degree in the PPI/DDI' if only_DDIs else 'Degree in the structural PPI'
        return network_sizes[organism][ppi_type],ppi_type

    if organism=='Human':
        if only_DDIs:
            n=52467

            #  every pathway degreee two sturctural PPIs and
            # 'Degree in the structural PPI' : degree in ppi annotated with DDI.DMI/PDB
            # 'Degree in the PPI/DDI' : degree in ppi annotated with DDI only
            ppi_type='Degree in the PPI/DDI'

        else:
            n=60235
            ppi_type='Degree in the structural PPI'

    elif organism=='Mouse':
        # The mouse pathway degrees and network sizes are computed from the bundled networks
        # with gene_sets.network_degrees('Mouse') (PPI/DDI: 13348 edges, structural PPI: 31144 edges)
        if only_DDIs:
            n=13348

            #  every pathway degreee two sturctural PPIs and
            # 'Degree in the structural PPI' : degree in ppi annotated with DDI.DMI/PDB
            # 'Degree in the PPI/DDI' : degree in ppi annotated with DDI only
            ppi_type='Degree in the PPI/DDI'

        else:
            n=31144
            ppi_type='Degree in the structural PPI'

    return n,ppi_type



def pathway_enrichment(g2edges,paths, mapping,organism,p_value_cutoff,only_DDIs,n=None,ppi_type=None):
    # General enrichment analysis
    # n, ppi_type: number of network edges and pathway degree column (default: built-in pathway tables)

    if n is None:
        n,ppi_type=builtin_network_size(organism,only_DDIs)

    # gene names, looked up once (not for every pathway)
    names={}

    # number of effected edges
    affected_edges=len([item for sublist in g2edges.values() for item in sublist])



    
    # for every path :
    path_name=[]
    path_id=[]
    source=[]
    genes=[]
    p_values=[]
    score=[]
    
    # iterate over rows (not names): some pathway names occur in several sources
    for path, path_source, path_external_id, path_genes, p in zip(paths['pathway'],
                                                                   paths['source'],
                                                                   paths['external_id'],
                                                                   paths['entrez_gene_ids'],
                                                                   paths[ppi_type]):
        # count of affected edges connected to the pathway
        # specific to that pathway list p
        # initialise the variable for every path
        connected=0
        genes_tmp=[]
        gene_count=0

        # path total degree "p"
        p=int(p)
        path_genes=set(path_genes)

        for gene in g2edges:
                # for every affected gene
                # count affected gene edges connected to the 
                # specific to the gene and to the pathway list
                tmp=len([x for x in g2edges[gene] if x in path_genes ])

                if tmp>0:

                    # increment for path edges
                    connected=connected+tmp
                    
                    # add gene to the gene list of the pathway
                    if gene not in names:
                        names[gene]=Entrez_to_name(gene,mapping,organism)
                    genes_tmp.append(names[gene] +" ("+str(tmp)+")")
                    

                    # gene specific test
                    _,p_gene=edge_enrich(tmp ,len(g2edges[gene])-tmp , p, n)

                    if p_gene<=0.05:
                        # gene with edges siginifically connected to the pathway
                    
                        gene_count=gene_count+1
                    
        
        
        #  affected edges not connected to tha pathway
        not_connected=affected_edges-connected
        
        # Join function is slow can be optimized later 
        if genes_tmp==[]: genes_tmp= ''
        else: genes_tmp=(", ").join(genes_tmp)
            
        # Run hypergeometric test on affected edges
        _,p_value_tmp=edge_enrich(connected ,not_connected , p, n)
        
    

        p_values.append(p_value_tmp)
        
        #compute combined score
        if p_value_tmp<p_value_cutoff:
            s= -(np.sqrt(gene_count) * np.log10(p_value_tmp))
        else: s=0
            
        score.append( s  )
        path_name.append(path)
        source.append(path_source)
        path_id.append(path_external_id)
        genes.append(genes_tmp)

        
    # save results
    Enrichment = pd.DataFrame(list(zip(path_id, path_name,source, genes,p_values,score)), 
                              columns= ['Pathway ID', 'Pathway name', 'Source', 'Spliced genes (number of interactions affecting the pathway)', "p_value",'Nease score'] )
    

    
    
    return Enrichment.sort_values(['p_value'], ascending=True, kind='stable')



def single_path_enrich(path_id,Pathways,g2edges,mapping,organism,only_DDIs,n=None,ppi_type=None):

        # n, ppi_type: number of network edges and pathway degree column (default: built-in pathway tables)
        if n is None:
            n,ppi_type=builtin_network_size(organism,only_DDIs)

        p=int(Pathways [Pathways['external_id']==path_id][ppi_type].iloc[0])
        path_genes=list(Pathways [Pathways['external_id']==path_id]['entrez_gene_ids'])[0]
        
        #collect:
        spliced_genes=[]
        spliced_genes_entrez=[]
        gene_association=[]
        num=[]
        affected_edges=[]
        affected_edges_entrez=[]
        p_val=[]
        
        # graph to save affected edges
        G = nx.Graph()
        
        for g in g2edges:
            
            # affected edges of the gene g
            affected=g2edges[g]
            
            # edges connected to the pathway
            edges=[x for x in affected if x in path_genes ]
            a=len(edges)
            
            # Not connected
            b=len(affected)-a

            if a!=0:
                # calculate gene specific p_value:
                _,p_value=edge_enrich(a,b,p,n)
                
                # Save results
                spliced_genes.append(Entrez_to_name(g,mapping,organism))
                spliced_genes_entrez.append(g)
                gene_association.append(g in path_genes)
                num.append(str(a)+'/'+str(a+b))
                affected_edges.append((',').join([ Entrez_to_name(x,mapping,organism) for x in edges]))
                affected_edges_entrez.append((',').join(edges))
                p_val.append(p_value)
                
                # save affected edges
                G.add_edges_from( [ (g,x) for x in edges ]  )
                
        Enrichment = pd.DataFrame(list(zip(spliced_genes,spliced_genes_entrez, gene_association, num,p_val, affected_edges,affected_edges_entrez)), 
                          columns= ['Spliced genes','NCBI gene ID','Gene is known to be in the pathway','Percentage of edges associated to the pathway', 'p_value', 'Affected binding (edges)','Affected binding (NCBI)'] )
        
        
        
        
        
        return Enrichment,G

    
# stats functions
# Statistical test
# Edge enrichment test

def edge_enrich(a,b,p,n):
    #function to calculate P value
    #test if affected edges by AS are significally enriched in a pathway
    # fisher exact test
    # a+b affected adges, with a the one linked to pathway p
    #p total degree of pathway p
    #n total edges in the ppi used.
    
    #linked to pathway but not affected
    c=p-a
    if c < 0:
        return np.nan, 1.0
    # background of test: not linked to p and not affected edges
    d=(2*n)-p-b
    
    # retun oddsratio and pvalue from fisher exact test
    return   stats.fisher_exact([[a, b], [c, d]], alternative='greater')







# Visualization Function

def extract_subnetwork(path_genes,
                       ppi,
                       affected_genes,
                       all_spliced_genes,
                       k,
                       mapping,
                       affected_graph,
                       significant,
                       organism='Human'):


        # Affected_genes: genes with lost/gained interaction in the pathway
        # all_spliced_genes: all genes affected with splicing
        all_spliced_genes=[ Ensemb_to_entrez(x,mapping,organism) for x in all_spliced_genes]
        

        # Extract the pathway module for the complete PPI
        # We would like to visualize the pathway with affected edges:

        G=ppi.subgraph(path_genes+affected_genes)
        G = nx.Graph(G)
        G.remove_nodes_from(list(nx.isolates(G)))

        
        
        #Position nodes using Fruchterman-Reingold force-directed algorithm.
        # (fixed seed: the same pathway is always drawn the same way)
        pos = nx.spring_layout(G, k=k, iterations=100, seed=0)
        #pos = nx.kamada_kawai_layout(G)
        

        
        # Prepare the visualization
        for n, p in pos.items():
            G.nodes[n]['pos'] = p
            
        
        # Collect nodes and edges in lists and create the plotly traces once at the end
        # (growing a trace with += re-validates all of its points every time, which is very slow for large pathways)
        node_x,node_y,node_text,node_color,node_size=[],[],[],[],[]
        
        for node in G.nodes():
            x, y = G.nodes[node]['pos']
            node_x.append(x)
            node_y.append(y)
            node_info=Entrez_to_name(node,mapping,organism)
            node_text.append(node_info)
            
            if not node in path_genes:
                #Node not in pathway
                color='orange'
                if node in significant:
                    size=50
                else:size=30

            elif int(node) in all_spliced_genes:
                # spliced and part of the pathway
                color='red'
                if node in significant:
                    size=50
                else:size=30
                
            else:
                # other pathway nodes
                color='#888'
                if node in significant:
                    size=50
                else:size=20
                
            node_color.append(color)
            node_size.append(size)
            
            
        edge_x,edge_y,colored_edge_x,colored_edge_y=[],[],[],[]
        
        for edge in G.edges():
            
                x0, y0 = G.nodes[edge[0]]['pos']
                x1, y1 = G.nodes[edge[1]]['pos']
                
                

                # Check if edge is affected
                if affected_graph.has_edge(*edge):
                        colored_edge_x += [x0, x1, None]
                        colored_edge_y += [y0, y1, None]
                    
                    
                else:
                        edge_x += [x0, x1, None]
                        edge_y += [y0, y1, None]


        # Define node and edges in plott
        node_trace = go.Scatter(x=node_x,
                                y=node_y,
                                text=node_text,
                                mode='markers+text',
                                hoverinfo='text',
                                textposition='top center',
                                #textfont_size = 22,
                                marker=dict(
                                reversescale=True,
                                    color=node_color,
                                    size=node_size,
                                    line=dict(width=0)))
        
        
        edge_trace = go.Scatter(x=edge_x,
                                y=edge_y,
                                mode='lines',
                                line=dict(width=0.7,
                                          color='#888'),
                                hoverinfo='none')
        
        colored_edge_trace = go.Scatter(x=colored_edge_x,
                                y=colored_edge_y,
                                mode='lines',
                                line=dict(width=4,
                                          color='red'),
                                hoverinfo='none')
        

        return [ colored_edge_trace,node_trace,edge_trace]
            
        
        

        
        
# plot domains stats

def stats_domains(affecting_percentage,
                  number_of_features,
                  domain_number,
                  elm_number,
                  pdb_number,
                  file_path):
    
    
        # from https://matplotlib.org/stable/gallery/pie_and_polar_charts/bar_of_pie.html#sphx-glr-gallery-pie-and-polar-charts-bar-of-pie-py
        # make figure and assign axis objects
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7))
        fig.subplots_adjust(wspace=0)

        # pie chart parameters
        ratios = [affecting_percentage, 1-affecting_percentage]
        labels = ['Affecting protein features','Not affecting any feature']
        explode = [0.1, 0,]
        # rotate so that first wedge is split by the x-axis
        angle = -180 * ratios[0]
        ax1.pie(ratios, autopct='%1.1f%%', startangle=angle,
                labels=labels, explode=explode, shadow=True,)
        ax1.set_title("Genes with AS affecting protein features")
        # bar chart parameters

        xpos = 0
        bottom = 0
        numbers = [elm_number, pdb_number, domain_number]
        ratios = [ round(x/number_of_features, 2) for x in numbers]
        width = .2
        colors = ['#F0D0C8','#B09880','#9B412B']
        features = ['Linear motifs', 'Residues','Domains']

        # only draw features that are present (e.g. no motifs/residues with only_DDIs=True)
        bars = []
        for number, height, color, feature in zip(numbers, ratios, colors, features):
            if number == 0:
                continue
            bar = ax2.bar(xpos, height, width, bottom=bottom, color=color, label=feature)
            bars.append(bar)
            ypos = bottom + height / 2
            bottom += height
            ax2.text(xpos, ypos, "%d%%" % (height * 100),
                     ha='center')

        ax2.set_title('Affected features')
        ax2.legend(handles=bars)
        ax2.axis('off')
        ax2.set_xlim(- 2.5 * width, 2.5 * width)

        # use ConnectionPatch to draw lines between the two plots
        # get the wedge data
        theta1, theta2 = ax1.patches[0].theta1, ax1.patches[0].theta2
        center, r = ax1.patches[0].center, ax1.patches[0].r
        bar_height = sum([item.get_height() for item in ax2.patches])

        # draw top connecting line
        x = r * np.cos(np.pi / 180 * theta2) + center[0]
        y = r * np.sin(np.pi / 180 * theta2) + center[1]
        con = ConnectionPatch(xyA=(-width / 2, bar_height), coordsA=ax2.transData,
                              xyB=(x, y), coordsB=ax1.transData)
        con.set_color([0, 0, 0])
        con.set_linewidth(4)
        ax2.add_artist(con)

        # draw bottom connecting line
        x = r * np.cos(np.pi / 180 * theta1) + center[0]
        y = r * np.sin(np.pi / 180 * theta1) + center[1]
        con = ConnectionPatch(xyA=(-width / 2, 0), coordsA=ax2.transData,
                              xyB=(x, y), coordsB=ax1.transData)
        con.set_color([0, 0, 0])
        ax2.add_artist(con)
        con.set_linewidth(4)
        
        file_path=os.path.join(os.path.dirname(file_path),'NEASE_stats.pdf')

        plt.savefig(file_path,format='pdf',bbox_inches='tight')
        plt.show()
        
        # Save the figure


            
        return