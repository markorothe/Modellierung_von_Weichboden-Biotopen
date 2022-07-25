# -*- coding: utf-8 -*-
"""
Created on Sun May 22 21:04:30 2022

@author: marko
"""

# PERFORM CLUSTERING ----------------------------------------------------------
def run_assignment_for_cluster_runs(character_stat_scores):
    import pandas as pd
    from src import data_filtering_functions
    
    def move_column_inplace(df, col, pos):
        col = df.pop(col)
        df.insert(pos, col.name, col)
        return df
    
    cluster_assignments = pd.read_pickle('LABELS.pkl')
    HAUL_DB = pd.read_pickle('HAUL_DB.pkl')                        

    for iassignment in range(0,len(cluster_assignments.index)):
        name_of_run =  list(cluster_assignments.index)[iassignment]
        
        labels = cluster_assignments['Cluster-IDs'].iloc[iassignment]
        haul_ids = cluster_assignments['Haul_IDs'].iloc[iassignment]
        labels = pd.DataFrame(data={'Haul_IDs':haul_ids, 'Labels':labels})
        print(name_of_run, len(labels))
        
        region = name_of_run.split('-')[0]
        haul_db = data_filtering_functions.clip_geomety(HAUL_DB, clipping_object = region)
        PIV = haul_db.pivot_table(index='Haul_ID',columns='ScientificName_Accepted', values='SPECIMENHAUL',aggfunc='sum', fill_value=0)

        ch_sub, ch_full = assign_labels(PIV, labels, character_stat_scores)
        # returns lists of length n_clusters containing analysis-dataframe per cluster
        
        cid = 0
        for idf in ch_sub:
            df = idf.copy()
            df['run'] = name_of_run
            df['CL_ID'] = cid
            df = move_column_inplace(df, 'run', 0)
            df = move_column_inplace(df, 'CL_ID', 1)
            df.index.name='\n Species'
            cid += 1
            df.to_csv('character_selection.csv',mode='a')
        
        cid = 0
        for idf in ch_full:
            df = idf.copy()
            df['run'] = name_of_run
            df['CL_ID'] = cid
            df = move_column_inplace(df, 'run', 0)
            df = move_column_inplace(df, 'CL_ID', 1)
            df.index.name='\n Species'
            cid += 1
            df.to_csv('character_full.csv',mode='a')

# -----------------------------------------------------------------------------

def species_total_presence(piv):
    total_sum = piv.sum(axis = 0)
    piv2 = piv.copy()
    piv2[piv2>0] = 1
    total_count = piv2.sum(axis = 0)
    return total_count, total_sum

def site_and_sum_cl(Ch_Cl,icluster):
    group = Ch_Cl.loc[Ch_Cl['Labels']==icluster]
    group = group.transpose()
    sum_cl = group.sum(axis=1)
    group2 = group.copy()
    group2[group2>0] = 1
    site_cl = group2.sum(axis=1)
    return site_cl, sum_cl

def simper(PIV,labels):
    import warnings
    import pandas as pd
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri#, default_converter
    # activate R/Py connection
    pandas2ri.activate()     
    vegan = importr('vegan')    
    ri_df_piv = pandas2ri.py2ri_pandasdataframe(PIV)
    ri_df_lab = pandas2ri.py2ri(labels)
    
    # run similarity percentages
    sim = vegan.simper(ri_df_piv, ri_df_lab)
    result_ids = pandas2ri.ri2py(sim[0].names)    
    species = pandas2ri.ri2py(sim[0][next((i for i, j in enumerate(list(result_ids=='species')) if j), None)])
    sd      = pandas2ri.ri2py(sim[0][next((i for i, j in enumerate(list(result_ids=='sd')) if j), None)])
    
    sim_df = pd.DataFrame({'species':species, 'simper_sd':sd})
    sim_df = sim_df.set_index('species')
    return sim_df

def rank_from_R(sd_column):
    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri#, default_converter
    base = importr('base')
    sd_column = -1* sd_column
    ri_df_sd = pandas2ri.py2ri(sd_column)
    rank_result = base.rank(ri_df_sd, na_last = True, ties_method = "min")
    rank_result = pandas2ri.ri2py(rank_result)
    return rank_result

def assign_labels(PIV, CLUSTERS, character_stat_scores):
    import pandas as pd
    import numpy as np
    n_clusters = len(np.unique(CLUSTERS['Labels']))

    CLUSTERS = CLUSTERS.set_index('Haul_IDs')
    Ch_Cl = PIV.join(CLUSTERS['Labels'])
    Ch_Cl = Ch_Cl.dropna(subset=["Labels"])
   
    list_of_subs=[]
    list_of_completes=[]
    
    for icluster in np.unique(CLUSTERS['Labels']):
        print('Cluster ',icluster)
        cl_count, cl_sum = site_and_sum_cl(Ch_Cl,icluster)
        cl_count = cl_count.drop(index='Labels')
        cl_sum = cl_sum.drop(index='Labels')
        total_count, total_sum = species_total_presence(Ch_Cl.loc[:, Ch_Cl.columns != 'Labels'])
        
        # Diss:
        group = Ch_Cl.copy()
        group["Labels"] += 1
        group.loc[group["Labels"] != (icluster+1), "Labels"] = 0
        # similarity percentages
        labels = group['Labels']
        sim_df = simper(Ch_Cl.loc[:, Ch_Cl.columns != 'Labels'],labels)
        # there area points in species names -> replace to spaces
        sim_df.index = total_count.index
        
        n_species_in_cluster = sum(group["Labels"]>0)
        
        character_df =  sim_df.join(pd.DataFrame({"total_count":total_count}))\
            .join(pd.DataFrame({"total_sum":total_sum}))\
            .join(pd.DataFrame({"cl_count":cl_count}))\
            .join(pd.DataFrame({"cl_sum":cl_sum}))
        
        # apply stats ----------
        # ND Numerische Dominanz: Abundanz einer Art / Gesamtabundanz im Cluster
        character_df["ND"] = character_df["cl_sum"] / sum(character_df["cl_sum"])

        # Abundanztreue AT: Individuenzahl einer Art in einem Cluster / gesamte Individuenzahl der Art im UG
        character_df["AT"] = character_df["cl_sum"].div(character_df["total_sum"])

        # PrÃ¤senz P: Anteil der Stationen innerhalb eines Clusters, an der die Art gefunden wurde.
        character_df["P"] = character_df["cl_count"] / n_species_in_cluster

        # PrÃ¤senztreue PT: Anzahl Stationen, an denen eine Art in einem Cluster vorkommt / Gesamtzahl der Stationen, an denen die Art im Untersuchungsgebiet vorkommt
        character_df["PT"] = character_df["cl_count"].div(character_df["total_count"])

        # Trennarten hohen Ranges nach DissimilaritÃ¤ts-Analyse (Clark 1993)
        character_df["T"] = rank_from_R(character_df["simper_sd"])

        # score stats ----------
        # Abfrage ND in einer neuen Spalte ND_Q (hier ND >= 3%); gibt alle wahren FÃ¤lle mit eines aus, die Ã¼brigen mit null.
        character_df["ND_Q"] = (character_df["ND"] >= character_stat_scores["ND"]).astype(int)
        
        # Abfrage AT in einer neuen Spalte AT_Q (hier ND >= 60%); gibt alle wahren FÃ¤lle mit eines aus, die Ã¼brigen mit null.
        character_df["AT_Q"] = (character_df["AT"] > character_stat_scores["AT"]).astype(int)

        # Abfrage P in einer neuen Spalte P_Q (hier P >= 60%); gibt alle wahren FÃ¤lle mit eines aus, die Ã¼brigen mit null.
        character_df["P_Q"] = (character_df["P"] > character_stat_scores["P"]).astype(int)
                
        # Abfrage PT in einer neuen Spalte PT_Q (hier PT >= 60%); gibt alle wahren FÃ¤lle mit eines aus, die Ã¼brigen mit null.
        character_df["PT_Q"] = (character_df["PT"] > character_stat_scores["PT"]).astype(int)

        # Abfrage T in einer neuen Spalte T_Q (hier T <= 8; sprich die hÃ¶chsten acht RÃ¤nge); gibt alle wahren FÃ¤lle mit eines aus, die Ã¼brigen mit null.
        character_df["T_Q"] = (character_df["T"] <= character_stat_scores["T"]).astype(int)

        # Bildung einer neuen Spalte, in der ND_Q, AT_Q, P_Q, PT_Q, T_Q aufsummiert werden. 
        character_df["Sum_Q"] = character_df['ND_Q'] +character_df['AT_Q'] + \
            character_df['P_Q'] + character_df['PT_Q'] + character_df['T_Q']

        list_of_completes.append(character_df)

        character_df_sub = character_df[character_df["Sum_Q"] >= character_stat_scores["Threshold_Sum"]]
        if(character_df_sub.empty):
            list_of_subs.append("character_df_sub empty for Cluster "+str(icluster))            
        else:
             list_of_subs.append(character_df_sub)
        
    return list_of_subs, list_of_completes