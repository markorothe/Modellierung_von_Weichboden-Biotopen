# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 13:19:05 2022

@author: marko
"""

# PERFORM CLUSTERING ----------------------------------------------------------
def run_specific_cluster_settings(settings_dict):
    import pandas as pd
    import geopandas as gpd
    from src import data_filtering_functions
    
    LABELS = []
    NAMES = []
    haul_ids = []
    
    for iregion in settings_dict.keys():
        print('running settings labelled as ',iregion)
        HAUL_DB = pd.read_pickle('HAUL_DB.pkl')
        #HAUL_DB = HAUL_DB.sample(frac=0.01) # ------------------------- <<<
        region = iregion.split('-')[0]
        HAUL_DB = data_filtering_functions.clip_geomety(HAUL_DB, clipping_object = region)
        PIV = HAUL_DB.pivot_table(index='Haul_ID',columns='ScientificName_Accepted', values='SPECIMENHAUL',aggfunc='sum', fill_value=0)
        PIV, GEOMETRY = data_filtering_functions.sampling_haul_data(PIV, HAUL_DB, 'first')
        
        iconsist = settings_dict[iregion]['consistency']
        imethod = settings_dict[iregion]['method']
        incluster = settings_dict[iregion]['ncluster']
        ifuzzy = settings_dict[iregion]['fuzzyfication']
        itransf = settings_dict[iregion]['transformation']
        
        NAMES.append(iregion)
        if itransf == 'HF_Euc':
            PIV_tr = data_filtering_functions.apply_hellinger_transformation(PIV)
            DIST = distance_euclidian(PIV_tr)
        elif itransf == '4th_BC':
            PIV_tr = data_filtering_functions.apply_fourth_root(PIV)
            DIST = distance_bray_curtis(PIV_tr)
        if imethod == 'AvgL':
            labels, _ = perform_avg_linkage(DIST, PIV_tr, incluster)
            show_clustering_in_map(GEOMETRY, labels, iregion+'_'+imethod+'_'+itransf+'_'+str(iconsist)+'_'+str(incluster)+'.png')
        elif imethod == 'FKM':                    
            labels, _, part_coef = perform_fkm(DIST, PIV_tr, incluster, ifuzzy)
            show_clustering_in_map(GEOMETRY, labels, iregion+'_'+imethod+'_'+itransf+'_'+str(iconsist)+'_'+str(incluster)+'_'+str(ifuzzy)+'.png')
        LABELS.append([x for x in labels['Labels']])
        haul_ids.append([x for x in GEOMETRY.index])
        GEOMETRY['cl_ID'] = [x for x in labels['Labels']]
        GEOMETRY.to_file('Geo'+iregion+'.json',driver='GeoJSON')
              
    LABELS = pd.DataFrame({'Calculation':NAMES, 'Haul_IDs':haul_ids, 'Cluster-IDs':LABELS})
    LABELS = LABELS.set_index('Calculation')
    LABELS.to_pickle('LABELS.pkl') 
        
#----------------------------------------------------------------------------------------------------------
'''
def large_loop_through_aois():
    import pandas as pd
    from src import data_filtering_functions
    
    for iregion in ['NW','NE','S','E']:
        #S: [633 rows x 458 columns] GEOMETRY.to_file('geo_S.json',driver='GeoJSON')
        #E: [702 rows x 520 columns]
        #NW:[671 rows x 661 columns]
        #NE:[1114 rows x 790 columns]
        HAUL_DB = pd.read_pickle('HAUL_DB.pkl')
        HAUL_DB = data_filtering_functions.replace_placeholders(HAUL_DB)

        sampling = 1
        #HAUL_DB = HAUL_DB.sample(frac=sampling)
        
        HAUL_DB = data_filtering_functions.clip_geomety(HAUL_DB, clipping_object = iregion)
        PIV = HAUL_DB.pivot_table(index='Haul_ID',columns='ScientificName_Accepted', values='SPECIMENHAUL',aggfunc='sum', fill_value=0)
        PIV, GEOMETRY = data_filtering_functions.sampling_haul_data(PIV, HAUL_DB, 'first')
        
        for iconsist in [0.1,0.05,0.02]:
            PIV = data_filtering_functions.consistency_criterium_spatial(PIV, GEOMETRY, iconsist, tile_dxy=5000)
            
            for ivariant in [0,1]:
                if ivariant == 0:
                    print(iregion+'--- Consistency:'+str(iconsist)+'--- AvL------------------------------')
                    PIV_tr = data_filtering_functions.apply_fourth_root(PIV)
                    DIST = distance_bray_curtis(PIV_tr)
                    test_for_n_cluster_avl(DIST, PIV_tr, GEOMETRY, n1=2, n2=20, dn=1, comment=iregion+'_Consist '+str(iconsist)+'_'+str(sampling*100)+'%, 4th, BC, AvgL')
                else:                    
                    print(iregion+'--- Consistency:'+str(iconsist)+'--- FKM------------------------------')
                    PIV_tr = data_filtering_functions.apply_hellinger_transformation(PIV)
                    DIST = distance_euclidian(PIV_tr)
                    test_for_n_cluster_fkm(DIST, PIV_tr, GEOMETRY, fuzzylevels=[1.1,1.2,1.3], clustersizes=[2,3,4,5], comment=iregion+'_Consist '+str(iconsist)+'_'+str(sampling*100)+'%, HF, Eucl, FKM')
                if ivariant == 0:
                    print(iregion+'--- Consistency:'+str(iconsist)+'--- AvL------------------------------')
                    PIV_tr = data_filtering_functions.apply_hellinger_transformation(PIV)
                    DIST = distance_euclidian(PIV_tr)
                    test_for_n_cluster_avl(DIST, PIV_tr, GEOMETRY, n1=2, n2=20, dn=1, comment=iregion+'_Consist '+str(iconsist)+'_'+str(sampling*100)+'%, HF, Eucl, AvgL')
                else:                    
                    print(iregion+'--- Consistency:'+str(iconsist)+'--- FKM------------------------------')
                    PIV_tr = data_filtering_functions.apply_fourth_root(PIV)
                    DIST = distance_bray_curtis(PIV_tr)
                    test_for_n_cluster_fkm(DIST, PIV_tr, GEOMETRY, fuzzylevels=[1.1,1.2,1.3], clustersizes=[2,3,4,5,6,7,8], comment=iregion+'_Consist '+str(iconsist)+'_'+str(sampling*100)+'%, 4th, BC, FKM')
'''

#----------------------------------------------------------------------------------------------------------
def save_results_from_spatial_consistency_criterium():
    import pandas as pd
    from src import data_filtering_functions
    
    for iregion in ['NW','NE','S','E']:

        HAUL_DB = pd.read_pickle('HAUL_DB.pkl')
        HAUL_DB = data_filtering_functions.replace_placeholders(HAUL_DB)
        HAUL_DB = data_filtering_functions.clip_geomety(HAUL_DB, clipping_object = iregion)
        PIV = HAUL_DB.pivot_table(index='Haul_ID',columns='ScientificName_Accepted', values='SPECIMENHAUL',aggfunc='sum', fill_value=0)
        PIV, GEOMETRY = data_filtering_functions.sampling_haul_data(PIV, HAUL_DB, 'first')
        PIV.to_pickle('tmp.pkl')
        GEOMETRY.to_pickle('tmp2.pkl')
        result = pd.DataFrame(
            data = {"region":[iregion],"consistency":[0],"n_spec":[len(PIV.columns.values)],"species":[PIV.columns.values],})
        result.to_csv("consistency_crit_tab.csv",mode='a', index=False, header=False)
        print(iregion+'   '+str(0)+'   '+str(len(PIV.columns.values)))
           
        grid_reduced, _ = data_filtering_functions.binning(GEOMETRY, 5000) 
        grid_reduced.to_file(iregion+'grid_binned_5000m.json',driver='GeoJSON')

        for iconsist in [0.5,0.4,0.3,0.2,0.1,0.05,0.02,0.01]:
            PIV = pd.read_pickle('tmp.pkl')
            GEOMETRY = pd.read_pickle('tmp2.pkl')
            PIV = data_filtering_functions.consistency_criterium_spatial(PIV, GEOMETRY, iconsist, tile_dxy=5000)
            print(iregion+'   '+str(iconsist)+'   '+str(len(PIV.columns.values)))
            result=pd.DataFrame(
                data = {"region":[iregion],"consistency":[0],"n_spec":[len(PIV.columns.values)],"species":[PIV.columns.values],})
            result.to_csv("consistency_crit_tab.csv",mode='a', index=False, header=False)
            

#----------------------------------------------------------------------------------------------------------

def test_for_n_cluster_fkm(dist_matrix, INPUT, GEOMETRY, fuzzylevels=[1.25,1.5,1.75], clustersizes=[2,3,4,5,6,7], comment=''):
    import pandas as pd
    for fuzzyness in fuzzylevels:
        for n_cluster in clustersizes:
            print("Fuzzyness:",fuzzyness,"N:",n_cluster)
            labels, silhuette_index, part_coef = perform_fkm(dist_matrix, INPUT, n_cluster, fuzzyness)
            result=pd.DataFrame(
                data={"Comment":[comment],
                      "NCluster":[n_cluster],
                      "Fuzzyfication":[fuzzyness],
                      "Silhuette-Index":[silhuette_index],
                      "Partition-Coeff":[part_coef],
                      #"Result":[[str(x).zfill(4) for x in labels]],
                      "Time":[pd.Timestamp.now()]           
                      })
            result.to_csv("large_tab.csv",mode='a', index=False, header=False)
            #show_clustering_in_map(GEOMETRY, labels, comment+'_fkm_N'+str(n_cluster)+'_fuzzy_'+str(fuzzyness)+'.png')
            

def test_for_n_cluster_avl(dist_matrix, INPUT, GEOMETRY, n1=2, n2=20, dn=1, comment=''):
    import pandas as pd
    for n_cluster in range(n1,n2+1,dn):
        labels, silhuette_index = perform_avg_linkage(dist_matrix, INPUT, n_cluster)
        result=pd.DataFrame(
            data={"Comment":[comment],
                  "NCluster":[n_cluster],
                  "Silhuette-Index":[silhuette_index],
                  #"Result":[[str(x).zfill(4) for x in labels['Labels']]],
                  "Time":[pd.Timestamp.now()]            
                  })
        result.to_csv("large_tab.csv",mode='a', index=False, header=False)
        #show_clustering_in_map(GEOMETRY, labels, 'avglin_N'+str(n_cluster)+'.png')

# ---------- 2nd approach: b-c-distances of 4th sq. roots

def distance_bray_curtis(df):
    # https://docs.scipy.org/doc//scipy-1.4.0/reference/generated/scipy.spatial.distance.cdist.html
    from scipy.spatial.distance import cdist
    import numpy as np
    distances = cdist(df,df,metric='braycurtis')
    #distances = np.tril(distances, -1) #reduce to lower triangle
    #distances[distances==0] = np.nan
    distances = np.nan_to_num(distances, nan=np.finfo(float).eps)
    return distances

def distance_euclidian(df):
    from sklearn.metrics.pairwise import euclidean_distances
    import numpy as np
    values = np.array(df.values)
    distances = euclidean_distances(values,values)
    #distances = np.tril(distances, -1) #reduce to lower triangle
    #distances[distances==0] = np.nan
    return distances

# ---------- 3rd approach.. edited KMeans/FuzzyKMeans functions

def part_coef(membershipmatrix_df):
    from rpy2.robjects.packages import importr
    from rpy2.robjects import numpy2ri, default_converter
    from rpy2.robjects.conversion import Converter, localconverter 
    numpy2ri.activate()     
    fclust = importr('fclust')
    U = membershipmatrix_df.values
    with localconverter(default_converter + numpy2ri.converter) as cv:          
        U_r = numpy2ri.py2ri(U)
    part_coef = float(numpy2ri.ri2py(fclust.PC(U_r)))
    return part_coef

def silhuette_index(dist_matrix, labels):
    from sklearn.metrics import silhouette_score 
    silhouette_avg = silhouette_score(dist_matrix, labels, metric='precomputed') 
    return silhouette_avg

def perform_fkm(dist_matrix, INPUT, n_cluster, fuzziness_factor):
    import numpy as np
    import pandas as pd
    from src import fuzzykmedoids as fkm
    
    # code doesn't perform when pariwise distances are 0
    dist_matrix[dist_matrix==0]=np.finfo(float).eps

    _, membership = fkm.FKM(dist_matrix,n_cluster,fuzziness_factor)
    membership=np.array(membership)

    medoid_ids = np.where(membership[0,]>0) # center of clusters

    # perform to better result form
    membership_dense=membership[np.ix_(~np.all(membership==False,axis=1),
                                       ~np.all(membership==False,axis=0))]

    membership_df = pd.DataFrame(membership_dense)
    membership_df.columns=[x for x in list(medoid_ids[0])]
    # calculate partition coefficient
    partition_coef = part_coef(membership_df) 

    membership_df['ordered_indices'] = membership_df.apply(lambda x: x.nlargest(2).index,axis=1)
    membership_df['1st_association'] = membership_df['ordered_indices'].apply(lambda x: x[0])
#    membership_df['2nd_association'] = membership_df['ordered_indices'].apply(lambda x: x[1])
#    share_1st=[]; share_2nd=[]
#    for irow in range(0,len(membership_df)):
#        share_1st.append(membership_df[membership_df['1st_association'].iloc[irow]].iloc[irow])
#        share_2nd.append(membership_df[membership_df['2nd_association'].iloc[irow]].iloc[irow])
#        membership_df['1st_association_share'] = share_1st
#        membership_df['2nd_association_share'] = share_2nd
    del membership_df['ordered_indices']
    
    # assign labels to haul_ids
    labels_df = pd.DataFrame(data={'Valid_Row':membership_df.index, 'MedoidID':membership_df['1st_association']})
    tmp = INPUT.iloc[membership_df.index]
    labels_df['Haul_ID'] = tmp.index
    labels_df = labels_df.set_index('Haul_ID')
    del labels_df['Valid_Row']

    tmp_num=np.arange(len(np.unique(labels_df['MedoidID'])))
    tmp_lab=np.unique(labels_df['MedoidID'])
    mapfile={}
    for ix in range(0,len(tmp_num)):
        mapfile[tmp_lab[ix]]=tmp_num[ix]
    labels_df['Labels'] = labels_df['MedoidID'].map(mapfile)

    #calulate silhuette_index
    silh_indx = silhuette_index(dist_matrix, labels_df['Labels'])
        
    return labels_df, silh_indx, partition_coef


# ---------- 3rd approach.. edited KMeans/FuzzyKMeans functions

def perform_avg_linkage(dist_matrix, INPUT, n_cluster):
    from sklearn.cluster import AgglomerativeClustering
    import numpy as np
    import pandas as pd
    
    DIST=dist_matrix
    mask_row = np.any(np.isnan(dist_matrix), axis=1)
    mask_column = np.any(np.isnan(dist_matrix), axis=0)
    dist_matrix = dist_matrix[:, ~mask_column]
    dist_matrix = dist_matrix[~mask_row, :]

    valid_rows = np.array([i for i in range(0,DIST.shape[1])])
    valid_rows = valid_rows[[int(i) for i in np.where(mask_row == False)[0]]]

    clustering = AgglomerativeClustering(n_clusters = n_cluster,
                                         affinity='precomputed',linkage='average')
    result=clustering.fit(dist_matrix)
    
    # assign labels to haul_ids
    labels_df = pd.DataFrame(data={'Valid_Row':valid_rows, 'Labels':result.labels_})
    tmp = INPUT.iloc[valid_rows]
    labels_df['Haul_ID'] = tmp.index
    labels_df = labels_df.set_index('Haul_ID')
    del labels_df['Valid_Row']
    
    #calulate silhuette_index
    silh_indx = silhuette_index(dist_matrix, labels_df['Labels'])
        
    return labels_df, silh_indx
    

def show_clustering_in_map(GEOMETRY, labels_df, pngname=''):
    from mycolorpy import colorlist as mcp
    import matplotlib.pyplot as plt
    #import contextily as cx
    import numpy as np

    GEOMETRY = GEOMETRY.join(labels_df)

    fig, (ax1,ax2) = plt.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]}, figsize = (13,10))
    df_wm = GEOMETRY.to_crs(epsg=3857)
    df_wm.plot(column='Labels',markersize=3,ax=ax1,cmap='nipy_spectral')
    colors=mcp.gen_color(cmap="nipy_spectral",n=np.max(labels_df['Labels'])+1)
    #cx.add_basemap(ax1,source=cx.providers.OpenSeaMap)    
    ax1.set_title('Cluster distribution')
    
    x=GEOMETRY.groupby(by='Labels').agg('sum')
    x['Count'].plot(kind="bar", color=colors, ax=ax2)
    ax2.set_title('Cluster histogram')
    
    if pngname == '':
        plt.show()
    else:
        plt.savefig(pngname)
        