# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 15:39:45 2022

"""

import pandas as pd
import numpy as np
import geopandas as gpd
pd.set_option('display.max_rows', 100)
   
# --------------------------------------------------------------------------
def apply_specific_filter_haul_db(haul_db):
    pd.options.mode.chained_assignment = None
    print('DB-entries: '+str(len(haul_db)))
#
    print('Ignore samples without geometry')
    haul_db = haul_db[haul_db['geometry'].is_empty == False]
    print('DB-entries: '+str(len(haul_db)))
#
    print('Ignore Phoronis entries')
    haul_db = haul_db[haul_db['ScientificName_Accepted'].str.contains('Phoronis')==0] #all species "Phoronis *" to be  ignored
    print('DB-entries: '+str(len(haul_db)))
#
    print('Ignore dublicates')
    haul_db = haul_db.drop_duplicates(subset=['ST_ID','Haul_ID','PARALLELNO','GEOLATSTARTHL','GEOLONSTARTHL','SAMPLINGINSTR','ScientificName_Accepted','APHIA_ID_accepted','SPECIMENHAUL','WETWPROBE'], keep='first')
    print('DB-entries: '+str(len(haul_db)))
#
    haul_db = haul_db.reset_index()
    return haul_db

# --------------------------------------------------------------------------
def look_for_placeholders(haul_db,col_list=['SPECIMENHAUL','WETWPROBE']):
    for cid in col_list:
        if cid not in ['Haul_ID','geometry','SAMPLINGINSTR','ScientificName_Accepted','ST_ID']:
            print(cid+':')
            unq, counts = np.unique(haul_db[cid].values,return_counts=True)
            df =pd.DataFrame({'value':unq,'count':counts})
            df = df.dropna()
            df = df.sort_values(by='value')
            print(df.head(30))
            print(df.tail(30))
            
# --------------------------------------------------------------------------
def look_for_dublicates(haul_db,col_list=['ST_ID','Haul_ID','PARALLELNO','GEOLATSTARTHL','GEOLONSTARTHL','SAMPLINGINSTR','ScientificName_Accepted']):
    x=haul_db[haul_db.duplicated(keep=False)]
    x=x.sort_values(by=col_list)
    x.to_csv('dubl.csv',index=False)
            

# --------------------------------------------------------------------------
def replace_placeholders(haul_db, dictionary={}):
    # replace certain placeholders in columns specimenhaul1 and wetwprobe1, 
    # values of -9 and 0 are replaced by NaN
    if dictionary=={}:
        dictionary={'SPECIMENHAUL':{-991.0:1, -333.0:1, -99.0:1, -9.0:np.nan, 0.0:1, 0.1:1, 0.2:1, 0.3:1, 0.4:1, 0.5:1, 0.7:1, 0.8:1},
                    'WETWPROBE':{-999000:np.nan, -99000:0.05, -9999.0:1, -9979.0:0.01, -666.0:np.nan,-9.0:np.nan, 0.0:np.nan}}
    colnames = [key for key in dictionary]
    for icol in colnames:
        haul_db = haul_db.replace({icol: dictionary[icol]})
    return haul_db

# --------------------------------------------------------------------------
def apply_hellinger_transformation(haul_db):
    from src import external_ecopy_transformations as ectr
    haul_db = ectr.transform(haul_db, axis=1, method='hellinger', breakNA=False)
    haul_db = haul_db.replace(np.nan,0)
    return haul_db

# --------------------------------------------------------------------------
def apply_fourth_root(piv):
    def fourthroot(x):
        return x**(1/4)
    piv = piv.apply(fourthroot)
    return piv

# --------------------------------------------------------------------------
def sampling_haul_data(PIV, HAUL_DB, choice='first'):
    import random
    # haul - station - connection
    tmp_haul_db = HAUL_DB[['Haul_ID','ST_ID','geometry']].dropna(subset=['Haul_ID'])
    # get count of repeated hauls
    tmp_stations = list(tmp_haul_db['ST_ID'])
    tmp_stations_counts = {item:tmp_stations.count(item) for item in tmp_stations}
    tmp_haul_db['Count'] = tmp_haul_db['ST_ID'].map(tmp_stations_counts)
    # get list of haul_ids for each station    
    tmp_samples_db = tmp_haul_db.groupby('ST_ID')['Haul_ID'].apply(list).reset_index(name='All_Hauls')
    #return custom coice of selection
    tmp_list = list(tmp_samples_db['All_Hauls'])
    if choice == 'random':
        tmp_samples_db['Selected_Haul'] = [random.sample(item,1)[0] for item in tmp_list]
    elif choice == 'first':
        tmp_samples_db['Selected_Haul'] = [item[0] for item in tmp_list]       
    
    #apply chosen hauls
    PIV = PIV[PIV.index.isin(tmp_samples_db['Selected_Haul'])]
    GEOMETRY = tmp_haul_db.set_index('Haul_ID')
    GEOMETRY = GEOMETRY[~GEOMETRY.index.duplicated(keep='first')]
    GEOMETRY = GEOMETRY.loc[PIV.index]
    return PIV, GEOMETRY

# --------------------------------------------------------------------------
# Für Arten: stetigkeitsprinzip: min 10% der Stationen
def consistency_criterium_percentage(PIV, perc_value=0.1):
    n_stations = len(PIV.index)
    precence_of_species = (PIV>0).sum(axis=0)
    keep_species = precence_of_species.index[precence_of_species >= perc_value*n_stations]
    PIV = PIV[keep_species.values.tolist()]
    return PIV

# --------------------------------------------------------------------------
def consistency_criterium_spatial(PIV, GEOMETRY, perc_value=0.1, tile_dxy=5e3):
    # spatial binning in tiles, 5km resolution
    _, assigned_tile_ids = binning(GEOMETRY, tile_dxy) 
    n_tiles = len(np.unique(assigned_tile_ids))
    
    help_piv = PIV.copy()
    help_piv['tile_id'] = assigned_tile_ids
    help_piv = help_piv.groupby(['tile_id'], as_index=True).agg(sum)
    precence_of_species = (help_piv>0).sum(axis=0)
    keep_species = precence_of_species[precence_of_species >= perc_value*n_tiles]

    PIV = PIV[list(keep_species.index)]
    return PIV


# spatial consistency criterium -----------------------------------------------
def binning(gdf, dxy):
    import math
    import numpy as np
    from shapely.geometry import Polygon
    
    # transform from geographic to UTM for square grid cells
    gdf=gdf.to_crs(25832)
    xmin, ymin, xmax, ymax = gdf.total_bounds
    xmin = math.floor(xmin/dxy)*dxy
    ymin = math.floor(ymin/dxy)*dxy
    cols = list(np.arange(xmin, xmax + dxy, dxy))
    rows = list(np.arange(ymin, ymax + dxy, dxy))

    tiles = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            tiles.append(Polygon([(x,y), (x+dxy, y), (x+dxy, y+dxy), (x, y+dxy)]))

    grid = gpd.GeoDataFrame({'geometry':tiles})
    grid = grid.set_crs(gdf.crs, allow_override=True)

    grid_assigned=gpd.sjoin(gdf,grid)
    assigned_tile_ids = grid_assigned['index_right'].values
    
    assigned_indices = np.unique(grid_assigned['index_right'].values)
    grid_reduced = grid.iloc[assigned_indices]
    
    return grid_reduced, assigned_tile_ids


# clip with geometry ---------------------------------------------------
def clip_geomety(gdf, clipping_object = 'NW'):
    # clipping_object = ['NW','NE','S','E']
    
    import geopandas as gpd
    inp_eez_regions = '../01_Daten/Geodaten/UGebiet_vierfach.shp'
        
    gdf = gdf.set_crs(4326, allow_override=True)

    gdf_clip = gpd.read_file(inp_eez_regions)
    gdf_clip = gdf_clip[gdf_clip['name']==clipping_object]
        
    gdf_clip = gdf_clip.to_crs("EPSG:4326") 
    clipped = gpd.clip(gdf, gdf_clip)
    
    return clipped


