# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 14:38:19 2022

"""

import pandas as pd
import geopandas as gpd

def multiple_excel_to_gdf(list_of_filenames):
    GDF = []
    for ifile in list_of_filenames:
        gdf = pd.read_excel(ifile,sheet_name='Messdaten')
        gdf = recude_columns_to_specified(gdf)
        if ifile == '../01_Daten/BioConsult_Jun22/DatenExport_AWI_critterbase_20220623_MK_TB.xlsx':
            gdf = generate_stid_awi_data(gdf)
        else:
            gdf['ST_ID'] = gdf['ST_ID'].apply(str)
        gdf = convert_to_numeric(gdf, ["GEOLONSTARTHL","GEOLATSTARTHL","WETWPROBE","SPECIMENHAUL"])
        if(gdf['PARALLELNO'].dtype != 'float64'):
            gdf['PARALLELNO'] = gdf['PARALLELNO'].astype('float64')
        gdf = pd_to_gpd(gdf)
        GDF.append(gdf)
    return GDF

def generate_stid_awi_data(df):
    df['ST_ID'] = df['Haul_ID'].apply(lambda x: "-".join(x.split('-')[:-1]))
    return df

def excel_to_gdf(filename,sheet_name='Messdaten'):
    df = pd.read_excel(filename,sheet_name=sheet_name)
    df = recude_columns_to_specified(df)
    gdf = pd_to_gpd(df)
    return gdf

def recude_columns_to_specified(df):
    df2 = df.copy()
    df2 = df2[['ST_ID','Haul_ID','PARALLELNO','GEOLATSTARTHL','GEOLONSTARTHL','SAMPLINGINSTR','ScientificName_Accepted','APHIA_ID_accepted','SPECIMENHAUL','WETWPROBE']]
    return df2

def convert_to_numeric(haul_db,col_list=['SPECIMENHAUL','WETWPROBE']):
    for icol in col_list:
        try:
            haul_db[icol] = pd.to_numeric(haul_db[icol], errors='coerce')
        except:
            print('all already numeric')
    return haul_db

def pd_to_gpd(df,collon='GEOLONSTARTHL',collat='GEOLATSTARTHL',crs='EPSG:4326'):
    gdf = gpd.GeoDataFrame(
        df, geometry=gpd.points_from_xy(df[collon], df[collat],crs=crs))
    return gdf
