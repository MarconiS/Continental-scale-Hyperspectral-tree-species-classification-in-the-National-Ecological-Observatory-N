#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 18:27:38 2021

@author: sergiomarconi
"""
import pandas as pd
from sklearn.preprocessing import normalize
from src.functions_brdf import *

def kld_transform(hsi_cube, kld_out):
    #brick = brick.values
    kld_groups = pd.read_csv(kld_out, header=None)
    kld_groups = kld_groups.rename(columns={0: "_kld_grp"})
    all_data = np.zeros([hsi_cube.shape[0],1])
    for jj in np.unique(kld_groups):
        which_bands = kld_groups._kld_grp == jj
        #min
        new_col = np.apply_along_axis(min, 1, hsi_cube[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
        #mean
        new_col = np.apply_along_axis(np.mean, 1, hsi_cube[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
        #max
        new_col = np.apply_along_axis(max, 1, hsi_cube[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
        #
    all_data = pd.DataFrame(all_data)
    all_data = all_data.drop([0], axis=1)
    return(all_data)


    # #shift back individualID to first row
    # cols = list(all_data.columns)
    # cols = [cols[-1]] + cols[:-1]
    # all_data = all_data[cols]
    # return all_data



# predict species labels for dbh dataset
hsi = pd.read_csv("/Volumes/Data/all_predictions_data.csv")
meta = pd.read_csv("/Volumes/Data/Spectral_library/Archive/all_data_metadata.csv")
#coeffs_brdf = pd.read_csv("./data/corrections/angle_all_all_brdf_coeffs_1920.csv")
#topo_coeff = pd.read_csv("./data/corrections/topo_all_all_brdf_coeffs_1920.csv")
#ave_coords = pd.read_csv("/Volumes/Data/Spectral_library/sites_coordinates.csv")
#site_effect = pd.read_csv("/Volumes/Data/Species_Classification/Outuputs_model/site_encoder.csv")
#bands_used = pd.read_csv("/Volumes/Data/Spectral_library/brdf_spectra_2103_all0.csv", nrows= 1).columns 
#transform into BRDF corrected tile
kld_out = "/Volumes/Data/Species_Classification/Outuputs_model/kld_grps_2104b.csv"

#hsi = pd.merge(hsi, ave_coords, left_on='site', right_on='siteID')
#hsi = pd.merge(hsi, site_effect, left_on=['latitude', 'longitude'], right_on=['latitude', 'longitude'])

coords = hsi[['latitude', 'longitude', 'siteID']]
elevation = hsi[['elevation']]
#brdf = calculate_brdf(hsi, coeffs_brdf, topo_coeff)
identifiers = hsi[['individualID']][~which_pixels_dropped]
hsi.columns
brdf = hsi.iloc[:,5:352].to_numpy()

del(hsi)


#filter for greennes and shadows
ndvi = (brdf[:,89]- brdf[:,57])/(brdf[:,57] + brdf[:,89]) <0.5
nir860 = (brdf[:,95] + brdf[:,96])/2 <0.1
brdf[ndvi,:]= np.NAN
brdf[nir860,:]= np.NAN

# apply KLD reduction
#brdf = brdf[:,10:357]

# map pixels with NAN
which_pixels_dropped = np.isnan(brdf).any(axis=1)
brdf = brdf[~np.isnan(brdf).any(axis=1)]
#remove bad pixels 
nbands = brdf.shape[1]
brdf = normalize(brdf)
brdf = kld_transform(brdf, kld_out)


elevation = elevation[~which_pixels_dropped]
coords = coords[~which_pixels_dropped]
# load model
import joblib
clf_bl = joblib.load('/Volumes/Data/Species_Classification/Outuputs_model/model_BRDF_April_model.pkl')


#make prediction at pixels
taxa_classes = clf_bl.meta_clf_.classes_
#brdf = pd.DataFrame(brdf[:,1:46])
brdf.reset_index(drop=True, inplace=True)
elevation.reset_index(drop=True, inplace=True)
coords.reset_index(drop=True, inplace=True)

features = pd.concat([elevation, brdf, coords], axis=1)
#brdf.columns = ['elevation',np.array(1:45),  'latitude', 'longitude']
predict_an_ = clf_bl.predict_proba(features)
predict_an_ = pd.DataFrame(predict_an_)
predict_an_.columns = taxa_classes

# append crownID and make prediction at crown level
identifiers = hsi[['individualID', 'year']][~which_pixels_dropped]
identifiers = identifiers
identifiers.reset_index(drop=True, inplace=True)
predict_an_.reset_index(drop=True, inplace=True)
pp2 = pd.concat([identifiers, predict_an_], axis=1)
pp2.to_csv("/Volumes/Data/Species_Classification/all_out.csv")
identifiers.to_csv("/Volumes/Data/Species_Classification/all_per_year.csv")
