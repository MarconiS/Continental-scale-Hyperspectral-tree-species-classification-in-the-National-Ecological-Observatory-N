import joblib
import sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
from pathlib import Path
import h5py
import sys
import ogr, os
import gdal, osr
import math
from math import pi
import random
import string
import warnings   

#import arguments
#from the tile path, get the info needed for BRDF correction and append them to a dataframe
full_path = "///orange/ewhite/NeonData/BART/DP1.30006.001/2017/FullSite/D01/2017_BART_3/L3/Spectrometer/Reflectance/NEON_D01_BART_DP3_314000_4881000_reflectance.h5"  
#mod = joblib.load('/blue/ewhite/s.marconi/NeonSpeciesClassification/outputs/model_BRDF_model.pkl')
epsg = 32610
year = 2018
site= "BART"	

model_pkl= sys.argv[1]
pt_kld = sys.argv[2]
year = sys.argv[3]
site= sys.argv[4]
full_path = sys.argv[5]

ras_dir = "/orange/idtrees-collab/species_classification/"


#global_brdf_correction functions
def calculate_geom_kernel(df, sns_zn = 0, ross = "thick", li = "dense"):
	relative_az = df["sns_az"] - df["sol_az"]
	solar_zn_ = np.arctan(10*np.tan(df["sol_zn"]))
	sensor_zn_ = np.arctan(10*np.tan(sns_zn))
	D = np.sqrt((np.tan(solar_zn_)**2) + (np.tan(sensor_zn_)**2) - \
	          2*np.tan(solar_zn_)* np.tan(sensor_zn_)*np.cos(relative_az))    
	# Eq 49. Wanner et al. JGRA 1995
	t_num = 2. * np.sqrt(D**2 + (np.tan(solar_zn_)*np.tan(sensor_zn_)* \
	                           np.sin(relative_az))**2) 
	t_denom = (1/np.cos(solar_zn_))  + (1/np.cos(sensor_zn_))
	t_ = np.minimum(1,np.maximum(t_num/t_denom, -1))
	t = np.arccos(t_)
	# Eq 33,48. Wanner et al. JGRA 1995
	O = (1/np.pi) * (t - np.sin(t)*np.cos(t)) * t_denom
	# Eq 51. Wanner et al. JGRA 1995
	cosPhase_ =  np.cos(solar_zn_)*np.cos(sensor_zn_) + \
	np.sin(solar_zn_)* np.sin(sensor_zn_)* np.cos(relative_az)
	#
	if(li == 'sparse'):
	# Eq 32. Wanner et al. JGRA 1995
	    k_geom = O - (1/np.cos(solar_zn_)) - (1/np.cos(sensor_zn_)) + \
	      0.5*(1+ cosPhase_) * (1/np.cos(sensor_zn_))
	elif(li == 'dense'):
	# Eq 47. Wanner et al. JGRA 1995
	    k_geom = (((1+cosPhase_) * (1/np.cos(sensor_zn_)))/ (t_denom - O)) - 2
	#
	return(k_geom)




def generate_volume_kernel(df, sns_zn = 0, ross = "thick", li = "dense"):
    relative_az = df["sns_az"] - df["sol_az"]
    #Ross kernels 
    ############
    # Eq 2. Schlapfer et al. IEEE-TGARS 2015
    phase = np.arccos(np.cos(df["sol_zn"])*np.cos(sns_zn) + \
               np.sin(df["sol_az"])*np.sin(sns_zn)* np.cos(relative_az))  
    if(ross == 'thick'):
    # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)*np.cos(phase) + \
                 np.sin(phase))/(np.cos(df["sns_zn"]) * np.cos(df["sol_zn"])) - np.pi/4
    elif(ross == 'thin'):
    # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)* np.cos(phase) + \
                 np.sin(phase))/(np.cos(df["sns_zn"])*np.cos(df["sol_zn"])) - np.pi/2
    return(k_vol)




def generate_topographic_coeffs(df):
    relative_az = df["aspect"] - df["sol_az"]
    cos_i = np.cos(df["sol_zn"]) * np.cos(df["slope"])+ np.sin(df["sol_zn"]) * \
      np.sin(df["slope"])*np.cos(relative_az)
    c1 = np.cos(df["sol_zn"])*np.cos(df["slope"])
    #Ross kernels 
    ############
    return(cos_i, c1)




#calculate brdf correction on the tile's DF
def calculate_brdf(hsi, coeffs_brdf, topo_coeff):
    # calculate scattering kernels 
    k_geom = calculate_geom_kernel(hsi, sns_zn = hsi["sns_zn"])
    k_vol = generate_volume_kernel(hsi, sns_zn = hsi["sns_zn"])
    # calculate scattering kernels (at NADIR)
    k_geom_nadir = calculate_geom_kernel(hsi)
    k_vol_nadir = generate_volume_kernel(hsi)
    #generate topographic coefficients
    topo_coeffs = generate_topographic_coeffs(hsi)
    #
   # k_geom = calculate_geom_kernel(hsi, sns_zn = hsi["sns_zn"])
   # k_vol = generate_volume_kernel(hsi, sns_zn = hsi["sns_zn"])        
    # calculate scattering kernels (at NADIR)
   # k_geom_nadir = calculate_geom_kernel(hsi)
   # k_vol_nadir = generate_volume_kernel(hsi)
    #generate topographic coefficients
    cos_i, c1 = generate_topographic_coeffs(hsi)
    #metadata = hsi %>% select(!contains("band"))
    #hsi = hsi %>% select(contains("band"))
    topo = cos_i
    #X = pd.concat([k_vol,k_geom, k_geom_nadir, k_vol_nadir, topo],axis=1)
    #if hsi.shape[1] > 367:
    hsi = hsi.rename(columns={"band_368": "CHM"})
        #
    #if hsi.shape[1] > 369:
    hsi = hsi.filter(regex='band')
        #
    refl = np.zeros(hsi.shape, dtype=np.float)
    for ii in range(refl.shape[1]):      
        y = hsi.iloc[:,ii]/10000
        #apply coefficients to perform correction
        #k_vol+k_geom
        brdf = coeffs_brdf.iloc[ii, 1] * k_vol + \
          coeffs_brdf.iloc[ii, 2] * k_geom  + \
          coeffs_brdf.iloc[ii, 0]
        brdf_nd = coeffs_brdf.iloc[ii, 1] * k_vol_nadir + \
          coeffs_brdf.iloc[ii, 2] * k_geom_nadir  + \
          coeffs_brdf.iloc[ii, 0]
        #calculate BRDF correction
        bdrf_cor = brdf_nd/brdf
        #apply coefficients to perform topographic correction
        topo_cor = (c1 * topo_coeff.iloc[ii, 1] + 
                      topo_coeff.iloc[ii, 0]) /  (topo * \
                      topo_coeff.iloc[ii, 1] + topo_coeff.iloc[ii, 0])   
        #bnd = bnd/10000
        refl[..., ii] = y * bdrf_cor * topo_cor
    #
    return(refl)




def tile_solar_angle(full_path):
	hdf5_file = h5py.File(full_path, 'r')
	file_attrs_string = str(list(hdf5_file.items()))
	file_attrs_string_split = file_attrs_string.split("'")
	sitename = file_attrs_string_split[1]    
	flight_paths = hdf5_file[sitename]["Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index"].attrs["Data_Files"]
	flight_paths=str(flight_paths).split(",")
	which_paths = np.unique(hdf5_file[sitename]["Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index"].value)
	solar_angle = []
	for pt in which_paths:
	    #if pt is negative, get any from the available to avoid error(the pixel is blank anyway)
	    if pt < 0:
	        flight = (flight_paths)[ which_paths[-1]].split("_")[5]
	    else:
	        flight = (flight_paths)[pt].split("_")[5]
	      #  
	    sol_az = hdf5_file[sitename]["Reflectance/Metadata/Logs/"][str(flight)]["Solar_Azimuth_Angle"].value
	    sol_zn = hdf5_file[sitename]["Reflectance/Metadata/Logs/"][str(flight)]["Solar_Zenith_Angle"].value
	    solar_angle.append([pt, sol_az, sol_zn])
	return(solar_angle)
    


def h5refl2array(full_path, epsg):
    #refl, refl_md, wavelengths, sol_az, sol_zn, sns_az, sns_zn, slope, aspect = h5refl2array(full_path, epsg = epsg)
    hdf5_file = h5py.File(full_path, 'r')
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("'")
    sitename = file_attrs_string_split[1]
    epsg = hdf5_file[sitename]["Reflectance/Metadata/Coordinate_System/EPSG Code"].value
    solar_angles = tile_solar_angle(full_path)
    #Extract the reflectance & wavelength datasets
    reflArray = hdf5_file[sitename]['Reflectance']
    refl =reflArray['Reflectance_Data'].value
    wavelengths = reflArray['Metadata']['Spectral_Data']['Wavelength'].value
    # Create dictionary containing relevant metadata information
    refl_md = {}
    refl_md['mapInfo'] = reflArray['Metadata']['Coordinate_System']['Map_Info'].value
    refl_md['wavelength'] = reflArray['Metadata']['Spectral_Data']['Wavelength'].value
    refl_md['shape'] = refl.shape
    #Extract no data value & scale factor
    refl_md['noDataVal'] = float(reflArray['Reflectance_Data'].attrs['Data_Ignore_Value'])
    refl_md['scaleFactor'] = float(reflArray['Reflectance_Data'].attrs['Scale_Factor'])
    #metadata['interleave'] = reflData.attrs['Interleave']
    refl_md['bad_band_window1'] = np.array([1340, 1445])
    refl_md['bad_band_window2'] = np.array([1790, 1955])
    refl_md['epsg'] = str(epsg).split("'")[1]
#    
    #get tiles for BRDF correction
    sns_az = hdf5_file[sitename]['Reflectance/Metadata/to-sensor_azimuth_angle']
    sns_zn = hdf5_file[sitename]['Reflectance/Metadata/to-sensor_zenith_angle']
    slope = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Slope']
    aspect = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Aspect']
    elevation = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Smooth_Surface_Elevation']
#    
    #get solar angles as array to leverage flightpaths mosaic
    flightpaths = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index'].value
    sol_zn = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index'].value
    sol_az = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index'].value
    for pt in range(len(solar_angles)):
        sol_az[flightpaths==solar_angles[pt][0]] = solar_angles[pt][1]
        sol_zn[flightpaths==solar_angles[pt][0]] = solar_angles[pt][2]
#
    mapInfo_string = str(refl_md['mapInfo']);
    mapInfo_split = mapInfo_string.split(",")
    mapInfo_split
#
    # Extract the resolution & convert to floating decimal number
    refl_md['res'] = {}
    refl_md['res']['pixelWidth'] = float(mapInfo_split[5])
    refl_md['res']['pixelHeight'] = float(mapInfo_split[6])
    # Extract the upper left-hand corner coordinates from mapInfo
    xMin = float(mapInfo_split[3])  # convert from string to floating point number
    yMax = float(mapInfo_split[4])
#
    # Calculate the xMax and yMin values from the dimensions
    xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
    yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
    refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
    refl_md['ext_dict'] = {}
    refl_md['ext_dict']['xMin'] = xMin
    refl_md['ext_dict']['xMax'] = xMax
    refl_md['ext_dict']['yMin'] = yMin
    refl_md['ext_dict']['yMax'] = yMax
    hdf5_file.close
#
    return refl, refl_md, sitename, wavelengths, sol_az, sol_zn, sns_az, sns_zn, slope, aspect, elevation


def stack_subset_bands(reflArray, reflArray_metadata, bands, clipIndex):
    subArray_rows = clipIndex['yMax'] - clipIndex['yMin']
    subArray_cols = clipIndex['xMax'] - clipIndex['xMin']
#
    stackedArray = np.zeros((subArray_rows, subArray_cols, len(bands)), dtype=np.int16)
    band_clean_dict = {}
    band_clean_names = []
#
    for i in range(len(bands)):
        band_clean_names.append("b" + str(bands[i]) + "_refl_clean")
        band_clean_dict[band_clean_names[i]] = subset_clean_band(reflArray, reflArray_metadata, clipIndex, bands[i])
        stackedArray[..., i] = band_clean_dict[band_clean_names[i]]
#
    return stackedArray



def subset_clean_band(reflArray, reflArray_metadata, clipIndex, bandIndex):
    bandCleaned = reflArray[clipIndex['yMin']:clipIndex['yMax'], clipIndex['xMin']:clipIndex['xMax'],
                  bandIndex - 1].astype(np.int16)
#
    return bandCleaned



def array2raster(newRaster, reflBandArray, reflArray_metadata, extent, ras_dir, epsg):
    NP2GDAL_CONVERSION = {
        "uint8": 1,
        "int8": 1,
        "uint16": 2,
        "int16": 3,
        "uint32": 4,
        "int32": 5,
        "float32": 6,
        "float64": 7,
        "complex64": 10,
        "complex128": 11,
    }
#
    pwd = os.getcwd()
    os.chdir(ras_dir)
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    bands = reflBandArray.shape[2]
    pixelWidth = float(reflArray_metadata['res']['pixelWidth'])
    pixelHeight = -float(reflArray_metadata['res']['pixelHeight'])
    originX = extent['xMin']
    originY = extent['yMax']
#
    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[reflBandArray.dtype.name]
    print(gdaltype)
    print(newRaster)
    print(cols, rows, bands)
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # outband = outRaster.GetRasterBand(1)
    # outband.WriteArray(reflBandArray[:,:,x])
    for band in range(bands):
        outRaster.GetRasterBand(band + 1).WriteArray(reflBandArray[:, :, band])
#
    outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg'])
    #outRasterSRS.ExportToWkt()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.FlushCache()
    os.chdir(pwd)


def calc_clip_index(clipExtent, h5Extent, xscale=1, yscale=1):
    h5rows = h5Extent['yMax'] - h5Extent['yMin']
    h5cols = h5Extent['xMax'] - h5Extent['xMin']
#
    ind_ext = {}
    ind_ext['xMin'] = round((clipExtent['xMin'] - h5Extent['xMin']) / xscale)
    ind_ext['xMax'] = round((clipExtent['xMax'] - h5Extent['xMin']) / xscale)
    ind_ext['yMax'] = round(h5rows - (clipExtent['yMin'] - h5Extent['yMin']) / yscale)
    ind_ext['yMin'] = round(h5rows - (clipExtent['yMax'] - h5Extent['yMin']) / yscale)
#
    return ind_ext




def extract_hsi_and_brdf_data(full_path, epsg, ras_dir, year,site, ross="thick", li="dense"):
	warnings.filterwarnings("ignore")
	refl, refl_md, sitename, wavelengths, sol_az, sol_zn, sns_az, sns_zn, slope, aspect, elevation = h5refl2array(full_path, epsg = "326")
	rgb = np.r_[0:425]
	rgb = np.delete(rgb, np.r_[419:426])
	rgb = np.delete(rgb, np.r_[281:313])
	rgb = np.delete(rgb, np.r_[191:211])
	xmin, xmax, ymin, ymax = refl_md['extent']
	#   
	clipExtent = {}
	clipExtent['xMin'] = xmin
	clipExtent['yMin'] = ymin
	clipExtent['yMax'] = ymax
	clipExtent['xMax'] = xmax
	print(clipExtent)
	subInd = calc_clip_index(clipExtent, refl_md['ext_dict'])
	subInd['xMax'] = int(subInd['xMax'])
	subInd['xMin'] = int(subInd['xMin'])
	subInd['yMax'] = int(subInd['yMax'])
	subInd['yMin'] = int(subInd['yMin'])
	refl = refl[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax'],rgb]
	sns_az = sns_az[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
	sns_zn = sns_zn[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
	sol_az = sol_az[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
	sol_zn = sol_zn[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
	slope = slope[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
	aspect = aspect[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
	elevation = elevation[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
	# mask away bad pixels
	ndvi = (refl[:, :,90] - refl[:,:,58])/(refl[:, :,58] +refl[:, :,90]) > 0.5
	nir860 = (refl[:, :,96] + refl[:, :,97])/20000 > 0.1
	mask = (sns_zn < 10000)  * (aspect < 10000) *  (slope < 10000) * (sns_az < 10000) * ndvi * nir860
	#
	# convert degrees in radiants
	slope = (slope * pi) / 180
	aspect  = (aspect * pi) / 180
	sns_az  = (sns_az * pi) / 180
	sns_zn  = (sns_zn * pi) / 180
	sol_az = (sol_az * pi) / 180
	sol_zn = (sol_zn * pi) / 180
	#
	subArray_rows = subInd['yMax'] - subInd['yMin']
	subArray_cols = subInd['xMax'] - subInd['xMin']
	hcp = np.zeros((subArray_rows, subArray_cols, len(rgb)), dtype=np.int16)
	#load info in multi-layer array
	band_clean_dict = {}
	band_clean_names = []
	for i in range(len(rgb)):
	    band_clean_names.append("b" + str([i]) + "_refl_clean")
	    band_clean_dict[band_clean_names[i]] = np.squeeze(refl[:, :, [i]].astype(np.int16))
	    hcp[..., i] = band_clean_dict[band_clean_names[i]]
	#
	del(refl, ndvi, nir860)
	tmp = elevation.reshape([elevation.shape[0],elevation.shape[1],1])
	del(elevation)
	hcp = np.concatenate([tmp, hcp], -1)
#
	external_dat = aspect.reshape([aspect.shape[0],aspect.shape[1],1])
	del(aspect)
	tmp = slope.reshape([slope.shape[0],slope.shape[1],1])
	external_dat = np.concatenate([external_dat, tmp], -1)
	del(slope)
	tmp = sns_zn.reshape([sns_zn.shape[0],sns_zn.shape[1],1])
	external_dat = np.concatenate([external_dat, tmp], -1)
	del(sns_zn)
	tmp = sns_az.reshape([sns_az.shape[0],sns_az.shape[1],1])
	external_dat = np.concatenate([external_dat, tmp], -1)
	del(sns_az)
	tmp = sol_zn.reshape([sol_zn.shape[0],sol_zn.shape[1],1])
	external_dat = np.concatenate([external_dat, tmp], -1)
	del(sol_zn)
	tmp = sol_az.reshape([sol_az.shape[0],sol_az.shape[1],1])
	external_dat = np.concatenate([external_dat, tmp], -1)
	del(sol_az, tmp)
	hcp = np.concatenate([hcp, external_dat],-1)
	del(external_dat)
	#save hcp into a tiff file [reflectance]
	#itc_id =  str(int(year)) + "_" +site+"_" + str(int(xmin)) + "_" + str(int(ymin)) 
	#ii = str(itc_id + ".tif")
	#array2raster(str(ii), hcp, refl_md, clipExtent, ras_dir = str(ras_dir), epsg = int(refl_md['epsg']))
	# output the dataframe
	hcp = hcp.reshape(-1,hcp.shape[2])
	hcp = pd.DataFrame(hcp)
	cl_nm = ["band_"+ str(i).zfill(1) for i in range(1,368)]
	hcp.columns = cl_nm + ["CHM", "aspect","slope","sns_zn","sns_az","sol_zn", "sol_az"]
	return(hcp, refl_md)




def kld_transform(brick, pt_kld):
    kld_groups = pd.read_csv(pt_kld, header=None)
    kld_groups = kld_groups.rename(columns={0: "_kld_grp"})
    all_data = np.zeros([brick.shape[0],1])
    for jj in np.unique(kld_groups):
        which_bands = kld_groups._kld_grp == jj
        #min
        new_col = np.apply_along_axis(min, 1, brick[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
        #mean
        new_col = np.apply_along_axis(np.mean, 1, brick[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
        #max
        new_col = np.apply_along_axis(max, 1, brick[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
#
    all_data = pd.DataFrame(all_data)
    all_data = all_data.drop([0], axis=1)
    return(all_data)


#main():
output_df, refl_md = extract_hsi_and_brdf_data(full_path, epsg, ras_dir, year,site)
#os.chdir('/blue/ewhite/s.marconi/NeonSpeciesClassification/')

coeffs_brdf = pd.read_csv("/blue/ewhite/s.marconi/NeonSpeciesClassification//data/corrections/angle_all_all_brdf_coeffs_1920.csv")
topo_coeff = pd.read_csv("/blue/ewhite/s.marconi/NeonSpeciesClassification//data/corrections/topo_all_all_brdf_coeffs_1920.csv")
ave_coords = pd.read_csv("/blue/ewhite/s.marconi/NeonSpeciesClassification//data/site_encoder.csv")

#transform into BRDF corrected tile
which_site = ave_coords.siteName == site
coords = ave_coords[which_site][['latitude', 'longitude','siteID']]
coords = pd.concat([coords]*output_df.shape[0])


brdf = calculate_brdf(output_df, coeffs_brdf, topo_coeff)
elevation = output_df[['CHM']]
#filter for greennes and shadows
ndvi = (brdf[:,89]- brdf[:,57])/(brdf[:,57] + brdf[:,89]) <0.5
nir860 = (brdf[:,95] + brdf[:,96])/2 <0.1
brdf[ndvi,:]= np.NAN
brdf[nir860,:]= np.NAN

# apply KLD reduction
brdf = brdf[:,10:357]

# map pixels with NAN
which_pixels_dropped = np.isnan(brdf).any(axis=1)
brdf = brdf[~np.isnan(brdf).any(axis=1)]
#remove bad pixels 
nbands = brdf.shape[1]
brdf[:,1:nbands] = normalize(brdf[:,1:nbands])
brdf = normalize(brdf)
brdf = kld_transform(brdf, pt_kld)
coords = coords[~which_pixels_dropped]
elevation = output_df.CHM[~which_pixels_dropped]

#brdf = pd.DataFrame(brdf[:,1:46])
brdf.reset_index(drop=True, inplace=True)
elevation.reset_index(drop=True, inplace=True)
coords.reset_index(drop=True, inplace=True)

brdf = pd.concat([elevation, brdf, coords], axis=1)

del(elevation, coords, ndvi, nir860, output_df)
# load model
import joblib
clf_bl = joblib.load(model_pkl)
#make prediction at pixels
taxa_classes = clf_bl.meta_clf_.classes_
#brdf = pd.DataFrame(brdf[:,1:46])
#brdf.reset_index(drop=True, inplace=True)
#elevation.reset_index(drop=True, inplace=True)
#coords.reset_index(drop=True, inplace=True)

#brdf = pd.concat([elevation, brdf, coords], axis=1)
predict_an = clf_bl.predict_proba(brdf)
predict_an = pd.DataFrame(predict_an)
predict_an.columns = taxa_classes

import os, psutil
process = psutil.Process(os.getpid())
print(process.memory_info().rss)  # in bytes 

outpred = np.empty([refl_md['shape'][0]*refl_md['shape'][1], taxa_classes.shape[0]])
outpred[:] = np.NaN
outpred[~which_pixels_dropped,:]=predict_an
outpred = np.reshape(outpred, [refl_md['shape'][0], refl_md['shape'][1], taxa_classes.shape[0]])
itc_id =  str(int(year)) + "_" +site+"_" + str(int(refl_md['extent'][0])) + "_" + str(int(refl_md['extent'][2])) 
ii = str(itc_id + ".tif")
clipExtent =  {}
clipExtent['xMin'] = refl_md['extent'][0]
clipExtent['yMin'] = refl_md['extent'][2]
clipExtent['yMax'] = refl_md['extent'][3]
clipExtent['xMax'] = refl_md['extent'][1]
array2raster(str(ii), outpred, refl_md, clipExtent, ras_dir = str(ras_dir)+'/predictions_raster/'+site+"/", epsg = int(refl_md['epsg']))

#plot transformed HSI
outpred = np.empty([refl_md['shape'][0]*refl_md['shape'][1], brdf.shape[1]])
outpred[:] = np.NaN
outpred[~which_pixels_dropped,:]=brdf
outpred = np.reshape(outpred, [refl_md['shape'][0], refl_md['shape'][1], brdf.shape[1]])
array2raster(str(ii), outpred, refl_md, clipExtent, ras_dir = str(ras_dir)+'/compiled_hsi/'+site+"/", epsg = int(refl_md['epsg']))












