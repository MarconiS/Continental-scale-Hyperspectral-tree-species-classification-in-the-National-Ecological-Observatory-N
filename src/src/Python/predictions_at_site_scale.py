os.chdir("/blue/ewhite/s.marconi/NeonSpeciesClassification/src/")
import os
site = "JERC"
year = "2019"
location = '/orange/ewhite/NeonData/'
md_pt =  "//blue/ewhite/s.marconi/NeonSpeciesClassification/outputs/model_ALL_final_model.pkl"
kld_pt = "//blue/ewhite/s.marconi/NeonSpeciesClassification/data/0411.csv"

# r=>root, d=>directories, f=>files
files_in_dir = []
for r, d, f in os.walk(location):
    for item in f:
        if site in item:
            if year in r:
                if '.h5' in item:
                    files_in_dir.append(os.path.join(r, item))


from random import sample
sub_sample = sample(files_in_dir,5)

#check the sampled tiles
import matplotlib.pyplot as plt
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize)
import rasterio
import h5py

for tile in sub_sample:
    hdf5_file = h5py.File(tile, 'r')
    reflArray = hdf5_file[site]['Reflectance']
    refl =reflArray['Reflectance_Data'][()]
    
    #pyplot.imshow(refl[:,:,35], cmap='viridis')
    norm = ImageNormalize(refl[:,:,[75]], interval=MinMaxInterval(),
                          stretch=SqrtStretch())
    # Display the image
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(refl[:,:,[75]], origin='lower', norm=norm)
    
md_pt
ras_dir = "/orange/idtrees-collab/species_classification/"

import os
import pandas as pd
#if not os.path.exists(str(ras_dir)+'/compiled_hsi/'+site+"/"):
#    os.makedirs(str(ras_dir)+'/compiled_hsi/'+site+"/")
    
if not os.path.exists(str(ras_dir)+'/predictions_raster/'+site+"/"):
    os.makedirs(str(ras_dir)+'/predictions_raster/'+site+"/")
    
    
! rm '/blue/ewhite/s.marconi/NeonSpeciesClassification/src/'{site}'launch_jobs.txt'
with open('/blue/ewhite/s.marconi/NeonSpeciesClassification/src/'+site+'launch_jobs.txt', 'w') as f:
    for tile in sub_sample:
        slurm_file = "/blue/ewhite/s.marconi/NeonSpeciesClassification/src/scale_classification.SLURM" 
        #{md_pt} {kld_pt} {year} {site} {tile} 
        f.write(slurm_file + " " + md_pt + " " + kld_pt + " " + year + " " + site + " " + tile +'\n' )

site+'launch_jobs.txt'
hsi = refl.reshape(-1,refl.shape[2])

import os
site = "OSBS"
year = "2018"
location = pt_preds+"/"+site+"/"
files_in_dir = []

for r, d, f in os.walk(location):
    for item in f:
        if site in item:
            if year in item:
                files_in_dir.append(item)
                
                
pattern = files_in_dir[0]
pattern[10:len(pattern)-4]

def get_box_center(x):
    lower_bnd = max(x.shape[0]-4, 0)
    upper_bnd = min(x.shape[0]+4, x.shape[0])
    return lower_bnd, upper_bnd
    

for pattern in files_in_dir:
    import fnmatch
    import os
    import rasterio
    import pandas as pd
    import geopandas as gpd
    import numpy as np
    def find_files(pt_boxes, pattern):
        return [n for n in fnmatch.filter(os.listdir(pt_boxes), pattern) if
            os.path.isfile(os.path.join(pt_boxes, n))]
    wcard = '*'+pattern[10:len(pattern)-4]+"_image.shp"
    bboxes_shp = find_files(pt_boxes, wcard)
    for fl in bboxes_shp:
        if year+"_"+site in fl:
            pt = fl
    #
    geodf = gpd.read_file(pt_boxes+pt)
    geodf['ID'] = geodf.index
    from geocube.api.core import make_geocube
    cube = make_geocube(
        geodf,
        resolution=(1, -1),
    )
    cube = cube.to_dataframe()
    src = rasterio.open(pt_preds+"/"+site+"/"+pattern)
    predictions = src.read()
    src.close()
    predictions = np.swapaxes(predictions,0,2)
    predictions = np.swapaxes(predictions,0,1)
    br_shape = predictions.shape
    predictions = np.reshape(predictions, (br_shape[0]*br_shape[1], br_shape[2]))
    predictions = pd.DataFrame(predictions)
    predictions = pd.concat([cube.reset_index(drop=True), predictions.reset_index(drop=True)], axis=1)
    predictions = predictions.groupby('ID').apply(lambda p: p.iloc[get_box_center(p)[0]:get_box_center(p)[1],:])
    predictions = predictions.reset_index(drop=True)
    predictions = predictions.groupby(['ID'], as_index=False).mean()
    eval_prob = predictions.drop(columns=["ID","left","bottom","right","top","height","area","spatial_ref"])
    label = eval_prob.idxmax(axis=1)
    score = eval_prob.max(axis=1)
    gdf = geodf.assign(score=score)
    gdf = gdf.assign(label=label)
    gdf.to_file(pt_results+pt)



