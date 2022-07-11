#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 16:25:56 2021

@author: sergiomarconi
"""

from yellowbrick.model_selection import FeatureImportances
from mlxtend.evaluate import feature_importance_permutation
import joblib
import pandas as pd
import numpy as np  
X_test = pd.read_csv("/blue/ewhite/s.marconi/NeonSpeciesClassification//data/X_test_14april.csv")
y_test = pd.read_csv("/blue/ewhite/s.marconi/NeonSpeciesClassification//data/y_test_14april.csv")
X_test = X_test.drop(columns = ['Unnamed: 0'])
model_pkl = "/blue/ewhite/s.marconi/NeonSpeciesClassification//outputs/model_ALL_final_model.pkl"
clf_bl3 = joblib.load(model_pkl)
indx_tst = y_test.groupby("individualID").sample(1).index
tst = X_test.iloc[indx_tst,:]
yst = y_test.iloc[indx_tst,:]

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
corr = spearmanr(tst).correlation

# Ensure the correlation matrix is symmetric
corr = (corr + corr.T) / 2
np.fill_diagonal(corr, 1)

# We convert the correlation matrix to a distance matrix before performing
# hierarchical clustering using Ward's linkage.
distance_matrix = 1 - np.abs(corr)
dist_linkage = hierarchy.ward(squareform(distance_matrix))
dendro = hierarchy.dendrogram(
    dist_linkage, labels=tst.columns.tolist(), ax=ax1, leaf_rotation=90
)
dendro_idx = np.arange(0, len(dendro["ivl"]))

ax2.imshow(corr[dendro["leaves"], :][:, dendro["leaves"]])
ax2.set_xticks(dendro_idx)
ax2.set_yticks(dendro_idx)
ax2.set_xticklabels(dendro["ivl"], rotation="vertical")
ax2.set_yticklabels(dendro["ivl"])
fig.tight_layout()
plt.show()


imp_vals, all_val = feature_importance_permutation(
    predict_method=clf_bl3.predict, 
    X=np.array(tst),
    y=np.array(yst.taxonID.ravel()),
    metric='accuracy',
    num_rounds=5,
    seed=1)

imp_vals.shape
all_val.shape
pd.Series(imp_vals).to_csv('./mods'+"_"+'_model2_importance_var.csv', index = False)
pd.DataFrame(all_val).to_csv('./mods'+"_"+'_model2_all_var.csv', index = False)


