#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 16:05:26 2021

@author: sergiomarconi
"""

import numpy as np
import pandas as pd
import glob

import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import BaggingClassifier
from mlxtend.classifier import StackingCVClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.experimental import enable_hist_gradient_boosting 
from sklearn.ensemble import HistGradientBoostingClassifier
from mlxtend.classifier import SoftmaxRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import normalize

import os

from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from imblearn.over_sampling import SMOTENC
from imblearn.over_sampling import ADASYN
from imblearn.under_sampling import TomekLinks 
from collections import Counter    
from hdr import *
from sklearn.preprocessing import normalize


from sklearn.calibration import CalibratedClassifierCV
from sklearn.svm import SVC
import category_encoders as ce

species_to_genus = {
        '2PLANT':"NA",
        'ABBA':"AB",
        'ABLAL':"AB",
        "ABLO":'AB',
        "ABMA":"AB",
        "ABAM":"AB",
        'ACNE2': "AC",
        'ACNEN': "AC",
        'ACPE':"AC",
        'ACRU': "AC",
        'ACSA3' : "AC",
        'ACSA2' :"AC",
        'ACFA':"AC2",
        'ACKO': "AC2",
        'ACGR':"AC2",
        'AIAL' : "AI",
        'ALRU2': "AL",
        'ALVI5':'AL',
        'AMLA' : "AM",
        'AMEL':'AM',
        'ARVIM':"AR",
        'BEAL2': "BE",
        'BEGL/BENA':"BE",
        'BEPA': "BE",
        'BELE': "BE",
        'BEPO': "BE",
        'BETUL':'BE',
        'BENE4' : "BE",
        'BUBU':"BU",
        'BUSI':"BU",
        'BOSU2':"BO",
        'CACA18':"CA1",
        'CADE27':"CA2",
        'CAGL8':"CA3",
        'CAOV2':"CA3",
        'CAOV3':"CA3",
        'CAAQ2':'CA',
        'CACO15': "CA3",
        'CATO6':"CA3",
        'CAIL2':"CA3",
        'CECA4':"CE1",
        'CELA':"CE2",
        'CEOC':"CE2",
        'CODR':"CO",
        'CODI8':"CO",
        'COFL2':"CO2",
        'DIVI5':"DI",
        'ELAN':"EL",
        'FAGR':"FA",
        'FRAM2':"FR",
        'FRAXI':'FR',
        'FRNI':'FR',
        'LARIX':'LA',
        'ILAN':'IL',
        'FRPE':"FR",
        'GYDI':"GY",
        'GUOF':"GU",
        'GUSA':"GU",
        'GLTR':"GL",
        'HALES':"HA",
        'JUNI':"JU1",
        'JUNIP':"JU2",
        'JUVI':"JU2",
        'JUOS':"JU2",
        'LIST2':"LI1",
        'LITU':"LI2",
        'MAPO':"MA",
        'MAFR':'MA',
        'MAGNO':'MA',
        'MORU2':"MO",
        'NYBI':"NY",
        'NYSY':"NY",
        'NYAQ2':'NY',
        'OXYDE':"OX",
        'OXAR':"OX",
        'OSVI':'OS',
        'PICEA':"PI1",
        'PIAL3':"PI2",
        'PIAC':"PI3",
        'PICO':"PI2",
        'PIEL':"PI2",
        'PIEN':"PI2",
        'PIEC2':"PI2",
        'PIFL2':"PI2",
        'PIGL':"PI2", 
        'PIMA':"PI2",
        'PINUS':'PI2',
        'PIPA2':"PI2",
        'PIPO':"PI2", 
        'PIRU':"PI2",
        'PIPOS':"PI2",
        'PIPU5':"PI2",
        'PIST':"PI2",
        'PITA':"PI2",
        'PIGL2':"PI2",
        'PIED':"PI",
        'PIJE':"PI",
        'PIRI':'PI',
        'PIVI2':'PI',
        'PINUS':"PI2",
        'PLOC':"PL",
        'POTR5':"PO",
        'POGR4':"PO",
        'PODE3':"PO",
        'PRVE':"PR",
        'PRVI':"PR",
        'PRAV':'PR',
        'PRSE2': "PR",
        'PRAN3':"PR",
        'PSME':"PS",
        'QUAL':"QU", 
        'QUCO2':"QU",
        'QUCH':"QU",
        'QUCH2':"QU",
        'QUHE2':'QU',
        'QUERC':"QU",
        'QUGE2':"QU", 
        'QUSH':"QU",
        'QULA2':'QU',
        "QUPH":"QU",
        'QULA3':"QU",
        'QUERCUS':"QU",
        'QULY':"QU", 
        'QUMA3':"QU",
        'QUMA13':"QU",
        'THUJA':"TU",
        'PISA2':"PI2",
        'TABR2':"TA",
        'QUDO':"QU",
        'MEPO5':'ME',
        'QUMI':"QU",
        'QUFA':"QU",
        'QUMO4':"QU",
        'QUMU':"QU",
        'QUNI':"QU",
        'QUKE':"QU",
        'QUVE':'QU',
        'QUWI2':"QU",
        'QUPA5':"QU",
        'QURU':"QU",
        'QUST':"QU",
        'RHGL':"RH",
        "ROPS":"RO",
        'SASSA':'SA',
        'SALIX':'SA',
        'SYOC':"SY",
        'SILA20':"SI",
        'SWMA2':"SW",
        'TRSE6':"TR",
        'TSCA':"TS",
        'TSHE':"TS",
        'TIAM':"TI",
        'TAHE':"TA",
        'ULAL':"UL",
        'ULAM':"UL",
        'ULMUS':"UL", 
        'ULCR':"UL",
        'ULRU':"UL",
        }

domainid = {
  "BART": "D01",
  "HARV": "D01",
  "STEI": "D05",
  "TREE": "D05",
  "UNDE": "D05",
  "SERC": "D02",
  "SCBI": "D02",
  "OSBS": "D03",
  "GRSM": "D07",
  "MLBS": "D07",
  "DELA": "D08",
  "TALL": "D08",
  "BLAN": "D02",
  "UKFS": "D06",
  "RMNP": "D10",
  "BONA": "D19",
  "MOAB": "D13",
  "DEJU": "D19",
  "LENO": "D08",
  "DSNY": "D03",
  "JERC": "D03",
  "GUAN": "D04",
  "HEAL": "D19",
  "KONZ": "D06",
  "CLBJ": "D11",
  "YELL": "D12",
  "NIWO": "D13",
  "ABBY": "D16",
  "WREF": "D16",
  "SJER": "D17",
  "PUUM": "D20"
}



# prepare input data
def prepare_inputs(X_train, X_test, cats = ['domainID', 'siteID']):
	X_train_enc, X_test_enc = list(), list()
	# label encode each column
	for i in  cats:
		le = LabelEncoder()
		le.fit(X_train[i])
		# encode
		train_enc = le.transform(X_train[i])
		test_enc = le.transform(X_test[i])
		# store
		X_train_enc.append(train_enc)
		X_test_enc.append(test_enc)
	return X_train_enc, X_test_enc


# dimensionality reduction
def kld_reduction(brick, kld_out):
    from sklearn import preprocessing
    refl = brick.drop(['individualID'], axis=1)
    scaler = preprocessing.StandardScaler().fit(refl)
    refl =  scaler.transform(refl)
    kld_groups = getClusters(refl, numBands = 15)
    np.savetxt(kld_out, kld_groups, delimiter=",")
    individualID=brick["individualID"]
#
    brick = brick.drop(columns=['individualID'])
    brick = brick.values
    all_data = np.zeros([brick.shape[0],1])
    for jj in np.unique(kld_groups):
        which_bands = kld_groups == jj
        #min
        new_col = np.apply_along_axis(min, 1, brick[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
        #mean
        new_col = np.apply_along_axis(np.mean, 1, brick[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
        #max
        new_col = np.apply_along_axis(max, 1, brick[:,which_bands])[...,None]
        all_data = np.append(all_data, new_col, 1)
    #reappend individualID on the all data dataframe
    all_data = pd.DataFrame(all_data)
    all_data["individualID"] = individualID
    all_data = all_data.drop([0], axis=1)
    #
    #shift back individualID to first row
    cols = list(all_data.columns)
    cols = [cols[-1]] + cols[:-1]
    all_data = all_data[cols]
    return all_data


feats_pt = glob.glob('./data/features*.csv')
brick = pd.read_csv(feats_pt[0])      #"./data/brdf_spectra_2104b.csv") 
metadata_pt = glob.glob('./data/metadata_*.csv')
metadata = pd.read_csv(metadata_pt[0])        #"./data/metadata_2104b.csv") 
metadata = metadata[["individualID", "groupID", "plotID","siteID","elevation","latitude", "longitude", "taxonID"]]  
kld_out="./data/0411.csv"#"./data/kld_grps_2104b.csv"
nbands = brick.shape[0]
brick.iloc[:,1:nbands] = normalize(brick.iloc[:,1:nbands])
brick = kld_reduction(brick, kld_out)
foo = brick.drop(columns=[ 'individualID'])

#metadata = metadata.join(elevation.set_index('individualID'), on='individualID')
metadata = metadata.dropna()

siteID =  "ALL" #"D01"
#foo = brick.drop(columns=['index', 'individualID'])
data = pd.concat([metadata, foo], axis=1)
#data = brick.set_index('individualID').join(metadata.set_index('individualID'))
data.reset_index(inplace=True)

is_bad_genus = ["MAGNO", "AMLA"]
is_bad_genus =  data['taxonID'].isin(is_bad_genus)
data = data[~is_bad_genus]

ave_coords= data[["latitude", "longitude", "siteID"]].groupby(['siteID'], as_index=False).mean()
data = data.drop(columns=["latitude", "longitude"]).set_index('siteID').join(ave_coords.set_index('siteID'))

data = data.dropna()
data = data.drop(columns=['index', 'plotID'])
species_id = data.taxonID.unique()

#splin into train and test by chosing columns
train = data.groupID == "train"
test = data.groupID == "test"

y_test =  data[['individualID','taxonID', 'groupID']][~train]
X_test = data.drop(columns=['individualID', 'taxonID', 'groupID'])[~train]
X_train = data.drop(columns=['individualID', 'taxonID', 'groupID'])[train]
y_train =  data[['taxonID']][train]
y = data[['taxonID']]



# oversample using SMOTENC in order not to loose the categorical effects
#get relatie frequency of each class
ratios_for_each = Counter(y_train.taxonID)
ratios_for_each = pd.DataFrame.from_dict(ratios_for_each, orient='index').reset_index()
#ratios_for_each.iloc[:,1] = ratios_for_each.iloc[:,1]

#remove mean of kld class (redundant feature)
cat_col = X_train.shape[1]
cols = np.arange(start = 2, stop = cat_col-2, step=3)
#X_train.drop(X_train.columns[cols],axis=1,inplace=True)
#X_test.drop(X_test.columns[cols],axis=1,inplace=True)
cat_col = X_train.shape[1]
X_train.columns

#undersample  
from imblearn.under_sampling import RandomUnderSampler
from imblearn.under_sampling import NeighbourhoodCleaningRule
from imblearn.under_sampling import TomekLinks 
from imblearn.combine import SMOTETomek 

min_class = Counter(y_train.taxonID)
min_class = min_class[min(min_class, key=min_class.get)]

print(min_class)
smotenc = SMOTENC(random_state=0, categorical_features = [0,cat_col-2,cat_col-1], k_neighbors=min_class-1)
smt = SMOTETomek(random_state=42, smote=smotenc)
X_res, y_res = smt.fit_resample(X_train, y_train)
#turn site in encoded effect
new_df = pd.merge(X_res, ave_coords,  how='left', left_on=['latitude', 'longitude'], right_on = ['latitude', 'longitude'])
df=new_df[~new_df['siteID'].isna()]
y_res = y_res[~new_df['siteID'].isna()]

def categorical_encoder(cats,y):
    import category_encoders as ce
    le = LabelEncoder()
    le.fit(y)
    le = le.transform(y)
    enc = ce.LeaveOneOutEncoder(cols=['siteID'])
    # enc = enc.fit(cats).transform(cats)
    train_enc = enc.fit_transform(cats,le)
    return(train_enc)



cat_encoder = categorical_encoder(df,y_res['taxonID'])
cat_encoder
labels_encoder = cat_encoder[['latitude', 'longitude','siteID']].drop_duplicates()
labels_encoder = labels_encoder.groupby(['latitude', 'longitude']).mean()




X_res = pd.merge(df.drop(columns=['siteID']), labels_encoder,  how='left', left_on=['latitude', 'longitude'], right_on = ['latitude', 'longitude']) 
X_test = pd.merge(X_test, labels_encoder,  how='left', left_on=['latitude', 'longitude'], right_on = ['latitude', 'longitude']) 
X_test.to_csv("./data/X_test_14april.csv")
y_test.to_csv("./data/y_test_14april.csv")
labels_encoder.to_csv("./data/site_14encoder.csv")


import numpy as np
import pandas as pd

import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import BaggingClassifier
from mlxtend.classifier import StackingClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.experimental import enable_hist_gradient_boosting 
from sklearn.ensemble import HistGradientBoostingClassifier
from mlxtend.classifier import SoftmaxRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier




from sklearn.calibration import CalibratedClassifierCV
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier

#define models
rf2 = make_pipeline(StandardScaler(),
                    RandomForestClassifier(random_state=0, oob_score = True,
        n_estimators = 300, max_features = 'sqrt', criterion = 'entropy'))


knn2 = make_pipeline(StandardScaler(),
                     KNeighborsClassifier(weights = 'distance', 
                                          p=1, n_neighbors=20))
gb2 = make_pipeline(StandardScaler(),
                    HistGradientBoostingClassifier(random_state=0, 
                    max_iter = 1000, learning_rate = 0.01, 
                max_depth = 25, loss = 'categorical_crossentropy', 
                l2_regularization = 0.5))
mlpc2 = make_pipeline(StandardScaler(), 
                      MLPClassifier(random_state=0, 
                                    beta_2=0.9, max_iter = 1200))

from sklearn.naive_bayes import GaussianNB
bayes2 = make_pipeline(StandardScaler(), GaussianNB())

bsvc2 =make_pipeline(StandardScaler(),
                     BaggingClassifier(
                         base_estimator=SVC(probability = True, C = 1000), 
                        n_jobs = 1, random_state=0))

logc = LogisticRegression(penalty = "elasticnet", solver = "saga", 
                          max_iter = 10000, n_jobs=3, l1_ratio = 0.5)


logc = LogisticRegression(penalty = "elasticnet", solver = "saga", 
                          max_iter = 10000, n_jobs=3, l1_ratio = 0.5)
from mlxtend.classifier import StackingCVClassifier
clf_bl2 = StackingCVClassifier(classifiers = [rf2, gb2, bsvc2, mlpc2, knn2],
            use_probas=True, cv = 3, n_jobs =1,
            meta_classifier= logc)

# grid.fit(X_res, y_res.taxonID.ravel())
clf_bl2.fit(X_res, y_res.taxonID.ravel())
print(clf_bl2.score(X_test, y_test['taxonID'].ravel()))

predict_an = clf_bl2.predict_proba(X_test)
predict_an = pd.DataFrame(predict_an)

from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
#from class_hierarchy import *
# format outputs for final evaluation
#knn = knn.fit(X_res, y_res.taxonID)
taxa_classes = clf_bl2.meta_clf_.classes_
colnames = np.append(['individualID', 'taxonID', 'groupID'], taxa_classes) #.tolist()
y_test.reset_index(drop=True, inplace=True)
predict_an.reset_index(drop=True, inplace=True)
eval_an = pd.concat([y_test, predict_an], axis=1)
eval_an.columns = colnames
#aggregate probabilities of each pixel to return a crown based evaluation. Using majority vote
eval_an = eval_an.groupby(['individualID', 'taxonID'], as_index=False).mean()

y_itc = eval_an['taxonID']
pi_itc = eval_an.drop(columns=['individualID', 'taxonID'])

# get the column name of max values in every row
pred_itc = pi_itc.idxmax(axis=1)
cm = confusion_matrix(y_itc, pred_itc, labels = taxa_classes)
cm = pd.DataFrame(cm, columns = taxa_classes, index = taxa_classes)
#mcm = multilabel_confusion_matrix(y_itc, pred_itc)
#f1_score(y_itc, pred_itc, average='macro')
#f1_score(y_itc, pred_itc, average='micro')
report = classification_report(y_itc, pred_itc, output_dict=True)
report = pd.DataFrame(report).transpose()
report

pd.set_option('display.max_rows', report.shape[0]+1)
report

result_check = pd.DataFrame(np.c_[eval_an[['individualID','taxonID']], pred_itc])
result_check = pd.DataFrame(np.c_[result_check.iloc[:,0].str.slice(start=13, stop=17),result_check])

from sklearn.metrics import f1_score
sites = result_check.iloc[:,0].unique()
result_check.columns =['siteID', 'individualID', 'obs', 'pred'] 
for st in sites:
    results_in_site = result_check[result_check.siteID.eq(st)]
    print(st + ": " + str(f1_score(results_in_site.obs, results_in_site.pred, average="micro")))
    print(st + ": " + str(f1_score(results_in_site.obs, results_in_site.pred, average="macro")))

obs = np.zeros(len(pred_itc)).astype(str)
pred = np.zeros(len(pred_itc)).astype(str)
for ii in range(len(pred_itc)): 
    obs[ii] = species_to_genus[y_itc[ii]]
    pred[ii] = species_to_genus[pred_itc[ii]]
    
rep_fam = classification_report(obs, pred, output_dict=True)
rep_fam = pd.DataFrame(rep_fam).transpose()
rep_fam 


#siteID = "BRDF_April"
test_families = pd.concat([pd.Series(obs), pd.Series(pred)], axis=1) 
test_families = test_families.rename(columns={0: "observed", 1: "predicted"})
test_families.to_csv("./outdir/"+siteID+"_final_"+"family_predictions.csv")   
cm_fam = confusion_matrix(obs, pred)
cm_fam = pd.DataFrame(cm_fam, columns = rep_fam.index[0:-3], index = rep_fam.index[0:-3])

eval_an.to_csv("./outdir/"+siteID+"_final_"+"kld_probabilities.csv")
pd.DataFrame(np.c_[eval_an[['individualID','taxonID']], pred_itc]).to_csv("./outdir/"+siteID+"_final_"+"kld_pairs.csv")
pd.DataFrame(taxa_classes).to_csv("./outdir/taxonID_dict"+"_final_"+siteID+".csv")


import joblib
mod_out_pt = "./outdir/model_"+siteID+'_final_model.pkl'
joblib.dump(clf_bl2, mod_out_pt)

