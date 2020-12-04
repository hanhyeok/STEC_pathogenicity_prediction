import os 
import pandas as pd
import numpy as np
import re

import sklearn
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.svm import SVC
#from sklearn.externals import joblib
import joblib

import pickle
import datetime

## gene classification
files = os.listdir('/data/classification')
tag_list = [f[:f.find('_assigned_locus')] for f in files]

with open('/STEC_prediction/Scoary_significant_gene_list.pickle', 'rb') as f:
    scoary = pickle.load(f)

scoary_pd = pd.DataFrame(columns = tag_list , index=scoary)

for file in files:
    tag = re.sub('_assigned_locus.tsv', '', file)
    read = pd.read_csv('/data/classification/{}'.format(file), sep='\t')    
    dicted_data = dict(zip(read['0'], read[tag]))
    scoary_pd[tag] = scoary_pd.apply(lambda x: dicted_data[x.name], axis=1 )

scoary_input = scoary_pd.T
scoary_input.columns = [col[col.find('__')+2:] for col in scoary_input.columns]


## Stx subtyping
col = ['stx1a_alpha', 'stx1a_beta', 'stx1c_alpha', 'stx1c_beta', 'stx1d_alpha', 'stx1d_beta', 'stx2a_alpha', 'stx2a_beta', 'stx2b_alpha',  'stx2b_beta', 
'stx2c_alpha', 'stx2c_beta', 'stx2d_alpha', 'stx2d_beta', 'stx2e_alpha', 'stx2e_beta', 'stx2f_alpha', 'stx2f_beta', 'stx2g_alpha', 'stx2g_beta'] 

col2 = ['stx1a', 'stx1c', 'stx1d', 'stx2a', 'stx2b', 'stx2c', 'stx2d', 'stx2e', 'stx2f', 'stx2g']

idx_list = [re.sub('_stx.blastp', '',i) for i in os.listdir('/data/blastp/blastp_stx/')]

def custom(col_value):
    beta = col_value[[1,3,5,7,9,11,13,15,17,19]]
    alpha = col_value[[0,2,4,6,8,10,12,14,16,18]]
    test_beta = [0 if val == 0 else 1 for val in beta]
    test_alpha = [0 if val == 0 else 1 for val in alpha]
    stx_list = []
    if sum(test_beta) != 0:
        beta_index = [i for i, e in enumerate(test_beta) if e == 1]
        stx_list1 = []
        stx_list1 = np.array(col2)[beta_index]
        stx_dict1 = dict(zip(stx_list1, beta[beta_index]))
        stx_dict1_locus = [int(re.findall('\d+', lc[lc.rfind('_')+1:])[0]) for lc in stx_dict1.values()]
        added_stx = []
        
        if sum(test_alpha) != 0:
            alpha_index = [i for i, e in enumerate(test_alpha) if e == 1]
            stx_list2 = []
            stx_list2 = np.array(col2)[alpha_index]
            stx_dict2 = dict(zip(stx_list2, alpha[alpha_index]))
            stx_dict2_locus = [int(re.findall('\d+', lc[lc.rfind('_')+1:])[0]) for lc in stx_dict2.values()]
            #for i in stx_dict1_locus:
            df = pd.DataFrame(columns = stx_dict2.keys())
            order = 0
            for i in stx_dict1_locus:
                order += 1
                df.loc[i] = abs(stx_dict2_locus-np.array(i)) > 10
            for clm in df.columns:
                if all(df[clm]):
                    added_stx.append(clm)
                    
        stx_list = list(stx_dict1.keys()) + added_stx
                    
    else:
        alpha_index = [i for i, e in enumerate(test_alpha) if e == 1]
        stx_list2 = []
        stx_list2 = np.array(col2)[alpha_index]
        
        stx_list = stx_list2
            
    return stx_list

stx_classification = pd.DataFrame(columns=col, index=idx_list)

for idx in stx_classification.index:
    try:
        data = pd.read_csv('/data/blastp/blastp_stx/{}_stx.blastp'.format(idx), sep = '\t', header=None)
        dicted = {}
        dicted = dict(zip([data.iloc[data[data[0]==i][2].index[0], 1] for i in set(data[0])], set(data[0])))
        stx_classification.loc[idx] = [dicted.get(c) for c in col]
    except:
        print('"{}" may not have Shiga toxin genes. Please check "Prediction_locus_tag.csv" file at result folder.'.format(idx))

stx_classification = stx_classification.fillna(0)

stx_df = pd.DataFrame(columns = col2)

for idx in stx_classification.index:
    stx_input = stx_classification.loc[idx]
    stx_exist = custom(stx_input)
    stx_df.loc[idx] = [1 if stx in stx_exist else 0 for stx in col2]


## Result data
scoary_summary = pd.concat([scoary_input, stx_classification.fillna(0)], join = 'inner', axis =1 )
ML_input = pd.concat([scoary_input.apply(pd.to_numeric, errors='coerce').fillna(1), stx_df], axis=1, join='inner')


## SVM model prediction
svm = joblib.load('/STEC_prediction/trained_SVM_model.joblib')

prediction_result = svm.predict(ML_input)
decision_function_value = svm.decision_function(ML_input)
result_df= pd.DataFrame(columns = ['Classification', 'Decision function value'],index = ML_input.index)
result_df['Classification'] = ['Clinical isolate' if i == 1 else 'Environmental isolate' for i in svm.predict(ML_input)]
result_df['Decision function value'] = np.round(svm.decision_function(ML_input), 3)

import sys
from tabulate import tabulate
from colorama import init
init(strip=not sys.stdout.isatty())

from termcolor import cprint
from pyfiglet import figlet_format
from datetime import date

cprint(figlet_format(' STEC              \n', font='standard')+' pathogenic potential prediction model', 'white', 'on_blue', attrs=['bold'])
#cprint('   Pathogenic Potential Prediction Model   \n\n', 'white', 'on_blue', attrs=['bold'])
pdtabulate = lambda df:tabulate(df, headers='keys', tablefmt='psql')

#cprint('   Pathogenic Potential Prediction Model   \n\n','white', 'on_blue', attrs=['bold'] )
print('\n'+'ver0.9       '+ datetime.datetime.now().strftime('%Y-%m-%d, %H:%M and %Ss')+ '\n')
print(pdtabulate(result_df))


## Result file 

if not os.path.exists('/data/result'):
	os.makedirs('/data/result')

scoary_summary.to_csv('/data/result/Prediction_locus_tag_{}.csv'.format(datetime.datetime.now().strftime('%Y%m%d%H%M%S')))
ML_input.to_csv('/data/result/Prediction_input_data_{}.csv'.format(datetime.datetime.now().strftime('%Y%m%d%H%M%S')))
result_df.to_csv('/data/result/Prediction_result_{}.csv'.format(datetime.datetime.now().strftime('%Y%m%d%H%M%S')))

