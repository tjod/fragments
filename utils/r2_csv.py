import pandas as pd #import the pandas module
import numpy as np
from sklearn.metrics import r2_score
import sys

file1 = 'roche_py.csv'
df1 = pd.read_csv (file1, sep=',')
print (df1.columns)
file2 = 'roche_R.csv'
df2 = pd.read_csv (file2, sep=',')
print (df2.columns)

a = df1['pred_logD_py']
b = df2['pred_logD_R']
print ( len(a), len(b) )
print ("%s:%s %.3f" % ('pred_logD_py', 'pred_logD_R', r2_score(a,b)))

# remove missing values
df1 = df1[(df1['pRED LogD_LOGD_AVG'] != 'None')]
df2 = df2[(df2['pRED LogD_LOGD_AVG'] != 'None')]

a = df1['pRED LogD_LOGD_AVG']
b = df1['pred_logD_py']
print ("%s:%s %.3f" % ('pRED LogD_LOGD_AVG', 'pred_logD_py', r2_score(a,b)))

a = df2['pRED LogD_LOGD_AVG']
b = df2['pred_logD_R']
print ("%s:%s %.3f" % ('pRED LogD_LOGD_AVG', 'pred_logD_R', r2_score(a,b)))

a = df1['pRED LogD_LOGD_AVG']
b = df2['pRED LogD_LOGD_AVG']
print ("%s:%s %.3f" % ('pRED LogD_LOGD_AVG', 'pRED LogD_LOGD_AVG', r2_score(a,b)))