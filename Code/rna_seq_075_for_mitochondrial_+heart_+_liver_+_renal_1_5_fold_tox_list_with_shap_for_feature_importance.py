# -*- coding: utf-8 -*-
"""RNA Seq 075 for Mitochondrial +Heart + Liver + Renal 1.5 fold tox list with SHAP for feature importance.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1c0insI_DGfE5qGpjUPQIBiOVb0MMdh4F
"""

# Commented out IPython magic to ensure Python compatibility.
!pip install shap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
import xgboost
import shap

from numpy import mean, std
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_score, recall_score, f1_score, make_scorer, roc_curve, auc
from sklearn.feature_selection import SelectKBest, chi2, RFE
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold,cross_val_predict
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import TomekLinks
from sklearn.pipeline import Pipeline, make_pipeline

# %matplotlib inline

! python --version

print("shap version:", shap.__version__)

from google.colab import drive
drive.mount('/content/drive')

df = pd.read_csv('/content/drive/MyDrive/Colab Notebooks/GSE152075_raw_counts_GEO _correct to CSV_for_COVID19 for Python ML.csv',dtype=object)
df.head()

# Transpose
data = df.values
index1 = list(df.keys())
data = list(map(list, zip(*data)))
data = pd.DataFrame(data, index=index1)
data.to_csv('TML.csv', header=0)

df = pd.read_csv('TML.csv',dtype=object)
df.head()

df = df.drop(df.columns[0], axis=1)
df.head()

df['target'] = df.target
df.target.value_counts()

# Mitochondrial +Heart + Liver + Renal tox list significantly expressed genes in GSE152075
df_removetarget = df.iloc[:, :-1]
df_removetarget = df_removetarget[['BAD',
'TLR4',
'CASP1',
'TLR2',
'CTNNB1',
'JUN',
'NDUFS6',
'NDUFA13',
'GPX1',
'GPX4',
'GSTP1',
'NDUFAB1',
'JUND',
'CCN1',
'SOCS3',
'USP18',
'CYBB',
'ACE2',
'PRKCB',
'SLC8A1',
'NTN1',
'FLT1',
'SLC2A1',
'TLR7',
'FOS',
'CD274',
'IER3',
'KMO',
'P2RX7',
'MIF',
'AIM2',
'TNFSF14',
'BIRC3',
'MT-ND6',
'NDUFV1',
'CYC1',
'NDUFB7',
'NDUFB10',
'COX5A',
'PRDX5',
'CLIC2',
'UQCR11',
'UQCRQ',
'NDUFA4',
'NDUFA11',
'ATP5F1E',
'COX7B',
'UQCRC1',
'FIS1',
'NDUFB2',
'SURF1',
'ACO2',
'COX5B',
'NDUFA2',
'TOMM7',
'COX4I1',
'COX6A1',
'ATP5ME',
'PIN1',
'NUB1',
'DVL1',
'PROX1',
'STEAP3',
'RRAD',
'CARD6',
'LMNA',
'CIB1',
'S100A6',
'MB',
'KLF15',
'LARP6',
'DAG1',
'KLF4',
'SERPINF1',
'NEXN',
'ABCC9',
'ALDH3A1',
'KRT8',
'CXCL10',
'ATG4B',
'TKT',
'FGL2',
'CCR5',
'EPHA2',
'CCL2',
'SIGIRR',
'GBP5',
'PHB2',
'PROS1',
'SELL',
'IRF8',
'BSG',
'GADD45B',
'KEAP1',
'CCL4',
'PTPRC',
'IL2RG',
'CHCHD2',
'TFF3',
'MEFV',
'APRT',
'TNFSF13B',
'STUB1',
'ALDH3B1',
'MLKL',
'RFXANK',
'C3AR1',
'IDO1',
'CLU',
'CX3CR1',
'CITED2',
'APOBEC3A',
'LRP5',
'OLR1',
'ZEB1',
'LAMB2',
'FCGR3A',
'BCL2L14',
'PRDX2',
'GNB2',
'FCGR3B',
'HSPA1B',
'CCR1',
'ERBB2',
'ZBP1',
'HSPA1A']]
print(df_removetarget)
df_removetarget.head()
#df_removetarget.to_csv('sig without target.csv')

X = df_removetarget
y = df.iloc[:, -1]
print(X)
print(y)

X = df_removetarget
y = df.iloc[:, -1]
X_number = X.values
y_index = y.values

# standardize the dataset
scaler = MinMaxScaler()
X_number = scaler.fit_transform(X_number)

print(X_number.shape)
print(y_index.shape)

# # for feature importance of XGBoost gene names
# X_number = pd.DataFrame(data=X_number, columns=['BAD',
# 'TLR4',
# 'CASP1',
# 'TLR2',
# 'CTNNB1',
# 'JUN',
# 'NDUFS6',
# 'NDUFA13',
# 'GPX1',
# 'GPX4',
# 'GSTP1',
# 'NDUFAB1',
# 'JUND',
# 'CCN1',
# 'SOCS3',
# 'USP18',
# 'CYBB',
# 'ACE2',
# 'PRKCB',
# 'SLC8A1',
# 'NTN1',
# 'FLT1',
# 'SLC2A1',
# 'TLR7',
# 'FOS',
# 'CD274',
# 'IER3',
# 'KMO',
# 'P2RX7',
# 'MIF',
# 'AIM2',
# 'TNFSF14',
# 'BIRC3',
# 'MT-ND6',
# 'NDUFV1',
# 'CYC1',
# 'NDUFB7',
# 'NDUFB10',
# 'COX5A',
# 'PRDX5',
# 'CLIC2',
# 'UQCR11',
# 'UQCRQ',
# 'NDUFA4',
# 'NDUFA11',
# 'ATP5F1E',
# 'COX7B',
# 'UQCRC1',
# 'FIS1',
# 'NDUFB2',
# 'SURF1',
# 'ACO2',
# 'COX5B',
# 'NDUFA2',
# 'TOMM7',
# 'COX4I1',
# 'COX6A1',
# 'ATP5ME',
# 'PIN1',
# 'NUB1',
# 'DVL1',
# 'PROX1',
# 'STEAP3',
# 'RRAD',
# 'CARD6',
# 'LMNA',
# 'CIB1',
# 'S100A6',
# 'MB',
# 'KLF15',
# 'LARP6',
# 'DAG1',
# 'KLF4',
# 'SERPINF1',
# 'NEXN',
# 'ABCC9',
# 'ALDH3A1',
# 'KRT8',
# 'CXCL10',
# 'ATG4B',
# 'TKT',
# 'FGL2',
# 'CCR5',
# 'EPHA2',
# 'CCL2',
# 'SIGIRR',
# 'GBP5',
# 'PHB2',
# 'PROS1',
# 'SELL',
# 'IRF8',
# 'BSG',
# 'GADD45B',
# 'KEAP1',
# 'CCL4',
# 'PTPRC',
# 'IL2RG',
# 'CHCHD2',
# 'TFF3',
# 'MEFV',
# 'APRT',
# 'TNFSF13B',
# 'STUB1',
# 'ALDH3B1',
# 'MLKL',
# 'RFXANK',
# 'C3AR1',
# 'IDO1',
# 'CLU',
# 'CX3CR1',
# 'CITED2',
# 'APOBEC3A',
# 'LRP5',
# 'OLR1',
# 'ZEB1',
# 'LAMB2',
# 'FCGR3A',
# 'BCL2L14',
# 'PRDX2',
# 'GNB2',
# 'FCGR3B',
# 'HSPA1B',
# 'CCR1',
# 'ERBB2',
# 'ZBP1',
# 'HSPA1A'])

# #XGBoost feature importance
# y_index = y_index.astype(int)
# classifier = xgboost.XGBClassifier(random_state=42).fit(X_number, y_index)
# xgboost.plot_importance(classifier, max_num_features=20)
# plt.title("xgboost.plot_importance(XGBClassifier)")
# plt.show()

# booster = classifier.get_booster()
# importance = booster.get_score(importance_type='weight')
# sorted_importance = sorted(importance.items(), key=lambda x: x[1], reverse=True)
# for feat, imp in sorted_importance:
#     print(f"{feat}: {imp}")

explainer = shap.Explainer(classifier)

shap_output = explainer(X_number)
shap_values = shap_output.values

if len(shap_values.shape) == 3:
    shap_values = shap_values[:, :, 0]

mean_abs_shap_values = abs(shap_values).mean(axis=0)

sorted_indices = mean_abs_shap_values.argsort()[::-1]
sorted_features = X_train.columns[sorted_indices]
sorted_importances = mean_abs_shap_values[sorted_indices]

for feat, imp in zip(sorted_features, sorted_importances):
    print(f"{feat}: {imp}")

df_shap_importance = pd.DataFrame({
    'Feature': sorted_features,
    'Importance': sorted_importances
})
df_shap_importance.to_excel("shap_importance.xlsx", index=False)

# 計算 "Importance" 列的平均值
importance_mean = df_shap_importance["Importance"].mean()

# 計算 "Importance" 列的中位數
importance_median = df_shap_importance["Importance"].median()

print("Mean Importance:", importance_mean)
print("Median Importance:", importance_median)

max_features = 50
shap.plots.bar(shap_output, max_display=max_features, show=False)
plt.ticklabel_format(style='plain', axis='x')  # Disable the scientific notation on x-axis
plt.gca().xaxis.set_major_formatter('{:.3f}'.format)  # Set the x-axis labels to 3 decimal places

plt.savefig("shap_importance_bar_plot_decimal.png", bbox_inches='tight')

plt.show()

import matplotlib.pyplot as plt
import numpy as np

mean_shap = np.mean(np.abs(shap_values), axis=0)

sorted_indices = np.argsort(mean_shap)[-50:]

remaining_indices = np.setdiff1d(np.arange(len(mean_shap)), sorted_indices)
sum_other_features = np.sum(mean_shap[remaining_indices])

num_other_features = len(mean_shap) - 50
sorted_features = ["Sum of " + str(num_other_features) + " other features"] + X_train.columns[sorted_indices].tolist()
sorted_values = [sum_other_features] + mean_shap[sorted_indices].tolist()

fig, ax = plt.subplots(figsize=(10, 8))
bars = ax.barh(sorted_features, sorted_values)

for bar in bars:
    width = bar.get_width()
    ax.text(width + 0.02, bar.get_y() + bar.get_height()/2,
            '{:.3f}'.format(width),
            va='center', ha='left', color='black')

ax.set_xlabel('Mean Absolute SHAP Value')
ax.xaxis.set_major_formatter('{:.3f}'.format)  # Set the x-axis labels to 3 decimal places

# Extend the x-axis limit to make room for the text
ax.set_xlim(0, max(sorted_values) + max(sorted_values) * 0.15)  # Extend by 15% of the maximum SHAP value

plt.tight_layout()
plt.savefig("shap_importance_bar_plot_mean_with_values_and_others_v2.png")

plt.show()