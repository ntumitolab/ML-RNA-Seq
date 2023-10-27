# -*- coding: utf-8 -*-
"""RNA Seq 075 for Mitochondrial +Heart + Liver + Renal 1.5 fold tox list with SHAP, CV & Train_test_split.ipynb

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

X_train, X_test, y_train, y_test = train_test_split(X_number, y_index, test_size=0.2, random_state=42)

# for feature importance gene names
X_train = pd.DataFrame(data=X_train, columns=['BAD',
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
'HSPA1A'])

X_test = pd.DataFrame(data=X_test, columns=['BAD',
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
'HSPA1A'])

X_number = pd.DataFrame(data=X_number, columns=['BAD',
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
'HSPA1A'])

print(X_train.shape)
print(y_train.shape)
print(X_test.shape)
print(y_test.shape)

# to avoid imbalance data - OverSampling SMOTE
smote = SMOTE(random_state=42)
X_train, y_train = smote.fit_resample(X_train, y_train)
X_test, y_test = smote.fit_resample(X_test, y_test)

X_number, y_index = smote.fit_resample(X_number, y_index)

# to avoid imbalance data - OverSampling Tomek Link
#X_number, y_index = TomekLinks().fit_resample(X_number, y_index)

# convert to int64
y_train = y_train.astype(int)
y_test = y_test.astype(int)

print(X_number)
print(y_index)
print(X_number.shape)
print(y_index.shape)
print(X_train.shape)
print(y_train.shape)
print(X_test.shape)
print(y_test.shape)

plt.figure( figsize=(10,5) )
y_train_series = pd.Series(y_index)
y_train_series.value_counts().plot( kind='pie', colors=['lightcoral','skyblue'], autopct='%1.2f%%' )
plt.title( 'Pos/Neg' )
plt.ylabel( '' )
plt.show()

# Cross Validation of mschine learning for 4 models results
models = {
    "XGBoost": XGBClassifier(colsample_bytree=0.9, learning_rate=0.1, max_depth=10, n_estimators=50, random_state=17),
    "RandomForestClassifier": RandomForestClassifier(max_depth=10, min_samples_split=5, n_estimators=100, random_state=1),
    "LogisticRegression": LogisticRegression(C=50, max_iter=5000),
    "SVC": SVC(kernel='rbf', C=100, gamma=0.01, probability=True)
}

y_index = y_index.astype(int)  # Ensure y_index is of type int

def confusion_metrics_and_report(conf_matrix):
    # Extract metrics from the confusion matrix
    TP = conf_matrix[1][1]
    TN = conf_matrix[0][0]
    FP = conf_matrix[0][1]
    FN = conf_matrix[1][0]

    print('True Positives:', TP)
    print('True Negatives:', TN)
    print('False Positives:', FP)
    print('False Negatives:', FN)

    # Calculate metrics
    accuracy = (TP + TN) / (TP + TN + FP + FN)
    mis_classification = 1 - accuracy
    sensitivity = TP / (TP + FN)
    specificity = TN / (TN + FP)
    precision = TP / (TP + FP)
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity)

    conf_accuracy = (float(TP+TN) / float(TP + TN + FP + FN))
    conf_misclassification = 1 - conf_accuracy
    conf_sensitivity = (TP / float(TP + FN))
    conf_specificity = (TN / float(TN + FP))
    conf_precision = (TP / float(TP + FP))
    conf_f1 = 2 * ((conf_precision * conf_sensitivity) / (conf_precision + conf_sensitivity))

   # Print metrics
    print('-'*50)
    print(f'Accuracy: {round(conf_accuracy,3)}')
    print(f'Mis-Classification: {round(conf_misclassification,3)}')
    print(f'Sensitivity: {round(conf_sensitivity,3)}')
    print(f'Specificity: {round(conf_specificity,3)}')
    print(f'Precision: {round(conf_precision,3)}')
    print(f'f_1 Score: {round(conf_f1,3)}')

    # plot confusion_metrics
    cm_df = pd.DataFrame(conf_matrix, columns=['Predicted Negative', 'Predicted Positive'], index=['Actual Negative', 'Actual Positive'])
    sns.heatmap(cm_df, annot=True, cmap='Blues', fmt='g')
    plt.title('Confusion Matrix')
    plt.show()

def evaluate_model(name, model, X, y):
    print(f"Evaluating model: {name}")
    print('-'*50)

    # Step 1: Perform cross-validation and get accuracy
    cv_scores = cross_val_score(model, X, y, cv=10, scoring='accuracy')
    print("Cross-validated Accuracy: %.3f +/- %.3f" % (cv_scores.mean(), cv_scores.std()))

    # Step 2: Get predictions from cross-validation
    cv_preds = cross_val_predict(model, X, y, cv=10)

    # Step 3: Compute the confusion matrix and related metrics
    conf_matrix = confusion_matrix(y, cv_preds)
    confusion_metrics_and_report(conf_matrix)
    print("\n\n")

for name, model in models.items():
    evaluate_model(name, model, X_number, y_index)

# Cross Validation of mschine learning for 4 models ROC curves

X_numbera = np.array(X_number)
y_indexa = np.array(y_index)

# Setting up plot for ROC Curves
plt.figure(figsize=(12, 8))
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', alpha=0.8)

cv = StratifiedKFold(n_splits=10)

for name, model in models.items():
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    for train, test in cv.split(X_numbera, y_indexa):
        model.fit(X_numbera[train], y_indexa[train])
        probs = model.predict_proba(X_numbera[test])
        predsroc = probs[:, 1]
        fpr, tpr, _ = roc_curve(y_indexa[test], predsroc)
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = np.mean(aucs)
    plt.plot(mean_fpr, mean_tpr, label='%s (AUC = %0.3f)' % (name, mean_auc), lw=2, alpha=0.8)

plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curves of Models')
plt.legend(loc="lower right")
plt.show()

# Train_test_split of mschine learning for 4 models results

# Convert to int for compatibility with the models
y_train = y_train.astype(int)
y_test = y_test.astype(int)



def evaluate_model(name, model):
    print(f"Evaluating model: {name}")
    print('-'*50)

    # Fit the model
    model.fit(X_train, y_train)

    # Predictions
    y_preds = model.predict(X_test)

    # Classification report
    print(classification_report(y_test, y_preds, digits=3))

    # Confusion matrix
    cm = confusion_matrix(y_test, y_preds)
    TP = cm[1][1]
    TN = cm[0][0]
    FP = cm[0][1]
    FN = cm[1][0]
    print('True Positives:', TP)
    print('True Negatives:', TN)
    print('False Positives:', FP)
    print('False Negatives:', FN)

    # Metrics
    accuracy = (TP + TN) / (TP + TN + FP + FN)
    mis_classification = 1 - accuracy
    sensitivity = TP / (TP + FN)
    specificity = TN / (TN + FP)
    precision = TP / (TP + FP)
    f1_score_val = 2 * (precision * sensitivity) / (precision + sensitivity)

    print('-'*50)
    print(f'Accuracy: {round(accuracy,3)}')
    print(f'Mis-Classification: {round(mis_classification,3)}')
    print(f'Sensitivity: {round(sensitivity,3)}')
    print(f'Specificity: {round(specificity,3)}')
    print(f'Precision: {round(precision,3)}')
    print(f'f_1 Score: {round(f1_score_val,3)}')

    # Train and Test Accuracy
    print('Training set accuracy:', model.score(X_train, y_train))
    print('Testing set accuracy:', model.score(X_test, y_test))

    # Plot confusion matrix
    cm_df = pd.DataFrame(cm, columns=['Predicted Negative', 'Predicted Positive'], index=['Actual Negative', 'Actual Positive'])
    sns.heatmap(cm_df, annot=True, cmap='Blues', fmt='g')
    plt.title('Confusion Matrix')
    plt.show()
    print('='*100)

for name, model in models.items():
    evaluate_model(name, model)

# Train_test_split of mschine learning for 4 models ROC curves
plt.figure(figsize=(10, 8))

for name, model in models.items():
    model.fit(X_train, y_train)
    probas_ = model.predict_proba(X_test)
    fpr, tpr, thresholds = roc_curve(y_test, probas_[:, 1])
    auc_score = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=f"{name} (AUC = {auc_score:.3f})")

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.show()

#XGBoost feature importance
y_index = y_index.astype(int)
classifier = xgboost.XGBClassifier(random_state=42).fit(X_number, y_index)
xgboost.plot_importance(classifier, max_num_features=20)
plt.title("xgboost.plot_importance(XGBClassifier)")
plt.show()

booster = classifier.get_booster()
importance = booster.get_score(importance_type='weight')
sorted_importance = sorted(importance.items(), key=lambda x: x[1], reverse=True)
for feat, imp in sorted_importance:
    print(f"{feat}: {imp}")

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