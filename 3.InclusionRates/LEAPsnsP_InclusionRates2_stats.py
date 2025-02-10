#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 17:07:01 2024

@author: riannehaartsen
"""
# This script reads in the *_IncllRates.csv files from each of the pipelines. 
# It then calculates the inclusion rates for connectivity and power (thresholds
# of 90 and 20 epochs respectively), for values calculated across all trials and
# for condition differences. For each measure, Cochran's Q-tests are applied among
# the high-density pipelines, and low-density pipelines separately, to test for 
# differences in inclusion rates between pipelines for different measures of interest. 

# This scripts was created by Rianne Haartsen, PhD.
# Birkbeck College, University of London

# This script is released under the GNU General Public License version 3.  



# Import libraries
import pandas as pd
import numpy as np

# Create dataframes in prep for stats: across all epochs ######################
# Functional connectivity; threshold = 90
# read in data
Manual_Ntrls = pd.read_csv('Manual_InclRates.csv', sep = ',', header=0) 
MADE_Ntrls = pd.read_csv('MADE_InclRates.csv', sep = ',', header=0) 
MADEBOND_Ntrls = pd.read_csv('MADEBOND_InclRates.csv', sep = ',', header=0) 
HAPPEv1_Ntrls = pd.read_csv('HAPPEv1_InclRates.csv', sep = ',', header=0) 
HAPPEv4_Ntrls = pd.read_csv('HAPPEv4_InclRates.csv', sep = ',', header=0) 
MADEBONDld_Ntrls = pd.read_csv('MADE-BONDld_InclRates.csv', sep = ',', header=0) 
HAPPILEE_Ntrls = pd.read_csv('HAPPILEE_InclRates.csv', sep = ',', header=0) 
miniMADE_Ntrls = pd.read_csv('miniMADE_InclRates.csv', sep = ',', header=0) 

# convert tot to incl vs excl if tot number of trials/epochs threshold
Threshold = 90
conditions = {
    0: (Manual_Ntrls['Ntrls_tot'] < Threshold),
    1: (Manual_Ntrls['Ntrls_tot'] >= Threshold),
}
Manual_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=Manual_Ntrls['Ntrls_tot'])
Manual_incl = Manual_Ntrls[['ID','Ntrls_tot']]
Manual_incl.rename(columns={'Ntrls_tot': 'P1_trls'}, inplace=True)

conditions = {
    0: (MADE_Ntrls['Ntrls_tot'] < Threshold),
    1: (MADE_Ntrls['Ntrls_tot'] >= Threshold),
}
MADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADE_Ntrls['Ntrls_tot'])
MADE_incl = MADE_Ntrls[['ID','Ntrls_tot']]
MADE_incl.rename(columns={'Ntrls_tot': 'P2_trls'}, inplace=True)

conditions = {
    0: (MADEBOND_Ntrls['Ntrls_tot'] < Threshold),
    1: (MADEBOND_Ntrls['Ntrls_tot'] >= Threshold),
}
MADEBOND_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBOND_Ntrls['Ntrls_tot'])
MADEBOND_incl = MADEBOND_Ntrls[['ID','Ntrls_tot']]
MADEBOND_incl.rename(columns={'Ntrls_tot': 'P3_trls'}, inplace=True)

conditions = {
    0: (HAPPEv1_Ntrls['Neps_tot'] < Threshold),
    1: (HAPPEv1_Ntrls['Neps_tot'] >= Threshold),
}
HAPPEv1_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv1_Ntrls['Neps_tot'])
HAPPEv1_incl = HAPPEv1_Ntrls[['ID','Neps_tot']]
HAPPEv1_incl.rename(columns={'Neps_tot': 'P4_trls'}, inplace=True)

conditions = {
    0: (HAPPEv4_Ntrls['Neps_tot'] < Threshold),
    1: (HAPPEv4_Ntrls['Neps_tot'] >= Threshold),
}
HAPPEv4_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv4_Ntrls['Neps_tot'])
HAPPEv4_incl = HAPPEv4_Ntrls[['ID','Neps_tot']]
HAPPEv4_incl.rename(columns={'Neps_tot': 'P5_trls'}, inplace=True)

conditions = {
    0: (MADEBONDld_Ntrls['Ntrls_tot'] < Threshold),
    1: (MADEBONDld_Ntrls['Ntrls_tot'] >= Threshold),
}
MADEBONDld_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBONDld_Ntrls['Ntrls_tot'])
MADEBONDld_incl = MADEBONDld_Ntrls[['ID','Ntrls_tot']]
MADEBONDld_incl.rename(columns={'Ntrls_tot': 'P6_trls'}, inplace=True)

conditions = {
    0: (HAPPILEE_Ntrls['Neps_tot'] < Threshold),
    1: (HAPPILEE_Ntrls['Neps_tot'] >= Threshold),
}
HAPPILEE_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPILEE_Ntrls['Neps_tot'])
HAPPILEE_incl = HAPPILEE_Ntrls[['ID','Neps_tot']]
HAPPILEE_incl.rename(columns={'Neps_tot': 'P7_trls'}, inplace=True)

conditions = {
    0: (miniMADE_Ntrls['Ntrls_tot'] < Threshold),
    1: (miniMADE_Ntrls['Ntrls_tot'] >= Threshold),
}
miniMADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=miniMADE_Ntrls['Ntrls_tot'])
miniMADE_incl = miniMADE_Ntrls[['ID','Ntrls_tot']]
miniMADE_incl.rename(columns={'Ntrls_tot': 'P8_trls'}, inplace=True)


Ntrls_all = Manual_incl.merge(MADE_incl, on='ID')\
    .merge(MADEBOND_incl, on='ID')\
        .merge(HAPPEv1_incl, on='ID')\
            .merge(HAPPEv4_incl, on='ID')\
                .merge(MADEBONDld_incl, on='ID')\
                    .merge(HAPPILEE_incl, on='ID')\
                        .merge(miniMADE_incl, on='ID')

# Check values
pd.value_counts(Ntrls_all['P1_trls'] == 1)
pd.value_counts(Ntrls_all['P2_trls'] == 1)
pd.value_counts(Ntrls_all['P3_trls'] == 1)
pd.value_counts(Ntrls_all['P4_trls'] == 1)
pd.value_counts(Ntrls_all['P5_trls'] == 1)
pd.value_counts(Ntrls_all['P6_trls'] == 1)
pd.value_counts(Ntrls_all['P7_trls'] == 1)
pd.value_counts(Ntrls_all['P8_trls'] == 1)

# Run stats test
data_1 = Ntrls_all.iloc[:,1:9]
data = data_1.to_numpy()

# Import functions
from statsmodels.stats.contingency_tables import cochrans_q
from statsmodels.stats.contingency_tables import mcnemar

# Function to create contingency table and run McNemar test
def run_mcnemar(data1, data2):
    # Create 2x2 contingency table
    table = np.zeros((2, 2))
    
    # Fill contingency table
    # [0,0]: both excluded (0,0)
    # [0,1]: excluded in 1st, included in 2nd (0,1)
    # [1,0]: included in 1st, excluded in 2nd (1,0)
    # [1,1]: both included (1,1)
    table[0,0] = np.sum((data1 == 0) & (data2 == 0))
    table[0,1] = np.sum((data1 == 0) & (data2 == 1))
    table[1,0] = np.sum((data1 == 1) & (data2 == 0))
    table[1,1] = np.sum((data1 == 1) & (data2 == 1))
    
    # Run McNemar test
    result = mcnemar(table, exact=True)  # exact=True for small frequencies
    return result.pvalue

# Function to run all pairwise comparisons
def pairwise_mcnemar(data):
    n_pipelines = data.shape[1]
    p_values = np.zeros((n_pipelines, n_pipelines))
    
    # Run McNemar test for each pair
    for i in range(n_pipelines):
        for j in range(i+1, n_pipelines):
            p_value = run_mcnemar(data[:, i], data[:, j])
            p_values[i,j] = p_value
            p_values[j,i] = p_value
    
    return p_values

# Across 8 pipelines
# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')


if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")


# High-density pipelines
data_1 = Ntrls_all.iloc[:,1:6]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")



# Low-density pipelines
data_1 = Ntrls_all.iloc[:,6:9]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+6}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")









# Spectral power; threshold = 20
# read in data
Manual_Ntrls = pd.read_csv('Manual_InclRates.csv', sep = ',', header=0) 
MADE_Ntrls = pd.read_csv('MADE_InclRates.csv', sep = ',', header=0) 
MADEBOND_Ntrls = pd.read_csv('MADEBOND_InclRates.csv', sep = ',', header=0) 
HAPPEv1_Ntrls = pd.read_csv('HAPPEv1_InclRates.csv', sep = ',', header=0) 
HAPPEv4_Ntrls = pd.read_csv('HAPPEv4_InclRates.csv', sep = ',', header=0) 
MADEBONDld_Ntrls = pd.read_csv('MADE-BONDld_InclRates.csv', sep = ',', header=0) 
HAPPILEE_Ntrls = pd.read_csv('HAPPILEE_InclRates.csv', sep = ',', header=0) 
miniMADE_Ntrls = pd.read_csv('miniMADE_InclRates.csv', sep = ',', header=0) 

# convert tot to incl vs excl if tot number of trials/epochs threshold
Threshold = 20
conditions = {
    0: (Manual_Ntrls['Ntrls_tot'] < Threshold),
    1: (Manual_Ntrls['Ntrls_tot'] >= Threshold),
}
Manual_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=Manual_Ntrls['Ntrls_tot'])
Manual_incl = Manual_Ntrls[['ID','Ntrls_tot']]
Manual_incl.rename(columns={'Ntrls_tot': 'P1_trls'}, inplace=True)

conditions = {
    0: (MADE_Ntrls['Ntrls_tot'] < Threshold),
    1: (MADE_Ntrls['Ntrls_tot'] >= Threshold),
}
MADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADE_Ntrls['Ntrls_tot'])
MADE_incl = MADE_Ntrls[['ID','Ntrls_tot']]
MADE_incl.rename(columns={'Ntrls_tot': 'P2_trls'}, inplace=True)

conditions = {
    0: (MADEBOND_Ntrls['Ntrls_tot'] < Threshold),
    1: (MADEBOND_Ntrls['Ntrls_tot'] >= Threshold),
}
MADEBOND_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBOND_Ntrls['Ntrls_tot'])
MADEBOND_incl = MADEBOND_Ntrls[['ID','Ntrls_tot']]
MADEBOND_incl.rename(columns={'Ntrls_tot': 'P3_trls'}, inplace=True)

conditions = {
    0: (HAPPEv1_Ntrls['Neps_tot'] < Threshold),
    1: (HAPPEv1_Ntrls['Neps_tot'] >= Threshold),
}
HAPPEv1_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv1_Ntrls['Neps_tot'])
HAPPEv1_incl = HAPPEv1_Ntrls[['ID','Neps_tot']]
HAPPEv1_incl.rename(columns={'Neps_tot': 'P4_trls'}, inplace=True)

conditions = {
    0: (HAPPEv4_Ntrls['Neps_tot'] < Threshold),
    1: (HAPPEv4_Ntrls['Neps_tot'] >= Threshold),
}
HAPPEv4_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv4_Ntrls['Neps_tot'])
HAPPEv4_incl = HAPPEv4_Ntrls[['ID','Neps_tot']]
HAPPEv4_incl.rename(columns={'Neps_tot': 'P5_trls'}, inplace=True)

conditions = {
    0: (MADEBONDld_Ntrls['Ntrls_tot'] < Threshold),
    1: (MADEBONDld_Ntrls['Ntrls_tot'] >= Threshold),
}
MADEBONDld_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBONDld_Ntrls['Ntrls_tot'])
MADEBONDld_incl = MADEBONDld_Ntrls[['ID','Ntrls_tot']]
MADEBONDld_incl.rename(columns={'Ntrls_tot': 'P6_trls'}, inplace=True)

conditions = {
    0: (HAPPILEE_Ntrls['Neps_tot'] < Threshold),
    1: (HAPPILEE_Ntrls['Neps_tot'] >= Threshold),
}
HAPPILEE_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPILEE_Ntrls['Neps_tot'])
HAPPILEE_incl = HAPPILEE_Ntrls[['ID','Neps_tot']]
HAPPILEE_incl.rename(columns={'Neps_tot': 'P7_trls'}, inplace=True)

conditions = {
    0: (miniMADE_Ntrls['Ntrls_tot'] < Threshold),
    1: (miniMADE_Ntrls['Ntrls_tot'] >= Threshold),
}
miniMADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=miniMADE_Ntrls['Ntrls_tot'])
miniMADE_incl = miniMADE_Ntrls[['ID','Ntrls_tot']]
miniMADE_incl.rename(columns={'Ntrls_tot': 'P8_trls'}, inplace=True)


Ntrls_all = Manual_incl.merge(MADE_incl, on='ID')\
    .merge(MADEBOND_incl, on='ID')\
        .merge(HAPPEv1_incl, on='ID')\
            .merge(HAPPEv4_incl, on='ID')\
                .merge(MADEBONDld_incl, on='ID')\
                    .merge(HAPPILEE_incl, on='ID')\
                        .merge(miniMADE_incl, on='ID')

# Check values
pd.value_counts(Ntrls_all['P1_trls'] == 1)
pd.value_counts(Ntrls_all['P2_trls'] == 1)
pd.value_counts(Ntrls_all['P3_trls'] == 1)
pd.value_counts(Ntrls_all['P4_trls'] == 1)
pd.value_counts(Ntrls_all['P5_trls'] == 1)
pd.value_counts(Ntrls_all['P6_trls'] == 1)
pd.value_counts(Ntrls_all['P7_trls'] == 1)
pd.value_counts(Ntrls_all['P8_trls'] == 1)

# Run stats test
data_1 = Ntrls_all.iloc[:,1:9]
data = data_1.to_numpy()

# Across 8 pipelines
# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')


if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")


# High-density pipelines
data_1 = Ntrls_all.iloc[:,1:6]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")



# Low-density pipelines
data_1 = Ntrls_all.iloc[:,6:9]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+6}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")












# Create dataframes in prep for stats: condition effects ######################
# read in data
Manual_Ntrls = pd.read_csv('Manual_InclRates.csv', sep = ',', header=0) 
MADE_Ntrls = pd.read_csv('MADE_InclRates.csv', sep = ',', header=0) 
MADEBOND_Ntrls = pd.read_csv('MADEBOND_InclRates.csv', sep = ',', header=0) 
HAPPEv1_Ntrls = pd.read_csv('HAPPEv1_InclRates.csv', sep = ',', header=0) 
HAPPEv4_Ntrls = pd.read_csv('HAPPEv4_InclRates.csv', sep = ',', header=0) 
MADEBONDld_Ntrls = pd.read_csv('MADE-BONDld_InclRates.csv', sep = ',', header=0) 
HAPPILEE_Ntrls = pd.read_csv('HAPPILEE_InclRates.csv', sep = ',', header=0) 
miniMADE_Ntrls = pd.read_csv('miniMADE_InclRates.csv', sep = ',', header=0) 

# convert tot to incl vs excl if both soc and toy reach threshold
Threshold = 90
conditions = {
    0: (Manual_Ntrls['Ntrls_soc'] < Threshold) | (Manual_Ntrls['Ntrls_toy'] < Threshold),
    1: (Manual_Ntrls['Ntrls_soc'] >= Threshold) & (Manual_Ntrls['Ntrls_toy'] >= Threshold),
}
Manual_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=Manual_Ntrls['Ntrls_tot'])
Manual_incl = Manual_Ntrls[['ID','Ntrls_tot']]
Manual_incl.rename(columns={'Ntrls_tot': 'P1_trls'}, inplace=True)

conditions = {
    0: (MADE_Ntrls['Ntrls_soc'] < Threshold) | (MADE_Ntrls['Ntrls_toy'] < Threshold),
    1: (MADE_Ntrls['Ntrls_soc'] >= Threshold) & (MADE_Ntrls['Ntrls_toy'] >= Threshold),
}
MADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADE_Ntrls['Ntrls_tot'])
MADE_incl = MADE_Ntrls[['ID','Ntrls_tot']]
MADE_incl.rename(columns={'Ntrls_tot': 'P2_trls'}, inplace=True)

conditions = {
    0: (MADEBOND_Ntrls['Ntrls_soc'] < Threshold) | (MADEBOND_Ntrls['Ntrls_toy'] < Threshold),
    1: (MADEBOND_Ntrls['Ntrls_soc'] >= Threshold) & (MADEBOND_Ntrls['Ntrls_toy'] >= Threshold),
}
MADEBOND_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBOND_Ntrls['Ntrls_tot'])
MADEBOND_incl = MADEBOND_Ntrls[['ID','Ntrls_tot']]
MADEBOND_incl.rename(columns={'Ntrls_tot': 'P3_trls'}, inplace=True)

conditions = {
    0: (HAPPEv1_Ntrls['Neps_soc'] < Threshold) | (HAPPEv1_Ntrls['Neps_toy'] < Threshold),
    1: (HAPPEv1_Ntrls['Neps_soc'] >= Threshold) & (HAPPEv1_Ntrls['Neps_toy'] >= Threshold),
}
HAPPEv1_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv1_Ntrls['Neps_tot'])
HAPPEv1_incl = HAPPEv1_Ntrls[['ID','Neps_tot']]
HAPPEv1_incl.rename(columns={'Neps_tot': 'P4_trls'}, inplace=True)

conditions = {
    0: (HAPPEv4_Ntrls['Neps_soc'] < Threshold) | (HAPPEv4_Ntrls['Neps_toy'] < Threshold),
    1: (HAPPEv4_Ntrls['Neps_soc'] >= Threshold) & (HAPPEv4_Ntrls['Neps_toy'] >= Threshold),
}
HAPPEv4_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv4_Ntrls['Neps_tot'])
HAPPEv4_incl = HAPPEv4_Ntrls[['ID','Neps_tot']]
HAPPEv4_incl.rename(columns={'Neps_tot': 'P5_trls'}, inplace=True)

conditions = {
    0: (MADEBONDld_Ntrls['Ntrls_soc'] < Threshold) | (MADEBONDld_Ntrls['Ntrls_toy'] < Threshold),
    1: (MADEBONDld_Ntrls['Ntrls_soc'] >= Threshold) & (MADEBONDld_Ntrls['Ntrls_toy'] >= Threshold),
}
MADEBONDld_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBONDld_Ntrls['Ntrls_tot'])
MADEBONDld_incl = MADEBONDld_Ntrls[['ID','Ntrls_tot']]
MADEBONDld_incl.rename(columns={'Ntrls_tot': 'P6_trls'}, inplace=True)

conditions = {
    0: (HAPPILEE_Ntrls['Neps_soc'] < Threshold) | (HAPPILEE_Ntrls['Neps_toy'] < Threshold),
    1: (HAPPILEE_Ntrls['Neps_soc'] >= Threshold) & (HAPPILEE_Ntrls['Neps_toy'] >= Threshold),
}
HAPPILEE_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPILEE_Ntrls['Neps_tot'])
HAPPILEE_incl = HAPPILEE_Ntrls[['ID','Neps_tot']]
HAPPILEE_incl.rename(columns={'Neps_tot': 'P7_trls'}, inplace=True)

conditions = {
    0: (miniMADE_Ntrls['Ntrls_soc'] < Threshold) | (miniMADE_Ntrls['Ntrls_toy'] < Threshold),
    1: (miniMADE_Ntrls['Ntrls_soc'] >= Threshold) & (miniMADE_Ntrls['Ntrls_toy'] >= Threshold),
}
miniMADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=miniMADE_Ntrls['Ntrls_tot'])
miniMADE_incl = miniMADE_Ntrls[['ID','Ntrls_tot']]
miniMADE_incl.rename(columns={'Ntrls_tot': 'P8_trls'}, inplace=True)


Ntrls_all = Manual_incl.merge(MADE_incl, on='ID')\
    .merge(MADEBOND_incl, on='ID')\
        .merge(HAPPEv1_incl, on='ID')\
            .merge(HAPPEv4_incl, on='ID')\
                .merge(MADEBONDld_incl, on='ID')\
                    .merge(HAPPILEE_incl, on='ID')\
                        .merge(miniMADE_incl, on='ID')

# Check values
pd.value_counts(Ntrls_all['P1_trls'] == 1)
pd.value_counts(Ntrls_all['P2_trls'] == 1)
pd.value_counts(Ntrls_all['P3_trls'] == 1)
pd.value_counts(Ntrls_all['P4_trls'] == 1)
pd.value_counts(Ntrls_all['P5_trls'] == 1)
pd.value_counts(Ntrls_all['P6_trls'] == 1)
pd.value_counts(Ntrls_all['P7_trls'] == 1)
pd.value_counts(Ntrls_all['P8_trls'] == 1)

# Run stats test
data_1 = Ntrls_all.iloc[:,1:9]
data = data_1.to_numpy()

# Across 8 pipelines
# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')


if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")


# High-density pipelines
data_1 = Ntrls_all.iloc[:,1:6]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")



# Low-density pipelines
data_1 = Ntrls_all.iloc[:,6:9]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+6}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")










# convert values for different pipelines based on Threshold 20 for spectral power
Threshold = 20
conditions = {
    0: (Manual_Ntrls['Ntrls_soc'] < Threshold) | (Manual_Ntrls['Ntrls_toy'] < Threshold),
    1: (Manual_Ntrls['Ntrls_soc'] >= Threshold) & (Manual_Ntrls['Ntrls_toy'] >= Threshold),
}
Manual_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=Manual_Ntrls['Ntrls_tot'])
Manual_incl = Manual_Ntrls[['ID','Ntrls_tot']]
Manual_incl.rename(columns={'Ntrls_tot': 'P1_trls'}, inplace=True)

conditions = {
    0: (MADE_Ntrls['Ntrls_soc'] < Threshold) | (MADE_Ntrls['Ntrls_toy'] < Threshold),
    1: (MADE_Ntrls['Ntrls_soc'] >= Threshold) & (MADE_Ntrls['Ntrls_toy'] >= Threshold),
}
MADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADE_Ntrls['Ntrls_tot'])
MADE_incl = MADE_Ntrls[['ID','Ntrls_tot']]
MADE_incl.rename(columns={'Ntrls_tot': 'P2_trls'}, inplace=True)

conditions = {
    0: (MADEBOND_Ntrls['Ntrls_soc'] < Threshold) | (MADEBOND_Ntrls['Ntrls_toy'] < Threshold),
    1: (MADEBOND_Ntrls['Ntrls_soc'] >= Threshold) & (MADEBOND_Ntrls['Ntrls_toy'] >= Threshold),
}
MADEBOND_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBOND_Ntrls['Ntrls_tot'])
MADEBOND_incl = MADEBOND_Ntrls[['ID','Ntrls_tot']]
MADEBOND_incl.rename(columns={'Ntrls_tot': 'P3_trls'}, inplace=True)

conditions = {
    0: (HAPPEv1_Ntrls['Neps_soc'] < Threshold) | (HAPPEv1_Ntrls['Neps_toy'] < Threshold),
    1: (HAPPEv1_Ntrls['Neps_soc'] >= Threshold) & (HAPPEv1_Ntrls['Neps_toy'] >= Threshold),
}
HAPPEv1_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv1_Ntrls['Neps_tot'])
HAPPEv1_incl = HAPPEv1_Ntrls[['ID','Neps_tot']]
HAPPEv1_incl.rename(columns={'Neps_tot': 'P4_trls'}, inplace=True)

conditions = {
    0: (HAPPEv4_Ntrls['Neps_soc'] < Threshold) | (HAPPEv4_Ntrls['Neps_toy'] < Threshold),
    1: (HAPPEv4_Ntrls['Neps_soc'] >= Threshold) & (HAPPEv4_Ntrls['Neps_toy'] >= Threshold),
}
HAPPEv4_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPEv4_Ntrls['Neps_tot'])
HAPPEv4_incl = HAPPEv4_Ntrls[['ID','Neps_tot']]
HAPPEv4_incl.rename(columns={'Neps_tot': 'P5_trls'}, inplace=True)

conditions = {
    0: (MADEBONDld_Ntrls['Ntrls_soc'] < Threshold) | (MADEBONDld_Ntrls['Ntrls_toy'] < Threshold),
    1: (MADEBONDld_Ntrls['Ntrls_soc'] >= Threshold) & (MADEBONDld_Ntrls['Ntrls_toy'] >= Threshold),
}
MADEBONDld_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=MADEBONDld_Ntrls['Ntrls_tot'])
MADEBONDld_incl = MADEBONDld_Ntrls[['ID','Ntrls_tot']]
MADEBONDld_incl.rename(columns={'Ntrls_tot': 'P6_trls'}, inplace=True)

conditions = {
    0: (HAPPILEE_Ntrls['Neps_soc'] < Threshold) | (HAPPILEE_Ntrls['Neps_toy'] < Threshold),
    1: (HAPPILEE_Ntrls['Neps_soc'] >= Threshold) & (HAPPILEE_Ntrls['Neps_toy'] >= Threshold),
}
HAPPILEE_Ntrls['Neps_tot'] = np.select(conditions.values(), conditions.keys(), default=HAPPILEE_Ntrls['Neps_tot'])
HAPPILEE_incl = HAPPILEE_Ntrls[['ID','Neps_tot']]
HAPPILEE_incl.rename(columns={'Neps_tot': 'P7_trls'}, inplace=True)

conditions = {
    0: (miniMADE_Ntrls['Ntrls_soc'] < Threshold) | (miniMADE_Ntrls['Ntrls_toy'] < Threshold),
    1: (miniMADE_Ntrls['Ntrls_soc'] >= Threshold) & (miniMADE_Ntrls['Ntrls_toy'] >= Threshold),
}
miniMADE_Ntrls['Ntrls_tot'] = np.select(conditions.values(), conditions.keys(), default=miniMADE_Ntrls['Ntrls_tot'])
miniMADE_incl = miniMADE_Ntrls[['ID','Ntrls_tot']]
miniMADE_incl.rename(columns={'Ntrls_tot': 'P8_trls'}, inplace=True)


Ntrls_all = Manual_incl.merge(MADE_incl, on='ID')\
    .merge(MADEBOND_incl, on='ID')\
        .merge(HAPPEv1_incl, on='ID')\
            .merge(HAPPEv4_incl, on='ID')\
                .merge(MADEBONDld_incl, on='ID')\
                    .merge(HAPPILEE_incl, on='ID')\
                        .merge(miniMADE_incl, on='ID')
                        
# Check values
pd.value_counts(Ntrls_all['P1_trls'] == 1)
pd.value_counts(Ntrls_all['P2_trls'] == 1)
pd.value_counts(Ntrls_all['P3_trls'] == 1)
pd.value_counts(Ntrls_all['P4_trls'] == 1)
pd.value_counts(Ntrls_all['P5_trls'] == 1)
pd.value_counts(Ntrls_all['P6_trls'] == 1)
pd.value_counts(Ntrls_all['P7_trls'] == 1)
pd.value_counts(Ntrls_all['P8_trls'] == 1)

# Run stats test
data_1 = Ntrls_all.iloc[:,1:9]
data = data_1.to_numpy()

# Across 8 pipelines
# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')


if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")


# High-density pipelines
data_1 = Ntrls_all.iloc[:,1:6]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+1}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")



# Low-density pipelines
data_1 = Ntrls_all.iloc[:,6:9]
data = data_1.to_numpy()

# First run Cochran's Q
result = cochrans_q(data)
print(f'Cochran Q Statistic: {result.statistic:.3f}')
print(f'p-value: {result.pvalue:.4f}')

if result.pvalue < 0.05:
    # Run pairwise comparisons
    p_values = pairwise_mcnemar(data)
    
    # Create DataFrame for better visualization
    pipeline_names = [f'Pipeline{i+6}' for i in range(data.shape[1])]
    p_value_df = pd.DataFrame(p_values, 
                             index=pipeline_names,
                             columns=pipeline_names)
    
    # Calculate Bonferroni-corrected alpha
    n_comparisons = (data.shape[1] * (data.shape[1] - 1)) / 2
    bonferroni_alpha = 0.05 / n_comparisons
    
    print(f"\nBonferroni-corrected significance level: {bonferroni_alpha:.4f}")
    print("\nPairwise McNemar test p-values:")
    print(p_value_df)
    
    # Optionally, show which comparisons are significant
    significant = p_value_df < bonferroni_alpha
    print("\nSignificant differences (after Bonferroni correction):")
    for i in range(len(pipeline_names)):
        for j in range(i+1, len(pipeline_names)):
            if p_value_df.iloc[i,j] < bonferroni_alpha:
                print(f"{pipeline_names[i]} vs {pipeline_names[j]}: p = {p_value_df.iloc[i,j]:.4f}")



