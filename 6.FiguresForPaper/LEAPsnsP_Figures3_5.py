#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 17:38:44 2024

@author: riannehaartsen
"""
# This script plots the correlation matrices from the between and within 
# pipeline comparisons: Figures 3 and 5 in the manuscript, resp. 

# The first part of the script creates the plots used in Figure 3 of the 
# manuscript. It reads in the ICC .csv files from each of the 
# pipelines from the between pipeline comparisons. It then plots the 
# correlation matrices using the heatmap function into 1 figure with the 16 
# matrices in separate subplots for the power and connectivity metrics 
# calculated across all epochs and for condition differences, for each 
# frequency band (delta, theta, alpha, beta).

# The second part of the script creates the plots used in Figure 5 of the 
# manuscript. It reads in the .csv files from the within 
# pipeline comparisons and plots the data using the heatmap function. 

# Further alterations to the Python figures were done in GIMP to improve 
# readability. 

# This scripts was created by Rianne Haartsen, PhD.
# Birkbeck College, University of London

# This script is released under the GNU General Public License version 3.  



# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

### Correlation matrices: comparison between pipelines ########################
# Power - all epochs
aP_delta = pd.read_csv('Pow_ICCs_alltrls_rvals_delta.csv', sep = ',', header=None) 
r_values1 = aP_delta.to_numpy()
r_values1[r_values1<0] = 0
mask1 = np.triu(np.ones_like(r_values1, dtype=bool)) # this masks the top part of the heatmap

aP_theta = pd.read_csv('Pow_ICCs_alltrls_rvals_theta.csv', sep = ',', header=None) 
r_values2 = aP_theta.to_numpy()
r_values2[r_values2<0] = 0
mask2 = np.triu(np.ones_like(r_values2, dtype=bool)) # this masks the top part of the heatmap

aP_alpha = pd.read_csv('Pow_ICCs_alltrls_rvals_alpha.csv', sep = ',', header=None) 
r_values3 = aP_alpha.to_numpy()
r_values3[r_values3<0] = 0
mask3 = np.triu(np.ones_like(r_values3, dtype=bool)) # this masks the top part of the heatmap

aP_beta = pd.read_csv('Pow_ICCs_alltrls_rvals_beta.csv', sep = ',', header=None) 
r_values4 = aP_beta.to_numpy()
r_values4[r_values4<0] = 0
mask4 = np.triu(np.ones_like(r_values4, dtype=bool)) # this masks the top part of the heatmap

# Power - condition differences
cP_delta = pd.read_csv('Pow_ICCs_cdiffs_rvals_delta.csv', sep = ',', header=None) 
r_values5 = cP_delta.to_numpy()
r_values5[r_values5<0] = 0
mask5 = np.triu(np.ones_like(r_values5, dtype=bool)) # this masks the top part of the heatmap

cP_theta = pd.read_csv('Pow_ICCs_cdiffs_rvals_theta.csv', sep = ',', header=None) 
r_values6 = cP_theta.to_numpy()
r_values6[r_values6<0] = 0
mask6 = np.triu(np.ones_like(r_values6, dtype=bool)) # this masks the top part of the heatmap

cP_alpha = pd.read_csv('Pow_ICCs_cdiffs_rvals_alpha.csv', sep = ',', header=None) 
r_values7 = cP_alpha.to_numpy()
r_values7[r_values7<0] = 0
mask7 = np.triu(np.ones_like(r_values7, dtype=bool)) # this masks the top part of the heatmap

cP_beta = pd.read_csv('Pow_ICCs_cdiffs_rvals_beta.csv', sep = ',', header=None) 
r_values8 = cP_beta.to_numpy()
r_values8[r_values4<0] = 0
mask8 = np.triu(np.ones_like(r_values8, dtype=bool)) # this masks the top part of the heatmap

# FC - all epochs
aF_delta = pd.read_csv('FC_ICCs_alltrls_rvals_delta.csv', sep = ',', header=None) 
r_values9 = aF_delta.to_numpy()
r_values9[r_values9<0] = 0
mask9 = np.triu(np.ones_like(r_values9, dtype=bool)) # this masks the top part of the heatmap

aF_theta = pd.read_csv('FC_ICCs_alltrls_rvals_theta.csv', sep = ',', header=None) 
r_values10 = aF_theta.to_numpy()
r_values10[r_values10<0] = 0
mask10 = np.triu(np.ones_like(r_values10, dtype=bool)) # this masks the top part of the heatmap

aF_alpha = pd.read_csv('FC_ICCs_alltrls_rvals_alpha.csv', sep = ',', header=None) 
r_values11 = aF_alpha.to_numpy()
r_values11[r_values3<0] = 0
mask11 = np.triu(np.ones_like(r_values11, dtype=bool)) # this masks the top part of the heatmap

aF_beta = pd.read_csv('FC_ICCs_alltrls_rvals_beta.csv', sep = ',', header=None) 
r_values12 = aF_beta.to_numpy()
r_values12[r_values4<0] = 0
mask12 = np.triu(np.ones_like(r_values12, dtype=bool)) # this masks the top part of the heatmap

# FC - condition differences
cF_delta = pd.read_csv('FC_ICCs_cdiffs_rvals_delta.csv', sep = ',', header=None) 
r_values13 = cF_delta.to_numpy()
r_values13[r_values5<0] = 0
mask13 = np.triu(np.ones_like(r_values13, dtype=bool)) # this masks the top part of the heatmap

cF_theta = pd.read_csv('FC_ICCs_cdiffs_rvals_theta.csv', sep = ',', header=None) 
r_values14 = cF_theta.to_numpy()
r_values14[r_values14<0] = 0
mask14 = np.triu(np.ones_like(r_values14, dtype=bool)) # this masks the top part of the heatmap

cF_alpha = pd.read_csv('FC_ICCs_cdiffs_rvals_alpha.csv', sep = ',', header=None) 
r_values15 = cF_alpha.to_numpy()
r_values15[r_values15<0] = 0
mask15 = np.triu(np.ones_like(r_values15, dtype=bool)) # this masks the top part of the heatmap

cF_beta = pd.read_csv('FC_ICCs_cdiffs_rvals_beta.csv', sep = ',', header=None) 
r_values16 = cF_beta.to_numpy()
r_values16[r_values16<0] = 0
mask16 = np.triu(np.ones_like(r_values16, dtype=bool)) # this masks the top part of the heatmap


# info for plots
Measures1 = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7','P8']
cmap2 = sns.color_palette('Spectral', as_cmap=True)
# generate figure with subplots
f,((ax1, ax2, ax3, ax4, axcb1),(ax5, ax6, ax7, ax8, axcb2), (ax9, ax10, ax11, ax12, axcb3), (ax13, ax14, ax15, ax16, axcb4)) = plt.subplots(4,5,gridspec_kw={'width_ratios':[1,1,1,1,0.08]})
# link axes
ax1.get_shared_y_axes().join(ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16)
# Plots the heatmaps in the subplots
# Power - all
g11 = sns.heatmap(np.round(r_values1, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask1, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8}, cbar=False,ax=ax1)
g11.set_title('Power all - \u03B4')

g12 = sns.heatmap(np.round(r_values2, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask2, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8}, cbar=False,ax=ax2)
g12.set_title('Power all - \u03B8')

g13 = sns.heatmap(np.round(r_values3, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask3, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},cbar=False,ax=ax3)
g13.set_title('Power all - \u03B1')

g14 = sns.heatmap(np.round(r_values4, decimals=2),cmap=cmap2, 
                  vmin=0, vmax=1.0, mask = mask4, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},ax=ax4, cbar_ax=axcb1)
g14.set_title('Power all - \u03B2')

# Power - condition differences
g21 = sns.heatmap(np.round(r_values5, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask5, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8}, cbar=False,ax=ax5)
g21.set_title('Power diffs - \u03B4')

g22 = sns.heatmap(np.round(r_values6, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask6, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8}, cbar=False,ax=ax6)
g22.set_title('Power diffs - \u03B8')

g23 = sns.heatmap(np.round(r_values7, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask7, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},cbar=False,ax=ax7)
g23.set_title('Power diffs - \u03B1')

g24 = sns.heatmap(np.round(r_values8, decimals=2),cmap=cmap2, 
                  vmin=0, vmax=1.0, mask = mask8, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},ax = ax8, cbar_ax=axcb2)
g24.set_title('Power diffs - \u03B2')

# FC - all
g31 = sns.heatmap(np.round(r_values9, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask9, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8}, cbar=False,ax=ax9)
g31.set_title('FC all - \u03B4')

g32 = sns.heatmap(np.round(r_values10, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask10, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},cbar=False,ax=ax10)
g32.set_title('FC all - \u03B8')

g33 = sns.heatmap(np.round(r_values11, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask11, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},cbar=False,ax=ax11)
g33.set_title('FC all - \u03B1')

g34 = sns.heatmap(np.round(r_values12, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask14, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8} ,ax=ax12, cbar_ax=axcb3)
g34.set_title('FC all - \u03B2')

# FC - condition differences
g41 = sns.heatmap(np.round(r_values13, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask13, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},cbar=False,ax=ax13)
g41.set_title('FC diffs - \u03B4')

g42 = sns.heatmap(np.round(r_values14, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask14, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},cbar=False,ax=ax14)
g42.set_title('FC diffs - \u03B8')

g43 = sns.heatmap(np.round(r_values15, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask15, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},cbar=False,ax=ax15)
g43.set_title('FC diffs - \u03B1')

g44 = sns.heatmap(np.round(r_values16, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, mask = mask16, square = True, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8},ax=ax16, cbar_ax=axcb4)
g44.set_title('FC diffs - \u03B2')

for ax in [g11,g12, g13, g14, g21, g22, g23, g24, g31, g32, g33, g34, g41, g42, g43, g44]:
    # tl = ax.get_xticklabels()
    # ax.set_xticklabels(tl, rotation=90)
    ax.set_xticks(np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]))
    ax.set_xticklabels(['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'])
    ax.set_yticks(np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]))
    ax.set_yticklabels(['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7','P8'])
    tly = ax.get_yticklabels()
    ax.set_yticklabels(tly, rotation=0)

plt.show()
plt.tight_layout()



### Correlation matrices: split-half reliability ##############################
# Power - all epochs
PowAll_corrs_r = pd.read_csv('SplitHalf_8P_all_rvals.csv', sep = ',', header=None) 
r_values_all = PowAll_corrs_r.to_numpy()

# Power - condition differences
PowCD_corrs_r = pd.read_csv('SplitHalf_8P_diff_rvals.csv', sep = ',', header=None) 
r_values_cd = PowCD_corrs_r.to_numpy()

cmap2 = sns.color_palette('Spectral', as_cmap=True)
Pipelines = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7','P8']
FreqBand = ['\u03B4', '\u03B8', '\u03B1', '\u03B2']

f,((ax1, ax2, axcb1)) = plt.subplots(1,3,gridspec_kw={'width_ratios':[1,1,0.08]})
# link axes
ax1.get_shared_y_axes().join(ax2)
# Plots the heatmaps in the subplots
# Power - all
g11 = sns.heatmap(np.round(r_values_all, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0,  center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8}, cbar=False,ax=ax1)
g11.set_xticklabels(['\u03B4', '\u03B8', '\u03B1', '\u03B2'])
g11.set_yticklabels(['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7','P8'])
g11.set_title('Across all trials')

g12 = sns.heatmap(np.round(r_values_cd, decimals=2),cmap=cmap2,
                  vmin=0, vmax=1.0, center = 0, 
                  annot=True, linewidths=.5, annot_kws={'fontsize': 8}, cbar=True,ax=ax2,cbar_ax=axcb1)
g12.set_xticklabels(['\u03B4', '\u03B8', '\u03B1', '\u03B2'])
g12.set_yticklabels(['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7','P8'])
g12.set_title('Condition differences')

for ax in [g11,g12]:
    # tl = ax.get_xticklabels()
    # ax.set_xticklabels(tl, rotation=90)
    tly = ax.get_yticklabels()
    ax.set_yticklabels(tly, rotation=0)

plt.tight_layout()


