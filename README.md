# LEAPsns-EEG-pipeline-comparisons
This repository contains the scripts used in the analyses reported in the manuscript entitled: ‘Comparison of EEG pre-processing pipelines for spontaneous EEG data in a large neurodivergent cohort’.

The repository consists of 6 folders for each of the different stages of the analyses:
1.	Preprocessing
2.	MetricsOfInterest
3.	InclusionRates
4.	BetweenPipelineComparisons
5.	WithinPipelineComparisons
6.	FiguresForPaper


The folder entitled ‘1.Preprocessing’ contains the folders with the scripts used for the different 8 pre-processing pipelines. The different pipelines and their steps are described in the manuscript. Dependencies and additional functions needed for the scripts are detailed in the first lines within each script. Briefly, the scripts for the following pipelines are included: 

Pipeline 2 in the manuscript: MADE
Maryland analysis of developmental EEG pipeline 
Reference: Debnath, R., Buzzell, G. A., Morales, S., Bowers, M. E., Leach, S. C., & Fox, N. A. (2020). The Maryland analysis of developmental EEG (MADE) pipeline. Psychophysiology, 57(6), e13580.
Main script used: LEAPsnsP_P2_PreprocMADE.m (in P2_MADE)

Pipeline 3 in the manuscript: MADE-BOND
the adapted MADE pipeline in the BOND lab 
Reference: Debnath, R., Buzzell, G. A., Morales, S., Bowers, M. E., Leach, S. C., & Fox, N. A. (2020). The Maryland analysis of developmental EEG (MADE) pipeline. Psychophysiology, 57(6), e13580.
Main script used: LEAPsnsP_P3_PreprocMADEBOND.m (in P3_MADEBOND)

Pipeline 4 in the manuscript: HAPPEv1
the Harvard Automated Processing Pipeline for EEG version 1 or legacy version 
Reference:  Gabard-Durnam, L. J., Mendez Leal, A. S., Wilkinson, C. L., & Levin, A. R. (2018). The Harvard Automated Processing Pipeline for Electroencephalography (HAPPE): standardized processing software for developmental and high-artifact data. Frontiers in neuroscience, 12, 97.
Main script used: LEAPsnsP_P4_PreprocHAPPEv1.m (in P4_HAPPEv1)

Pipeline 5 in the manuscript: HAPPEv4
the Harvard Automated Processing Pipeline for EEG version 4 
Reference: Lopez, K. L., Monachino, A. D., Morales, S., Leach, S. C., Bowers, M. E., & Gabard-Durnam, L. J. (2022). HAPPILEE: HAPPE In Low Electrode Electroencephalography, a standardized pre-processing software for lower density recordings. NeuroImage, 260, 119390.
Main script used: LEAPsnsP_P4_HAPPEv4_P7_HAPPILEE.m (in P5_HAPPEv4_AND_P7_HAPPILEE > HAPPE-master_March2024 > 1.pre-process, for the specific input parameters see P5_HAPPEv4_AND_P7_HAPPILEE > HAPPEv4_input_parameters_used > inputParameters.mat)

Pipeline 6 in the manuscript: MADE-BOND ld
the adapted MADE-BOND pipeline tested on a low-density layout  
Reference: Debnath, R., Buzzell, G. A., Morales, S., Bowers, M. E., Leach, S. C., & Fox, N. A. (2020). The Maryland analysis of developmental EEG (MADE) pipeline. Psychophysiology, 57(6), e13580.
Main script used: LEAPsnsP_P6_PreprocMADEBONDld.m (in P6_MADEBONDld)

Pipeline 7 in the manuscript: HAPPILEE
the HAPPE In Low Electroencephalography 
Reference: Lopez, K. L., Monachino, A. D., Morales, S., Leach, S. C., Bowers, M. E., & Gabard-Durnam, L. J. (2022). HAPPILEE: HAPPE In Low Electrode Electroencephalography, a standardized pre-processing software for lower density recordings. NeuroImage, 260, 119390.
Main script used: LEAPsnsP_P4_HAPPEv4_P7_HAPPILEE.m (in P5_HAPPEv4_AND_P7_HAPPILEE > HAPPE-master_March2024 > 1.pre-process, for the specific input parameters see P5_HAPPEv4_AND_P7_HAPPILEE > HAPPILEE_input_parameters_used > inputParameters.mat)

Pipeline 8 in the manuscript: miniMADE
The version of the MADE pipeline specifically adapted for EEG data recorded with low-density layouts
Reference: Troller‐Renfree, S. V., Morales, S., Leach, S. C., Bowers, M. E., Debnath, R., Fifer, W. P., ... & Noble, K. G. (2021). Feasibility of assessing brain activity using mobile, in‐home collection of electroencephalography: Methods and analysis. Developmental psychobiology, 63(6), e22128.
Main script used: LEAPsnsP_P8_PreprocminiMADE.m (in P8_miniMADE)

Note: For information on pipeline 1 in the manuscript: Manual
Details on the preprocessing steps applied by this pipeline are described in the following publication: Garcés, P., Baumeister, S., Mason, L., Chatham, C. H., Holiga, S., Dukart, J., ... & Hipp, J. F. (2022). Resting state EEG power spectrum and functional connectivity in autism: a cross-sectional analysis. Molecular autism, 13(1), 22.
The code uses a GUI to manually process the data which is currently still being developed and is therefore not being released in the public at this moment. However, the publication by Garcés and colleagues provides sufficient information for researchers to replicate the analyses performed here. 

For each of the pre-processing pipelines, the output is cleaned EEG data. A reporting table with the IDs of the datasets and tracking of the different data preprocessing steps is also included in the output. The cleaned data are saved in a ‘processed_data’ or ‘4 – processed’ data folder with the ‘Preproc_xxxPipelinexxx’ folder (where ‘xxxPipelinexx’ stands for the different pipelines, e.g. ‘Preproc_manual’ contains the output for the manual pipeline). The reporting table is also saved in the ‘Preproc_xxxPipelinexxx’ folder.  


The folder ‘2.MetricsOfInterest’ contains scripts to calculate the metrics of interest for each of the pipelines. For the current study, this involves spectral power and EEG functional connectivity across all epochs, and for the social and toy condition separately. Based on the pre-processing report table in the folder ‘Preproc_xxxPipelinexxx’, the scripts in this folder loop though the IDs of the dataset and read in the cleaned EEG data. The metrics of interest are calculated across all electrodes using the fieldtrip toolbox and in-house scripts. Power data are saved into the ‘A_Power_data’ subfolder within the ‘Preproc_xxxPipelinexxx’ folder, while connectivity data are saved into the ‘B_FC_data’ subfolder within the ‘Preproc_xxxPipelinexxx’ folder for each individual. In addition, power and connectivity spectra across all epochs and for the separate conditions, clean epoch numbers, and paths to saved datafiles for each individual are saved into the Metrics of Interest table and saved for further analysis. 


Folders 3, 4, and 5 contain scripts for the statistical analyses reported in the manuscript. The folder ‘3.InclusionRates’ contains the 2 scripts that were used for the analysis of the inclusion rates. ‘LEAPsnsP_InclusionRates1.m’ reads in the Metrics of Interest table for each of the pipelines and extracts the number of clean epochs (all epochs, and for separate epochs). It saves the epoch numbers for each of the pipeline in .csv files into the ‘DataForComparisons/data_csv’ folder. The Python script ‘LEAPsnsP_InclusionRates2_stats.py’ reads in the .csv files containing the inclusion rates and applies the statistical tests to the data. 

For the comparisons between pipelines reported in the manuscript, the scripts in the folder ‘4.BetweenPipelineComparisons’ were used. ‘LEAPsnsP_Between1_ICCs.m’ calculates the intra-class correlation between the metrics of interest from the different pipelines. Data are saved as ‘DATA_8pipelines.mat’ into the folder ‘DataForComparisons’ and .csv files are exported for later visualisation and saved into the ‘DataForComparisons/data_csv’ folder. The other script in the folder is ‘LEAPsnsP_Between2_ComparingICCs.m’. This code applies statistical tests to the data to test for differences in ICCs between pipelines. 

The 5th folder in the repository contains the script used for the within pipeline comparisons: ‘5.WithinPipelineComparisons’. For each of the pipelines, the Metrics of Interest table is read in from the ‘Preproc_xxxPipelinexxx’. Spectral power data is loaded. The epochs are split into datasets with odd and even epoch numbers (datasets A and B, resp.) for comparisons across all conditions, and within separate conditions to calculate condition differences. Power spectra are averaged across the selected epochs and saved into a table. This table is saved into a ‘xxxPipelinexxx_SplitHalf_reliability.mat’ file. After the calculation of the metrics for each pipeline, the SplitHalf_reliability tables are read in and the datasets with sufficient data are selected for further analysis (>= 20 epochs for power, >= 90 epochs for connectivity). The collated data are saved into a ‘DataForComparisons’ folder. Finally, the power values are averaged across frequency bands, and correlations between datasets A and B within each pipeline are calculated and saved into the ‘DataForComparisons’ folder. Correlations are also saved into .csv files in a ‘/data_csv/’ subfolder for visualisation. 

The final folder ‘6.FiguresForPaper’ contains the scripts use to generate the figures in the manuscript. The script ‘LEAPsnsP_Figures2_4_prep3.m’ plots the power and connectivity spectra, and the raincloud figures from Figure 2 in the manuscript based on the ‘DATA_8pipelines.mat’ file from the between pipeline comparisons. Next, Figure 4 from the paper is generated by plotting the statistical values from the between pipeline comparisons. Finally, the ICC values from the ‘DATA_8pipelines.mat’ file are saved into .csv files in the DataForComparisons/data_csv folder for plotting in Python. The other script in this folder is the Python script ‘LEAPsnsP_Figures3_5.py’. This script first plots the intra-class correlation matrices (while masking the upper part of the heatmap). There 16 matrices in total with 1 for each metric (power, connectivity), condition combination (across all condition, for condition differences), and for each of the 4 frequency bands of interest (delta, theta, alpha, beta). 

(written by Rianne Haartsen, PhD., 10th of February,2025)

