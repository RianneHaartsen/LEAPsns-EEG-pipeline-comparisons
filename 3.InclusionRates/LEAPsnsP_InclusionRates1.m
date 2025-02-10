%% Analysis of social and non-social videos in LEAP: Inclusion rates

% This script counts the number of participants with sufficent epochs for
% power (>= 20), and functional connectivity (>= 90) for comparisons across
% all conditions, and condition differences, for each of the 8 pipelines. 

% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')

% Created by Rianne Haartsen, PhD.; 08-2024 
% Birkbeck College, University of London

% This script is released under the GNU General Public License version 3.


%% Check inclusion rates

% load data from the Measures of Interest tables
% P1
load xxx/Preproc_Manual/Manual_MeasuresOfInterest_report.mat 
% P2
load xxx/Preproc_MADE/MADE_MeasuresOfInterest_report.mat 
% P3
load xxx/Preproc_MADEBOND/MADEBOND_MeasuresOfInterest_report.mat 
%P4
load xxx/Preproc_HAPPEv1/HAPPEv1_MeasuresOfInterest_report.mat 
%P5
load xxx/Preproc_HAPPEv4/HAPPEv4_MeasuresOfInterest_report.mat 
%P6
load xxx/Preproc_MADEBONDld/MADEBOND_ld_MeasuresOfInterest_report.mat 
%P7
load xxx/Preproc_HAPPILEE/HAPPILEE_MeasuresOfInterest_report.mat 
%P8
load xxx/Preproc_miniMADE/miniMADE_MeasuresOfInterest_report.mat 

% substitute empty cells with a value of 0 for those without any clean
% epochs (This was specific to the current study.)
Manual_MOI_table.Ntrls_tot{104} = 0;
Manual_MOI_table.Ntrls_soc{104} = 0;
Manual_MOI_table.Ntrls_toy{104} = 0;
Manual_MOI_table.Ntrls_tot{14} = 0;
Manual_MOI_table.Ntrls_soc{14} = 0;
Manual_MOI_table.Ntrls_toy{14} = 0;

MADE_MOI_table.Ntrls_tot{14} = 0;
MADE_MOI_table.Ntrls_soc{14} = 0;
MADE_MOI_table.Ntrls_toy{14} = 0;
MADE_MOI_table.Ntrls_tot{35} = 0;
MADE_MOI_table.Ntrls_soc{35} = 0;
MADE_MOI_table.Ntrls_toy{35} = 0;


%% Inclusion rates: counts and percentages
% Across all trials
Ntrls_tot = [cell2mat(Manual_MOI_table.Ntrls_tot) ...
    cell2mat(MADE_MOI_table.Ntrls_tot) ...
    cell2mat(MADEBOND_MOI_table.Ntrls_tot) ...
    cell2mat(HAPPEv1_MOI_table.Neps_tot) ...
    cell2mat(HAPPEv4_MOI_table.Neps_tot) ...
    cell2mat(MADEBOND_ld_MOI_table.Ntrls_tot) ...
    cell2mat(HAPPILEE_MOI_table.Neps_tot) ...
    cell2mat(miniMADE_MOI_table.Ntrls_tot)];
% for pow and fc
Ind_incl_pow = Ntrls_tot >= 90;
PerP = sum(Ind_incl_pow,1); 
PerP_perc = round(PerP/131*100,0); 
[PerP' PerP_perc']
% for pow
Ind_incl_pow = Ntrls_tot >= 20;
PerP = sum(Ind_incl_pow,1); 
PerP_perc = round(PerP/131*100,0); 
[PerP' PerP_perc']
% not enough
Ind_incl_pow = Ntrls_tot < 20 & Ntrls_tot > 0;
PerP = sum(Ind_incl_pow,1); 
PerP_perc = round(PerP/131*100,0); 
[PerP' PerP_perc']
% no epochs
Ind_incl_pow = Ntrls_tot == 0;
PerP = sum(Ind_incl_pow,1); 
PerP_perc = round(PerP/131*100,0); 
[PerP' PerP_perc']

% Condition differences
Ntrls_cdiffs = [cell2mat(Manual_MOI_table.Ntrls_soc) cell2mat(Manual_MOI_table.Ntrls_toy)...
    cell2mat(MADE_MOI_table.Ntrls_soc) cell2mat(MADE_MOI_table.Ntrls_toy) ...
    cell2mat(MADEBOND_MOI_table.Ntrls_soc) cell2mat(MADEBOND_MOI_table.Ntrls_toy) ...
    cell2mat(HAPPEv1_MOI_table.Neps_soc) cell2mat(HAPPEv1_MOI_table.Neps_toy)...
    cell2mat(HAPPEv4_MOI_table.Neps_soc) cell2mat(HAPPEv4_MOI_table.Neps_toy) ...
    cell2mat(MADEBOND_ld_MOI_table.Ntrls_soc) cell2mat(MADEBOND_ld_MOI_table.Ntrls_toy) ...
    cell2mat(HAPPILEE_MOI_table.Neps_soc) cell2mat(HAPPILEE_MOI_table.Neps_toy) ...
    cell2mat(miniMADE_MOI_table.Ntrls_soc) cell2mat(miniMADE_MOI_table.Ntrls_toy)];

% for pow and fc
Ind_incl_pow = Ntrls_cdiffs >= 90;
Ind_incl_pow_cd = [Ind_incl_pow(:,1)+Ind_incl_pow(:,2) Ind_incl_pow(:,3)+Ind_incl_pow(:,4)...
    Ind_incl_pow(:,5)+Ind_incl_pow(:,6) Ind_incl_pow(:,7)+Ind_incl_pow(:,8)...
    Ind_incl_pow(:,9)+Ind_incl_pow(:,10) Ind_incl_pow(:,11)+Ind_incl_pow(:,12)...
    Ind_incl_pow(:,13)+Ind_incl_pow(:,14) Ind_incl_pow(:,15)+Ind_incl_pow(:,16)];
PerP = sum(Ind_incl_pow_cd == 2,1); 
PerP_perc = round(PerP/131*100,0); 
[PerP' PerP_perc']
% for pow
Ind_incl_pow = Ntrls_cdiffs >= 20;
Ind_incl_pow_cd = [Ind_incl_pow(:,1)+Ind_incl_pow(:,2) Ind_incl_pow(:,3)+Ind_incl_pow(:,4)...
    Ind_incl_pow(:,5)+Ind_incl_pow(:,6) Ind_incl_pow(:,7)+Ind_incl_pow(:,8)...
    Ind_incl_pow(:,9)+Ind_incl_pow(:,10) Ind_incl_pow(:,11)+Ind_incl_pow(:,12)...
    Ind_incl_pow(:,13)+Ind_incl_pow(:,14) Ind_incl_pow(:,15)+Ind_incl_pow(:,16)];
PerP = sum(Ind_incl_pow_cd == 2,1); 
PerP_perc = round(PerP/131*100,0); 
[PerP' PerP_perc']
% no epochs in both
Ind_incl_pow = Ntrls_cdiffs == 0;
Ind_incl_pow_cd = [Ind_incl_pow(:,1)+Ind_incl_pow(:,2) Ind_incl_pow(:,3)+Ind_incl_pow(:,4)...
    Ind_incl_pow(:,5)+Ind_incl_pow(:,6) Ind_incl_pow(:,7)+Ind_incl_pow(:,8)...
    Ind_incl_pow(:,9)+Ind_incl_pow(:,10) Ind_incl_pow(:,11)+Ind_incl_pow(:,12)...
    Ind_incl_pow(:,13)+Ind_incl_pow(:,14) Ind_incl_pow(:,15)+Ind_incl_pow(:,16)];
PerP = sum(Ind_incl_pow_cd == 2,1); 
PerP_perc = round(PerP/131*100,0); 
[PerP' PerP_perc']
% the rest: aka not enough
AllP = ones(8,1)*131; 
    % more than 20 for both
    Ind_incl_pow = Ntrls_cdiffs >= 20;
    Ind_incl_pow_cd = [Ind_incl_pow(:,1)+Ind_incl_pow(:,2) Ind_incl_pow(:,3)+Ind_incl_pow(:,4)...
        Ind_incl_pow(:,5)+Ind_incl_pow(:,6) Ind_incl_pow(:,7)+Ind_incl_pow(:,8)...
        Ind_incl_pow(:,9)+Ind_incl_pow(:,10) Ind_incl_pow(:,11)+Ind_incl_pow(:,12)...
        Ind_incl_pow(:,13)+Ind_incl_pow(:,14) Ind_incl_pow(:,15)+Ind_incl_pow(:,16)];
    PerP = sum(Ind_incl_pow_cd == 2,1); 
    % 0 for both
    Ind_incl_pow2 = Ntrls_cdiffs == 0;
    Ind_incl_pow2_cd = [Ind_incl_pow2(:,1)+Ind_incl_pow2(:,2) Ind_incl_pow2(:,3)+Ind_incl_pow2(:,4)...
        Ind_incl_pow2(:,5)+Ind_incl_pow2(:,6) Ind_incl_pow2(:,7)+Ind_incl_pow2(:,8)...
        Ind_incl_pow2(:,9)+Ind_incl_pow2(:,10) Ind_incl_pow2(:,11)+Ind_incl_pow2(:,12)...
        Ind_incl_pow2(:,13)+Ind_incl_pow2(:,14) Ind_incl_pow2(:,15)+Ind_incl_pow2(:,16)];
    PerP2 = sum(Ind_incl_pow2_cd == 2,1); 
PerP_left = AllP - PerP' - PerP2';
PerP_left_perc = round(PerP_left /131*100,0);
[PerP_left PerP_left_perc]


% Select IDs with sufficient data
Ntrls = [cell2mat(Manual_MOI_table.Ntrls_soc) cell2mat(Manual_MOI_table.Ntrls_toy)...
    cell2mat(MADE_MOI_table.Ntrls_soc) cell2mat(MADE_MOI_table.Ntrls_toy) ...
    cell2mat(MADEBOND_MOI_table.Ntrls_soc) cell2mat(MADEBOND_MOI_table.Ntrls_toy) ...
    cell2mat(HAPPEv1_MOI_table.Neps_soc) cell2mat(HAPPEv1_MOI_table.Neps_toy)...
    cell2mat(HAPPEv4_MOI_table.Neps_soc) cell2mat(HAPPEv4_MOI_table.Neps_toy) ...
    cell2mat(MADEBOND_ld_MOI_table.Ntrls_soc) cell2mat(MADEBOND_ld_MOI_table.Ntrls_toy) ...
    cell2mat(HAPPILEE_MOI_table.Neps_soc) cell2mat(HAPPILEE_MOI_table.Neps_toy) ...
    cell2mat(miniMADE_MOI_table.Ntrls_soc) cell2mat(miniMADE_MOI_table.Ntrls_toy)];
Ind_incl_pow = find(sum(Ntrls >= 20, 2) ==16);
Ind_incl_fc = find(sum(Ntrls >= 90,2) ==16);


%% Save inclusion rates for each pipeline for later statistical analyses in Python

% P1 - Manual
load xxx/Preproc_Manual/Manual_MeasuresOfInterest_report.mat
Incl_Manual = Manual_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_Manual)
    if isempty(Incl_Manual.Ntrls_tot{ii}) || isnan(Incl_Manual.Ntrls_tot{ii})
        Incl_Manual.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_Manual.Ntrls_soc{ii}) || isnan(Incl_Manual.Ntrls_soc{ii})
        Incl_Manual.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_Manual.Ntrls_toy{ii}) || isnan(Incl_Manual.Ntrls_toy{ii})
        Incl_Manual.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_Manual,'Manual_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_Manual Manual_MOI_table

% P2 - MADE
load xxx/Preproc_MADE/MADE_MeasuresOfInterest_report.mat
Incl_MADE = MADE_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_MADE)
    if isempty(Incl_MADE.Ntrls_tot{ii}) || isnan(Incl_MADE.Ntrls_tot{ii})
        Incl_MADE.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_MADE.Ntrls_soc{ii}) || isnan(Incl_MADE.Ntrls_soc{ii})
        Incl_MADE.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_MADE.Ntrls_toy{ii}) || isnan(Incl_MADE.Ntrls_toy{ii})
        Incl_MADE.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_MADE,'MADE_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_MADE MADE_MOI_table

% P3 - MADE-BOND
load xxx/Preproc_MADEBOND/MADEBOND_MeasuresOfInterest_report.mat
Incl_MADEBOND = MADEBOND_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_MADEBOND)
    if isempty(Incl_MADEBOND.Ntrls_tot{ii}) || isnan(Incl_MADEBOND.Ntrls_tot{ii})
        Incl_MADEBOND.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_MADEBOND.Ntrls_soc{ii}) || isnan(Incl_MADEBOND.Ntrls_soc{ii})
        Incl_MADEBOND.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_MADEBOND.Ntrls_toy{ii}) || isnan(Incl_MADEBOND.Ntrls_toy{ii})
        Incl_MADEBOND.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_MADEBOND,'MADEBOND_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_MADEBOND MADEBOND_MOI_table

% P4 - HAPPEv1
load xxx/Preproc_HAPPEv1/HAPPEv1_MeasuresOfInterest_report.mat
Incl_HAPPEv1 = HAPPEv1_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_HAPPEv1)
    if isempty(Incl_HAPPEv1.Neps_tot{ii}) || isnan(Incl_HAPPEv1.Neps_tot{ii})
        Incl_HAPPEv1.Neps_tot{ii} = 0;
    end
    if isempty(Incl_HAPPEv1.Neps_soc{ii}) || isnan(Incl_HAPPEv1.Neps_soc{ii})
        Incl_HAPPEv1.Neps_soc{ii} = 0;
    end
    if isempty(Incl_HAPPEv1.Neps_toy{ii}) || isnan(Incl_HAPPEv1.Neps_toy{ii})
        Incl_HAPPEv1.Neps_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_HAPPEv1,'HAPPEv1_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_HAPPEv1 HAPPEv1_MOI_table

% P5 - HAPPEv4
load xxx/Preproc_HAPPEv4/HAPPE_MeasuresOfInterest_report.mat
Incl_HAPPEv4 = HAPPEv4_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_HAPPEv4)
    if isempty(Incl_HAPPEv4.Neps_tot{ii}) || isnan(Incl_HAPPEv4.Neps_tot{ii})
        Incl_HAPPEv4.Neps_tot{ii} = 0;
    end
    if isempty(Incl_HAPPEv4.Neps_soc{ii}) || isnan(Incl_HAPPEv4.Neps_soc{ii})
        Incl_HAPPEv4.Neps_soc{ii} = 0;
    end
    if isempty(Incl_HAPPEv4.Neps_toy{ii}) || isnan(Incl_HAPPEv4.Neps_toy{ii})
        Incl_HAPPEv4.Neps_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_HAPPEv4,'HAPPEv4_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_HAPPEv4 HAPPEv4_MOI_table

% P6 - MADE-BONDld
load xxx/Preproc_MADEBOND/MADEBOND_ld_MeasuresOfInterest_report.mat
Incl_MADEBOND = MADEBOND_ld_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_MADEBOND)
    if isempty(Incl_MADEBOND.Ntrls_tot{ii}) || isnan(Incl_MADEBOND.Ntrls_tot{ii})
        Incl_MADEBOND.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_MADEBOND.Ntrls_soc{ii}) || isnan(Incl_MADEBOND.Ntrls_soc{ii})
        Incl_MADEBOND.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_MADEBOND.Ntrls_toy{ii}) || isnan(Incl_MADEBOND.Ntrls_toy{ii})
        Incl_MADEBOND.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_MADEBOND,'MADEBONDld_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_MADEBOND MADEBOND_ld_MOI_table

% P7 - HAPPILEE
load xxx/Preproc_HAPPILEE/HAPPILEE_MeasuresOfInterest_report.mat
Incl_HAPPILEE = HAPPILEE_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_HAPPILEE)
    if isempty(Incl_HAPPILEE.Neps_tot{ii}) || isnan(Incl_HAPPILEE.Neps_tot{ii})
        Incl_HAPPILEE.Neps_tot{ii} = 0;
    end
    if isempty(Incl_HAPPILEE.Neps_soc{ii}) || isnan(Incl_HAPPILEE.Neps_soc{ii})
        Incl_HAPPILEE.Neps_soc{ii} = 0;
    end
    if isempty(Incl_HAPPILEE.Neps_toy{ii}) || isnan(Incl_HAPPILEE.Neps_toy{ii})
        Incl_HAPPILEE.Neps_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_HAPPILEE,'HAPPILEE_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_HAPPILEE HAPPILEE_MOI_table

% P8 - miniMADE
load xxx/Preproc_miniMADE/miniMADE_MeasuresOfInterest_report.mat
Incl_miniMADE = miniMADE_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_miniMADE)
    if isempty(Incl_miniMADE.Ntrls_tot{ii}) || isnan(Incl_miniMADE.Ntrls_tot{ii})
        Incl_miniMADE.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_miniMADE.Ntrls_soc{ii}) || isnan(Incl_miniMADE.Ntrls_soc{ii})
        Incl_miniMADE.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_miniMADE.Ntrls_toy{ii}) || isnan(Incl_miniMADE.Ntrls_toy{ii})
        Incl_miniMADE.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/DataForComparisons/data_csv
writetable(Incl_miniMADE,'miniMADE_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_miniMADE miniMADE_MOI_table

