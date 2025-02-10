%% Analysis of social and non-social videos in LEAP: Comparison between pipelines

% This script extracts the EEG metrics of interest for each of the 8 
% pipelines: Manual, MADE, BOND-MADE, HAPPEv1, HAPPEv4, MADE-BONDld, 
% HAPPILEE, and miniMADE.

% The following metrics are extracted:

% 1) Whole brain power values: for delta, theta, alpha, and beta frequencies
% 2) Whole brain connectivity values: for delta, theta, alpha, and beta frequencies

% The current script calls to the ICC function which can be downloaded
% from:
% https://uk.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc

% Data are saved into the DataForComparison folder, and are exported to csv 
% files for further analysis and visualisation in Python

% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')

% Created by Rianne Haartsen, PhD.; 05-2024 
% Birkbeck College, University of London

% This script is released under the GNU General Public License version 3.

%% Extract spectra for spectral power and connectivity for each of the pipelines
% load data from the Measures of Interest tables
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
% epochs
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

% construct structure
DATA_8pipelines = struct();
DATA_8pipelines.pipelines_names = {'Manual','MADE','MADE-BOND','HAPPEv1', ...
    'HAPPEv4','MADE-BONDld','HAPPILEE','miniMADE'};
DATA_8pipelines.pipelines_MoI_tables = {Manual_MOI_table, MADE_MOI_table, MADEBOND_MOI_table, HAPPEv1_MOI_table, ...
    HAPPEv4_MOI_table, MADEBOND_ld_MOI_table, HAPPILEE_MOI_table, miniMADE_MOI_table};

% Select those with sufficient data
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

    % Get data - all trials
    for ii = 1:length(DATA_8pipelines.pipelines_names)
        Sel_powerspectra = cell2mat(table2array(DATA_8pipelines.pipelines_MoI_tables{1,ii}(Ind_incl_pow,10)));
        DATA_8pipelines.PowerSpectraA_data{ii} = Sel_powerspectra; 
        clear Sel_powerspectra
        Sel_fcspectra = cell2mat(table2array(DATA_8pipelines.pipelines_MoI_tables{1,ii}(Ind_incl_fc,14)));
        DATA_8pipelines.FunConSpectraA_data{ii} = Sel_fcspectra; 
        clear Sel_fcspectra
    end
    % Get data - soc trials
    for ii = 1:length(DATA_8pipelines.pipelines_names)
        Sel_powerspectra = cell2mat(table2array(DATA_8pipelines.pipelines_MoI_tables{1,ii}(Ind_incl_pow,11)));
        DATA_8pipelines.PowerSpectraS_data{ii} = Sel_powerspectra; 
        clear Sel_powerspectra
        Sel_fcspectra = cell2mat(table2array(DATA_8pipelines.pipelines_MoI_tables{1,ii}(Ind_incl_fc,15)));
        DATA_8pipelines.FunConSpectraS_data{ii} = Sel_fcspectra; 
        clear Sel_fcspectra
    end
    % Get data - toy trials
    for ii = 1:length(DATA_8pipelines.pipelines_names)
        Sel_powerspectra = cell2mat(table2array(DATA_8pipelines.pipelines_MoI_tables{1,ii}(Ind_incl_pow,12)));
        DATA_8pipelines.PowerSpectraT_data{ii} = Sel_powerspectra; 
        clear Sel_powerspectra
        Sel_fcspectra = cell2mat(table2array(DATA_8pipelines.pipelines_MoI_tables{1,ii}(Ind_incl_fc,16)));
        DATA_8pipelines.FunConSpectraT_data{ii} = Sel_fcspectra; 
        clear Sel_fcspectra
    end
    % Get IDs
    DATA_8pipelines.PowerSpectra_IDxFreq = {table2array(DATA_8pipelines.pipelines_MoI_tables{1,1}(Ind_incl_pow,1)), 1:1:32};
    DATA_8pipelines.FunConSpectra_IDxFreq = {table2array(DATA_8pipelines.pipelines_MoI_tables{1,1}(Ind_incl_fc,1)), 1:1:32};

% save data structure
cd xxx/DataForComparisons
for ii = 1:length(DATA_8pipelines.pipelines_names)
    % power
    writematrix(DATA_8pipelines.PowerSpectraA_data{ii}, strcat(DATA_8pipelines.pipelines_names{1,ii},'_pow_spectra_alltrls.csv'))
    % FC
    writematrix(DATA_8pipelines.FunConSpectraA_data{ii}, strcat(DATA_8pipelines.pipelines_names{1,ii},'_fc_spectra_alltrls.csv'))
end
% other
save('Incl_indices.mat','Ind_incl_pow', 'Ind_incl_fc')
save('DATA_8pipelines.mat','DATA_8pipelines')


%% ICC between pipelines: across all trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /Users/riannehaartsen/Documents/MATLAB/ICC
% Code for calculation of the ICC can be found on: 
% https://uk.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc

Freqs = 1:1:32;
Delta_ind = [find(Freqs == 2),find(Freqs == 3)];
Theta_ind = [find(Freqs == 4),find(Freqs == 6)];
Alpha_ind = [find(Freqs == 7),find(Freqs == 12)];
Beta_ind = [find(Freqs == 13),find(Freqs == 30)];

Freqs_inds = [Delta_ind; Theta_ind; Alpha_ind; Beta_ind];

% Power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Powerband_values_all = zeros(size(DATA_8pipelines.PowerSpectra_IDxFreq{1,1},1), length(DATA_8pipelines.pipelines_names), size(Freqs_inds,1));
    for ff = 1:size(Freqs_inds,1)
        for pp = 1:length(DATA_8pipelines.pipelines_names)
            Powerband_values_all(:,pp,ff) = mean(DATA_8pipelines.PowerSpectraA_data{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
        end
    end
    clear pp ff
    
    Powerband_forICC = reshape(Powerband_values_all, [size(Powerband_values_all,1), size(Powerband_values_all,2)*size(Powerband_values_all,3)]);
    
    % run ICC for all combinations
    r_vals = zeros(size(Powerband_forICC,2));
    LB_vals = zeros(size(Powerband_forICC,2));
    UB_vals = zeros(size(Powerband_forICC,2));
    p_vals = zeros(size(Powerband_forICC,2));
    for rr = 1:size(Powerband_forICC,2)
        for cc = 1:size(Powerband_forICC,2)
            Mat = [Powerband_forICC(:,rr) Powerband_forICC(:,cc)];
            [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
            r_vals(rr, cc) = r;
            LB_vals(rr, cc) = LB;
            UB_vals(rr, cc) = UB;
            p_vals(rr,cc) = p;
            clear Mar r p
        end
    end
    
    % visual checks
    % ICCs
    figure; imagesc(r_vals)
    c = colorbar; c.Label.String = 'ICC'; %c.Limits = [0 1];
    title('ICCs across log power for different frequency bands and pipelines')
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    yticks(1:8:32)
    yticklabels({'Delta','Theta','Alpha','Beta'})
    % power band values
    figure; imagesc(Powerband_forICC)
    c = colorbar; c.Label.String = 'Log power'; 
    ylabel('Subjects');
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    title('Power within bands for each subject')
    
    % r_vals(r_vals < 0) = 0;
    
    % save data
    cd xxx/DataForComparisons/data_csv
    % power
    writematrix(r_vals, 'Pow_ICCs_alltrls_rvals_masked.csv')
    writematrix(LB_vals, 'Pow_ICCs_alltrls_LBvals.csv')
    writematrix(UB_vals, 'Pow_ICCs_alltrls_UBvals.csv')
    writematrix(p_vals, 'Pow_ICCs_alltrls_pvals.csv')
    % data
    DATA_8pipelines.Powerband_all.Values = Powerband_values_all;
    DATA_8pipelines.Powerband_all.ForICC = Powerband_forICC;
    DATA_8pipelines.Powerband_all.Freq_inds = Freqs_inds;
    DATA_8pipelines.Powerband_all.ICCs_r_LB_UB_pvals = {r_vals, LB_vals, UB_vals, p_vals};
    cd xxx/DataForComparisons/
    save('DATA_8pipelines.mat','DATA_8pipelines')

    clear r_vals LB_vals UB_vals p_vals
    clear Powerband_values_all Powerband_forICC

% Connectivity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FunConband_values_all = zeros(size(DATA_8pipelines.FunConSpectra_IDxFreq{1,1},1), length(DATA_8pipelines.pipelines_names), size(Freqs_inds,1));
        for ff = 1:size(Freqs_inds,1)
            for pp = 1:length(DATA_8pipelines.pipelines_names)
                FunConband_values_all(:,pp,ff) = mean(DATA_8pipelines.FunConSpectraA_data{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
            end
        end
        clear pp ff
        FunConband_forICC = reshape(FunConband_values_all, [size(FunConband_values_all,1), size(FunConband_values_all,2)*size(FunConband_values_all,3)]);
        
    % run ICC for all combinations
    r_vals = zeros(size(FunConband_forICC,2));
    LB_vals = zeros(size(FunConband_forICC,2));
    UB_vals = zeros(size(FunConband_forICC,2));
    p_vals = zeros(size(FunConband_forICC,2));
    for rr = 1:size(FunConband_forICC,2)
        for cc = 1:size(FunConband_forICC,2)
            Mat = [FunConband_forICC(:,rr) FunConband_forICC(:,cc)];
            [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
            r_vals(rr, cc) = r;
            LB_vals(rr, cc) = LB;
            UB_vals(rr, cc) = UB;
            p_vals(rr,cc) = p;
            clear Mar r p
        end
    end
    clear rr cc

    % r_vals(r_vals < 0) = 0;

    % visual checks
    % ICCs
    figure; imagesc(r_vals)
    c = colorbar; c.Label.String = 'ICC'; 
    title('ICCs across dbWPLI for different frequency bands and pipelines')
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    yticks(1:8:32)
    yticklabels({'Delta','Theta','Alpha','Beta'})
    % fc band values
    figure; imagesc(FunConband_forICC)
    c = colorbar; c.Label.String = 'dbWPLI'; 
    ylabel('Subjects');
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    title('FC within bands for each subject')
    

    % save data
    cd xxx/DataForComparisons/data_csv
    % functional connectivity
    writematrix(r_vals, 'FC_ICCs_alltrls_rvals_masked.csv')
    writematrix(LB_vals, 'FC_ICCs_alltrls_LBvals.csv')
    writematrix(UB_vals, 'FC_ICCs_alltrls_UBvals.csv')
    writematrix(p_vals, 'FC_ICCs_alltrls_pvals.csv')
    % data
    DATA_8pipelines.FunConband_all.Values = FunConband_values_all;
    DATA_8pipelines.FunConband_all.ForICC = FunConband_forICC;
    DATA_8pipelines.FunConband_all.Freq_inds = Freqs_inds;
    DATA_8pipelines.FunConband_all.ICCs_r_LB_UB_pvals = {r_vals, LB_vals, UB_vals, p_vals};
    cd xxx/DataForComparisons/
    save('DATA_8pipelines.mat','DATA_8pipelines')
    clear r_vals LB_vals UB_vals p_vals
    clear FunConband_values_all FunConband_forICC

%% ICC between pipelines: Condition differences %%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/Users/riannehaartsen/Documents/MATLAB/ICC')

Freqs_inds = DATA_8pipelines.Powerband_all.Freq_inds;

% For power differences
    Powerband_values_diffs = zeros(size(DATA_8pipelines.PowerSpectra_IDxFreq{1,1},1), length(DATA_8pipelines.pipelines_names), size(Freqs_inds,1));
        for ff = 1:size(Freqs_inds,1)
            for pp = 1:length(DATA_8pipelines.pipelines_names)
                Powerband_values_diffs(:,pp,ff) = [mean(DATA_8pipelines.PowerSpectraS_data{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2) - ...
                    mean(DATA_8pipelines.PowerSpectraT_data{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2)];
            end
        end
        clear ii ff
 
    Powerband_forICC = reshape(Powerband_values_diffs, [size(Powerband_values_diffs,1), size(Powerband_values_diffs,2)*size(Powerband_values_diffs,3)]);

    % run ICC for all combinations
    r_vals = zeros(size(Powerband_forICC,2));
    LB_vals = zeros(size(Powerband_forICC,2));
    UB_vals = zeros(size(Powerband_forICC,2));
    p_vals = zeros(size(Powerband_forICC,2));
    for rr = 1:size(Powerband_forICC,2)
        for cc = 1:size(Powerband_forICC,2)
            Mat = [Powerband_forICC(:,rr) Powerband_forICC(:,cc)];
            [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
            r_vals(rr, cc) = r;
            LB_vals(rr, cc) = LB;
            UB_vals(rr, cc) = UB;
            p_vals(rr,cc) = p;
            clear Mar r p
        end
    end
    
    % r_vals(r_vals < 0) = 0;

    % visual checks
    % ICCs
    figure; imagesc(r_vals)
    c = colorbar; c.Label.String = 'ICC'; %c.Limits = [0 1];
    title('ICCs across log power condition differences for different frequency bands and pipelines')
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    yticks(1:8:32)
    yticklabels({'Delta','Theta','Alpha','Beta'})
    % power band values
    figure; imagesc(Powerband_forICC)
    c = colorbar; c.Label.String = 'Log power'; 
    ylabel('Subjects');
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    title('Power condition differences within bands for each subject')
    
    % save data
    cd xxx/DataForComparisons/data_csv
    % power
    writematrix(r_vals, 'Pow_ICCs_diffs_rvals_masked.csv')
    writematrix(LB_vals, 'Pow_ICCs_diffs_LBvals.csv')
    writematrix(UB_vals, 'Pow_ICCs_diffs_UBvals.csv')
    writematrix(p_vals, 'Pow_ICCs_diffs_pvals.csv')
    % data
    DATA_8pipelines.Powerband_cdiffs.Values = Powerband_values_diffs;
    DATA_8pipelines.Powerband_cdiffs.ForICC = Powerband_forICC;
    DATA_8pipelines.Powerband_cdiffs.Freq_inds = Freqs_inds;
    DATA_8pipelines.Powerband_cdiffs.ICCs_r_LB_UB_pvals = {r_vals, LB_vals, UB_vals, p_vals};
    cd xxx/DataForComparisons/
    save('DATA_8pipelines.mat','DATA_8pipelines')
    
    clear r_vals LB_vals UB_vals p_vals
    clear Powerband_values_diffs Powerband_forICC


% For connectivity
    FunConband_values_diffs = zeros(size(DATA_8pipelines.FunConSpectra_IDxFreq{1,1},1), length(DATA_8pipelines.pipelines_names), size(Freqs_inds,1));
        for ff = 1:size(Freqs_inds,1)
            for pp = 1:length(DATA_8pipelines.pipelines_names)
                FunConband_values_diffs(:,pp,ff) = [mean(DATA_8pipelines.FunConSpectraS_data{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2) - ...
                    mean(DATA_8pipelines.FunConSpectraT_data{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2)];
            end
        end
        clear ii ff
        
        FunConband_forICC = reshape(FunConband_values_diffs, [size(FunConband_values_diffs,1), size(FunConband_values_diffs,2)*size(FunConband_values_diffs,3)]);
    
    % run ICC for all combinations
    r_vals = zeros(size(FunConband_forICC,2));
    LB_vals = zeros(size(FunConband_forICC,2));
    UB_vals = zeros(size(FunConband_forICC,2));
    p_vals = zeros(size(FunConband_forICC,2));
    for rr = 1:size(FunConband_forICC,2)
        for cc = 1:size(FunConband_forICC,2)
            Mat = [FunConband_forICC(:,rr) FunConband_forICC(:,cc)];
            [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
            r_vals(rr, cc) = r;
            LB_vals(rr, cc) = LB;
            UB_vals(rr, cc) = UB;
            p_vals(rr,cc) = p;
            clear Mar r p
        end
    end
    
    % r_vals(r_vals < 0) = 0;

    % visual checks
    % ICCs
    figure; imagesc(r_vals)
    c = colorbar; c.Label.String = 'ICC'; %c.Limits = [0 1];
    title('ICCs across dbWPLI condition differences for different frequency bands and pipelines')
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    yticks(1:8:32)
    yticklabels({'Delta','Theta','Alpha','Beta'})
    % fc band values
    figure; imagesc(FunConband_forICC)
    c = colorbar; c.Label.String = 'dbWPLI'; 
    ylabel('Subjects');
    xticks(1:8:32)
    xticklabels({'Delta','Theta','Alpha','Beta'})
    title('FC condition differenwithin bands for each subject')
    
    % save data
    cd xxx/DataForComparisons/data_csv
    % power
    writematrix(r_vals, 'FC_ICCs_diff_rvals_masked.csv')
    writematrix(LB_vals, 'FC_ICCs_diff_LBvals.csv')
    writematrix(UB_vals, 'FC_ICCs_diff_UBvals.csv')
    writematrix(p_vals, 'FC_ICCs_diff_pvals.csv')
    % data
    DATA_8pipelines.FunConband_cdiffs.Values = FunConband_values_diffs;
    DATA_8pipelines.FunConband_cdiffs.ForICC = FunConband_forICC;
    DATA_8pipelines.FunConband_cdiffs.Freq_inds = Freqs_inds;
    DATA_8pipelines.FunConband_cdiffs.ICCs_r_LB_UB_pvals = {r_vals, LB_vals, UB_vals, p_vals};
    cd xxx/DataForComparisons/
    save('DATA_8pipelines.mat','DATA_8pipelines')
    clear r_vals LB_vals UB_vals p_vals
    clear FunConband_values_diffs FunConband_forICC
    

    

