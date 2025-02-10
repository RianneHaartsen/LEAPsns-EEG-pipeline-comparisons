%% LEAP soc-nsoc videos Pipelines: get power and connectivity values from the data preprocessed with the manual pipeline 
% P1 - pipeline 1 in the manuscript

% This script takes the data processed with the Manual pipeline and converges
% it to a fieldtrip structure. Then power is calculated
% for 1 to 32 Hz in steps of 1Hz for each trial. Finally, data are split
% into conditions and connectivity and power are calculated for the social 
% condition, non-social condition, and across all conditions. 

% Processing is looped through participants with paths from the
% Manual_preprocessing_report.mat file. Preprocessed data are saved. 
% The measures of interest are stored into a reporting table named
% Manual_MOI_table and saved as Manual_MeasuresOfInterest_report.mat. 

% This script call to functions from fieldtrip. Fieldtrip can be downloaded
% from: https://www.fieldtriptoolbox.org/

% Also cite the following in the manuscript using these functions:

% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. 
% FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and 
% Invasive Electrophysiological Data. Computational Intelligence and 
% Neuroscience, 2011; 2011:156869.


% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')
% additional comments by Rianne Haartsen are noted with '%RH note;'

% Adapted by Rianne Haartsen, PhD.; 06-2022 
% Birkbeck College, University of London

% This script is released under the GNU General Public License version 3.

%%

clear % clear matlab workspace
clc % clear matlab command window
addpath('xxx');% enter the path of the EEGLAB folder in this line
eeglab % call eeglab to set up the plugin
addpath('xxx');  % add fieldtrip path
ft_defaults
addpath(genpath('xxx/Preproc_Manual'))

load xxx/Preproc_Manual/Manual_preprocessing_report.mat

Pow_location = 'xxx/Preproc_Manual/A_Power_data';
FC_location = 'xxx/Preproc_Manual/B_FC_data';
output_location = 'xxx/Preproc_Manual';
Fig_location = 'xxx/Preproc_Manual/Report_figs';

for subject = 1 :height(Manual_report_table)
    
    CurrID = Manual_report_table.ID{subject};
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, CurrID);
    
    try
        % load preprocessed data
        filename_path_preproc = Manual_report_table.datafile_names_preproc{subject};
        load(filename_path_preproc);
        ft_EEG = data;
        clear data

        %% Split into conditions
        Ntrls_tot = size(ft_EEG.trialinfo,1);
        % social
        if sum(ft_EEG.trialinfo == 910) > 0
            cfg                 = [];
            cfg.trials          = ft_EEG.trialinfo == 910;
            ft_EEG_soc = ft_selectdata(cfg, ft_EEG);
        end
        Ntrls_soc = sum(ft_EEG.trialinfo == 910);
        % toy
        if sum(ft_EEG.trialinfo == 911) > 0
            cfg                 = [];
            cfg.trials          = ft_EEG.trialinfo == 911;
            ft_EEG_toy = ft_selectdata(cfg, ft_EEG);
        end
        Ntrls_toy = sum(ft_EEG.trialinfo == 911);
           
        %% Parameters across all analyses
        % for power analysis
            cfgpow                 = [];
            cfgpow.output          = 'pow';
            cfgpow.method          = 'mtmfft';
            cfgpow.taper           = 'hanning';
            cfgpow.foi             = 1:1:32;  
            cfgpow.keeptrials      = 'yes';
        % for FFT in connectivity analysis
            cfgfft             = [];
            cfgfft.method      = 'mtmfft';
            cfgfft.taper       =  'hanning'; 
            cfgfft.output      = 'fourier'; 
            cfgfft.tapsmofrq   = 1;
            cfgfft.foi             = 1:1:32; 
            cfgfft.keeptrials  = 'yes'; 
            % frequencies (in points) to examine
            F = 1:32;
            
        %% For all trials
        if Ntrls_tot > 1
            % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Power_all = ft_freqanalysis(cfgpow, ft_EEG);
                Power_all.logpowspctrm = log(Power_all.powspctrm);
                Power_all.logpowspctrm_avg = squeeze(mean(Power_all.logpowspctrm,1));
                Power_all.logpowspctrm_avg_global = squeeze(mean(Power_all.logpowspctrm_avg,1));
            % dbWPLI analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find segments with overlap
                if size(ft_EEG.trialinfo,1) > 1
                    OVERLAP_all = zeros(size(ft_EEG.trialinfo,1), size(ft_EEG.trialinfo,1));
                    for i=1: size(ft_EEG.trialinfo,1)-1
                        if ft_EEG.sampleinfo(i,2) > ft_EEG.sampleinfo((i+1),1) % if overlaps
                           OVERLAP_all (i, i+1)=1;  OVERLAP_all (i+1, i)=1; 
                        end
                    end
                end
                OVERLAP_all(find(eye(size(ft_EEG.trialinfo,1))==1)) = 1;
                % FFT
                FFT_all = ft_freqanalysis(cfgfft, ft_EEG);
                % dbWPLI
                PLI = []; WPLI = []; ubPLI = []; dbWPLI = [];
                [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(FFT_all.fourierspctrm, size(FFT_all.fourierspctrm,1), F, OVERLAP_all); %%dPLIbasic
                FC_all = struct();
                FC_all.PLI = PLI;
                FC_all.WPLI = WPLI;
                FC_all.ubPLI = ubPLI;
                FC_all.dbWPLI = dbWPLI;
                FC_all.Freqs_Hz = FFT_all.freq(F);
                FC_all.Chan_label = FFT_all.label;
                FC_all.overlap = OVERLAP_all;
                % global dbWPLI
                GlobdbWPLI_all = zeros(size(FC_all.dbWPLI,1),1);
                ff = [];
                for ff = 1:size(FC_all.dbWPLI,1)
                    FCmat = squeeze(FC_all.dbWPLI(ff,:,:));
                    Mat1 = tril(FCmat,-1);
                    Mat1(Mat1 == 0) = NaN;
                    GlobdbWPLI_all(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                    Mat1 =[]; FCmat = [];
                end
                FC_all.dbWPLI_glob = GlobdbWPLI_all';  
        end
        
        
        %% For social trials
        if Ntrls_soc > 1
            % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Power_soc = ft_freqanalysis(cfgpow, ft_EEG_soc);
                Power_soc.logpowspctrm = log(Power_soc.powspctrm);
                Power_soc.logpowspctrm_avg = squeeze(mean(Power_soc.logpowspctrm,1));
                Power_soc.logpowspctrm_avg_global = squeeze(mean(Power_soc.logpowspctrm_avg,1));
            % dbWPLI analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find segments with overlap
                if size(ft_EEG_soc.trialinfo,1) > 1
                OVERLAP_soc = zeros(size(ft_EEG_soc.trialinfo,1), size(ft_EEG_soc.trialinfo,1));
                    for i=1: size(ft_EEG_soc.trialinfo,1)-1
                        if ft_EEG_soc.sampleinfo(i,2) > ft_EEG_soc.sampleinfo((i+1),1) % if overlaps
                           OVERLAP_soc (i, i+1)=1;  OVERLAP_soc (i+1, i)=1; 
                        end
                    end
                end
                OVERLAP_soc(find(eye(size(ft_EEG_soc.trialinfo,1))==1)) = 1;
                % FFT
                FFT_soc = ft_freqanalysis(cfgfft, ft_EEG_soc);
                % dbWPLI
                PLI = []; WPLI = []; ubPLI = []; dbWPLI = [];
                [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(FFT_soc.fourierspctrm, size(FFT_soc.fourierspctrm,1), F, OVERLAP_soc); %%dPLIbasic
                FC_soc = struct();
                FC_soc.PLI = PLI;
                FC_soc.WPLI = WPLI;
                FC_soc.ubPLI = ubPLI;
                FC_soc.dbWPLI = dbWPLI;
                FC_soc.Freqs_Hz = FFT_soc.freq(F);
                FC_soc.Chan_label = FFT_soc.label;
                FC_soc.overlap = OVERLAP_soc;
                % global dbWPLI
                GlobdbWPLI_soc = zeros(size(FC_soc.dbWPLI,1),1);
                ff = [];
                for ff = 1:size(FC_soc.dbWPLI,1)
                    FCmat = squeeze(FC_soc.dbWPLI(ff,:,:));
                    Mat1 = tril(FCmat,-1);
                    Mat1(Mat1 == 0) = NaN;
                    GlobdbWPLI_soc(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                    Mat1 = []; FCmat = [];
                end
                FC_soc.dbWPLI_glob = GlobdbWPLI_soc';    
        end
        
        %% For toy trials
        if Ntrls_toy > 1
            % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Power_toy = ft_freqanalysis(cfgpow, ft_EEG_toy);
                Power_toy.logpowspctrm = log(Power_toy.powspctrm);
                Power_toy.logpowspctrm_avg = squeeze(mean(Power_toy.logpowspctrm,1));
                Power_toy.logpowspctrm_avg_global = squeeze(mean(Power_toy.logpowspctrm_avg,1));
            % dbWPLI analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find segments with overlap
                if size(ft_EEG_toy.trialinfo,1) > 1
                OVERLAP_toy = zeros(size(ft_EEG_toy.trialinfo,1), size(ft_EEG_toy.trialinfo,1));
                    for i=1: size(ft_EEG_toy.trialinfo,1)-1
                        if ft_EEG_toy.sampleinfo(i,2) > ft_EEG_toy.sampleinfo((i+1),1) % if overlaps
                           OVERLAP_toy (i, i+1)=1;  OVERLAP_toy (i+1, i)=1; 
                        end
                    end
                end
            OVERLAP_toy(find(eye(size(ft_EEG_toy.trialinfo,1))==1)) = 1;
            % FFT
                FFT_toy = ft_freqanalysis(cfgfft, ft_EEG_toy);
                % dbWPLI
                PLI = []; WPLI = []; ubPLI = []; dbWPLI = [];
                [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(FFT_toy.fourierspctrm, size(FFT_toy.fourierspctrm,1), F, OVERLAP_toy); %%dPLIbasic
                FC_toy = struct();
                FC_toy.PLI = PLI;
                FC_toy.WPLI = WPLI;
                FC_toy.ubPLI = ubPLI;
                FC_toy.dbWPLI = dbWPLI;
                FC_toy.Freqs_Hz = FFT_toy.freq(F);
                FC_toy.Chan_label = FFT_toy.label;
                FC_toy.overlap = OVERLAP_toy;
                % global dbWPLI
                GlobdbWPLI_toy = zeros(size(FC_toy.dbWPLI,1),1);
                ff = [];
                for ff = 1:size(FC_toy.dbWPLI,1)
                    FCmat = squeeze(FC_toy.dbWPLI(ff,:,:));
                    Mat1 = tril(FCmat,-1);
                    Mat1(Mat1 == 0) = NaN;
                    GlobdbWPLI_toy(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                    Mat1 = []; FCmat = [];
                end
                FC_toy.dbWPLI_glob = GlobdbWPLI_toy';
        end
        
    catch
        warning('Unable to calculate measures of interest')
    end
    
    %% Check data and create empty variables if not created
        % fieldtrip data
        if ~exist('ft_EEG','var')
            ft_EEG = [];
        end
        if ~exist('ft_EEG_soc','var')
            ft_EEG_soc = [];
        end
        if ~exist('ft_EEG_toy','var')
            ft_EEG_toy = [];
        end
        % power data
        if ~exist('Power_all','var')
            Power_all = [];
        end
        if ~exist('Power_soc','var')
            Power_soc = [];
        end
        if ~exist('Power_toy','var')
            Power_toy = [];
        end
        % fc data
        if ~exist('FFT_all','var')
            FFT_all = [];
        end
        if ~exist('FC_all','var')
            FC_all = [];
        end
        if ~exist('FFT_soc','var')
            FFT_soc = [];
        end
        if ~exist('FC_soc','var')
            FC_soc = [];
        end
        if ~exist('FFT_toy','var')
            FFT_toy = [];
        end
        if ~exist('FC_toy','var')
            FC_toy = [];
        end
        % number of trials
        if ~exist('Ntrls_tot','var')
            Ntrls_tot = [];
        end
        if ~exist('Ntrls_soc','var')
            Ntrls_soc = [];
        end
        if ~exist('Ntrls_toy','var')
            Ntrls_toy = [];
        end

    
    %% save data and add to current subject to report table 
    % keep track of saved data location and subject preprocessed
    % Info
    preproc_folder = 'xxx';
    Info_file = strcat(preproc_folder, filesep, 'Preproc_', CurrID,'_social_nonsocial_videos.mat');
    load(Info_file)
    Site{subject} = data.cfg.euaims.site;
    clear data preproc_folder Info_file
    DateCur{subject} = datestr(now,'dd-mm-yyyy');
    Ft_file{subject} = strcat(output_location, filesep, 'processed_data', filesep, CurrID, '_fieldtrip_data.mat');
    save(Ft_file{subject},'ft_EEG','ft_EEG_soc','ft_EEG_toy')
    % Power
    Power_file{subject} = strcat(Pow_location, filesep, CurrID, '_Power_data.mat');
    save(Power_file{subject}, 'Power_all','Power_soc','Power_toy')
    if ~isempty(Power_all)
        Power_global_all = Power_all.logpowspctrm_avg_global;
        Freqs = Power_all.freq;
    else
        Power_global_all = NaN(1,32);
        Freqs = NaN(1,32);
    end
    if ~isempty(Power_soc)
        Power_global_soc = Power_soc.logpowspctrm_avg_global;
    else
        Power_global_soc = NaN(1,32);
    end
    if ~isempty(Power_toy)
        Power_global_toy = Power_toy.logpowspctrm_avg_global;
    else
        Power_global_toy = NaN(1,32);
    end
    % FFT and FC
    FC_file{subject} = strcat(FC_location, filesep, CurrID,  '_FC_data.mat');
    save(FC_file{subject}, 'FFT_all','FFT_soc','FFT_toy','FC_all','FC_soc','FC_toy')
    if ~isempty(FC_all)
        dbWPLI_global_all = FC_all.dbWPLI_glob;
    else
        dbWPLI_global_all = NaN(1,32);
    end
    if ~isempty(FC_soc)
        dbWPLI_global_soc = FC_soc.dbWPLI_glob;
    else
        dbWPLI_global_soc = NaN(1,32);
    end
    if ~isempty(FC_toy)
        dbWPLI_global_toy = FC_toy.dbWPLI_glob;
    else
        dbWPLI_global_toy = NaN(1,32);
    end
    
    % add to tracking table
    cd(output_location)
    load Manual_MeasuresOfInterest_report.mat
        data_table_newrow = table({CurrID}, Site(1,subject), DateCur(1,subject), Ft_file(1,subject), {Ntrls_tot}, {Ntrls_soc}, {Ntrls_toy}, ...
        {Freqs}, Power_file(1,subject), {Power_global_all},{Power_global_soc},{Power_global_toy},...
        FC_file(1,subject), {dbWPLI_global_all},{dbWPLI_global_soc},{dbWPLI_global_toy});
    data_table_newrow.Properties.VariableNames={'ID','Site', 'Date','Fieldtrip_file','Ntrls_tot','Ntrls_soc','Ntrls_toy', ...
        'Freqs','Power_path','GlobPow_all','GlobPow_soc','GlobPow_toy',...
        'FC_path','GlobFC_all','GlobFC_soc','GlobFC_toy'};
    Manual_MOI_table(subject,:) = data_table_newrow;
    save('Manual_MeasuresOfInterest_report.mat','Manual_MOI_table');
    clear Manual_MOI_table
    
    % clear up variables
    clear ft_EEG ft_EEG_soc ft_EEG_toy 
    clear Power_all Power_soc Power_toy 
    clear FFT_all FFT_soc FFT_toy FC_all FC_soc FC_toy 
    clear Ntrls_tot Ntrls_soc Ntrls_toy
    clear OVERLAP_all OVERLAP_soc OVERLAP_toy
    clear GlobdbWPLI_all GlobdbWPLI_soc GlobdbWPLI_toy
    clear dbWPLI_global_all dbWPLI_global_soc dbWPLI_global_toy
    clear Power_global_all Power_global_soc Power_global_toy
    clear Site FC_file Pow_file data_table_newrow
    clear cfg cfgfft cfgpow CurrID DateCur
    clear EEG
    clear PLI ubPLI WPLI dbWPLI
    
 
end





 