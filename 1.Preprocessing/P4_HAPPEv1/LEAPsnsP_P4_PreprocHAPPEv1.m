%% Analysis of social and non-social videos in LEAP: Pipeline 3 - HAPPEv1
% The HAPPE version 1 pipeline was adjusted by Philipp Bomatter and Pilar Garces at Roche for 
% pre-processing of the resting state task in the LEAP study. Adjustments
% included making the pipeline compatible with EEG data formats from the
% LEAP study. 

% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')
% additional comments by Rianne Haartsen are noted with '%RH note;'

% Adapted by Rianne Haartsen, PhD.; 12-2022
% Birkbeck College, University of London

% This adapted pipeline is released under the GNU General Public License version 3.




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The Harvard Automated Pre-Processing Pipeline for EEG (HAPPE)
% Version 1.0
%
% Developed at Boston Children's Hospital Labs of Cognitive Neuroscience
%
% For a detailed description of the pipeline and user options, please see the following manuscript:
% Gabard-Durnam, et al., (submitted).To be updated upon acceptance for publication.
% Contributors to HAPPE:
%
% Laurel Joy Gabard-Durnam (laurel.gabarddurnam@gmail.com)
% Adriana S. Mendez Leal (asmendezleal@gmail.com)
% Carol L. Wilkinson (carol.wilkinson@childrens.harvard.edu)
% April R. Levin (april.levin@childrens.harvard.edu)
%
% HAPPE includes code that is dependent on the following third-party software.
% Please reference this third-party software in any manuscripts making use of HAPP-E as below:
%
% EEGLab A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
% analysis of single-trial EEG dynamics. Journal of Neuroscience Methods 134:9-21
%
% Cleanline by Tim Mullen (as an EEGlab plug-in): Mullen, T. (2012). NITRC:
% CleanLine: Tool/Resource Info. Available online at: http://www.nitrc.org/projects/cleanline.
%
% MARA by Irene Winkler (as an EEGlab plug-in): Winkler, et al. Automatic
% Classification of Artifactual ICA-Components for Artifact Removal in
% EEG Signals. Behavioral and Brain Functions 7:30 (2011).
%
% FASTER segment-level channel interpolation code: Nolan*, H., Whelan*, R.,
% & Reilly, R.B. (2010). FASTER: Fully Automated Statistical Thresholding
% for EEG artifact 	Rejection. Journal of Neuroscience Methods, 192,152-162.
%
% I have modified Matlab central code for the wICA function originally posted by
% Jordan Sorokin (2015) https://www.mathworks.com/matlabcentral/fileexchange/55413-wica-data-varargin-
%
% Any code that is not part of the third-party dependencies is released
% under the GNU General Public License version 3.

% This software is being distributed with the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See GNU General
% Public License for more details.
%
% In no event shall Boston Children�s Hospital (BCH), the BCH Division of Developmental Medicine, the
% Laboratories of Cognitive Neuroscience (LCN), or pipeline contributors to BEAPP be
% liable to any party for direct, indirect, special, incidental, or
% consequential damages, including lost profits, arising out of the use of
% this software and its documentation, even if Boston Children�s Hospital,
% the Lab of Cognitive Neuroscience, and software contributors have been
% advised of the possibility of such damage. Software and documentation is
% provided �as is.� Boston Children�s Hospital, the Lab of Cognitive
% Neuroscience, and software contributors are under no obligation to
% provide maintenance, support, updates, enhancements, or modifications.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License (version 3) as
% published by the Free Software Foundation.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%   10 USER INPUTS TO EDIT:
clear all
cd ' xxx path for folder with scripts '

% 1. enter path to the folder that has the datasets you want to analyze
src_folder_name = 'xxx/RawEEG_data/';
out_data_folder_name = 'xxx/Preproc_HAPPEv1/processed_data';
out_folder_name = 'xxx/Preproc_HAPPEv1';

% 2. Which acquisition layout would you like to use?
% RH note: We take care of this for the LEAP data.

% Line Noise Frequency
line_noise_freq = 50;

% Option to resume: will check the given output folder for files that have
% already been processed and skip them.
resume = true;

% 3. list channels of interest, including the 10-20 channels. User defined
% channels occur at the end of the sequence e.g. 'E39' 'E40' the 18 "10-20"
% channels that NEED to be in the chan_IDs: 'FP1' 'FP2' 'F3' 'F4' 'F7' 'F8'
% 'C3' 'C4' 'T3' 'T4' 'PZ' 'O1' 'O2' 'T5' 'T6' 'P3' 'P4' 'Fz' 
% Note: the HAPPE script uses the old nomenclature where
% T3/T4 = T7/T8 and T5/T6 = P7/P8.

channels_of_interest = {'FP1' 'FP2' 'AF7' 'AF3' 'AF4' 'AF8' ...
    'F7' 'F5' 'F3' 'F1' 'Fz' 'F2' 'F4' 'F6' 'F8' ...
    'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'FC6' 'FT8' ...
    'T3' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T4' ...
    'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' ...
    'T5' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'T6' ...
    'PO7' 'PO3' 'POz' 'PO4' 'PO8' ...
    'O1' 'Oz' 'O2'};

full_channel_locs = load('channel_locs.mat').channel_locs;


% 4. run HAPPE in semi-automated setting with visualizations (=1) or fully-automated, no visualizations setting ( = 0)
pipeline_visualizations_semiautomated = 0;

% if semi-automated, what is the minimum and maximum frequency you want to visualize in the power spectrum figure for each file?
vis_freq_min = 2;
vis_freq_max = 57;

% if semi-automated, which frequencies do you want to generate spatial topoplots for within the figure?
freq_to_plot = [6 10 20 30 55];

% 5. for resting-state EEG, set task_EEG_processing = 0
% for task-related EEG, set task_EEG_processing = 1
% for social/non-social videos, set task_EEG_processing = 2;
task_EEG_processing = 2;

% if task-related EEG:
task_conditions = {'near', 'devr'}; % enter the stimulus condition tags

% 6. do you want to segment your data? yes (=1) or no (=0)
segment_data = 1;

% if you are segmenting your task-related EEG:
% parameters to segment the data for each stimulus, in seconds:
task_segment_start = -0.5;
task_segment_end = 1.5;

% if you are segmenting your resting-state EEG:
% how long do you want your segments to be in seconds? here, 2 seconds is the default
segment_length = 1;

% 7. do you want to interpolate the specific channels' data determined to be artifact/bad within each segment?
% yes = 1, no = 0.
% This is segment-level channel interpolation from the FASTER EEGlab plug-in.
segment_interpolation = 1;  % 0 to prevent discontinuities within the same channel due to segment interpolation

% 8. do you want to do segment rejection (using amplitude and joint probability criteria)?
% yes = 1, no = 0.
segment_rejection = 1;  % set = 0 to only mark bad segments without rejecting them

% if you are rejecting segments, what minimum/maximum signal amplitude do you want to use as the artifact threshold?
reject_min_amp = -40;
reject_max_amp = 40;

% do you want to do segment rejection using all user-specified channels above ( = 0) or a subset of channels in an ROI ( = 1)?
ROI_channels_only = 0;

% if you want to do ROI segment rejection, which channels should be used in the ROI?
% ex ROI_channels = {'E27','E20'};
ROI_channels = {'E27','E20'};

% 9. Select the type of re-referencing you want. Average re-reference (=1)
% or re-referencing to another channel/subset of channels (=0)
average_rereference = 1;

% if you are referencing to another channel/subset of channels, what are they?
% make sure to use the channel name given above
% ex ROI_channels = {'E57','E100'};
NO_AVERAGE_REREF_channel_subset = {'E57','E100'};

% 10. Select the format to save your processed data at the end of HAPPE!
%save_as_format = 1 will save the processed data as a .txt file.(electrodes as columns, time as rows)
%save_as_format = 2 will save the processed data as a .mat file (matlab format)
%save_as_format = 3 will save the processed data as a .set file (EEGlab format)
save_as_format = 3;

%~~~~~~~~~~~~~~~~~~~~~~ no need to edit beyond this point ~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% create output folder:
if ~exist(out_data_folder_name, 'dir')
    mkdir(out_data_folder_name)
end

%% add relevant folders to path

% add HAPPE script path
% happe_directory_path = fileparts(which('HAPPE_pipeline_v1_0.m'));
happe_directory_path = 'xxx';

% will eventually allow users to set own eeglab path -- for now, assume
% using eeglab14_0_0b included in HAPPE 
eeglab_path = [happe_directory_path filesep 'Packages' filesep 'eeglab14_0_0b'];
% eeglab_path = '/Users/riannehaartsen/Documents/MATLAB/eeglab2022.0';

% add HAPPE subfolders and EEGLAB plugin folders to path
addpath([happe_directory_path filesep 'scripts'], eeglab_path, genpath([eeglab_path filesep 'functions']));
rmpath(genpath([eeglab_path filesep 'functions' filesep 'octavefunc']));
eeglab nogui;
pop_editoptions('option_storedisk', 0, 'option_savetwofiles', 0, 'option_saveversion6', 1, ...
            'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, ...
            'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, ...
            'option_donotusetoolboxes', 0, 'option_checkversion', 0, 'option_chat', 0);

plugin_directories = dir([eeglab_path filesep 'plugins']);
plugin_directories = strcat(eeglab_path,filesep,'plugins',filesep,{plugin_directories.name},';');
addpath([plugin_directories{:}]);

% add cleanline path
if exist('cleanline','file')
    cleanline_path = which('eegplugin_cleanline.m');
    cleanline_path = cleanline_path(1:findstr(cleanline_path,'eegplugin_cleanline.m')-1);
    addpath(genpath(cleanline_path));
else
    error('Please make sure cleanline is on your path');
end

%% pull file names to feed script and initialize the arrays to store file specific data quality metrics
src_file_ext = '.set';
FileNames = dir([src_folder_name filesep '*' src_file_ext]);
FileNames = {FileNames.name};

% If enabled, discard already processed files
% if resume
%     out_files = {dir([out_folder_name filesep '*' '_processed.set']).name};
%     processed_files = arrayfun(@(f) ismember(strrep(erase(f, 'Unfiltered-'), src_file_ext, '_processed.set'), out_files), FileNames);
%     FileNames = FileNames(~processed_files);
% end

n_files = length(FileNames);

% intialize report metrics
Number_ICs_Rejected = zeros(1,n_files);
Number_Good_Channels_Selected = zeros(1,n_files);
Interpolated_Channel_IDs = cell(1,n_files);
Percent_ICs_Rejected = zeros(1,n_files);
Percent_Variance_Kept_of_Post_Waveleted_Data = zeros(1,n_files);
File_Length_In_Secs = zeros(1,n_files);
File_Length_afterRSsegementing_In_Secs = zeros(1,n_files);
Number_Channels_User_Selected = zeros(1,n_files);
Percent_Good_Channels_Selected = zeros(1,n_files);
Median_Artifact_Probability_of_Kept_ICs = zeros(1,n_files);
Mean_Artifact_Probability_of_Kept_ICs = zeros(1,n_files);
Range_Artifact_Probability_of_Kept_ICs = zeros(1,n_files);
Min_Artifact_Probability_of_Kept_ICs = zeros(1,n_files);
Max_Artifact_Probability_of_Kept_ICs = zeros(1,n_files);
Number_Segments_Post_Segment_Rejection = zeros(1,n_files);
Length_ok4ICA = ones(1,n_files);

% Error log
error_log = cell(1,n_files);

%% iterate the following preprocessing pipeline over all your data files:
for current_file = 125:n_files
    try
%         tic

        %% Load file
        EEGloaded = pop_loadset(FileNames{current_file}, src_folder_name);

        if task_EEG_processing == 1 || task_EEG_processing == 2
            events=EEGloaded.event;
            complete_event_info=EEGloaded.urevent;
        end

        % Rename some electrode labels for compatibility with MARA
        EEGloaded = rename_channels_Roche2HAPPE(EEGloaded);
        % save chan locs common 59 chs
        chan_IDs_EEGloaded = intersect(channels_of_interest, {EEGloaded.chanlocs.labels}); 
        EEG_chans = pop_select(EEGloaded, 'channel', chan_IDs_EEGloaded);  
        EEG_chans = eeg_checkset(EEG_chans);
        full_selected_channels = EEG_chans.chanlocs;
        clear EEG_chans chan_IDs_EEGloaded
        
        % Deal with flat channels; remove if channel is 0 for 95% of the time
        EEGloaded = clean_flatlines(EEGloaded, EEGloaded.xmax * 0.95, 20);

        % Check if the required 10-20 Electrodes are available for the file.
        assert(all(ismember({'FP1', 'FP2', 'F3', 'F4', 'F7', 'F8', 'C3', 'C4', 'T3', 'T4', ...
            'PZ', 'O1', 'O2', 'T5', 'T6', 'P3', 'P4', 'Fz'}, {EEGloaded.chanlocs.labels})), ...
            'File %i does not contain all needed 10-20 channels.', current_file);

        srate=double(EEGloaded.srate);
        
        EEGloaded.setname='rawEEG';
        EEG = eeg_checkset(EEGloaded);
        File_Length_In_Secs(current_file) = EEG.xmax;

        chan_IDs = intersect(channels_of_interest, {EEG.chanlocs.labels});  % intersect because of missing channels
        EEG = eeg_checkset(EEG);

        %% filter the data with 1hz highpass (for srate 250), bandpass 1hz-249hz (for srate 500, ICA doesn't reliably work well with frequencies above 250hz)
        if srate<500
            EEG = pop_eegfiltnew(EEG, [],1,[],1,[],0);
        elseif srate >= 500
            EEG = pop_eegfiltnew(EEG, 1,249,[],0,[],0);
        end
        EEG.setname='rawEEG_f';
        EEG = eeg_checkset(EEG);

        %% select EEG channels of interest for analyses and 10-20 channels from the list you specified at the top of the script
        EEG = pop_select(EEG, 'channel', chan_IDs);  
        EEG.setname = 'rawEEG_f_cs';
        EEG = eeg_checkset(EEG);
%         full_selected_channels = EEG.chanlocs;

        %% reduce line noise in the data (note: may not completely eliminate, re-referencing helps at the end as well)
        disp('Applying Cleanline for line noise removal.')
        EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan,'computepower',1,'linefreqs',...
            [line_noise_freq (2 * line_noise_freq)] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
            'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
        EEG.setname='rawEEG_f_cs_ln';
        EEG = eeg_checkset(EEG);

        % close window if visualizations are turned off
        if pipeline_visualizations_semiautomated == 0
            close all;
        end

        %% Segment dynamic video data into trials; remove non-videos data (breaks between trials)
        % RH note: this section is specific for the current task in the
        % LEAP study.
        
            % LEAP site
            SiteCur = EEG.euaims.site;
            % video length
            rest_video_length = 60;
            subject = current_file;
            FlashDatapath = 'xxx'; % this folder contains data to check the log of the events from Task Engine and is specific to the LEAP EEG data

        
            EEG.condition_raw_aligned = zeros(1, EEG.pnts);  % To record sample annotations for later analysis. (1 x n_samples) vector with 0 = rejected/breaks, 210 = social, 211 = non-social
            EEG_eoi.event = EEG.event(~strcmp({EEG.event.type}, 'boundary'));
                % correct events for UMCU site
                if strcmpi(SiteCur,'UMCU') %@Rianne: this was needed when I did the resting preprocessing because their data had different format, remove if no longer the case
                    error('Check event types for site UMCU!')
                end

            segment_boundaries = zeros(length(EEG_eoi.event), 2);
            condition_oi = zeros(length(EEG_eoi.event), 1);
            % define boundaries between videos    
            samples_per_video = rest_video_length * EEG.srate;
                for i = 1:length(EEG_eoi.event)
                    segment_boundaries(i,:) = [EEG_eoi.event(i).latency, (EEG_eoi.event(i).latency + samples_per_video - 1)];
                    condition_oi(i) = str2double(EEG_eoi.event(i).type);
                end
                clear i


             % Check whether the end flash is after the time end of last video
             Subj = extractBefore(FileNames{subject},'_social_nonsocial_videos.set');
                Flash_name = fullfile(FlashDatapath,sprintf('Tablesocialnonsocial - %i.mat',str2double(Subj)));
                if exist(Flash_name,'file') == 2
                    load(Flash_name);
                    VidLength = restingvideostable.TrialOffsetTime - restingvideostable.TrialOnsetTime;
                    % check if resting videos table and trl are same length
                    if isequal(size(VidLength,1),size(segment_boundaries,1))
                        FlashCheck = 1;
                        for tt = 1:size(VidLength,1)
                            if VidLength(tt,1) < rest_video_length % check if flash before 60sec after video onset
                                disp('Flash before end of video; correcting now')
                                Newtrialnsamples = VidLength(tt,1) * EEG.srate;
                                % adjust in trl
                                segment_boundaries(tt,2) = segment_boundaries(tt,1) + Newtrialnsamples-1;
                                FlashCheck = 2;
                            end
                        end
                        clear tt VidLength restingvideostable Newtrialnsamples
                    else
                        warning('Mismatch nr videos in trl and resting videos table')
                        FlashCheck = 999;
                    end
                else
                    warning('No flash information found, skipping this step')
                    FlashCheck = 0;
                end

            % Check whether any trials overlap 
            if any((segment_boundaries(2:end,1)-segment_boundaries(1:end-1,2))<=0)
                warning('trials overlap: adjusting boundaries now');
                for i = 1:length(segment_boundaries)-1
                    if segment_boundaries(i,2) > segment_boundaries(i+1,1)
                        segment_boundaries(i,2) = segment_boundaries(i+1,1) - EEG.srate; %cut previous video to 1 second before start of next video
                    end
                end
                clear i
            end   

            % adjust EEG.condition_raw_aligned 
            for i = 1:length(segment_boundaries)
                if isequal(condition_oi(i,1),10) || isequal(condition_oi(i,1),12) || isequal(condition_oi(i,1),13) || isequal(condition_oi(i,1),14)
                    %social condition in different languages
                    CondCur = 10;
                elseif isequal(condition_oi(i,1),11) %toy condition
                    CondCur = 11;
                end
                EEG.condition_raw_aligned(floor(segment_boundaries(i,1)):floor(segment_boundaries(i,2)))= CondCur;
                clear CondCur
            end
            
            % adjust EEG events for social different languages
            for i = 1:length(EEG.event)
                if strcmp(EEG.event(i).type,'12') || strcmp(EEG.event(i).type,'13') || strcmp(EEG.event(i).type,'14') 
                    %change event to 10
                    EEG.event(i).type = '10';
                end
            end

            EEG = pop_select(EEG, 'point', segment_boundaries);

            % remove boundary events to prevent unnecessary epoch rejection later on (we have ensured that introduced boundaries will coincide with epoch boundaries).
            EEG.event = EEG.event(~strcmp({EEG.event.type}, 'boundary'));

        % Save event info for later
        EEG_afterRSsegmenting = EEG;
        File_Length_afterRSsegementing_In_Secs(current_file) = EEG.xmax;
        
        
        %% crude bad channel detection using spectrum criteria and 3SDeviations as channel outlier threshold, done twice
        EEG = pop_rejchan(EEG, 'elec', 1:EEG.nbchan, 'threshold', [-3 3], 'norm', 'on', 'measure', 'spec', 'freqrange', [1 125]);
        EEG.setname='rawEEG_f_cs_ln_badc';
        EEG = eeg_checkset(EEG);

        EEG = pop_rejchan(EEG, 'elec', 1:EEG.nbchan, 'threshold', [-3 3], 'norm', 'on', 'measure', 'spec', 'freqrange', [1 125]);
        EEG.setname='rawEEG_f_cs_ln_badc2x';
        EEG = eeg_checkset(EEG);
        selected_channel_locations = EEG.chanlocs;

        %save the names of the rejected channels for output table after the pipeline finishes
        selected_channel_labels={selected_channel_locations.labels};
        bad_channels_removed = setdiff(chan_IDs, selected_channel_labels);
        [~,ROI_indices_in_selected_chanlocs] = intersect(selected_channel_labels,ROI_channels);

        % Raise warning if sample size is smaller than recommended 
        if 30 * (EEG.nbchan) ^ 2 > EEG.pnts * EEG.trials
            warning('Sample size for file %i (%i samples) is smaller than recommended for robust, stable ICA decomposition (%i samples).', ...
                current_file, EEG.pnts * EEG.trials, 30 * (EEG.nbchan) ^ 2)
            Length_ok4ICA(current_file) = 0;
        end

        %% Run wavelet-ICA (ICA first for clustering the data, then wavelet thresholding on the ICs)
        %  Uses a soft, global threshold for the wavelets, wavelet family is
        %  coiflet (level 5), threshold multiplier .75 to remove more high
        %  frequency noise. For details, see wICA.m function.

        try 
            if pipeline_visualizations_semiautomated == 0
                [wIC, A, W, IC] = wICA(EEG,'runica', 1, 0, [], 5);
            elseif pipeline_visualizations_semiautomated == 1
                [wIC, A, W, IC] = wICA(EEG,'runica', 1, 1, srate, 5);
            end
        catch wica_err
            if strcmp ('Output argument "wIC" (and maybe others) not assigned during call to "wICA".',wica_err.message)
                error('Error during wICA, most likely due to memory settings. Please confirm your EEGLAB memory settings are set according to the description in the HAPPE ReadMe')
            else
                rethrow(wica_err)
            end
        end

        %reconstruct artifact signal as channelsxsamples format from the wavelet coefficients
        artifacts = A*wIC;

        %reshape EEG signal from EEGlab format to channelsxsamples format
        EEG2D = reshape(EEG.data, size(EEG.data,1), []);

        %subtract out wavelet artifact signal from EEG signal
        wavcleanEEG = EEG2D-artifacts;

        EEG = pop_importdata('dataformat', 'array', 'nbchan', 0, 'data', wavcleanEEG, 'srate', srate, 'pnts', 0, 'xmin', 0, 'chanlocs', selected_channel_locations);
        EEG.setname = 'wavcleanedEEG';
        EEG = eeg_checkset(EEG);

        %% run ICA to evaluate components this time
        EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on');
        EEG = eeg_checkset(EEG);

        %% use MARA to flag artifactual IComponents automatically if artifact probability > .5
        [~,EEG,~] = processMARA(EEG, EEG, EEG, [0, 0, pipeline_visualizations_semiautomated,...
            pipeline_visualizations_semiautomated , pipeline_visualizations_semiautomated]);

        EEG.reject.gcompreject = zeros(size(EEG.reject.gcompreject));
        EEG.reject.gcompreject(EEG.reject.MARAinfo.posterior_artefactprob > 0.5) = 1;
        EEG.setname='wavcleanedEEG_ICA_MARA';
        EEG = eeg_checkset( EEG );

        % store MARA related variables to assess ICA/data quality
        index_ICs_kept=(EEG.reject.MARAinfo.posterior_artefactprob < 0.5);
        median_artif_prob_good_ICs = median(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        mean_artif_prob_good_ICs = mean(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        range_artif_prob_good_ICs = range(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        min_artif_prob_good_ICs = min(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        max_artif_prob_good_ICs = max(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));

        %store IC variables and calculate variance of data that will be kept after IC rejection:
        ICs_to_keep = find(EEG.reject.gcompreject == 0);
        ICA_act = EEG.icaact;
        ICA_winv = EEG.icawinv;

        %variance of wavelet-cleaned data to be kept = varianceWav:
        [projWav, varianceWav] = compvar(EEG.data, ICA_act, ICA_winv, ICs_to_keep);

        %% reject the ICs that MARA flagged as artifact
        artifact_ICs = find(EEG.reject.gcompreject == 1);
        EEG = pop_subcomp(EEG, artifact_ICs, 0);
        EEG.setname = 'wavcleanedEEG_ICA_MARA_rej';
        EEG = eeg_checkset(EEG);

        %% segment data according to data type
        EEG.xmax = EEG.pnts / EEG.srate; % computation of xmax as done in eeg_checkset by (EEG.pnts-1)/EEG.srate+EEG.xmin leads to loss of the last trial

        if segment_data
            
             % import event tags 
            EEG.event = EEG_afterRSsegmenting.event;
            EEG.urevent = EEG_afterRSsegmenting.urevent;
        
            if task_EEG_processing == 0
                EEG = eeg_regepochs(EEG, 'recurrence', segment_length, 'limits', [0 segment_length], 'rmbase', [NaN]);
            elseif task_EEG_processing == 1
                EEG = pop_epoch(EEG, task_conditions, [task_segment_start task_segment_end], 'verbose', 'no', 'epochinfo', 'yes');
            elseif task_EEG_processing == 2


                % RH note: For the intended analyses, we want to segment
                % the social and non-social videos into 1-second
                % overlapping epochs, consistent with a resting state paradigm. 
                % For consistency across pipelines, we based the following code on 
                % the MADE code for continous paradigms and adapted this to be
                % compatible with the LEAP EEG data format. 
        
                % From the MADE developers:
                % Often resting state EEG is collected in eyes close and eyes open conditions
                % This script inserts event at specific time interval for resting state EEG
                % and creates separate events for eyes close and eyes open resting EEG
                % data. See MADE paper for details.   
 

                % Create segment info for epoching
                % 1. Name of eyes close and eyes open markers
                rest_event_markers = {'10','11'}; %RH: video markers
                % 2. Name of new eyes close and eyes open markers
                new_rest_markers={'910','911'}; % enter markers for epochs
                % 3. Does the data have a trial end marker?
                trial_end_marker=0; % 0=NO (no trial end marker in data), 1=YES (data have trial end marker)
                trial_end_marker_name=('xxx'); % enter trial end marker name
                % 4. Length of epochs in seconds
                rest_epoch_length=1;
                % 5. Do you want to create overlapping epoch?
                overlap_epoch = 1; % 0 = NO (do not create overlapping epoch), 1 = YES (50% overlapping epoch)

                % Insert markers
                if overlap_epoch==1
                    time_samples = (rest_epoch_length/2)*EEG.srate; % convert time window into samples or data points
                else
                    time_samples = rest_epoch_length*EEG.srate;
                end
                
                EEG.urevent_EEGses = EEG.urevent;
                EEG.urevent_vidseg = EEG.event;
                EEG.urevent = EEG.event;

                tm=1;
                for ue=1:length(EEG.urevent)
                    for rm=1:length(rest_event_markers)
                        if strcmp(EEG.urevent(ue).type, rest_event_markers{rm})==1
                            if trial_end_marker == 1
                                for te=ue:length(EEG.urevent)
                                    if strcmp(EEG.urevent(te).type, trial_end_marker_name)==1
                                        trial_end_latency=EEG.urevent(te).latency;
                                        break;
                                    end
                                end
                            else
                                if ue < length(EEG.urevent)
                                    trial_end_latency=EEG.urevent(ue+1).latency-rest_epoch_length-time_samples;
                                else
                                    trial_end_latency=length(EEG.times)-rest_epoch_length-time_samples;
                                end
                            end
                            event_times = EEG.urevent(ue).latency;
                            while event_times <trial_end_latency
                                EEG.event(tm).type=new_rest_markers{rm};
                                EEG.event(tm).latency=event_times;
                                event_times = event_times+time_samples;
                                tm=tm+1;
                            end
                        end
                    end
                end
                
                % adjust the new EEG.events.value and duration
                for eps = 1:length(EEG.event)
                        EEG.event(eps).value = 'EEGdummy';
                        EEG.event(eps).duration = 1;
                end
                % adjust the latencies to round numbers
                for eps = 1:length(EEG.event)
                        EEG.event(eps).latency = round(EEG.event(eps).latency);
                end
                % restore and adjust the EEG.urevent for each site
                switch SiteCur
                    case {'KCL','CIMH','RUNMC'}
                        Nlasturevent = length(EEG.urevent_EEGses);
                        NewUrevents = EEG.urevent_EEGses;
                        NewEvents = rmfield(EEG.event,'urevent');
                        for eps = 1:length(EEG.event)
                            % add event into EEG.urevent
                            NewUrevents(Nlasturevent + eps) = NewEvents(eps);
                            % add in index urevent_orig into event field
                            EEG.event(eps).urevent = Nlasturevent + eps;
                        end
                        EEG.urevent = NewUrevents;
                    case 'UCBM'
                        NewEvents = rmfield(EEG.event,'urevent');
                        % create new urevent structure 
                        NewUrevents = NewEvents(1);
                        % add original event in new format
                        for old_ev = 1:length(EEG.urevent_EEGses)
                            NewUrevents(old_ev).type = num2str(EEG.urevent_EEGses(old_ev).type);
                            NewUrevents(old_ev).value = 'OrigStimMarker';
                            NewUrevents(old_ev).latency = EEG.urevent_EEGses(old_ev).latency;
                            NewUrevents(old_ev).duration = 1;
                        end
                        % add new dummies to urevent structure
                        Nlasturevent = length(EEG.urevent_EEGses);
                        for eps = 1:length(EEG.event)
                            % add event into EEG.urevent
                            NewUrevents(Nlasturevent + eps) = NewEvents(eps);
                            % add in index urevent_orig into event field
                            EEG.event(eps).urevent = Nlasturevent + eps;
                        end
                        EEG.urevent = NewUrevents; 
                    otherwise
                        error('Unrecognised LEAP site')
                end
        
                % create epoch
                EEG = eeg_checkset(EEG);
                EEG = pop_epoch(EEG, new_rest_markers, [0 rest_epoch_length], 'epochinfo', 'yes');
                EEG = eeg_checkset(EEG);
        
                EEG_afterEpoching = EEG;
        
            end
        end
        
        
        

        %% if selected option, interpolate bad data within segments from "good channels" only:
        if segment_interpolation

            % use only the good channels to evaluate data:
            eeg_chans = [1:length(selected_channel_locations)];

            % evaluate the channels for each segment and interpolate channels
            % with bad data for that segment using the FASTER program,
            % interpolating channels scoring above/below z threshold of 3 for
            % a segment:  
            ext_chans = [];
            o = struct();
            o.epoch_interp_options.rejection_options.measure = [1 1 1 1];
            o.epoch_interp_options.rejection_options.z = [3 3 3 3];

            if  length(size(EEG.data)) > 2
                status = '';
                lengths_ep = cell(1,size(EEG.data,3));
                for v=1:size(EEG.data,3)
                    list_properties = single_epoch_channel_properties(EEG,v,eeg_chans);
                    lengths_ep{v} = eeg_chans(logical(min_z(list_properties,o.epoch_interp_options.rejection_options))); % indices of channels to be interpolated in epoch v
                    status = [status sprintf('%d: ',v) sprintf('%d ',lengths_ep{v}) newline];
                end
                EEG=h_epoch_interp_spl(EEG,lengths_ep,ext_chans);
                EEG.saved = 'no';

                %add the info about which channels were interpolated for each segment to the EEG file
                EEG.etc.epoch_interp_info = [status];
            end
        end

        %% rejection of bad segments using amplitude-based and joint probability artifact detection
        if ROI_channels_only == 0
            EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan] ,[reject_min_amp],[reject_max_amp],[EEG.xmin],[EEG.xmax],2,0);
            EEG = pop_jointprob(EEG,1,[1:EEG.nbchan],3,3,pipeline_visualizations_semiautomated,...
                0,pipeline_visualizations_semiautomated,[],pipeline_visualizations_semiautomated);
        else
            EEG = pop_eegthresh(EEG,1,[ROI_indices_in_selected_chanlocs]',[reject_min_amp],[reject_max_amp],[EEG.xmin],[EEG.xmax],2,0);
            EEG = pop_jointprob(EEG,1,[ROI_indices_in_selected_chanlocs]',3,3,pipeline_visualizations_semiautomated,...
                0,pipeline_visualizations_semiautomated,[],pipeline_visualizations_semiautomated);
        end

        EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal] ,0); 
        EEG = eeg_checkset(EEG);

        %% interpolate the channels that were flagged as bad earlier:
%         EEG = pop_interp(EEG, full_channel_locs, 'spherical');
        EEG = pop_interp(EEG, full_selected_channels, 'spherical');
        EEG.setname='wavcleanedEEG_ICA_MARA_rej_chan_int';
        EEG = eeg_checkset(EEG);

        %% re-reference the data: average reference used here
        if average_rereference == 1
            EEG = pop_reref(EEG, []);
            EEG.setname='wavcleanedEEG_ICA_MARA_rej_chan_int_avgreref';
            EEG = eeg_checkset(EEG);
        else
            [~,ref_chan_indices_in_full_selected_chanlocs] = intersect({full_selected_channels.labels},NO_AVERAGE_REREF_channel_subset);
            EEG = pop_reref(EEG, ref_chan_indices_in_full_selected_chanlocs);
            EEG.setname='wavcleanedEEG_ICA_MARA_rej_chan_int_chansubsetreref';
            EEG = eeg_checkset(EEG);
        end

        %% store outputs and report metrics
        Number_Channels_User_Selected(current_file) = size(chan_IDs, 2);
        Number_ICs_Rejected(current_file) = length(artifact_ICs);
        Number_Good_Channels_Selected(current_file) = size(selected_channel_locations,2);
        Percent_Good_Channels_Selected(current_file) = Number_Good_Channels_Selected(current_file)/Number_Channels_User_Selected(current_file)* 100;
        Percent_ICs_Rejected(current_file) = Number_ICs_Rejected(current_file)/Number_Good_Channels_Selected(current_file)* 100;
        Percent_Variance_Kept_of_Post_Waveleted_Data(current_file) = varianceWav;
        if isempty(bad_channels_removed)
            Interpolated_Channel_IDs{current_file} = 'none';
        else
            Interpolated_Channel_IDs{current_file} = [sprintf('%s ',bad_channels_removed{1:end-1}),bad_channels_removed{end}];
        end
        Median_Artifact_Probability_of_Kept_ICs(current_file) = median_artif_prob_good_ICs;
        Mean_Artifact_Probability_of_Kept_ICs(current_file) = mean_artif_prob_good_ICs;
        Range_Artifact_Probability_of_Kept_ICs(current_file) = range_artif_prob_good_ICs;
        Min_Artifact_Probability_of_Kept_ICs(current_file) = min_artif_prob_good_ICs;
        Max_Artifact_Probability_of_Kept_ICs(current_file) = max_artif_prob_good_ICs;
        Number_Segments_Post_Segment_Rejection(current_file) = EEG.trials;


        %% save preprocessed dataset with subject ID as either txt file (user specified) or eeglab .set file 
        EEG = rename_channels_HAPPE2Roche(EEG);  % Rename channels such that they match the Roche format

        switch save_as_format
            case 1 % txt file
                pop_export(EEG,[out_data_folder_name filesep strrep(erase(FileNames{current_file}, 'Unfiltered-'), src_file_ext,'_processed.txt')],'transpose','on','precision',8);
            case 2 % .mat file
                save([out_data_folder_name filesep strrep(erase(FileNames{current_file}, 'Unfiltered-'), src_file_ext,'_processed.mat')], 'EEG');
            case 3 % .set file
                EEG = pop_saveset(EEG, 'filename', strrep(erase(FileNames{current_file}, 'Unfiltered-'), src_file_ext, '_processed.set'), 'filepath', out_data_folder_name);
        end

        %% generate power spectrum and topoplot visualization if user requested:
        if pipeline_visualizations_semiautomated == 1
            %%
            EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal] ,0); % remove bad segments for computation of powerspectrum
            EEG = eeg_checkset(EEG );
            figure; pop_spectopo(EEG, 1, [], 'EEG' , 'freq', [[freq_to_plot]], 'freqrange',[[vis_freq_min] [vis_freq_max]],'electrodes','off');
            saveas (gcf, [out_data_folder_name filesep strrep(erase(FileNames{current_file}, 'Unfiltered-'), src_file_ext,'_processedspectrum.jpg')]);
        end

       %% Display progress
       fprintf('File %i / %i processed.\n', current_file, n_files);
%        toc
       
    catch ME
        warning('File %s was skipped because an error occured during processing: %s', FileNames{current_file}, ME.message);
        error_log{current_file} = {FileNames{current_file}, ME};
    end
    
    %% HAPPE output: generate output table in the "preprocessed" subfolder listing the subject file name and relevant variables for assesssing how good/bad that datafile was and how well the pipeline worked
    cd(out_folder_name)
    if exist('HAPPE_preprocessing_report.mat','file') == 2
        load HAPPE_preprocessing_report.mat
    end
    outputtable=table(FileNames',File_Length_In_Secs',File_Length_afterRSsegementing_In_Secs', Number_Channels_User_Selected',Number_Segments_Post_Segment_Rejection',...
        Number_Good_Channels_Selected', Percent_Good_Channels_Selected', Interpolated_Channel_IDs',Number_ICs_Rejected',...
        Percent_ICs_Rejected', Percent_Variance_Kept_of_Post_Waveleted_Data',Median_Artifact_Probability_of_Kept_ICs',...
        Mean_Artifact_Probability_of_Kept_ICs',Range_Artifact_Probability_of_Kept_ICs',Min_Artifact_Probability_of_Kept_ICs',...
        Max_Artifact_Probability_of_Kept_ICs', Length_ok4ICA');
    outputtable.Properties.VariableNames ={'FileNames','File_Length_In_Secs','File_Length_afterRSsegementing_In_Secs','Number_Channels_User_Selected','Number_Segments_Post_Segment_Rejection',...
        'Number_Good_Channels_Selected', 'Percent_Good_Channels_Selected', 'Interpolated_Channel_IDs','Number_ICs_Rejected',...
        'Percent_ICs_Rejected', 'Percent_Variance_Kept_of_Post_Waveleted_Data','Median_Artifact_Probability_of_Kept_ICs',...
        'Mean_Artifact_Probability_of_Kept_ICs','Range_Artifact_Probability_of_Kept_ICs','Min_Artifact_Probability_of_Kept_ICs',...
        'Max_Artifact_Probability_of_Kept_ICs','Length_ok4ICA'};
    
    HAPPE_report_table(current_file,:) = outputtable(current_file,:);
    
    save('HAPPE_preprocessing_report.mat','HAPPE_report_table');

    if ~isempty(error_log{current_file})
        disp('Errors have occured during the processing of files. Saving error log with file names and corresponding MException objects to output directory.')
        if exist('HAPPE_preprocessing_error_log.mat','file')
                load HAPPE_preprocessing_error_log.mat
        end
        HAPPEpreproc_error_log{current_file} = error_log{current_file};
        save('HAPPEpreproc_error_log.mat','HAPPEpreproc_error_log')
    end
    
end
