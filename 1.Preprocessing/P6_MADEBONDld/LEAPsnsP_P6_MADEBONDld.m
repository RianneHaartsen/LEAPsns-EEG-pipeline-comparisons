%% Analysis of social and non-social videos in LEAP: Pipeline 6 - BOND-MADEld
% using the BOND-MADE pipeline which includes adaptations to the original 
% MADE pipeline on a lower-density layout: 
% 
% Adjustments made relative to original MADE pipeline (P2):
% - adjustment for flat channels
% - tracking n channels for each stage of preprocessing
% Adjustments made relative to MADE-BOND pipeline (P3):
% - channels of interest set to 19 channels present in other low-density
% studies

% and adaptations compatible with the LEAP EEG
% data formats. 

% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')
% additional comments by Rianne Haartsen are noted with '%RH note;'

% Adapted by Rianne Haartsen, PhD.; 03-2024 
% Birkbeck College, University of London

% This adapted pipeline is released under the GNU General Public License version 3.



% ************************************************************************
% The Maryland Analysis of Developmental EEG (UMADE) Pipeline
% Version 1.0
% Developed at the Child Development Lab, University of Maryland, College Park

% Contributors to MADE pipeline:
% Ranjan Debnath (rdebnath@umd.edu)
% George A. Buzzell (gbuzzell@umd.edu)
% Santiago Morales Pamplona (moraless@umd.edu)
% Stephanie Leach (sleach12@umd.edu)
% Maureen Elizabeth Bowers (mbowers1@umd.edu)
% Nathan A. Fox (fox@umd.edu)

% MADE uses EEGLAB toolbox and some of its plugins. Before running the pipeline, you have to install the following:
% EEGLab:  https://sccn.ucsd.edu/eeglab/downloadtoolbox.php/download.php

% You also need to download the following plugins/extensions from here: https://sccn.ucsd.edu/wiki/EEGLAB_Extensions

% Specifically, download:
% MFFMatlabIO: https://github.com/arnodelorme/mffmatlabio/blob/master/README.txt
% FASTER: https://sourceforge.net/projects/faster/
% ADJUST: https://www.nitrc.org/projects/adjust/
% Adjusted ADJUST (included in this pipeline):  https://github.com/ChildDevLab/MADE-EEG-preprocessing-pipeline

% After downloading these plugins (as zip files), you need to place it in the eeglab/plugins folder.
% For instance, for FASTER, you uncompress the downloaded extension file (e.g., 'FASTER.zip') and place it in the main EEGLAB "plugins" sub-directory/sub-folder.
% After placing all the required plugins, add the EEGLAB folder to your path by using the following code:

% addpath(genpath(('...')) % Enter the path of the EEGLAB folder in this line

% Please cite the following references for in any manuscripts produced utilizing MADE pipeline:

% EEGLAB: A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
% analysis of single-trial EEG dynamics. Journal of Neuroscience Methods, 134, 9?21.

% firfilt (filter plugin): developed by Andreas Widmann (https://home.uni-leipzig.de/biocog/content/de/mitarbeiter/widmann/eeglab-plugins/)

% FASTER: Nolan, H., Whelan, R., Reilly, R.B., 2010. FASTER: Fully Automated Statistical
% Thresholding for EEG artifact Rejection. Journal of Neuroscience Methods, 192, 152?162.

% ADJUST: Mognon, A., Jovicich, J., Bruzzone, L., Buiatti, M., 2011. ADJUST: An automatic EEG
% artifact detector based on the joint use of spatial and temporal features. Psychophysiology, 48, 229?240.
% Our group has modified ADJUST plugin to improve selection of ICA components containing artifacts

% This pipeline is released under the GNU General Public License version 3.

% ************************************************************************

%% User input: user provide relevant information to be used for data processing
% Preprocessing of EEG data involves using some common parameters for
% every subject. This part of the script initializes the common parameters.

clear % clear matlab workspace
clc % clear matlab command window

addpath(' xxx ');% enter the path of the EEGLAB folder in this line
eeglab nogui; % call eeglab to set up the plugin
addpath(genpath(' xxxx path for folder with scripts '))

% 1. Enter the path of the folder that has the raw data to be analyzed
rawdata_location = 'xxx/RawEEG_data/';
FlashDatapath = 'xxx'; % RH note; this variable is specific to the LEAP social/non-social videos dataset only

% 2. Enter the path of the folder where you want to save the processed data
output_location = 'xxx/Preproc_MADEBONDld';

% 3. Enter the path of the channel location file
channel_locations = ['xxx' filesep 'euaimschan.mat'];

% 4. Do your data need correction for anti-aliasing filter and/or task related time offset?
adjust_time_offset = 0; % 0 = NO (no correction), 1 = YES (correct time offset)

% If your data need correction for time offset, initialize the offset time (in milliseconds)
% YOU NEED TO CHANGE THE "xx" TO A REAL NUMBER WITHOUT QUOTATION MARKS!
filter_timeoffset = 0;     % anti-aliasing time offset (in milliseconds). 0 = No time offset
stimulus_timeoffset   = 0; % stimulus related time offset (in milliseconds). 0 = No time offset
response_timeoffset = 0;    % response related time offset (in milliseconds). 0 = No time offset
stimulus_markers = {'xxx', 'xxx'};      % enter the stimulus makers that need to be adjusted for time offset
respose_markers = {'xxx', 'xxx'};       % enter the response makers that need to be adjusted for time offset

% 5. Do you want to down sample the data?
down_sample = 0; % 0 = NO (no down sampling), 1 = YES (down sampling)
sampling_rate = 1000; % set sampling rate (in Hz), if you want to down sample

% 6. Do you want to delete the outer layer of the channels? (Rationale has been described in MADE manuscript)
% RH: for LEAP, we select the channels common across all sites immediately
% after reading in the data
    %    This function can also be used to down sample electrodes. For example, if EEG was recorded with 128 channels but you would
    %    like to analyse only 64 channels, you can assign the list of channnels to be excluded in the 'outerlayer_channel' variable.    
    delete_outerlayer = 0; % 0 = NO (do not delete outer layer), 1 = YES (delete outerlayer);
    % If you want to delete outer layer, make a list of channels to be deleted
    outerlayer_channel = {'list of channels'}; % list of channels
    % recommended list for EGI 128 chanenl net: {'E17' 'E38' 'E43' 'E44' 'E48' 'E49' 'E113' 'E114' 'E119' 'E120' 'E121' 'E125' 'E126' 'E127' 'E128' 'E56' 'E63' 'E68' 'E73' 'E81' 'E88' 'E94' 'E99' 'E107'}
    % RH: for LEAP, some channels were not included in some sites. We will
    % select the common electrodes across all sites only to keep pre-processing
    % consistent. 

% 7. Initialize the filters
highpass = 1; % High-pass frequency
lowpass  = 32; % Low-pass frequency. We recommend low-pass filter at/below line noise frequency (see manuscript for detail)

% 8. Are you processing task-related or resting-state EEG data?
task_eeg = 2; %RH:  0 = resting, 1 = task, 2 = resting state with 2 conditions
task_event_markers = {'10', '11'}; % enter all the event/condition markers

% 9. Do you want to epoch/segment your data?
% RH Note: this information needs to be given at step 12. See later

% epoch_data = 1; % 0 = NO (do not epoch), 1 = YES (epoch data)
% task_epoch_length = [0 0]; % epoch length in second
% rest_video_length = 60; %RH: video length in seconds
% rest_epoch_length = 1; % for resting EEG continuous data will be segmented into consecutive epochs of a specified length (here 1 second) by adding dummy events
% overlap_epoch = 1;     % 0 = NO (do not create overlapping epoch), 1 = YES (50% overlapping epoch)
% dummy_events ={'910','911'}; % enter dummy events name

% 10. Do you want to remove/correct baseline?
remove_baseline = 1; % 0 = NO (no baseline correction), 1 = YES (baseline correction)
baseline_window = []; % baseline period in milliseconds (MS) [] = entire epoch

% 11. Do you want to remove artifact laden epoch based on voltage threshold?
voltthres_rejection = 1; % 0 = NO, 1 = YES
volt_threshold = [-100 100]; % lower and upper threshold (in uV)

% 12. Do you want to perform epoch level channel interpolation for artifact laden epoch? (see manuscript for detail)
interp_epoch = 1; % 0 = NO, 1 = YES.
frontal_channels = {'AF8', 'AF7'}; % If you set interp_epoch = 1, enter the list of frontal channels to check (see manuscript for detail)
% recommended list for EGI 128 channel net: {'E1', 'E8', 'E14', 'E21', 'E25', 'E32', 'E17'}
% equivalent in LEAP EEG layout: {'F8', 'AF8', ' Fp2', 'Fp1', 'AF7', 'F7'}, but for the low-density layout, only AF7/8 are available 

%13. Do you want to interpolate the bad channels that were removed from data?
interp_channels = 1; % 0 = NO (Do not interpolate), 1 = YES (interpolate missing channels)

% 14. Do you want to rereference your data?
rerefer_data = 1; % 0 = NO, 1 = YES
reref=[]; % Enter electrode name/s or number/s to be used for rereferencing
% For channel name/s enter, reref = {'channel_name', 'channel_name'};
% For channel number/s enter, reref = [channel_number, channel_number];
% For average rereference enter, reref = []; default is average rereference

% 15. Do you want to save interim results?
save_interim_result = 1; % 0 = NO (Do not save) 1 = YES (save interim results)

% 16. How do you want to save your data? .set or .mat
output_format = 1; % 1 = .set (EEGLAB data structure), 2 = .mat (Matlab data structure)

% ********* no need to edit beyond this point for EGI .mff data **********
% ********* for non-.mff data format edit data import function ***********
% ********* below using relevant data import plugin from EEGLAB **********

%% Read files to analyses
datafile_names=dir(rawdata_location);
datafile_names=datafile_names(~ismember({datafile_names.name},{'.', '..', '.DS_Store'}));
datafile_names=datafile_names(~contains({datafile_names.name},{'.fdt'}));
datafile_names={datafile_names.name};
[filepath,name,ext] = fileparts(char(datafile_names{1}));

%% Check whether EEGLAB and all necessary plugins are in Matlab path.
if exist('eeglab','file')==0
    error(['Please make sure EEGLAB is on your Matlab path. Please see EEGLAB' ...
        'wiki page for download and instalation instructions']);
end

if strcmp(ext, '.mff')==1
    if exist('mff_import', 'file')==0
        error(['Please make sure "mffmatlabio" plugin is in EEGLAB plugin folder and on Matlab path.' ...
            ' Please see EEGLAB wiki page for download and instalation instructions of plugins.' ...
            ' If you are not analysing EGI .mff data, edit the data import function below.']);
    end
else
    warning('Your data are not EGI .mff files. Make sure you edit data import function before using this script');
end

if exist('pop_firws', 'file')==0
    error(['Please make sure  "firfilt" plugin is in EEGLAB plugin folder and on Matlab path.' ...
        ' Please see EEGLAB wiki page for download and instalation instructions of plugins.']);
end

if exist('channel_properties', 'file')==0
    error(['Please make sure "FASTER" plugin is in EEGLAB plugin folder and on Matlab path.' ...
        ' Please see EEGLAB wiki page for download and instalation instructions of plugins.']);
end

if exist('ADJUST', 'file')==0
    error(['Please make sure you download modified "ADJUST" plugin from GitHub (link is in MADE manuscript)' ...
        ' and ADJUST is in EEGLAB plugin folder and on Matlab path.']);
end

%% Create output folders to save data
if save_interim_result ==1
    if exist([output_location filesep 'filtered_data'], 'dir') == 0
        mkdir([output_location filesep 'filtered_data'])
    end
    if exist([output_location filesep 'ica_data'], 'dir') == 0
        mkdir([output_location filesep 'ica_data'])
    end
end
if exist([output_location filesep 'processed_data'], 'dir') == 0
    mkdir([output_location filesep 'processed_data'])
end

%% Initialize output variables
reference_used_for_faster=cell(1, length(datafile_names)); % reference channel used for running faster to identify bad channel/s
% bad channel identification
n_all_channels=cell(1, length(datafile_names)); %number of channels from the start
n_flat_bad_channels=cell(1, length(datafile_names)); %number of bad channel/s identified as flat
flat_bad_channels=cell(1, length(datafile_names)); %labels of bad channel/s identified as flat
n_faster_bad_channels=cell(1, length(datafile_names)); %number of bad channel/s identified by faster
faster_bad_channels=cell(1, length(datafile_names)); %labels of bad channel/s identified by faster
% ica
ica_prep_input_n_channels=cell(1, length(datafile_names)); %number channels going into ica prep
ica_input_n_channels=cell(1, length(datafile_names)); %number channels going into ica
n_ica_preparation_bad_channels=cell(1, length(datafile_names)); %number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
ica_preparation_bad_channels=cell(1, length(datafile_names)); %labels of bad channel/s due to channel/s exceeding xx% of artifacted epochs
length_ica_data=[]; % length of data (in second) fed into ICA decomposition
total_ICs=[]; % total independent components (ICs)
n_ICs_removed=cell(1, length(datafile_names)); % number of artifacted ICs
ICs_removed=cell(1, length(datafile_names)); % artifacted ICs index
% epochs rejection
total_epochs_before_artifact_rejection=[];
total_epochs_after_artifact_rejection=[];
% final report
n_total_channels_interpolated=[]; % total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels > for the whole EEG task
total_channels_interpreted=cell(1, length(datafile_names)); %labels of bad channel/s interpolated in total
datafile_name_preproc_done={}; % name of preprocessed data
preproc_path = {};
date_preproc = {};
ExpComments = {};
FlashChecks = [];
error_log = cell(1,length(datafile_names));


%% Loop over all data files
for subject = 1:length(datafile_names)
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});

    try
        EEG=[];

        %% STEP 1: Import EGI data file and relevant information
        EEG = pop_loadset([rawdata_location filesep datafile_names{subject}]);
        EEG = eeg_checkset(EEG);

        % select common channels only
        % RH note; here we only select the 19 channels that correspond to
        % the low-density layout in a study within a clinical context
        Arbaclofen_19channels = {'FCz', 'AF7', 'AF8', 'F3', 'Fz', 'F4', 'T7', 'C3', 'Cz', 'C4', 'T8', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO7', 'PO8', 'Oz'};
        EEG = pop_select(EEG, 'nochannel', find(~ismember(lower({EEG.chanlocs.labels}),lower(Arbaclofen_19channels)))); % keep only channels that correspond to the Arbaclofen channels 
        EEG = eeg_checkset(EEG);

        % LEAP site
        SiteCur = EEG.euaims.site;
        % comments
        ExpComments = EEG.euaims.comments;
        
        % Edit this data import function and use appropriate plugin from EEGLAB
        % for non-.mff data. For example, to import biosemi data, use biosig plugin.
        % The example codes for 64 channels biosemi data:
    %     EEG = pop_biosig([rawdata_location, filesep, datafile_names{subject}]);
    %     EEG = eeg_checkset(EEG);
    %     EEG = pop_select( EEG,'nochannel', 65:72); % delete redundant channels

        %% STEP 1.5: Delete discontinuous data from the raw data file (OPTIONAL, but necessary for most EGI files)
        % Note: code below may need to be modified to select the appropriate markers (depends on the import function)
        % remove discontinous data at the start of the file
    %    disconMarkers = find(strcmp({EEG.event.type}, 'boundary')); % boundary markers often indicate discontinuity
    %    EEG = eeg_eegrej( EEG, [1 EEG.event(disconMarkers(1)).latency] ); % remove discontinuous chunk... if not EGI, MODIFY BEFORE USING THIS SECTION
    %    EEG = eeg_checkset( EEG );
        % remove data after last trsp (OPTIONAL for EGI files... useful when file has noisy data at the end)
    %    trsp_flags = find(strcmp({EEG.event.type},'TRSP')); % find indices of TRSP flags
    %    EEG = eeg_eegrej( EEG, [(EEG.event(trsp_flags(end)).latency+(1.5*EEG.srate)) EEG.pnts] ); % remove everything 1.5 seconds after the last TRSP
    %    EEG = eeg_checkset( EEG );

        %% STEP 2: Import channel locations
        % RH; LEAP comes with own channel locations in data, correct to the
        % euaimschan.mat info here
        load xxx/euaimschan.mat
        for ch_row = 1:size(EEG.chanlocs,2)
            Ch_curr = EEG.chanlocs(ch_row).labels;
            Index = zeros(size(euaimschan.chanlocs,2),1);
            for cc = 1:length(Index)
                if strcmp(Ch_curr, euaimschan.chanlocs(cc).labels)
                    Index(cc) = 1;
                end
            end
            clear cc
            EEG.chanlocs(ch_row).theta = euaimschan.chanlocs(Index == 1).theta;
            EEG.chanlocs(ch_row).radius = euaimschan.chanlocs(Index == 1).radius;
            EEG.chanlocs(ch_row).X = euaimschan.chanlocs(Index == 1).X;
            EEG.chanlocs(ch_row).Y = euaimschan.chanlocs(Index == 1).Y;
            EEG.chanlocs(ch_row).Z = euaimschan.chanlocs(Index == 1).Z;
            EEG.chanlocs(ch_row).sph_theta = euaimschan.chanlocs(Index == 1).sph_theta;
            EEG.chanlocs(ch_row).sph_radius = euaimschan.chanlocs(Index == 1).sph_radius;
        end
        clear ch_row

        % Check whether the channel locations were properly imported. The EEG signals and channel numbers should be same.
        if size(EEG.data, 1) ~= length(EEG.chanlocs)
            error('The size of the data does not match with channel numbers.');
        end

        %% STEP 2.5: Label the task (OPTIONAL)
        % insert the call to a labeling script below (script will need to be on your path)
    %    edit_event_markers_example() % an example of a labeling script from the appendix folder that can be called

        %% STEP 3: Adjust anti-aliasing and task related time offset
        if adjust_time_offset==1
            % adjust anti-aliasing filter time offset
            if filter_timeoffset~=0
                for aafto=1:length(EEG.event)
                    EEG.event(aafto).latency=EEG.event(aafto).latency+(filter_timeoffset/1000)*EEG.srate;
                end
            end
            % adjust stimulus time offset
            if stimulus_timeoffset~=0
                for sto=1:length(EEG.event)
                    for sm=1:length(stimulus_markers)
                        if strcmp(EEG.event(sto).type, stimulus_markers{sm})
                            EEG.event(sto).latency=EEG.event(sto).latency+(stimulus_timeoffset/1000)*EEG.srate;
                        end
                    end
                end
            end
            % adjust response time offset
            if response_timeoffset~=0
                for rto=1:length(EEG.event)
                    for rm=1:length(response_markers)
                        if strcmp(EEG.event(rto).type, response_markers{rm})
                            EEG.event(rto).latency=EEG.event(rto).latency-(response_timeoffset/1000)*EEG.srate;
                        end
                    end
                end
            end
        end

        %% STEP 4: Change sampling rate
        % RH note: resample to 1000Hz if srate is not 1000Hz to make sure all are
        % preprocessed with same sampling rate
        if ~isequal(EEG.srate, 1000)
            EEG = pop_resample( EEG, 1000);
            EEG = eeg_checkset( EEG );
        end

        %% STEP 5: Delete outer layer of channels
        chans_labels=cell(1,EEG.nbchan);
        for i=1:EEG.nbchan
            chans_labels{i}= EEG.chanlocs(i).labels;
        end
        [chans,chansidx] = ismember(outerlayer_channel, chans_labels);
        outerlayer_channel_idx = chansidx(chansidx ~= 0);
        if delete_outerlayer==1
            if isempty(outerlayer_channel_idx)==1
                error(['None of the outer layer channels present in channel locations of data.'...
                    ' Make sure outer layer channels are present in channel labels of data (EEG.chanlocs.labels).']);
            else
                EEG = pop_select( EEG,'nochannel', outerlayer_channel_idx);
                EEG = eeg_checkset( EEG );
            end
        end

        %% STEP 6: Filter data
        % Calculate filter order using the formula: m = dF / (df / fs), where m = filter order,
        % df = transition band width, dF = normalized transition width, fs = sampling rate
        % dF is specific for the window type. Hamming window dF = 3.3

        high_transband = highpass; % high pass transition band
        low_transband = 10; % low pass transition band

        hp_fl_order = 3.3 / (high_transband / EEG.srate);
        lp_fl_order = 3.3 / (low_transband / EEG.srate);

        % Round filter order to next higher even integer. Filter order is always even integer.
        if mod(floor(hp_fl_order),2) == 0
            hp_fl_order=floor(hp_fl_order);
        elseif mod(floor(hp_fl_order),2) == 1
            hp_fl_order=floor(hp_fl_order)+1;
        end

        if mod(floor(lp_fl_order),2) == 0
            lp_fl_order=floor(lp_fl_order)+2;
        elseif mod(floor(lp_fl_order),2) == 1
            lp_fl_order=floor(lp_fl_order)+1;
        end

        % Calculate cutoff frequency
        high_cutoff = highpass/2;
        low_cutoff = lowpass + (low_transband/2);

        % Performing high pass filtering
        EEG = eeg_checkset( EEG );
        EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
        EEG = eeg_checkset( EEG );

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        % pop_firws() - filter window type hamming ('wtype', 'hamming')
        % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)

        % Performing low pass filtering
        EEG = eeg_checkset( EEG );
        EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
        EEG = eeg_checkset( EEG );

        % pop_firws() - transition band width: 10 Hz
        % pop_firws() - filter window type hamming ('wtype', 'hamming')
        % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)

        %% Step 6 1/2: Remove non-SNSvid data in LEAP recordings and check the flash recordings
        % RH note; This step was added to ensure the timing in the EEG is
        % consistent with the timing externally recorded by the Task Engine tracker. 

            EEG.condition_raw_aligned = zeros(1, EEG.pnts);  % To record sample annotations for later analysis. (1 x n_samples) vector with 0 = rejected/breaks, 210 = social, 211 = non-social
            EEG_eoi.event = EEG.event(~strcmp({EEG.event.type}, 'boundary'));
                % correct events for UMCU site
                if strcmpi(SiteCur,'UMCU') %@Rianne: this was needed when I did the resting preprocessing because their data had different format, remove if no longer the case
                    % take out non-events
                    idbad = ismember({EEG_eoi.event.value},{'Epoch';'CM_in_range'});
                    EEG_eoi.event = EEG_eoi.event(~idbad);
                    % take out events of interest only
                    idgood = ismember({EEG_eoi.event.type},{'11';'12'});
                    EEG_eoi.event = EEG_eoi.event(idgood);
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
             Subj = extractBefore(datafile_names{subject},'_social_nonsocial_videos.set');
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
                elseif strcmpi(SiteCur,'UMCU') %compare with the video offset event
                    % find latency for offset
                    FlashCheck = 1;
                    idoffset = ismember({EEG.event.type},{'15'});
                    Latency_offsetvid = EEG.event(idoffset).latency;
                    if Latency_offsetvid < segment_boundaries(end,2)
                        disp('Flash before end of video; correcting now')
                        Newtrialnsamples = VidLength(tt,1) * EEG.srate;
                        % adjust in trl
                        segment_boundaries(tt,2) = segment_boundaries(tt,1) + Newtrialnsamples-1;
                        FlashCheck = 2;
                    end   
                else
                    warning('No flash information found, skipping this step')
                    FlashCheck = 0;
                end
                 FlashChecks(subject) = FlashCheck; clear FlashCheck

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

            % clean up
            clear condition_oi EEG_eoi

        %% STEP 7: Find flat channels and run faster to find bad channels
        
        % First check whether reference channel (i.e. zeroed channels) is present in data
        % reference channel is needed to run faster
        %RH: in LEAP the REF channel is the same for all datasets
        ref_chan=[]; FlatbadChans=[]; all_chan_bad_Flat=0; FASTbadChans=[]; all_chan_bad_FAST=0;
        ref_chan=find(strcmp({EEG.chanlocs.labels}, 'FCz')); % find FCz channel index
        n_all_channels{subject}= EEG.nbchan; 
        %%% Step 7.1: Find flat channels and remove from data %%%%%%%%%%%%%
        FlatbadIdx = all(abs(EEG.data) < .0001,2);
        FlatbadChans=find(FlatbadIdx==1);
        FlatbadChans=FlatbadChans(FlatbadChans~=ref_chan);
        EEG = eeg_checkset(EEG);
        channels_analysed=EEG.chanlocs;
        N_channels_analysed = size(channels_analysed,2);
        
        % Keep track of rejected channels
        % flat
        EEG.preproc.flat_channels.number = length(FlatbadChans);
        EEG.preproc.flat_channels.labels = {EEG.chanlocs(FlatbadChans).labels};
        % all
        EEG.preproc.rejected_channels.number = length(FlatbadChans);
        EEG.preproc.rejected_channels.labels = {EEG.chanlocs(FlatbadChans).labels};
        
        % If all or more than 10% of channels are identified as flat channels, save the dataset
        % at this stage and ignore the remaining of the preprocessing.
        if numel(FlatbadChans)==EEG.nbchan || numel(FlatbadChans)+1==EEG.nbchan ...
                || numel(FlatbadChans) >= round((N_channels_analysed-1)/100*10)
            all_chan_bad_Flat=1;
            warning(['No usable data for datafile ', datafile_names{subject}, ': too many flat channels']);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_flat_channels'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_flat_channels.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_flat_channels.set');
                DataFileLocation = [output_location filesep 'processed_data'];
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_flat_channels.mat')], 'EEG'); % save .mat format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_flat_channels.mat');
                DataFileLocation = [output_location filesep 'processed_data'];
            end
        else
            % Reject channels that are flat
            EEG = pop_select( EEG,'nochannel', FlatbadChans);
            EEG = eeg_checkset(EEG);
        end
        
        %keep track for report table
        n_flat_bad_channels{subject}= EEG.preproc.flat_channels.number;
        if numel(FlatbadChans)==0
            flat_bad_channels{subject}='0';
        else
            flat_bad_channels{subject}=EEG.preproc.flat_channels.labels;
        end

        if all_chan_bad_Flat==1 % set rest of report variables to 0
            % faster bad channel identification
            reference_used_for_faster{subject}={'NaN'};
            n_faster_bad_channels{subject}=0;
            faster_bad_channels{subject}='0';
            % ica
            ica_prep_input_n_channels{subject}=0;
            ica_input_n_channels{subject}=0;
            n_ica_preparation_bad_channels{subject}=0;
            ica_preparation_bad_channels{subject}='0';
            ica_input_n_channels{subject}=0;
            length_ica_data(subject)=0;
            total_ICs(subject)=0;
            n_ICs_removed{subject}=0;
            ICs_removed{subject}='0';
            % epoch rejection
            total_epochs_before_artifact_rejection(subject)=0;
            total_epochs_after_artifact_rejection(subject)=0;
            %final report
            n_total_channels_interpolated(subject)=0;
            total_channels_interpreted{subject} = '0';
            % info
            datafile_name_preproc_done{subject} = DataFileName;
            preproc_path{subject} = DataFileLocation;
            date_preproc{subject} = datetime('now');
            
            % add to current subject to report table 
            cd('/Users/riannehaartsen/Documents/02a_LEAP_EEG/LEAP_SocNSocVids/01_EEGdata/Preproc_MADE')
            ID = extractBefore(datafile_names{subject},'_social');
            % report table
            if exist('MADE_preprocessing_report.mat','file')
                load MADE_preprocessing_report.mat
            end
            report_table_newrow = table({ID}, datafile_name_preproc_done(1,subject), preproc_path(1,subject), date_preproc(1,subject), ...
                {SiteCur}, {ExpComments}, {FlashChecks(1,subject)}, ...
                n_all_channels(1,subject), n_flat_bad_channels(1,subject), flat_bad_channels(1,subject),...
                reference_used_for_faster(1,subject), n_faster_bad_channels(1,subject), faster_bad_channels(1,subject),...
                ica_prep_input_n_channels(1,subject), ica_input_n_channels(1,subject), ...
                n_ica_preparation_bad_channels(1,subject), ica_preparation_bad_channels(1,subject), ...
                length_ica_data(subject), total_ICs(1,subject), n_ICs_removed(1,subject), ICs_removed(1,subject), ...
                total_epochs_before_artifact_rejection(subject), total_epochs_after_artifact_rejection(subject), ...
                n_total_channels_interpolated(subject), {total_channels_interpreted{subject}});
            report_table_newrow.Properties.VariableNames={'ID','datafile_names_preproc', 'path_preproc','date_preproc',...
                'site','exp_comments','flash_check',...
                'n_allChs', 'n_flatbadChs', 'label_flatbadChs', ...
                'reference_used_for_faster', 'n_FASTERbadChs','label_FASTERbadChs', ...
                'n_pre_prepica_Chs','n_post_prepica_Chs',...
                'n_prepicabadChs', 'label_prepicabadChs',...
                'length_ica_data', 'total_ICs', 'n_ICs_removed', 'label_ICs_removed',...
                'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
                'n_total_channels_interpolated', 'label_total_channels_interpolated'};
            MADE_report_table(subject,:) = report_table_newrow;
            save('MADE_preprocessing_report.mat','MADE_report_table');
            
            clear DataFileName DataFileLocation
            
            continue % ignore rest of the processing and go to next subject
        end
        
        % RH note: after this we go back to the approach of the original
        % MADE pipeline:

        %%% Step 7.2: Run faster and remove from data %%%%%%%%%%%%%%%%%%%%%
        ref_chan=[]; ref_chan=find(strcmp({EEG.chanlocs.labels}, 'FCz')); % find FCz channel index
        % run faster
        list_properties = channel_properties(EEG, 1:EEG.nbchan, ref_chan); % run faster
        FASTbadIdx=min_z(list_properties);
        FASTbadChans=find(FASTbadIdx==1);
        FASTbadChans=FASTbadChans(FASTbadChans~=ref_chan);
        reference_used_for_faster{subject}={EEG.chanlocs(ref_chan).labels};
        EEG = eeg_checkset(EEG);
        
        % Keep track of rejected channels
        % faster
        EEG.preproc.faster_channels.number = length(FASTbadChans);
        EEG.preproc.faster_channels.labels = {EEG.chanlocs(FASTbadChans).labels};
        % all (flat + faster)
        EEG.preproc.rejected_channels.number = EEG.preproc.flat_channels.number + length(FASTbadChans);
        EEG.preproc.rejected_channels.labels = union(EEG.preproc.rejected_channels.labels, {EEG.chanlocs(FASTbadChans).labels});
        % for consistency (for some reason, union transposes cell if one of the inputs is empty)
        if size(EEG.preproc.rejected_channels.labels, 1) > size(EEG.preproc.rejected_channels.labels, 2)
            EEG.preproc.rejected_channels.labels = EEG.preproc.rejected_channels.labels';
        end

        
        % If FASTER identifies all channels as bad channels, save the dataset
        % at this stage and ignore the remaining of the preprocessing.
        if numel(FASTbadChans)==EEG.nbchan || numel(FASTbadChans)+1==EEG.nbchan || ...
                numel(FASTbadChans)==N_channels_analysed || numel(FASTbadChans)+1==N_channels_analysed || ...
                numel(FASTbadChans) >= round((N_channels_analysed-1)/100*10)
            all_chan_bad_FAST=1;
            warning(['No usable data for datafile', datafile_names{subject}, ': too many bad FASTER channels']);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_faster_bad_channels'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_faster_bad_channels.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_faster_bad_channels.set');
                DataFileLocation = [output_location filesep 'processed_data'];
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_faster_bad_channels.mat')], 'EEG'); % save .mat format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_faster_bad_channels.set');
                DataFileLocation = [output_location filesep 'processed_data'];
            end
        else
            % Reject channels that are bad as identified by Faster
            EEG = pop_select( EEG,'nochannel', FASTbadChans);
            EEG = eeg_checkset(EEG);
            ref_chan=[]; ref_chan=find(strcmp({EEG.chanlocs.labels}, 'FCz'));   
            EEG = pop_select( EEG,'nochannel', ref_chan); % remove reference channel
        end

        %keep track for report table
        n_faster_bad_channels{subject}=EEG.preproc.faster_channels.number;
        if numel(FASTbadChans)==0
            faster_bad_channels{subject}='0';
        else
            faster_bad_channels{subject}=EEG.preproc.faster_channels.labels;
        end

        if all_chan_bad_FAST==1 % set rest of report variables to 0
            % ica
            ica_prep_input_n_channels{subject}=0;
            ica_input_n_channels{subject}=0;
            n_ica_preparation_bad_channels{subject}=0;
            ica_preparation_bad_channels{subject}='0';
            ica_input_n_channels{subject}=0;
            length_ica_data(subject)=0;
            total_ICs(subject)=0;
            n_ICs_removed{subject}=0;
            ICs_removed{subject}='0';
            % epoch rejection
            total_epochs_before_artifact_rejection(subject)=0;
            total_epochs_after_artifact_rejection(subject)=0;
            %final report
            n_total_channels_interpolated(subject)=0;
            total_channels_interpreted{subject} = '0';
            % info
            datafile_name_preproc_done{subject} = DataFileName;
            preproc_path{subject} = DataFileLocation;
            date_preproc{subject} = datetime('now');
            
            % add to current subject to report table 
            cd('/Users/riannehaartsen/Documents/02a_LEAP_EEG/LEAP_SocNSocVids/01_EEGdata/Preproc_MADE')
            ID = extractBefore(datafile_names{subject},'_social');
            % report table
            if exist('MADE_preprocessing_report.mat','file')
                load MADE_preprocessing_report.mat
            end
            report_table_newrow = table({ID}, datafile_name_preproc_done(1,subject), preproc_path(1,subject), date_preproc(1,subject), ...
                {SiteCur}, {ExpComments}, {FlashChecks(1,subject)}, ...
                n_all_channels(1,subject), n_flat_bad_channels(1,subject), flat_bad_channels(1,subject),...
                reference_used_for_faster(1,subject), n_faster_bad_channels(1,subject), faster_bad_channels(1,subject),...
                ica_prep_input_n_channels(1,subject), ica_input_n_channels(1,subject), ...
                n_ica_preparation_bad_channels(1,subject), ica_preparation_bad_channels(1,subject), ...
                length_ica_data(subject), total_ICs(1,subject), n_ICs_removed(1,subject), ICs_removed(1,subject), ...
                total_epochs_before_artifact_rejection(subject), total_epochs_after_artifact_rejection(subject), ...
                n_total_channels_interpolated(subject), {total_channels_interpreted{subject}});
            report_table_newrow.Properties.VariableNames={'ID','datafile_names_preproc', 'path_preproc','date_preproc',...
                'site','exp_comments','flash_check',...
                'n_allChs', 'n_flatbadChs', 'label_flatbadChs', ...
                'reference_used_for_faster', 'n_FASTERbadChs','label_FASTERbadChs', ...
                'n_pre_prepica_Chs','n_post_prepica_Chs',...
                'n_prepicabadChs', 'label_prepicabadChs',...
                'length_ica_data', 'total_ICs', 'n_ICs_removed', 'label_ICs_removed',...
                'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
                'n_total_channels_interpolated', 'label_total_channels_interpolated'};
            MADE_report_table(subject,:) = report_table_newrow;
            save('MADE_preprocessing_report.mat','MADE_report_table');
            
            clear DataFileName DataFileLocation
            
            continue % ignore rest of the processing and go to next subject
        end

        %% Save data after running filter and FASTER function, if saving interim results was preferred
        if save_interim_result ==1
            if output_format==1
                EEG = eeg_checkset( EEG );
                EEG = pop_editset(EEG, 'setname', strrep(datafile_names{subject}, ext, '_filtered_data'));
                EEG = pop_saveset( EEG,'filename',strrep(datafile_names{subject}, ext, '_filtered_data.set'),'filepath', [output_location filesep 'filtered_data' filesep]); % save .set format
            elseif output_format==2
                save([[output_location filesep 'filtered_data' filesep ] strrep(datafile_names{subject}, ext, '_filtered_data.mat')], 'EEG'); % save .mat format
            end
        end

        %% STEP 8: Prepare data for ICA
        EEG_copy=[];
        EEG_copy=EEG; % make a copy of the dataset
        EEG_copy = eeg_checkset(EEG_copy);
        % note down number of channels used for ica prep
        ica_prep_input_n_channels{subject} = EEG_copy.nbchan;

        % Perform 1Hz high pass filter on copied dataset
        transband = 1;
        fl_cutoff = transband/2;
        fl_order = 3.3 / (transband / EEG.srate);

        if mod(floor(fl_order),2) == 0
            fl_order=floor(fl_order);
        elseif mod(floor(fl_order),2) == 1
            fl_order=floor(fl_order)+1;
        end

        EEG_copy = pop_firws(EEG_copy, 'fcutoff', fl_cutoff, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', fl_order, 'minphase', 0);
        EEG_copy = eeg_checkset(EEG_copy);

        % Create 1 second epoch
        EEG_copy=eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1], 'rmbase', [NaN], 'eventtype', '999'); % insert temporary marker 1 second apart and create epochs
        EEG_copy = eeg_checkset(EEG_copy);

        % Find bad epochs and delete them from dataset
        vol_thrs = [-1000 1000]; % [lower upper] threshold limit(s) in uV.
        emg_thrs = [-100 30]; % [lower upper] threshold limit(s) in dB.
        emg_freqs_limit = [20 40]; % [lower upper] frequency limit(s) in Hz.

        % Find channel/s with xx% of artifacted 1-second epochs and delete them
        chanCounter = 1; ica_prep_badChans = [];
        numEpochs =EEG_copy.trials; % find the number of epochs
        all_bad_channels_ica=0;

        for ch=1:EEG_copy.nbchan
            % Find artifacted epochs by detecting outlier voltage
            EEG_copy = pop_eegthresh(EEG_copy,1, ch, vol_thrs(1), vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax, 0, 0);
            EEG_copy = eeg_checkset( EEG_copy );

            % 1         : data type (1: electrode, 0: component)
            % 0         : display with previously marked rejections? (0: no, 1: yes)
            % 0         : reject marked trials? (0: no (but store the  marks), 1:yes)

            % Find artifaceted epochs by using thresholding of frequencies in the data.
            % this method mainly rejects muscle movement (EMG) artifacts
            EEG_copy = pop_rejspec( EEG_copy, 1,'elecrange',ch ,'method','fft','threshold', emg_thrs, 'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);

            % method                : method to compute spectrum (fft)
            % threshold             : [lower upper] threshold limit(s) in dB.
            % freqlimits            : [lower upper] frequency limit(s) in Hz.
            % eegplotplotallrej     : 0 = Do not superpose rejection marks on previous marks stored in the dataset.
            % eegplotreject         : 0 = Do not reject marked trials (but store the  marks).

            % Find number of artifacted epochs
            EEG_copy = eeg_checkset( EEG_copy );
            EEG_copy = eeg_rejsuperpose( EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
            artifacted_epochs=EEG_copy.reject.rejglobal;

            % Find bad channel / channel with more than 20% artifacted epochs
            if sum(artifacted_epochs) > (numEpochs*20/100)
                ica_prep_badChans(chanCounter) = ch;
                chanCounter=chanCounter+1;
            end
        end

        % If all channels are bad, save the dataset at this stage and ignore the remaining of the preprocessing.
        if numel(ica_prep_badChans)==EEG.nbchan || numel(ica_prep_badChans)+1==EEG.nbchan || ...
                numel(ica_prep_badChans)==N_channels_analysed || numel(ica_prep_badChans)+1==N_channels_analysed || ...
                numel(ica_prep_badChans) >= round((N_channels_analysed-1)/100*10)
            all_bad_channels_ica=1;
            warning(['No usable data for datafile', datafile_names{subject}, ': too many artifacted channels from ICA']);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels_ica_prep'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels_ica_prep.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels_ica_prep.set');
                DataFileLocation = [output_location filesep 'processed_data'];
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels_ica_prep.mat')], 'EEG'); % save .mat format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels_ica_prep.mat');
                DataFileLocation = [output_location filesep 'processed_data'];
            end

        else
            % Reject bad channel - channel with more than xx% artifacted epochs
            EEG_copy = pop_select( EEG_copy,'nochannel', ica_prep_badChans);
            EEG_copy = eeg_checkset(EEG_copy);
        end
        
        %keep track for report table
        n_ica_preparation_bad_channels{subject}=length(ica_prep_badChans);
        if numel(ica_prep_badChans)==0
            ica_preparation_bad_channels{subject}='0';
        else
            ica_preparation_bad_channels{subject}={EEG.chanlocs(ica_prep_badChans).labels};
        end

        if all_bad_channels_ica == 1 % set rest of report variables to 0
            % rest of ica variables
            ica_prep_input_n_channels{subject}=0;
            n_ica_preparation_bad_channels{subject}=0;
            ica_preparation_bad_channels{subject}='0';
            ica_input_n_channels{subject}=0;
            length_ica_data(subject)=0;
            total_ICs(subject)=0;
            n_ICs_removed{subject}=0;
            ICs_removed{subject}='0';
            % epoch rejection
            total_epochs_before_artifact_rejection(subject)=0;
            total_epochs_after_artifact_rejection(subject)=0;
            %final report
            n_total_channels_interpolated(subject)=0;
            total_channels_interpreted{subject} = '0';
            % info
            datafile_name_preproc_done{subject} = DataFileName;
            preproc_path{subject} = DataFileLocation;
            date_preproc{subject} = datestr('now');
            
            % add to current subject to report table 
            cd('xxx')
            ID = extractBefore(datafile_names{subject},'_social');
            % report table
            if exist('MADE_preprocessing_report.mat','file')
                load MADE_preprocessing_report.mat
            end
            report_table_newrow = table({ID}, datafile_name_preproc_done(1,subject), preproc_path(1,subject), date_preproc(1,subject), ...
                {SiteCur}, {ExpComments}, {FlashChecks(1,subject)}, ...
                n_all_channels(1,subject), n_flat_bad_channels(1,subject), flat_bad_channels(1,subject),...
                reference_used_for_faster(1,subject), n_faster_bad_channels(1,subject), faster_bad_channels(1,subject),...
                ica_prep_input_n_channels(1,subject), ica_input_n_channels(1,subject), ...
                n_ica_preparation_bad_channels(1,subject), ica_preparation_bad_channels(1,subject), ...
                length_ica_data(subject), total_ICs(1,subject), n_ICs_removed(1,subject), ICs_removed(1,subject), ...
                total_epochs_before_artifact_rejection(subject), total_epochs_after_artifact_rejection(subject), ...
                n_total_channels_interpolated(subject), {total_channels_interpreted{subject}});
            report_table_newrow.Properties.VariableNames={'ID','datafile_names_preproc', 'path_preproc','date_preproc',...
                'site','exp_comments','flash_check',...
                'n_allChs', 'n_flatbadChs', 'label_flatbadChs', ...
                'reference_used_for_faster', 'n_FASTERbadChs','label_FASTERbadChs', ...
                'n_pre_prepica_Chs','n_post_prepica_Chs',...
                'n_prepicabadChs', 'label_prepicabadChs',...
                'length_ica_data', 'total_ICs', 'n_ICs_removed', 'label_ICs_removed',...
                'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
                'n_total_channels_interpolated', 'label_total_channels_interpolated'};
            MADE_report_table(subject,:) = report_table_newrow;
            save('MADE_preprocessing_report.mat','MADE_report_table');
            
            clear DataFileName DataFileLocation
            
            continue % ignore rest of the processing and go to next datafile
        end

        % Find the artifacted epochs across all channels and reject them before doing ICA.
        EEG_copy = pop_eegthresh(EEG_copy,1, 1:EEG_copy.nbchan, vol_thrs(1), vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax,0,0);
        EEG_copy = eeg_checkset(EEG_copy);

        % 1         : data type (1: electrode, 0: component)
        % 0         : display with previously marked rejections? (0: no, 1: yes)
        % 0         : reject marked trials? (0: no (but store the  marks), 1:yes)

        % Find artifaceted epochs by using power threshold in 20-40Hz frequency band.
        % This method mainly rejects muscle movement (EMG) artifacts.
        EEG_copy = pop_rejspec(EEG_copy, 1,'elecrange', 1:EEG_copy.nbchan, 'method', 'fft', 'threshold', emg_thrs ,'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);

        % method                : method to compute spectrum (fft)
        % threshold             : [lower upper] threshold limit(s) in dB.
        % freqlimits            : [lower upper] frequency limit(s) in Hz.
        % eegplotplotallrej     : 0 = Do not superpose rejection marks on previous marks stored in the dataset.
        % eegplotreject         : 0 = Do not reject marked trials (but store the  marks).

        % Find the number of artifacted epochs and reject them
        EEG_copy = eeg_checkset(EEG_copy);
        EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
        reject_artifacted_epochs=EEG_copy.reject.rejglobal;
        EEG_copy = pop_rejepoch(EEG_copy, reject_artifacted_epochs, 0);

        %% STEP 9: Run ICA
        ica_input_n_channels{subject} = EEG_copy.nbchan; % number of channels used for ica 
        length_ica_data(subject)=EEG_copy.trials; % length of data (in second) fed into ICA
        EEG_copy = eeg_checkset(EEG_copy);
        EEG_copy = pop_runica(EEG_copy, 'icatype', 'runica', 'extended', 1, 'stop', 1E-7, 'interupt','off');

        % Find the ICA weights that would be transferred to the original dataset
        ICA_WINV=EEG_copy.icawinv;
        ICA_SPHERE=EEG_copy.icasphere;
        ICA_WEIGHTS=EEG_copy.icaweights;
        ICA_CHANSIND=EEG_copy.icachansind;

        % If channels were removed from copied dataset during preparation of ica, then remove
        % those channels from original dataset as well before transferring ica weights.
        
        % keep track of additionally rejected channels
        % ica prep
        EEG.preproc.ica_prep_channels.number = length(ica_prep_badChans);
        EEG.preproc.ica_prep_channels.labels = {EEG.chanlocs(ica_prep_badChans).labels};
        % all 
        EEG.preproc.rejected_channels.number = EEG.preproc.rejected_channels.number + length(ica_prep_badChans);
        EEG.preproc.rejected_channels.labels = union(EEG.preproc.rejected_channels.labels, {EEG.chanlocs(ica_prep_badChans).labels});
        
        % for consistency (for some reason, union transposes cell if one of the inputs is empty)
        if size(EEG.preproc.rejected_channels.labels, 1) > size(EEG.preproc.rejected_channels.labels, 2)
            EEG.preproc.rejected_channels.labels = EEG.preproc.rejected_channels.labels';
        end
        
        EEG = eeg_checkset(EEG);
        EEG = pop_select(EEG,'nochannel', ica_prep_badChans);

        % Transfer the ICA weights of the copied dataset to the original dataset
        EEG.icawinv=ICA_WINV;
        EEG.icasphere=ICA_SPHERE;
        EEG.icaweights=ICA_WEIGHTS;
        EEG.icachansind=ICA_CHANSIND;
        
        EEG.preproc.ica.W = ICA_WEIGHTS * ICA_SPHERE;
        EEG.preproc.ica.A = ICA_WINV;
        EEG.preproc.ica.ica_channels = {EEG_copy.chanlocs(ICA_CHANSIND).labels};
        
        EEG = eeg_checkset(EEG);

        %% STEP 10: Run adjust to find artifacted ICA components
        badICs=[]; EEG_copy =[];
        EEG_copy = EEG;
        EEG_copy =eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1], 'rmbase', [NaN], 'eventtype', '999'); % insert temporary marker 1 second apart and create epochs
        EEG_copy = eeg_checkset(EEG_copy);
        
        
        cd([output_location filesep 'ica_data'])
        if size(EEG_copy.icaweights,1) == size(EEG_copy.icaweights,2)
            if save_interim_result==1
                badICs = adjusted_ADJUST(EEG_copy, [[output_location filesep 'ica_data' filesep] strrep(datafile_names{subject}, ext, '_adjust_report')]);
            else
                badICs = adjusted_ADJUST(EEG_copy, [[output_location filesep 'processed_data' filesep] strrep(datafile_names{subject}, ext, '_adjust_report')]);
            end
            close all;
        else % if rank is less than the number of electrodes, throw a warning message
            warning('The rank is less than the number of electrodes. ADJUST will be skipped. Artefacted ICs will have to be manually rejected for this participant');
        end
        
        % track bad ics
        EEG.preproc.ica.bad_components = zeros(size(ICA_WEIGHTS, 1), 1);
        EEG.preproc.ica.bad_components(badICs) = 1;
        EEG.preproc.ica.bad_components = logical(EEG.preproc.ica.bad_components);

        % Mark the bad ICs found by ADJUST
        for ic=1:length(badICs)
            EEG.reject.gcompreject(1, badICs(ic))=1;
            EEG = eeg_checkset(EEG);
        end
        total_ICs(subject)=size(EEG.icasphere, 1);
        n_ICs_removed{subject} = numel(badICs);
        if numel(badICs)==0
            ICs_removed{subject}='0';
        else
            ICs_removed{subject}=num2str(double(badICs));
        end

        %% Save dataset after ICA, if saving interim results was preferred
        if save_interim_result==1
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_ica_data'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_ica_data.set'),'filepath', [output_location filesep 'ica_data' filesep ]); % save .set format
            elseif output_format==2
                save([[output_location filesep 'ica_data' filesep ] strrep(datafile_names{subject}, ext, '_ica_data.mat')], 'EEG'); % save .mat format
            end
        end

        %% STEP 11: Remove artifacted ICA components from data
        all_bad_ICs=0;
        ICs2remove=find(EEG.reject.gcompreject); % find ICs to remove

        % If all ICs and bad, save data at this stage and ignore rest of the preprocessing for this subject.
        if numel(ICs2remove)==total_ICs(subject)
            all_bad_ICs=1;
            warning(['No usable data for datafile', datafile_names{subject}, ': too many artifacted ICs from ICA']);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.set');
                DataFileLocation = [output_location filesep 'processed_data'];
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.mat')], 'EEG'); % save .mat format
                DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.mat');
                DataFileLocation = [output_location filesep 'processed_data'];
            end
        else
            EEG = eeg_checkset( EEG );
            EEG = pop_subcomp( EEG, ICs2remove, 0); % remove ICs from dataset
        end

        if all_bad_ICs==1
            total_epochs_before_artifact_rejection(subject)=0;
            total_epochs_after_artifact_rejection(subject)=0;
            n_total_channels_interpolated(subject)=0;
            total_channels_interpreted{subject} = '0';
            % info
            datafile_name_preproc_done{subject} = DataFileName;
            preproc_path{subject} = DataFileLocation;
            date_preproc{subject} = datetime('now');
            
            
            % add to current subject to report table 
            cd('xxx')
            ID = extractBefore(datafile_names{subject},'_social');
            % report table
            if exist('MADE_preprocessing_report.mat','file')
                load MADE_preprocessing_report.mat
            end
            report_table_newrow = table({ID}, datafile_name_preproc_done(1,subject), preproc_path(1,subject), date_preproc(1,subject), ...
                {SiteCur}, {ExpComments}, {FlashChecks(1,subject)}, ...
                n_all_channels(1,subject), n_flat_bad_channels(1,subject), flat_bad_channels(1,subject),...
                reference_used_for_faster(1,subject), n_faster_bad_channels(1,subject), faster_bad_channels(1,subject),...
                ica_prep_input_n_channels(1,subject), ica_input_n_channels(1,subject), ...
                n_ica_preparation_bad_channels(1,subject), ica_preparation_bad_channels(1,subject), ...
                length_ica_data(subject), total_ICs(1,subject), n_ICs_removed(1,subject), ICs_removed(1,subject), ...
                total_epochs_before_artifact_rejection(subject), total_epochs_after_artifact_rejection(subject), ...
                n_total_channels_interpolated(subject), {total_channels_interpreted{subject}});
            report_table_newrow.Properties.VariableNames={'ID','datafile_names_preproc', 'path_preproc','date_preproc',...
                'site','exp_comments','flash_check',...
                'n_allChs', 'n_flatbadChs', 'label_flatbadChs', ...
                'reference_used_for_faster', 'n_FASTERbadChs','label_FASTERbadChs', ...
                'n_pre_prepica_Chs','n_post_prepica_Chs',...
                'n_prepicabadChs', 'label_prepicabadChs',...
                'length_ica_data', 'total_ICs', 'n_ICs_removed', 'label_ICs_removed',...
                'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
                'n_total_channels_interpolated', 'label_total_channels_interpolated'};
            MADE_report_table(subject,:) = report_table_newrow;
            save('MADE_preprocessing_report.mat','MADE_report_table');
            
            clear DataFileName DataFileLocatino
            
            continue % ignore rest of the processing and go to next datafile
        end

        %% STEP 12: Segment data into fixed length epochs
        if epoch_data==1
            if task_eeg==1 % task eeg
                EEG = eeg_checkset(EEG);
                EEG = pop_epoch(EEG, task_event_markers, task_epoch_length, 'epochinfo', 'yes');
            elseif task_eeg==0 % resting eeg - 1 condition
                if overlap_epoch==1
                    EEG=eeg_regepochs(EEG,'recurrence',(rest_epoch_length/2),'limits',[0 rest_epoch_length], 'rmbase', [NaN], 'eventtype', char(dummy_events));
                    EEG = eeg_checkset(EEG);
                else
                    EEG=eeg_regepochs(EEG,'recurrence',rest_epoch_length,'limits',[0 rest_epoch_length], 'rmbase', [NaN], 'eventtype', char(dummy_events));
                    EEG = eeg_checkset(EEG);
                end
            elseif task_eeg==2 % resting eeg - 2 condition

                % RH note: For the intended analyses, we want to segment
                % the social and non-social videos into 1-second
                % overlapping epochs, consistent with a resting state paradigm. 
                % The following code is based on the MADE code
                % for continous paradigms and has been adapted to be
                % compatible with the LEAP EEG data format. 

                % From the MADE developers:
                % Often resting state EEG is collected in eyes close and eyes open conditions
                % This script inserts event at specific time interval for resting state EEG
                % and creates separate events for eyes close and eyes open resting EEG
                % data. See MADE paper for details.

                % Initialize variables

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
                    case {'KCL','CIMH','RUNMC','UMCU'}
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
                end
            end
                % create epoch
                EEG = eeg_checkset(EEG);
                EEG = pop_epoch(EEG, new_rest_markers, [0 rest_epoch_length], 'epochinfo', 'yes');
                EEG = eeg_checkset(EEG);

        end

        total_epochs_before_artifact_rejection(subject)=EEG.trials;

        %% STEP 13: Remove baseline
        if remove_baseline==1
            EEG = eeg_checkset( EEG );
            EEG = pop_rmbase( EEG, baseline_window);
        end

        %% STEP 14: Artifact rejection
        all_bad_epochs=0;
        if voltthres_rejection==1 % check voltage threshold rejection
            if interp_epoch==1 % check epoch level channel interpolation
                chans=[]; chansidx=[];chans_labels2=[];
                chans_labels2=cell(1,EEG.nbchan);
                for i=1:EEG.nbchan
                    chans_labels2{i}= EEG.chanlocs(i).labels;
                end
                [chans,chansidx] = ismember(frontal_channels, chans_labels2);
                frontal_channels_idx = chansidx(chansidx ~= 0);
                badChans = zeros(EEG.nbchan, EEG.trials);
                badepoch=zeros(1, EEG.trials);
                if isempty(frontal_channels_idx)==1 % check whether there is any frontal channel in dataset to check
                    warning('No frontal channels from the list present in the data. Only epoch interpolation will be performed.');
                else
                    % find artifaceted epochs by detecting outlier voltage in the specified channels list and remove epoch if artifacted in those channels
                    for ch =1:length(frontal_channels_idx)
                        EEG = pop_eegthresh(EEG,1, frontal_channels_idx(ch), volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
                        EEG = eeg_checkset( EEG );
                        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                        badChans(ch,:) = EEG.reject.rejglobal;
                    end
                    for ii=1:size(badChans, 2)
                        badepoch(ii)=sum(badChans(:,ii));
                    end
                    badepoch=logical(badepoch);
                end

                % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
                if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                    all_bad_epochs=1;
                    warning(['No usable data for datafile', datafile_names{subject}]);
                    if output_format==1
                        EEG = eeg_checkset(EEG);
                        EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epoch'));
                        EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epoch.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
                    elseif output_format==2
                        save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.mat')], 'EEG'); % save .mat format
                    end
                else
                    EEG = pop_rejepoch( EEG, badepoch, 0);
                    EEG = eeg_checkset(EEG);
                end

                if all_bad_epochs==1
                    warning(['No usable data for datafile', datafile_names{subject}]);
                else
                    % Interpolate artifacted data for all reaming channels
                    badChans = zeros(EEG.nbchan, EEG.trials);
                    % Find artifacted epochs by detecting outlier voltage but don't remove
                    for ch=1:EEG.nbchan
                        EEG = pop_eegthresh(EEG,1, ch, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
                        EEG = eeg_checkset(EEG);
                        EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                        badChans(ch,:) = EEG.reject.rejglobal;
                    end
                    tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
                    for e = 1:EEG.trials
                        % Initialize variables EEGe and EEGe_interp;
                        EEGe = []; EEGe_interp = []; badChanNum = [];
                        % Select only this epoch (e)
                        EEGe = pop_selectevent( EEG, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
                        badChanNum = find(badChans(:,e)==1); % find which channels are bad for this epoch
                        EEGe_interp = eeg_interp(EEGe,badChanNum); %interpolate the bad channels for this epoch
                        tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
                    end
                    EEG.data = tmpData; % now that all of the epochs have been interpolated, write the data back to the main file

                    % If more than 10% of channels in an epoch were interpolated, reject that epoch
                    badepoch=zeros(1, EEG.trials);
                    for ei=1:EEG.trials
                        NumbadChan = badChans(:,ei); % find how many channels are bad in an epoch
                        if sum(NumbadChan) > round((10/100)*(N_channels_analysed-1))% check if more than 10% are bad
                            badepoch (ei)= sum(NumbadChan);
                        end
                    end
                    badepoch=logical(badepoch);
                end
                % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
                if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                    all_bad_epochs=1;
                    warning(['No usable data for datafile', datafile_names{subject}]);
                    if output_format==1
                        EEG = eeg_checkset(EEG);
                        EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs'));
                        EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
                    elseif output_format==2
                        save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.mat')], 'EEG'); % save .mat format
                    end
                else
                    EEG = pop_rejepoch(EEG, badepoch, 0);
                    EEG = eeg_checkset(EEG);
                end
            else % if no epoch level channel interpolation
                EEG = pop_eegthresh(EEG, 1, (1:EEG.nbchan), volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);
                EEG = eeg_checkset(EEG);
                EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                 
            end % end of epoch level channel interpolation if statement

            % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
            if sum(EEG.reject.rejthresh)==EEG.trials || sum(EEG.reject.rejthresh)+1==EEG.trials
                all_bad_epochs=1;
                warning(['No usable data for datafile', datafile_names{subject}]);
                if output_format==1
                    EEG = eeg_checkset(EEG);
                    EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs'));
                    EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
                    DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.set');
                    DataFileLocation = [output_location filesep 'processed_data'];
                elseif output_format==2
                    save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.mat')], 'EEG'); % save .mat format
                    DataFileName = strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.mat');
                    DataFileLocation = [output_location filesep 'processed_data'];
                end
            else
                EEG = pop_rejepoch(EEG,(EEG.reject.rejthresh), 0);
                EEG = eeg_checkset(EEG);
            end
        end % end of voltage threshold rejection if statement

        % if all epochs are found bad during artifact rejection
        if all_bad_epochs==1
            total_epochs_after_artifact_rejection(subject)=0;
            n_total_channels_interpolated(subject)=0;
            total_channels_interpreted{subject}='0';
            % info
            datafile_name_preproc_done{subject} = DataFileName;
            preproc_path{subject} = DataFileLocation;
            date_preproc{subject} = datetime('now');
            
            
            % add to current subject to report table 
            cd('xxx')
            ID = extractBefore(datafile_names{subject},'_social');
            % report table
            if exist('MADE_preprocessing_report.mat','file')
                load MADE_preprocessing_report.mat
            end
            report_table_newrow = table({ID}, datafile_name_preproc_done(1,subject), preproc_path(1,subject), date_preproc(1,subject), ...
                {SiteCur}, {ExpComments}, {FlashChecks}, ...
                n_all_channels(1,subject), n_flat_bad_channels(1,subject), flat_bad_channels(1,subject),...
                reference_used_for_faster(1,subject), n_faster_bad_channels(1,subject), faster_bad_channels(1,subject),...
                ica_prep_input_n_channels(1,subject), ica_input_n_channels(1,subject), ...
                n_ica_preparation_bad_channels(1,subject), ica_preparation_bad_channels(1,subject), ...
                length_ica_data(subject), total_ICs(1,subject), n_ICs_removed(1,subject), ICs_removed(1,subject), ...
                total_epochs_before_artifact_rejection(subject), total_epochs_after_artifact_rejection(subject), ...
                n_total_channels_interpolated(subject), {total_channels_interpreted{subject}});
            report_table_newrow.Properties.VariableNames={'ID','datafile_names_preproc', 'path_preproc','date_preproc',...
                'site','exp_comments','flash_check',...
                'n_allChs', 'n_flatbadChs', 'label_flatbadChs', ...
                'reference_used_for_faster', 'n_FASTERbadChs','label_FASTERbadChs', ...
                'n_pre_prepica_Chs','n_post_prepica_Chs',...
                'n_prepicabadChs', 'label_prepicabadChs',...
                'length_ica_data', 'total_ICs', 'n_ICs_removed', 'label_ICs_removed',...
                'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
                'n_total_channels_interpolated', 'label_total_channels_interpolated'};
            MADE_report_table(subject,:) = report_table_newrow;
            save('MADE_preprocessing_report.mat','MADE_report_table');
            
            clear DataFileName DataFileLocation
            
            continue % ignore rest of the processing and go to next datafile
        else
            total_epochs_after_artifact_rejection(subject)=EEG.trials;
        end

        %% STEP 15: Interpolate deleted channels
        if interp_channels==1
            EEG = eeg_interp(EEG, channels_analysed);
            EEG = eeg_checkset(EEG);
        end
        if numel(FlatbadChans)==0 && numel(FASTbadChans)==0 && numel(ica_prep_badChans)==0
            n_total_channels_interpolated(subject)=0;
        else
            n_total_channels_interpolated(subject)=numel(FlatbadChans) + numel(FASTbadChans)+ numel(ica_prep_badChans);
        end
        
        % additional check
        if isequal(n_total_channels_interpolated(subject), EEG.preproc.rejected_channels.number)
            total_channels_interpreted{subject}=EEG.preproc.rejected_channels.labels;
        else
            warning('Mismatch in number of channels rejected in EEG.preproc.rejected_channels and script')
            total_channels_interpreted{subject} = 'mismatch';
        end
        
        

        %% STEP 16: Rereference data
        if rerefer_data==1
            if iscell(reref)==1
                reref_idx=zeros(1, length(reref));
                for rr=1:length(reref)
                    reref_idx(rr)=find(strcmp({EEG.chanlocs.labels}, reref{rr}));
                end
                EEG = eeg_checkset(EEG);
                EEG = pop_reref( EEG, reref_idx);
            else
                EEG = eeg_checkset(EEG);
                EEG = pop_reref(EEG, reref);
            end
        end

        %% Save processed data
        if output_format==1
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_processed_data'));
            EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_processed_data.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_processed_data.mat')], 'EEG'); % save .mat format
        end
        % keep track of name of preprocessed file
        datafile_name_preproc_done{subject} = strrep(datafile_names{subject}, ext, '_processed_data.set');
        preproc_path{subject} = strcat(output_location, filesep, 'processed_data', filesep, strrep(datafile_names{subject}, ext, '_processed_data.set'));
        date_preproc{subject} = datetime('now');
        
    catch ME
        warning('File %s was skipped because an error occured during processing: %s', datafile_names{subject}, ME.message);
        % update error log
        cd('xxx')
        if exist('error_log.mat','file')
            load error_log.mat
        end
        error_log{subject} = {datafile_names{subject}, ME};
         % error log
        error_log = error_log(~cellfun(@isempty, error_log));
        if ~isempty(error_log)
            disp('Errors have occured during the processing of files. Saving error log with file names and corresponding MException objects to output directory.')
            if exist('MADE_preprocessing_error_log.mat','file')
                load MADE_preprocessing_error_log.mat
            end
            MADEpreproc_error_log{subject} = error_log;
            save('MADEpreproc_error_log.mat','MADEpreproc_error_log')
        end
        % update other tracking variables
        datafile_name_preproc_done{subject} = 'NA';
        preproc_path{subject} = 'NA';
        date_preproc{subject} = datetime('now');
        if length(FlashChecks) < subject || isempty(FlashChecks(subject))
            FlashChecks(subject) = NaN;
        end
        if isempty(reference_used_for_faster{subject})
            reference_used_for_faster{subject} = 'NA';
        end
        if isempty(n_all_channels{subject})
            n_all_channels{subject} = 'NA';
        end
        if isempty(n_flat_bad_channels{subject})
            n_flat_bad_channels{subject} = 'NA';
        end
        if isempty(flat_bad_channels{subject})
            flat_bad_channels{subject} = 'NA';
        end
        if isempty(n_faster_bad_channels{subject})
            n_faster_bad_channels{subject} = 'NA';
        end
        if isempty(faster_bad_channels{subject})
            faster_bad_channels{subject} = 'NA';
        end
        if isempty(ica_prep_input_n_channels{subject})
            ica_prep_input_n_channels{subject} = 'NA';
        end
        if isempty(ica_input_n_channels{subject})
            ica_input_n_channels{subject} = 'NA';
        end
        if isempty(ica_preparation_bad_channels{subject})
            ica_preparation_bad_channels{subject} = 'NA';
        end
        if isempty(n_ica_preparation_bad_channels{subject})
            n_ica_preparation_bad_channels{subject} = 'NA';
        end
        if isempty(ica_preparation_bad_channels{subject})
            ica_preparation_bad_channels{subject} = 'NA';
        end
        if length(length_ica_data) < subject || isempty(length_ica_data(subject))  
            length_ica_data(subject) = NaN;
        end
        if length(total_ICs) < subject || isempty(total_ICs(subject))
            total_ICs(subject) = NaN;
        end
        if isempty(n_ICs_removed)
            n_ICs_removed{subject} = 'NA';
        end
        if isempty(ICs_removed{subject})
            ICs_removed{subject} = 'NA';
        end 
        if length(total_epochs_before_artifact_rejection) < subject || isempty(total_epochs_before_artifact_rejection(subject))
            total_epochs_before_artifact_rejection(subject) = NaN;
        end 
        if length(total_epochs_after_artifact_rejection) < subject || isempty(total_epochs_after_artifact_rejection(subject))
            total_epochs_after_artifact_rejection(subject) = NaN;
        end 
        if length(n_total_channels_interpolated) < subject || isempty(n_total_channels_interpolated(subject))
            n_total_channels_interpolated(subject) = NaN;
        end 
        if isempty(total_channels_interpreted{subject})
            total_channels_interpreted{subject} = 'NA';
        end
    end
    
    %% add to current subject to report table 
    cd('xxx')
    ID = extractBefore(datafile_names{subject},'_social');
    % report table
    if exist('MADEBOND_ld_preprocessing_report.mat','file')
        load MADEBOND_ld_preprocessing_report.mat
    end
    report_table_newrow = table({ID}, datafile_name_preproc_done(1,subject), preproc_path(1,subject), date_preproc(1,subject), ...
        {SiteCur}, {ExpComments}, {FlashChecks(1,subject)}, ...
        n_all_channels(1,subject), n_flat_bad_channels(1,subject), flat_bad_channels(1,subject),...
        reference_used_for_faster(1,subject), n_faster_bad_channels(1,subject), faster_bad_channels(1,subject),...
        ica_prep_input_n_channels(1,subject), ica_input_n_channels(1,subject), ...
        n_ica_preparation_bad_channels(1,subject), ica_preparation_bad_channels(1,subject), ...
        length_ica_data(subject), total_ICs(1,subject), n_ICs_removed(1,subject), ICs_removed(1,subject), ...
        total_epochs_before_artifact_rejection(subject), total_epochs_after_artifact_rejection(subject), ...
        n_total_channels_interpolated(subject), {total_channels_interpreted{subject}});
    report_table_newrow.Properties.VariableNames={'ID','datafile_names_preproc', 'path_preproc','date_preproc',...
        'site','exp_comments','flash_check',...
        'n_allChs', 'n_flatbadChs', 'label_flatbadChs', ...
        'reference_used_for_faster', 'n_FASTERbadChs','label_FASTERbadChs', ...
        'n_pre_prepica_Chs','n_post_prepica_Chs',...
        'n_prepicabadChs', 'label_prepicabadChs',...
        'length_ica_data', 'total_ICs', 'n_ICs_removed', 'label_ICs_removed',...
        'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
        'n_total_channels_interpolated', 'label_total_channels_interpolated'};
    MADEBOND_ld_report_table(subject,:) = report_table_newrow;
    save('MADEBOND_ld_preprocessing_report.mat','MADEBOND_ld_report_table');

end % end of subject loop



%% Correct for missing data

for ii = 1:height(MADEBOND_ld_report_table)
    if isempty(MADEBOND_ld_report_table.ID{ii})
        subject = ii;
        ID = extractBefore(datafile_names{subject},'_social');
        report_table_newrow = table({ID}, datafile_name_preproc_done(1,subject), preproc_path(1,subject), date_preproc(1,subject), ...
            {SiteCur}, {ExpComments}, {FlashChecks(1,subject)}, ...
            n_all_channels(1,subject), n_flat_bad_channels(1,subject), flat_bad_channels(1,subject),...
            reference_used_for_faster(1,subject), n_faster_bad_channels(1,subject), faster_bad_channels(1,subject),...
            ica_prep_input_n_channels(1,subject), ica_input_n_channels(1,subject), ...
            n_ica_preparation_bad_channels(1,subject), ica_preparation_bad_channels(1,subject), ...
            length_ica_data(subject), total_ICs(1,subject), n_ICs_removed(1,subject), ICs_removed(1,subject), ...
            total_epochs_before_artifact_rejection(subject), total_epochs_after_artifact_rejection(subject), ...
            n_total_channels_interpolated(subject), {total_channels_interpreted{subject}});
        report_table_newrow.Properties.VariableNames={'ID','datafile_names_preproc', 'path_preproc','date_preproc',...
            'site','exp_comments','flash_check',...
            'n_allChs', 'n_flatbadChs', 'label_flatbadChs', ...
            'reference_used_for_faster', 'n_FASTERbadChs','label_FASTERbadChs', ...
            'n_pre_prepica_Chs','n_post_prepica_Chs',...
            'n_prepicabadChs', 'label_prepicabadChs',...
            'length_ica_data', 'total_ICs', 'n_ICs_removed', 'label_ICs_removed',...
            'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
            'n_total_channels_interpolated', 'label_total_channels_interpolated'};
        MADEBOND_ld_report_table(subject,:) = report_table_newrow;
    end
end


