%% Analysis of social and non-social videos in LEAP: Pipeline 7 - miniMADE
% using the miniMADE pipeline developed to analyse at-home EEG recordings
% with low-density lay-outs. Adaptations to the pipeline were made for
% compatibility with the EEG data formats of the LEAP study, and segmenting
% approach. No further adaptations were made.

% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')
% additional comments by Rianne Haartsen are noted with '%RH note;'

% Adapted by Rianne Haartsen, PhD.; 03-2024 
% Birkbeck College, University of London

% This adapted pipeline is released under the GNU General Public License version 3.



% ************************************************************************
% The Maryland Analysis of Developmental EEG (MADE) Pipeline
% Version 1.1
% Developed at the Child Development Lab, University of Maryland, College Park

% Contributors to MADE pipeline:
% Ranjan Debnath (ranjan.ju@gmail.com)
% George A. Buzzell (georgebuzzell@gmail.com)
% Santiago Morales Pamplona (moraless@umd.edu)
% Stephanie C. Leach (sleach12@umd.edu)
% Maureen Elizabeth Bowers (mbowers1@umd.edu)
% Nathan A. Fox (fox@umd.edu)

% ----------------------------------------------------------------------- %
% Versions Log
%
% v1_1: Sept. 23, 2020 - modified by Stephanie C. Leach
%       -Fixed bugs arising from changes in EEGLab and Matlab functions
%           * pop_rmbase() no longer works with [] as the second argument.
%             The baseline removal code has been updated to change [] to 
%             an accepted argument in the baseline removal section
%       -Fixed a bug with the FASTER bad channel section
%           * removing the reference channel before removing the channels
%             marked as bad by FASTER causes an error if the reference
%             channel isn't the last row in the data matrix. We've now 
%             moved the reference channel removal to after bad channel 
%             removal (reference channel identified again)
%       -Added a flat channel check to the ICA prep (cleaning) code
%           * helps prevent ICA decompositions with less ICs than electrodes
%       -Changed ICA prep so that participants who lost more than 20% of
%        electrodes are saved as not having any usable data
%           * prevents crashing during automated IC rejection
%           * also prevents datasets with >20% channels interpolated
%       -Added a rank check before ADJUST or Adjusted-ADJUST
%           * avoids crashing when rank is < total number of channels
%       -Added option to run a version of MADE optimized for low-density EEG
%        systems (<32 channels)... called miniMADE
%           * skips the FASTER and ICA
%           * skipping interpolation is optional, but recommended
%           * modifies artifact rejection steps to include checks for flat
%             channels and voltage jumps
%           * advanced users have the option to replace bad channels with NaNs
%             instead of removing them from the code (NOTE: this will require
%             special code at later analysis steps to remove the NaNs before 
%             doing any more preprocessing steps or analyses)
%           * advanced users also have the option to use the user entered frontal
%             electrodes for additional epoch rejection before NaN replacement
%       -Added BIDS as a 3rd formatting option for saving data
%           * saves raw AND preprocessed data in BIDS format
%           * requires the user to fill out a few addition fields
% ----------------------------------------------------------------------- %

% MADE uses EEGLAB toolbox and some of its plugins. Before running the pipeline, you have to install the following:
% EEGLab:  https://sccn.ucsd.edu/eeglab/downloadtoolbox.php/download.php

% You also need to download the following plugins/extensions from here: https://sccn.ucsd.edu/wiki/EEGLAB_Extensions
% Specifically, download:
% MFFMatlabIO: https://github.com/arnodelorme/mffmatlabio/blob/master/README.txt
% FASTER: https://sourceforge.net/projects/faster/
% ADJUST: https://www.nitrc.org/projects/adjust/ AND Adjusted-ADJUST: included with this pipeline

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
%   Our group has modified ADJUST plugin to improve selection of ICA components containing artifacts.
%   If using our modified version, please cite the following reference.
%   Adjusted-ADJUST: Leach, S.C., Morales, S., Bowers, M. E., Buzzell, G. A., Debnath, R., Beall, D., 
%   Fox, N. A., 2020. Adjusting ADJUST: Optimizing the ADJUST Algorithm for Pediatric Data Using Geodesic 
%   Nets. Psychophysiology, 57(8), e13566.

% This pipeline is released under the GNU General Public License version 3.

% ************************************************************************

%% User input: user provide relevant information to be used for data processing
% Preprocessing of EEG data involves using some common parameters for
% every subject. This part of the script initializes the common parameters.

clear % clear matlab workspace
clc % clear matlab command window

% add path to MADE pipeline
addpath(genpath('xxx/MADE-EEG-preprocessing-pipeline'))
% for example: addpath(genpath('C:\Users\Berger\Documents\MADE-EEG-preprocessing-pipeline-master'));

addpath('xxx');% enter the path of the EEGLAB folder in this line
% for example: addpath(genpath('C:\Users\Berger\Documents\eeglab13_4_4b'));
eeglab nogui;% open eeglab


% Do you want to use miniMADE (recommended for low density (<32 channels) systems)
run_miniMADE = 1; % 0 = NO (run full MADE pipeline),  = YES (run MADE pipeline with minimal preprocessing steps)
% Note: Running miniMADE will skip the FASTER and ICA steps. Epoch level interpolation can still be performed, but is not recommended
% miniMADE also skips interim saving regardless of user selection

% 1. Enter the path of the folder that has the raw data to be analyzed
rawdata_location = 'xxx/RawEEG_data/';
FlashDatapath = 'xxx'; % RH note; this variable is specific to the LEAP social/non-social videos dataset only

% 2. Enter the path of the folder where you want to save the processed data
output_location = 'xxx/Preproc_miniMADE';

% 3. Enter the path of the channel location file
channel_locations = 'xxx/euaimschan.mat';

% 4. Do your data need correction for anti-aliasing filter and/or task related time offset?
adjust_time_offset = 0; % 0 = NO (no correction), 1 = YES (correct time offset)
% If your data need correction for time offset, initialize the offset time (in milliseconds)
filter_timeoffset   = 0; % anti-aliasing time offset (in milliseconds). 0 = No time offset
stimulus_timeoffset = 0; % stimulus related time offset (in milliseconds). 0 = No time offset
response_timeoffset = 0; % response related time offset (in milliseconds). 0 = No time offset
stimulus_markers = {'xxx', 'xxx'}; % enter the stimulus makers that need to be adjusted for time offset
respose_markers  = {'xxx', 'xxx'}; % enter the response makers that need to be adjusted for time offset

% 5. Do you want to down sample the data?
down_sample = 1; % 0 = NO (no down sampling), 1 = YES (down sampling)
sampling_rate = 1000; % set sampling rate (in Hz), if you want to down sample

% 6. Do you want to delete the outer layer of the channels? (Rationale has been described in MADE manuscript)
%    This function can also be used to down sample electrodes. For example, if EEG was recorded with 128 channels but you would
%    like to analyse only 64 channels, you can assign the list of channnels to be excluded in the 'outerlayer_channel' variable.    
delete_outerlayer = 0; % 0 = NO (do not delete outer layer), 1 = YES (delete outerlayer);
% If you want to delete outer layer, make a list of channels to be deleted
outerlayer_channel = {'list of channels'}; % list of channels
% recommended list for EGI 128 channel net: {'E17' 'E38' 'E43' 'E44' 'E48' 'E49' 'E113' 'E114' 'E119' 'E120' 'E121' 'E125' 'E126' 'E127' 'E128' 'E56' 'E63' 'E68' 'E73' 'E81' 'E88' 'E94' 'E99' 'E107'}

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
baseline_window = []; % baseline period in milliseconds (MS), [] = entire epoch

% 11. Do you want to remove artifact laden epoch based on voltage threshold?
voltthres_rejection = 1; % 0 = NO, 1 = YES
volt_threshold = [-100 100]; % lower and upper threshold (in uV)

% 12. Do you want to perform epoch level channel interpolation for artifact laden epoch? (see manuscript for detail)
% Note: interpolation is not recommended for systems with less than 20 channels
interp_epoch = 1; % 0 = NO, 1 = YES.
frontal_channels = {'AF8', 'AF7'}; % If you set interp_epoch = 1, enter the list of frontal channels to check (see manuscript for detail)
% recommended list for EGI 128 channel net: {'E1', 'E8', 'E14', 'E21', 'E25', 'E32', 'E17'}
% recommended list for EGI 64 channel net: {'E1', 'E5', 'E10', 'E17'}
% RH note; for the low-density layout, only AF7/8 are available 


%13. Do you want to interpolate the bad channels that were removed from data?
% Note: because miniMADE automatically skips FASTER and ICA, this field will not affect miniMADE preprocessing
interp_channels = 0; % 0 = NO (Do not interpolate), 1 = YES (interpolate missing channels)

% 14. Do you want to rereference your data?
rerefer_data = 1; % 0 = NO, 1 = YES
reref=[]; % Enter electrode name/s or number/s to be used for rereferencing
% For channel name/s enter, reref = {'channel_name', 'channel_name'};
% For channel number/s enter, reref = [channel_number, channel_number];
% For average rereference enter, reref = []; default is average rereference

% 15. Do you want to save interim results?
save_interim_result = 0; % 0 = NO (Do not save) 1 = YES (save interim results)

% 16. How do you want to save your data? .set or .mat
output_format = 1; % 1 = .set (EEGLAB data structure), 2 = .mat (Matlab data structure), 3 = BIDS format
% If you chose BIDS format, specify subject number location in the file name and the task name
subject_number_loc = [1 2]; % should enter as [start stop] locations (e.g., par001_eeg.mff would be entered as [4 6])
task_name = 'task_name'; % should enter the eeg task name you want included in the file name

% ---------------- ADVANCED OPTIONS ---------------- %
% 17. Do you want to allow missing channels in epochs?
% Note: matlab matrices do not allow for missing rows (channels for data matrix). As such, channels removed from epochs will be replaced with
%       NaNs, which will need to be removed when averaging epochs (in the case of ERPs) or calculating other metrics (e.g., spectral power)
allow_missing_chans = 0; % this will replace bad channels with NaNs, 0 = NO and 1 = YES
% If allow_missing_chans = 1 (YES), volt_threshold values will be used to determine which channels are bad and will be replaced with NaN
% If allow_missing_chans = 1 (YES), interp_epoch & interp_channels CANNOT also = 1 (YES)
% If allow_missing_chans = 1 (YES), an average rereference cannot be used (reref cannot = [])
blink_check = 0; % this will check for blinks using the frontal_channels (defined in #12) and remove epochs containing them before replacing bad channels with NaNs
% This field will only be considered if allow_missing_chans = 1 (YES)
% WARNING: Make sure frontal_channels contains a list of frontal channels to check... If this variable is not properly defined the code will crash
chan_thresh = 0; % acceptable values are 0-1 and represent the percent of channels that must be good to keep an epoch
% for example: chan_thresh = 0.8 would remove epochs where greater than 80% of channels were replaced by NaNs

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

%% Check that ADVANCED pipeline selections are compatible with other preprocessing selections
if allow_missing_chans == 1 && (interp_epoch == 1 || interp_channels == 1)
    error(['The allow_missing_chans option (ADVANCED) cannot be turned on if channel interpolation is on...' ...
        ' allow_missing_chans does not allow for channel interpolation. Please make sure interp_epoch and interp_channels are off']);
end

if allow_missing_chans == 1 && rerefer_data == 1 && isempty(reref)
    error(['An average rereference cannot be used if the allow_missing_chans option (ADVANCED) is on...' ...
        ' allow_missing_chans does not allow for average reference. Please ensure only a subset of channels are used for rereferencing']);
end

if allow_missing_chans == 1 && voltthres_rejection == 0
    warning('voltage threshold rejection thresholds will still be used to select bad channels and replace them with NaNs');
end

%% Create output folders to save data
if output_format < 3 % if not BIDS format
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
elseif output_format == 3 % if BIDS format
    % check if derivatives folder already exists and create it if not (where we will save preprocessed data)
    if exist([ output_location filesep 'derivatives']) == 0
        mkdir([ output_location filesep 'derivatives'])
    end
    % check if eegpreprocess folder already exists and create it if not (where we will save preprocessed data)
    if exist([ output_location filesep 'derivatives' filesep 'eegpreprocess']) == 0
        mkdir([ output_location filesep 'derivatives' filesep 'eegpreprocess'])
    end
    % create file path for derivatives folder
    output_location_derivatives = [output_location filesep 'derivatives'];
end

%% Initialize output variables
reference_used_for_faster={}; % reference channel used for running faster to identify bad channel/s
faster_bad_channels=[]; % number of bad channel/s identified by faster
faster_bad_channels_labels=[]; % number of bad channel/s identified by faster
ica_preparation_bad_channels={}; % number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
ica_preparation_bad_channels_labels={}; % number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
length_ica_data=[]; % length of data (in second) fed into ICA decomposition
total_ICs=[]; % total independent components (ICs)
ICs_removed=[]; % number of artifacted ICs
total_epochs_before_artifact_rejection=[];
total_epochs_after_artifact_rejection=[];
total_channels_interpolated=[]; % total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels
curdate=datetime('now'); % set current date here so that table won't crash if date changes

%% Loop over all data files
for subject=3:length(datafile_names)
    EEG=[];
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});
    
    %% STEP 1: Import EGI data file and relevant information
    % EEG = mff_import([rawdata_location filesep datafile_names{subject}]);
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
    % for EEGLAB formatted data, use the following line
%     EEG = pop_loadset('filename',datafile_names{subject},'filepath',rawdata_location);
    % for non-.mff raw data. For example, to import biosemi data, use biosig plugin.
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
%    % remove data after last trsp (OPTIONAL for EGI files... useful when file has noisy data at the end)
%    trsp_flags = find(strcmp({EEG.event.type},'TRSP')); % find indices of TRSP flags
%    EEG = eeg_eegrej( EEG, [(EEG.event(trsp_flags(end)).latency+(1.5*EEG.srate)) EEG.pnts] ); % remove everything 1.5 seconds after the last TRSP
%    EEG = eeg_checkset( EEG );
    
    %% if BIDS output format requested
    if output_format == 3
        current_subject = ['sub-' datafile_names{subject}(subject_number_loc(1):subject_number_loc(2))];
        % create folder to save BIDS formatted raw data (no preprocessing)
        if exist([ output_location filesep current_subject]) == 0
            mkdir([ output_location filesep current_subject]) % create subject folder
            mkdir([ output_location filesep current_subject filesep 'eeg']) % create eeg folder in subject folder
        end
        % save raw data in BIDS format
        EEG = eeg_checkset(EEG);
        EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01_eeg']);
        EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01_eeg.set'],'filepath', [output_location filesep current_subject filesep 'eeg']); % save BIDS format
        % save other variables in BIDS format
        cd([output_location filesep current_subject filesep 'eeg'])
        events_to_tsv(EEG); % save event info
        channelloc_to_tsv(EEG); % save EEG.chanlocs in .tsv format
        electrodes_to_tsv(EEG); % save electrode info
        eeg_descript = jsonread( [fileparts(which('template_eeg.json')) filesep 'template_eeg.json'] );
        try eeg_descript.TaskName = task_name; eeg_descript.EEGChannelCount = EEG.nbchan; eeg_descript.RecordingDuration = EEG.xmax; eeg_descript.SamplingFrequency = EEG.srate; catch; warning('Task name, channel count, recording duration, and sampling rate may not have auto populated in eeg.json file'); end
        jsonwrite([output_location filesep current_subject filesep 'eeg' filesep current_subject '_task-' task_name '_run-01_eeg.json'], eeg_descript,struct('indent','  '));
        % create subject folder within derivatives folder for preprocessed data
        if exist([ output_location_derivatives filesep 'eegpreprocess' filesep current_subject]) == 0
            mkdir([ output_location_derivatives filesep 'eegpreprocess' filesep current_subject])
        end
    end
    
    
    %% STEP 2: Import channel locations
    % EEG=pop_chanedit(EEG, 'load',{channel_locations 'filetype' 'autodetect'});
    % EEG = eeg_checkset( EEG );
    
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
    % Throw a warning if system is low density and user has selected to run the full MADE pipeline or interpolation
    if length(EEG.chanlocs) < 32 % if a low density system
        if run_miniMADE == 0
            warning('Running MADE is not recommended for low-density systems. To run miniMADE instead set run_miniMADE equal to 0 at the top of the script');
        end
        if interp_channels == 1 && length(EEG.chanlocs) < 20
            warning('Channel interpolation is not recommended for low-density systems with fewer than 20 channels');
        end
    end
    % if reref is a vector (type double) instead of a cell array, grab the electrode name(s)
    if ~iscell(reref); reref = {EEG.chanlocs(reref).labels}; end % this will prevent using the wrong channel in cases where we remove channels
    
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


    if down_sample==1
        if floor(sampling_rate) > EEG.srate
            error ('Sampling rate cannot be higher than recorded sampling rate');
        elseif floor(sampling_rate) ~= EEG.srate
            EEG = pop_resample( EEG, sampling_rate);
            EEG = eeg_checkset( EEG );
        end
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
    

    %% RH: Step 6 1/2: Remove non-SNSvid data in LEAP recordings and check the flash recordings
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




    %% STEP 7: Run faster to find bad channels
    if run_miniMADE == 0
        % First check whether reference channel (i.e. zeroed channels) is present in data
        % reference channel is needed to run faster
        ref_chan=[]; FASTbadChans=[]; all_chan_bad_FAST=0;
        ref_chan=find(any(EEG.data, 2)==0);
        if numel(ref_chan)>1
            error(['There are more than 1 zeroed channel (i.e. zero value throughout recording) in data.'...
                ' Only reference channel should be zeroed channel. Delete the zeroed channel/s which is not reference channel.']);
        elseif numel(ref_chan)==1
            list_properties = channel_properties(EEG, 1:EEG.nbchan, ref_chan); % run faster
            FASTbadIdx=min_z(list_properties);
            FASTbadChans=find(FASTbadIdx==1);
            FASTbadChans=FASTbadChans(FASTbadChans~=ref_chan);
            reference_used_for_faster{subject}={EEG.chanlocs(ref_chan).labels};
            EEG = eeg_checkset(EEG);
            channels_analysed=EEG.chanlocs; % keep full channel locations to use later for interpolation of bad channels
        elseif numel(ref_chan)==0
            warning('Reference channel is not present in data. Cz channel will be used as reference channel.');
            ref_chan=find(strcmp({EEG.chanlocs.labels}, 'Cz')); % find Cz channel index
            EEG_copy=[];
            EEG_copy=EEG; % make a copy of the dataset
            EEG_copy = pop_reref( EEG_copy, ref_chan,'keepref','on'); % rerefer to Cz in copied dataset
            EEG_copy = eeg_checkset(EEG_copy);
            list_properties = channel_properties(EEG_copy, 1:EEG_copy.nbchan, ref_chan); % run faster on copied dataset
            FASTbadIdx=min_z(list_properties);
            FASTbadChans=find(FASTbadIdx==1);
            channels_analysed=EEG.chanlocs;
            reference_used_for_faster{subject}={EEG.chanlocs(ref_chan).labels};
        end

        % If FASTER identifies all channels as bad channels, save the dataset
        % at this stage and ignore the remaining of the preprocessing.
        if numel(FASTbadChans)==EEG.nbchan || numel(FASTbadChans)+1==EEG.nbchan
            all_chan_bad_FAST=1;
            warning(['No usable data for datafile', datafile_names{subject}]);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.mat')], 'EEG'); % save .mat format
            elseif output_format==3
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_channels']);
                EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_channels.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
            end
        else
            % Reject channels that are bad as identified by Faster
            EEG = pop_select( EEG,'nochannel', FASTbadChans);
            EEG = eeg_checkset(EEG);
            if numel(ref_chan)==1
                ref_chan=find(any(EEG.data, 2)==0);
                EEG = pop_select( EEG,'nochannel', ref_chan); % remove reference channel
            end
        end

        if numel(FASTbadChans)==0
            faster_bad_channels{subject}='0';
            faster_bad_channels_labels{subject}='0';
        else
            faster_bad_channels{subject}=num2str(FASTbadChans');
            faster_bad_channels_labels{subject}=strjoin({channels_analysed(FASTbadChans).labels});
        end

        if all_chan_bad_FAST==1
            faster_bad_channels{subject}='0';
            faster_bad_channels_labels{subject}='0';
            ica_preparation_bad_channels{subject}='0';
            ica_preparation_bad_channels_labels{subject}='0';
            length_ica_data(subject)=0;
            total_ICs(subject)=0;
            ICs_removed{subject}='0';
            total_epochs_before_artifact_rejection(subject)=0;
            total_epochs_after_artifact_rejection(subject)=0;
            total_channels_interpolated(subject)=0;
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
            elseif output_format==3
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01', '_desc-filtered_eeg']);
                EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01', '_desc-filtered_eeg.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
            end
        end
    
    %% STEP 8: Prepare data for ICA
        EEG_copy=[];
        EEG_copy=EEG; % make a copy of the dataset
        EEG_copy = eeg_checkset(EEG_copy);

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
        all_bad_channels=0;

        for ch=1:EEG_copy.nbchan
            % Find artifaceted epochs by detecting outlier voltage
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

            % Find flat channels
            flatChans = find(range(squeeze(EEG_copy.data(ch,:,:)),1) < 1);
            if ~isempty(flatChans)
                artifacted_epochs(flatChans) = 1;
            end
            
            % Find bad channel / channel with more than 20% artifacted epochs
            if sum(artifacted_epochs) > (numEpochs*20/100)
                ica_prep_badChans(chanCounter) = ch;
                chanCounter=chanCounter+1;
            end
        end

        % If all channels are bad, save the dataset at this stage and ignore the remaining of the preprocessing.
        if numel(ica_prep_badChans)==EEG.nbchan || numel(ica_prep_badChans)+1==EEG.nbchan
            all_bad_channels=1;
            warning(['No usable data for datafile', datafile_names{subject}]);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.mat')], 'EEG'); % save .mat format
            elseif output_format==3
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_channels']);
                EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_channels.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
            end
        % If > 20%, but not all of the channels are bad, save the dataset at this stage and ignore the remaining of the preprocessing.
        elseif numel(ica_prep_badChans) > 0.2*EEG.nbchan
            all_bad_channels=1;
            warning(['No usable data for datafile', datafile_names{subject}]);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_too_many_bad_channels'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_too_many_bad_channels.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_too_many_bad_channels.mat')], 'EEG'); % save .mat format
            elseif output_format==3
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_too_many_bad_channels']);
                EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_too_many_bad_channels.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
            end
        else
            % Reject bad channel - channel with more than xx% artifacted epochs
            ica_preparation_bad_channels_labels{subject}=strjoin({EEG.chanlocs(ica_prep_badChans).labels});
            EEG_copy = pop_select( EEG_copy,'nochannel', ica_prep_badChans);
            EEG_copy = eeg_checkset(EEG_copy);
        end

        if numel(ica_prep_badChans)==0
            ica_preparation_bad_channels{subject}='0';
            ica_preparation_bad_channels_labels{subject}='0';
        else
            ica_preparation_bad_channels{subject}=num2str(ica_prep_badChans);
        end

        if all_bad_channels == 1
            length_ica_data(subject)=0;
            total_ICs(subject)=0;
            ICs_removed{subject}='0';
            total_epochs_before_artifact_rejection(subject)=0;
            total_epochs_after_artifact_rejection(subject)=0;
            total_channels_interpolated(subject)=0;
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
        EEG = eeg_checkset(EEG);
        EEG = pop_select(EEG,'nochannel', ica_prep_badChans);

        % Transfer the ICA weights of the copied dataset to the original dataset
        EEG.icawinv=ICA_WINV;
        EEG.icasphere=ICA_SPHERE;
        EEG.icaweights=ICA_WEIGHTS;
        EEG.icachansind=ICA_CHANSIND;
        EEG = eeg_checkset(EEG);
    
    %% STEP 10: Run adjust to find artifacted ICA components
        badICs=[]; EEG_copy =[];
        EEG_copy = EEG;
        EEG_copy =eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1], 'rmbase', [NaN], 'eventtype', '999'); % insert temporary marker 1 second apart and create epochs
        EEG_copy = eeg_checkset(EEG_copy);
        
        if size(EEG_copy.icaweights,1) == size(EEG_copy.icaweights,2)
            if output_format == 3 % if BIDS format
                cd([output_location_derivatives filesep current_subject]) % move to this folder so that bump jpeg saves in derivatives folder
                badICs = adjusted_ADJUST(EEG_copy, [output_location_derivatives filesep current_subject '_adjust_report']); % save in raw data folder
            else % if not BIDS format
                if save_interim_result==1
                    badICs = adjusted_ADJUST(EEG_copy, [[output_location filesep 'ica_data' filesep] strrep(datafile_names{subject}, ext, '_adjust_report')]);
                else
                    badICs = adjusted_ADJUST(EEG_copy, [[output_location filesep 'processed_data' filesep] strrep(datafile_names{subject}, ext, '_adjust_report')]);
                end
            end
            close all;
        else % if rank is less than the number of electrodes, throw a warning message
            warning('The rank is less than the number of electrodes. ADJUST will be skipped. Artefacted ICs will have to be manually rejected for this participant');
        end
        
        % Mark the bad ICs found by ADJUST
        for ic=1:length(badICs)
            EEG.reject.gcompreject(1, badICs(ic))=1;
            EEG = eeg_checkset(EEG);
        end
        total_ICs(subject)=size(EEG.icasphere, 1);
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
            elseif output_format==3
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01', '_desc-ica_eeg']);
                EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01', '_desc-ica_eeg.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
            end
        end
    
    %% STEP 11: Remove artifacted ICA components from data
        all_bad_ICs=0;
        ICs2remove=find(EEG.reject.gcompreject); % find ICs to remove

        % If all ICs and bad, save data at this stage and ignore rest of the preprocessing for this subject.
        if numel(ICs2remove)==total_ICs(subject)
            all_bad_ICs=1;
            warning(['No usable data for datafile', datafile_names{subject}]);
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.mat')], 'EEG'); % save .mat format
            elseif output_format==3
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_ICs']);
                EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_ICs.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
            end
        else
            EEG = eeg_checkset( EEG );
            EEG = pop_subcomp( EEG, ICs2remove, 0); % remove ICs from dataset
        end

        if all_bad_ICs==1
            total_epochs_before_artifact_rejection(subject)=0;
            total_epochs_after_artifact_rejection(subject)=0;
            total_channels_interpolated(subject)=0;
            continue % ignore rest of the processing and go to next datafile
        end
    end % end if miniMADE check
    
    %% STEP 12: Segment data into fixed length epochs
    if epoch_data==1
        if task_eeg ==1 % task eeg
            EEG = eeg_checkset(EEG);
            EEG = pop_epoch(EEG, task_event_markers, task_epoch_length, 'epochinfo', 'yes');
        elseif task_eeg==0 % resting eeg
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

                    % create epoch
                    EEG = eeg_checkset(EEG);
                    EEG = pop_epoch(EEG, new_rest_markers, [0 rest_epoch_length], 'epochinfo', 'yes');
                    EEG = eeg_checkset(EEG);

        end
                
    end
    
    total_epochs_before_artifact_rejection(subject)=EEG.trials;
    
    %% STEP 13: Remove baseline
    if remove_baseline==1
        if isempty(baseline_window) % set up for entire epoch
            baseline_window = [EEG.times(1) EEG.times(end)]; % set start and stop times for the baseline window
        end
        EEG = eeg_checkset( EEG );
        EEG = pop_rmbase( EEG, baseline_window);
    end
    
    %% STEP 14: Artifact rejection
    all_bad_epochs=0;
    if allow_missing_chans == 0 
        if voltthres_rejection==1 % check voltage threshold rejection
            if interp_epoch==1 % check epoch level channel interpolation
                % first pass: loop through frontal channels and reject bad epochs
                %   - removes epochs with (possible) blinks
                %   - same steps for MADE and miniMADE (if miniMADE uses interpolation)
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
                else
                    EEG = pop_rejepoch( EEG, badepoch, 0);
                    EEG = eeg_checkset(EEG);
                end

                % second pass: loop through all channels and interpolate remaining bad channels at the epoch level
                %   - miniMADE has extra artifact checks at this step
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
                    if run_miniMADE == 0
                        for e = 1:EEG.trials
                            % Initialize variables EEGe and EEGe_interp;
                            EEGe = []; EEGe_interp = []; badChanNum = [];
                            % Select only this epoch (e)
                            EEGe = pop_selectevent( EEG, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
                            badChanNum = find(badChans(:,e)==1); % find which channels are bad for this epoch
                            EEGe_interp = eeg_interp(EEGe,badChanNum); %interpolate the bad channels for this epoch
                            tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
                        end
                    elseif run_miniMADE == 1
                        for e = 1:EEG.trials
                            EEGe = []; EEGe_interp = []; badChanNum = []; % Initialize variables EEGe and EEGe_interp;
                            %select only this epoch (e)
                            EEGe = pop_selectevent( EEG, 'epoch',e,'deleteevents','off','deleteepochs','on','invertepochs','off');
                            badChanNum = find(badChans(:,e)==1); %find which channels are bad for this epoch
                            % find and add flat chans to the bad chans list
                            flatChanNum = find(range(EEGe.data,2) < 1);
                            badChanNum  = unique([badChanNum; flatChanNum]);
                            % find chans with a large jump/deflection (also bad chans) by taking 1st derivative
                            [jump_chans, ~] = find( abs(diff(EEGe.data,1,2) ./ repmat(diff(1:EEGe.pnts),EEGe.nbchan,1)) > 50);
                            badChanNum = unique([badChanNum; unique(jump_chans)]);
                            % interpolate using bad channel list with extra checks
                            if length(badChanNum) < EEGe.nbchan - 1% script will crash if we try to interpolate with 0 or 1 channels left 
                                EEGe_interp = eeg_interp(EEGe,badChanNum); %interpolate the bad channels for this epoch
                                tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
                            end
                            badChans(badChanNum,e) = 1; % modify
                            % keep track of flat and jump channel information
                            %flat_mat(e) = ~isempty(flatChanNum);
                            %jump_mat(e) = length(unique(jump_chans));
                        end
                    end
                    EEG.data = tmpData; % now that all of the epochs have been interpolated, write the data back to the main file

                    % If more than 10% of channels in an epoch were interpolated, reject that epoch
                    badepoch=zeros(1, EEG.trials);
                    for ei=1:EEG.trials
                        NumbadChan = badChans(:,ei); % find how many channels are bad in an epoch
                        if sum(NumbadChan) > round((10/100)*EEG.nbchan)% check if more than 10% are bad
                            badepoch (ei)= sum(NumbadChan);
                        end
                    end
                    badepoch=logical(badepoch);
                end
                % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
                if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                    all_bad_epochs=1;
                    warning(['No usable data for datafile', datafile_names{subject}]);
                else
                    EEG = pop_rejepoch(EEG, badepoch, 0);
                    EEG = eeg_checkset(EEG);
                end
            else % if no epoch level channel interpolation
                if run_miniMADE == 0
                    EEG = pop_eegthresh(EEG, 1, (1:EEG.nbchan), volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);
                    EEG = eeg_checkset(EEG);
                    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                elseif run_miniMADE == 1
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
                        EEGe = []; badChanNum = []; % Initialize variables EEGe and badChanNum;
                        %select only this epoch (e)
                        EEGe = pop_selectevent( EEG, 'epoch',e,'deleteevents','off','deleteepochs','on','invertepochs','off');
                        badChanNum = find(badChans(:,e)==1); %find which channels are bad for this epoch
                        % find and add flat chans to the bad chans list
                        flatChanNum = find(range(EEGe.data,2) < 1);
                        badChanNum  = unique([badChanNum; flatChanNum]);
                        % find chans with a large jump/deflection (also bad chans) by taking 1st derivative
                        [jump_chans, ~] = find( abs(diff(EEGe.data,1,2) ./ repmat(diff(1:EEGe.pnts),EEGe.nbchan,1)) > 50);
                        badChanNum = unique([badChanNum; unique(jump_chans)]);
                        % add any new bad channels back to the bad chan list
                        badChans(badChanNum,e) = 1; % modify
                        % keep track of flat and jump channel information
                        %flat_mat(e) = ~isempty(flatChanNum);
                        %jump_mat(e) = length(unique(jump_chans));
                    end
                    badepoch=zeros(1, EEG.trials);
                    for ei=1:EEG.trials
                        if sum(badChans(:,ei)) > 0 % check if there are any bad chans
                            badepoch (ei)= 1;
                        end
                    end
                    badepoch=logical(badepoch);
                end
                % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
                if sum(EEG.reject.rejthresh)==EEG.trials || sum(EEG.reject.rejthresh)+1==EEG.trials || sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                    all_bad_epochs=1;
                    warning(['No usable data for datafile', datafile_names{subject}]);
                else
                    if run_miniMADE == 0
                        EEG = pop_rejepoch(EEG,(EEG.reject.rejthresh), 0);
                        EEG = eeg_checkset(EEG);
                    elseif run_miniMADE == 1
                        EEG = pop_rejepoch(EEG, badepoch, 0);
                        EEG = eeg_checkset(EEG);
                    end
                end
            end % end of epoch level channel interpolation if statement
        end % end of voltage threshold rejection if statement
        
    elseif allow_missing_chans == 1 % If advanced option to replace bad channels with NaNs is selected
        if rerefer_data==1
            % grab channels for rereferencing
            if iscell(reref)==1
                reref_idx=zeros(1, length(reref));
                for rr=1:length(reref)
                    reref_idx(rr)=find(strcmp({EEG.chanlocs.labels}, reref{rr}));
                end
                reref_chans = reref_idx;
            else
                reref_chans = reref; 
            end
            % perform traditional artifact rejection for rereference channels
            EEG = pop_eegthresh(EEG, 1, reref_chans, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);
            EEG = eeg_checkset( EEG );
            EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
            
            if length(find(EEG.reject.rejthresh)) == EEG.trials || length(find(EEG.reject.rejthresh))+1 == EEG.trials % if all epochs marked for rejection
                all_bad_epochs = 1; % set at 1 so we don't save the file below
            else
                % reject epochs where rereference channels are bad
                EEG = pop_rejepoch( EEG, (EEG.reject.rejthresh), 0);
                EEG = eeg_checkset(EEG);
                % re-reference to new reference channels
                EEG = eeg_checkset(EEG);
                EEG = pop_reref( EEG, reref_chans);
            end
        end
        
        % grab new channel labels in case some channels have been removed
        chans=[]; chans_labels2=[];
        chans_labels2=cell(1,EEG.nbchan);
        for i=1:EEG.nbchan
            chans_labels2{i}= EEG.chanlocs(i).labels;
        end
            
        if all_bad_epochs == 0 && blink_check == 1 % look at user entered frontal channels and remove artefacted epochs
            frontal_channels_idx=zeros(1, length(frontal_channels));
            for rr=1:length(frontal_channels)
                frontal_channels_idx(rr)=find(strcmp({EEG.chanlocs.labels}, frontal_channels{rr}));
            end
            % perform traditional artifact rejection for rereference channels
            EEG = pop_eegthresh(EEG, 1, frontal_channels_idx, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);
            EEG = eeg_checkset( EEG );
            EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
            % only mark epochs for rejection if all frontal channels are bad
            blink_epochs = find(sum(EEG.reject.rejthreshE(frontal_channels_idx,:))==length(frontal_channels_idx));
            if length(blink_epochs) == EEG.trials || length(blink_epochs)+1 == EEG.trials  % if all epochs marked for rejection
                all_bad_epochs = 1; % set at 1 so we don't save the file below
            else
                % reject epochs where rereference channels are bad
                EEG = pop_rejepoch( EEG, blink_epochs, 0);
                EEG = eeg_checkset(EEG);
            end
        end
        
        if all_bad_epochs == 0
            % replace remaining bad channels with NaNs
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
                EEGe = []; badChanNum = []; % Initialize variables EEGe and EEGe_interp;
                %select only this epoch (e)
                EEGe = pop_selectevent( EEG, 'epoch',e,'deleteevents','off','deleteepochs','on','invertepochs','off');
                badChanNum = find(badChans(:,e)==1); %find which channels are bad for this epoch
                % find and add flat chans to the bad chans list
                flatChanNum = [find(range(EEGe.data(:,1:(EEG.pnts/2)),2) < 1); find(range(EEGe.data(:,(EEG.pnts/2):EEG.pnts),2) < 1)];
                badChanNum  = unique([badChanNum; unique(flatChanNum)]);
                % find chans with a large jump/deflection (also bad chans) by taking 1st derivative
                [jump_chans, ~] = find( abs(diff(EEGe.data,1,2) ./ repmat(diff(1:EEGe.pnts),EEGe.nbchan,1)) > 50);
                badChanNum = unique([badChanNum; unique(jump_chans)]);
                % add any new bad channels back to the bad chan list
                badChans(badChanNum,e) = 1; % modify badChan list
                % remove the bad chans for this epoch (replace with NaN)
                EEGe = eeg_checkset( EEGe );
                EEGe.data(badChanNum,:) = NaN;
                tmpData(:,:,e) = EEGe.data; % store NaN replaced data into matrix
                % keep track of flat and jump channel information
                %flat_mat(e) = ~isempty(flatChanNum);
                %jump_mat(e) = length(unique(jump_chans));
            end
            EEG.data = tmpData; % save data containing NaNs back to EEG struct
            
            % reject epochs where NaN channel numbers are above the user entered threshold
            badepoch=zeros(1, EEG.trials);
            for ei=1:EEG.trials
                if sum(badChans(:,ei)) > chan_thresh*EEG.nbchan % check if there are any bad chans
                    badepoch(ei)= 1;
                end
            end
            badepoch=logical(badepoch);
            if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                all_bad_epochs=1;
                warning(['No usable data for datafile', datafile_names{subject}]);
            else
                EEG = pop_rejepoch(EEG, badepoch, 0);
                EEG = eeg_checkset(EEG);
                badChans = badChans(:,~badepoch);
            end
                
            % reformat badChans for output table
            total_epochs_by_chan_after_artifact_rejection = {'';[]}; badChansSum = [];
            total_epochs_by_chan_after_artifact_rejection(1,1:EEG.nbchan) = strcat(chans_labels2,'_total_epochs_after_artifact_rejection');
            badChansSum = sum((badChans-1)*-1,2)'; 
            for cch = 1:EEG.nbchan
                total_epochs_by_chan_after_artifact_rejection{2,cch} = badChansSum(cch);
            end
        else
            total_epochs_by_chan_after_artifact_rejection = {'';[]}; badChansSum = [];
            total_epochs_by_chan_after_artifact_rejection(1,1:EEG.nbchan) = strcat(chans_labels2,'_total_epochs_after_artifact_rejection');
            for cch = 1:EEG.nbchan
                total_epochs_by_chan_after_artifact_rejection{2,cch} = 0;
            end
        end
        total_epochs_by_chan_after_artifact_rejection{1,EEG.nbchan+1} = 'datafile_names';
        total_epochs_by_chan_after_artifact_rejection{2,EEG.nbchan+1} =  datafile_names{subject};
    end % end advanced option to use NaNs
    
    % if all epochs are found bad during artifact rejection
    if all_bad_epochs==1
        total_epochs_after_artifact_rejection(subject)=0;
        total_channels_interpolated(subject)=0;
        if output_format==1
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs'));
            EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.mat')], 'EEG'); % save .mat format
        elseif output_format==3
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_epochs']);
            EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_epochs.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
        end
        continue % ignore rest of the processing and go to next datafile
    else
        total_epochs_after_artifact_rejection(subject)=EEG.trials;
    end
    
    %% STEP 15: Interpolate deleted channels
    if run_miniMADE == 0
        if interp_channels==1 && allow_missing_chans == 0
            EEG = eeg_interp(EEG, channels_analysed);
            EEG = eeg_checkset(EEG);
        end
        if numel(FASTbadChans)==0 && numel(ica_prep_badChans)==0
            total_channels_interpolated(subject)=0;
        else
            total_channels_interpolated(subject)=numel(FASTbadChans)+ numel(ica_prep_badChans);
        end
    end
    
    %% STEP 16: Rereference data
    if allow_missing_chans == 0
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
    end
    
    %% Save processed data
    if output_format==1
        EEG = eeg_checkset(EEG);
        EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_processed_data'));
        EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_processed_data.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
    elseif output_format==2
        save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_processed_data.mat')], 'EEG'); % save .mat format
    elseif output_format==3
        EEG = eeg_checkset(EEG);
        EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01', '_desc-processed_eeg']);
        EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01', '_desc-processed_eeg.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
    end
    
    %% Create the report table for all the data files with relevant preprocessing outputs.
    if output_format < 3
        cd(output_location)
    else
        cd(output_location_derivatives)
    end
    if run_miniMADE == 0
        % check if file already exists and add to instead of overwriting if it does
        if exist([output_location filesep 'MADE_preprocessing_report.csv'], 'file') ~= 0
            % load existing table (do every iteration in case of crash)
            report_table = readtable([output_location filesep 'MADE_preprocessing_report.csv'], 'Format', '%s%s%s%s%s%s%d%d%s%d%d%d');
            % make table for current subject
            report_table2=table(datafile_names(subject)', reference_used_for_faster{subject}', faster_bad_channels(subject)', faster_bad_channels_labels(subject)', ica_preparation_bad_channels(subject)', ica_preparation_bad_channels_labels(subject)', length_ica_data(subject)', ...
                total_ICs(subject)', ICs_removed(subject)', total_epochs_before_artifact_rejection(subject)', total_epochs_after_artifact_rejection(subject)',total_channels_interpolated(subject)');
            report_table2.Properties.VariableNames={'datafile_names', 'reference_used_for_faster', 'faster_bad_channels', 'faster_bad_channels_labels', ...
                'ica_preparation_bad_channels', 'ica_preparation_bad_channels_labels', 'length_ica_data', 'total_ICs', 'ICs_removed', 'total_epochs_before_artifact_rejection', ...
                'total_epochs_after_artifact_rejection', 'total_channels_interpolated'};
            % add current subject to existing output table
            report_table = outerjoin(report_table,report_table2,'MergeKeys', true);
        else % create the file for the first iteration
            report_table=table(datafile_names(subject)', reference_used_for_faster{subject}', faster_bad_channels', faster_bad_channels_labels', ica_preparation_bad_channels', ica_preparation_bad_channels_labels',length_ica_data', ...
                total_ICs', ICs_removed', total_epochs_before_artifact_rejection', total_epochs_after_artifact_rejection',total_channels_interpolated');
            report_table.Properties.VariableNames={'datafile_names', 'reference_used_for_faster', 'faster_bad_channels', 'faster_bad_channels_labels', ...
                'ica_preparation_bad_channels', 'ica_preparation_bad_channels_labels', 'length_ica_data', 'total_ICs', 'ICs_removed', 'total_epochs_before_artifact_rejection', ...
                'total_epochs_after_artifact_rejection', 'total_channels_interpolated'};
        end
        writetable(report_table, ['MADE_preprocessing_report.csv']);
    elseif run_miniMADE == 1
        % check if file already exists and add to instead of overwriting if it does
        if exist([output_location filesep 'miniMADE_preprocessing_report.csv'], 'file') ~= 0
            % load existing table (do every iteration in case of crash)
            report_table = readtable([output_location filesep 'miniMADE_preprocessing_report.csv']);
            % make table for current subject
            report_table2=table(datafile_names(subject)', total_epochs_before_artifact_rejection(subject)', total_epochs_after_artifact_rejection(subject)');
            report_table2.Properties.VariableNames={'datafile_names', 'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection'};
            if allow_missing_chans == 1 % add epochs by channel to output
                epochT = cell2table(total_epochs_by_chan_after_artifact_rejection(2,:), 'VariableNames', total_epochs_by_chan_after_artifact_rejection(1,:));
                report_table2 = outerjoin(report_table2,epochT,'MergeKeys', true);
            end
            % add current subject to existing output table
            report_table = outerjoin(report_table,report_table2,'MergeKeys', true);
        else % create the file for the first iteration
            report_table=table(datafile_names(subject)', total_epochs_before_artifact_rejection', total_epochs_after_artifact_rejection');
            report_table.Properties.VariableNames={'datafile_names', 'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection'};
            if allow_missing_chans == 1 % add epochs by channel to output
                epochT = cell2table(total_epochs_by_chan_after_artifact_rejection(2,:), 'VariableNames', total_epochs_by_chan_after_artifact_rejection(1,:));
                report_table = outerjoin(report_table,epochT,'MergeKeys', true);
            end
        end
        writetable(report_table, ['miniMADE_preprocessing_report.csv']);
    end
end % end of subject loop
