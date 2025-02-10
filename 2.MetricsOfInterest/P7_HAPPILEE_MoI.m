%% LEAP soc-nsoc videos Pipelines: get power and connectivity values from the low-density data preprocessed with the HAPPILEE pipeline 
% P7 - pipeline 7 in the manuscript

% This script takes the data processed with the HAPPILEE (=HAPPEv4 for low 
% density) pipeline and converges it to a fieldtrip structure. Then power 
% is calculated for 1 to 32 Hz in steps of 1Hz for each trial. Finally, 
% data are split into conditions and connectivity and power are calculated 
% for the social condition, non-social condition, and across all conditions. 

% Optional: This script also creates a reporting figure that includes plots of the
% global power and connectivity spectra for metrics across all and separate
% conditions, and topoplots for theta and alpha power and connectivity for
% the separate conditions, and condition differences. Finally, the number
% of epochs per condition, and notes on cleaning are plotted. 

% Processing is looped through participants using the IDs from the raw data 
% folder. Preprocessed data are saved. 
% The measures of interest are stored into a reporting table named
% HAPPILEE_MOI_table and saved as HAPPILEE_MeasuresOfInterest_report.mat. 


% This script call to functions from EEGLAB and fieldtrip toolboxes.
% Before running this script, you have to install the following:
% EEGLab:  https://sccn.ucsd.edu/eeglab/downloadtoolbox.php/download.php
% Fieldtrip: https://www.fieldtriptoolbox.org/

% Also cite the following in the manuscript using these functions:

% EEGLAB: A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
% analysis of single-trial EEG dynamics. Journal of Neuroscience Methods, 134, 9?21.

% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. 
% FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and 
% Invasive Electrophysiological Data. Computational Intelligence and 
% Neuroscience, 2011; 2011:156869.


% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')
% additional comments by Rianne Haartsen are noted with '%RH note;'

% Adapted by Rianne Haartsen, PhD.; 03-2023 
% Birkbeck College, University of London

% This script is released under the GNU General Public License version 3.




%%

clear % clear matlab workspace
clc % clear matlab command window
addpath(genpath('xxx/eeglab2022.0/'));% enter the path of the EEGLAB folder in this line
eeglab nogui;% call eeglab to set up the plugin
addpath('xxx/fieldtrip-20180925');  % add fieldtrip path
ft_defaults
addpath(genpath('xxx/Preproc_HAPPILEE')) % any additional folders needed. 

InfoPath = 'xxx/Preproc_HAPPILEE/';

% create new folders
if exist(strcat(InfoPath,'/A_Power_data'),'dir') == 0
    mkdir(strcat(InfoPath,'/A_Power_data'))
end
if exist(strcat(InfoPath,'/B_FC_data'),'dir') == 0
    mkdir(strcat(InfoPath,'/B_FC_data'))
end
if exist(strcat(InfoPath,'/Report_figs'),'dir') == 0
    mkdir(strcat(InfoPath,'/Report_figs'))
end
if exist(strcat(InfoPath,'/processed_data'),'dir') == 0
    mkdir(strcat(InfoPath,'/processed_data'))
end

Pow_location = strcat(InfoPath,'/A_Power_data');
FC_location = strcat(InfoPath,'/B_FC_data');
output_location = InfoPath;

% get IDs for input
FilesSet = dir(fullfile(strcat('xxx/RawEEG_data/'), '*.set')); % folder with the raw data (to check the files that should've been preprocessed)
% plot data and topoplots; 1 = yes, 0 = no
PlotDataNow = 0; 
if PlotDataNow == 1
    Fig_location = strcat(InfoPath,'/Report_figs');
end


for subject = 1:height(FilesSet)
    
    CurrFileName = FilesSet(subject).name;
    ID = CurrFileName(1:(strfind(CurrFileName,'_')-1));
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, ID);

    ProcessedFileName = strcat(ID, '_social_nonsocial_videos_processed.set');
    
    if exist(strcat(InfoPath, '/4 - processed/', ProcessedFileName), 'file') == 2
    
        try
            %% Step 1: load preprocessed data
            Processedfullfilepath = strcat(InfoPath, '/4 - processed/', ProcessedFileName);
            EEG = pop_loadset(Processedfullfilepath);
            EEG = eeg_checkset(EEG);

            %% Step 2: Convert preprocessed eeglab data to fieldtrip data
            ft_EEG = eeglab2fieldtrip(EEG,'raw');

            % create trialinfo and sampleinfo fields in accordance with fieldtrip
                ft_EEG.trialinfo_eeglab = ft_EEG.trialinfo;
                % go through epochs
                trialinfo_cur = zeros(size(ft_EEG.trialinfo_eeglab,1),1);
                sampleinfo_cur = zeros(size(ft_EEG.trialinfo_eeglab,1),2);
                eps_start_urevent = zeros(size(ft_EEG.trialinfo_eeglab,1),1);
                for eps = 1:size(ft_EEG.trialinfo_eeglab,1)
                    trialinfo_cur(eps) = str2double(ft_EEG.trialinfo_eeglab.type{eps});
                    % check for consistency of event marker
                    eps_start_urevent(eps) = EEG.epoch(eps).eventurevent{1};
                    if ~isequal(trialinfo_cur(eps),str2double(EEG.urevent(eps_start_urevent(eps)).type))
                        error('Inconsistency ft and eeglab event types')
                    end
                    % get sampleinfo
                    begsample = EEG.urevent(eps_start_urevent(eps)).latency;
                    endsample = begsample + size(ft_EEG.time{eps},2)-1;
                    sampleinfo_cur(eps,:) = [begsample endsample];
                    clear ind_event begsample endsample
                end
                ft_EEG.trialinfo = trialinfo_cur;
                ft_EEG.sampleinfo = sampleinfo_cur;
                clear trialinfo_cur sampleinfo_cur eps eps_start_urevent

            % create time in accordance with fieldtrip
                ft_EEG.time_eeglab = ft_EEG.time;
                time = cell(1,size(ft_EEG.time_eeglab,2));
                % time vector
                time_end = ft_EEG.sampleinfo(end,2)*(1/ft_EEG.fsample);
                time_full = 0:(1/ft_EEG.fsample):time_end;
                for eps = 1:size(ft_EEG.time_eeglab,2)
                    time_cur = time_full(ft_EEG.sampleinfo(eps,1):ft_EEG.sampleinfo(eps,2));
                    time{1,eps} = time_cur;
                    clear time_cur
                end
                ft_EEG.time = time;
                clear time eps time_end time_full

            %% Step 3: Split into conditions
            if exist('ft_EEG','var')
                Neps_tot = size(ft_EEG.trialinfo,1);
            else
                Neps_tot = 0;
            end
                % social
                if sum(ft_EEG.trialinfo == 910) > 0
                    cfg                 = [];
                    cfg.trials          = ft_EEG.trialinfo == 910;
                    ft_EEG_soc = ft_selectdata(cfg, ft_EEG);
                end
                if exist('ft_EEG_soc','var')
                    Neps_soc = size(ft_EEG_soc.trial,2);
                else
                    Neps_soc = 0;
                end
                % toy
                if sum(ft_EEG.trialinfo == 911) > 0
                    cfg                 = [];
                    cfg.trials          = ft_EEG.trialinfo == 911;
                    ft_EEG_toy = ft_selectdata(cfg, ft_EEG);
                end
                if exist('ft_EEG_toy','var')
                    Neps_toy = size(ft_EEG_toy.trial,2);
                else
                    Neps_toy = 0;
                end

            %% Step 4: calculate power and dbWPLI if the numbers of epochs are equal between conditions
                % Parameters across all analyses
                % for power analysis
                    cfgpow                 = [];
                    cfgpow.output          = 'pow';
                    cfgpow.method          = 'mtmfft';
                    cfgpow.taper           = 'hanning';
                    cfgfft.tapsmofrq       = 1;
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
                if Neps_tot > 1
                    % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Power_all = ft_freqanalysis(cfgpow, ft_EEG);
                        Power_all.logpowspctrm = log(Power_all.powspctrm);
                        Power_all.logpowspctrm_avg = squeeze(mean(Power_all.logpowspctrm,1));
                        Power_all.logpowspctrm_avg_global = squeeze(mean(Power_all.logpowspctrm_avg,1));
                        Power_all.freq = round(Power_all.freq);
                    % dbWPLI analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % find segments with overlap
                        if size(ft_EEG.trialinfo,1) > 1
                            OVERLAP_all = zeros(size(ft_EEG.trialinfo,1), size(ft_EEG.trialinfo,1));
                            for i=1: size(ft_EEG.trialinfo,1)-1
                                if ft_EEG.sampleinfo(i,2) > ft_EEG.sampleinfo((i+1),1) % if overlaps
                                   OVERLAP_all (i, i+1)=1;  OVERLAP_all (i+1, i)=1; 
                                end
                            end
                            clear i
                        end
                        OVERLAP_all(find(eye(size(ft_EEG.trialinfo,1))==1)) = 1;
                        % FFT
                        FFT_all = ft_freqanalysis(cfgfft, ft_EEG);
                        % dbWPLI
                        [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(FFT_all.fourierspctrm, size(FFT_all.fourierspctrm,1), F, OVERLAP_all); %%dPLIbasic
                        FC_all = struct();
                        FC_all.PLI = PLI;
                        FC_all.WPLI = WPLI;
                        FC_all.ubPLI = ubPLI;
                        FC_all.dbWPLI = dbWPLI;
                        FC_all.Freqs_Hz = round(FFT_all.freq(F));
                        FC_all.Chan_label = FFT_all.label;
                        FC_all.overlap = OVERLAP_all;
                        clear PLI WPLI ubPLI dbWPLI 
                        % global dbWPLI
                        GlobdbWPLI_all = zeros(size(FC_all.dbWPLI,1),1);
                        for ff = 1:size(FC_all.dbWPLI,1)
                            FCmat = squeeze(FC_all.dbWPLI(ff,:,:));
                            Mat1 = tril(FCmat,-1);
                            Mat1(Mat1 == 0) = NaN;
                            GlobdbWPLI_all(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                            clear Mat1 FCmat
                        end
                        FC_all.dbWPLI_glob = GlobdbWPLI_all';  
                        clear ff GlobaldbWPLI_soc
                end


                %% For social trials
                if Neps_soc > 1
                    % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Power_soc = ft_freqanalysis(cfgpow, ft_EEG_soc);
                        Power_soc.logpowspctrm = log(Power_soc.powspctrm);
                        Power_soc.logpowspctrm_avg = squeeze(mean(Power_soc.logpowspctrm,1));
                        Power_soc.logpowspctrm_avg_global = squeeze(mean(Power_soc.logpowspctrm_avg,1));
                        Power_soc.freq = round(Power_soc.freq);
                    % dbWPLI analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % find segments with overlap
                        if size(ft_EEG_soc.trialinfo,1) > 1
                        OVERLAP_soc = zeros(size(ft_EEG_soc.trialinfo,1), size(ft_EEG_soc.trialinfo,1));
                            for i=1: size(ft_EEG_soc.trialinfo,1)-1
                                if ft_EEG_soc.sampleinfo(i,2) > ft_EEG_soc.sampleinfo((i+1),1) % if overlaps
                                   OVERLAP_soc (i, i+1)=1;  OVERLAP_soc (i+1, i)=1; 
                                end
                            end
                            clear i
                        end
                        OVERLAP_soc(find(eye(size(ft_EEG_soc.trialinfo,1))==1)) = 1;
                        % FFT
                        FFT_soc = ft_freqanalysis(cfgfft, ft_EEG_soc);
                        % dbWPLI
                        [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(FFT_soc.fourierspctrm, size(FFT_soc.fourierspctrm,1), F, OVERLAP_soc); %%dPLIbasic
                        FC_soc = struct();
                        FC_soc.PLI = PLI;
                        FC_soc.WPLI = WPLI;
                        FC_soc.ubPLI = ubPLI;
                        FC_soc.dbWPLI = dbWPLI;
                        FC_soc.Freqs_Hz = round(FFT_soc.freq(F));
                        FC_soc.Chan_label = FFT_soc.label;
                        FC_soc.overlap = OVERLAP_soc;
                        clear PLI WPLI ubPLI dbWPLI 
                        % global dbWPLI
                        GlobdbWPLI_soc = zeros(size(FC_soc.dbWPLI,1),1);
                        for ff = 1:size(FC_soc.dbWPLI,1)
                            FCmat = squeeze(FC_soc.dbWPLI(ff,:,:));
                            Mat1 = tril(FCmat,-1);
                            Mat1(Mat1 == 0) = NaN;
                            GlobdbWPLI_soc(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                            clear Mat1 FCmat
                        end
                        FC_soc.dbWPLI_glob = GlobdbWPLI_soc';    
                        clear ff GlobaldbWPLI_soc
                end

                %% For toy trials
                if Neps_toy > 1
                    % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Power_toy = ft_freqanalysis(cfgpow, ft_EEG_toy);
                        Power_toy.logpowspctrm = log(Power_toy.powspctrm);
                        Power_toy.logpowspctrm_avg = squeeze(mean(Power_toy.logpowspctrm,1));
                        Power_toy.logpowspctrm_avg_global = squeeze(mean(Power_toy.logpowspctrm_avg,1));
                        Power_toy.freq = round(Power_toy.freq);
                    % dbWPLI analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % find segments with overlap
                        if size(ft_EEG_toy.trialinfo,1) > 1
                        OVERLAP_toy = zeros(size(ft_EEG_toy.trialinfo,1), size(ft_EEG_toy.trialinfo,1));
                            for i=1: size(ft_EEG_toy.trialinfo,1)-1
                                if ft_EEG_toy.sampleinfo(i,2) > ft_EEG_toy.sampleinfo((i+1),1) % if overlaps
                                   OVERLAP_toy (i, i+1)=1;  OVERLAP_toy (i+1, i)=1; 
                                end
                            end
                            clear i
                        end
                    OVERLAP_toy(find(eye(size(ft_EEG_toy.trialinfo,1))==1)) = 1;
                    % FFT
                        FFT_toy = ft_freqanalysis(cfgfft, ft_EEG_toy);
                        % dbWPLI
                        [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(FFT_toy.fourierspctrm, size(FFT_toy.fourierspctrm,1), F, OVERLAP_toy); %%dPLIbasic
                        FC_toy = struct();
                        FC_toy.PLI = PLI;
                        FC_toy.WPLI = WPLI;
                        FC_toy.ubPLI = ubPLI;
                        FC_toy.dbWPLI = dbWPLI;
                        FC_toy.Freqs_Hz = round(FFT_toy.freq(F));
                        FC_toy.Chan_label = FFT_toy.label;
                        FC_toy.overlap = OVERLAP_toy;
                        clear PLI WPLI ubPLI dbWPLI 
                        % global dbWPLI
                        GlobdbWPLI_toy = zeros(size(FC_toy.dbWPLI,1),1);
                        for ff = 1:size(FC_toy.dbWPLI,1)
                            FCmat = squeeze(FC_toy.dbWPLI(ff,:,:));
                            Mat1 = tril(FCmat,-1);
                            Mat1(Mat1 == 0) = NaN;
                            GlobdbWPLI_toy(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                            clear Mat1 FCmat
                        end
                        FC_toy.dbWPLI_glob = GlobdbWPLI_toy';
                        clear ff GlobaldbWPLI_toy
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
            if ~exist('ft_EEG_orig','var')
                ft_EEG_orig = [];
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
            if ~exist('Neps_tot','var')
                Neps_tot = 0;
            end
            if ~exist('Neps_soc','var')
                Neps_soc = 0;
            end
            if ~exist('Neps_toy','var')
                Neps_toy = 0;
            end

        %% Plot the data
        if PlotDataNow == 1

            if Neps_tot > 0
                Spectra_Fig = figure('Position',[78 286 1437 669]);

                % Power spectrum across all channels, and plot topoplots for alpha and
                % theta log power %%%
                figure(Spectra_Fig)
                subplot(4,5,[1 2]) % power spectrum across all channels
                    h1 = plot(Power_all.freq, Power_all.logpowspctrm_avg_global,'Color',[0 0 .5],'LineWidth',2); %all
                    hold on 
                    if Neps_soc > 0
                        h2 = plot(Power_soc.freq, Power_soc.logpowspctrm_avg_global,'Color',[.4 0 .6],'LineWidth',2,'LineStyle',':'); %soc
                    end
                    if Neps_toy > 0
                        h3 = plot(Power_toy.freq, Power_toy.logpowspctrm_avg_global,'Color',[.9 .5 .1],'LineWidth',2,'LineStyle','-.'); %toy
                    end
                    xlim([1 32]);
                    props = gcf; 
                    if Neps_soc == 0
                        h2 = plot(Power_all.freq,ones(1,size(Power_all.freq,2))* props.CurrentAxes.YLim(1,1),'Color',[.4 0 .6],'LineWidth',2,'LineStyle',':'); %soc
                    end
                    if Neps_toy == 0
                        h3 = plot(Power_all.freq,ones(1,size(Power_all.freq,2))* props.CurrentAxes.YLim(1,1),'Color',[.9 .5 .1],'LineWidth',2,'LineStyle','-.'); %toy
                    end
                    % add lines for fois
                    h4 = plot([3.8 3.8],props.CurrentAxes.YLim,'k-');
                    h5 = plot([6.2 6.2],props.CurrentAxes.YLim,'k-');
                    h6 = plot([6.8 6.8],props.CurrentAxes.YLim,'k:');
                    h7 = plot([12.2 12.2],props.CurrentAxes.YLim,'k:');
                    title('FFT power spectrumm (all Chs)');
                    legend([h1 h2 h3],{'All','Soc', 'Toy'});
                    xlabel('Frequency (Hz)'); ylabel('Log power (\muV^2/Hz)');

                    clear h1 h2 h3 h4 h5 h6 h7

                % theta powerplots
                if Neps_soc > 0 && Neps_toy > 0
                    Th_ind_low =  find(Power_soc.freq == 4);
                    Th_ind_high = find(Power_soc.freq == 6);
                    Power_soc.logpow_theta = mean(Power_soc.logpowspctrm_avg(:,Th_ind_low:Th_ind_high),2);
                    clear Th_ind_low Th_ind_high
                    Th_ind_low =  find(Power_toy.freq == 4);
                    Th_ind_high = find(Power_toy.freq == 6);
                    Power_toy.logpow_theta = mean(Power_toy.logpowspctrm_avg(:,Th_ind_low:Th_ind_high),2);
                    clear Th_ind_low Th_ind_high
                minPow = min([Power_soc.logpow_theta; Power_toy.logpow_theta]);
                maxPow = max([Power_soc.logpow_theta; Power_toy.logpow_theta]);
                zLimUB = ceil(maxPow*10)/10; zLimLB = floor(minPow*10)/10;
                figure(Spectra_Fig)
                subplot(4,5,3) 
                        ttt.avg=Power_soc.logpow_theta; % each row is 1 channel
                        ttt.label=Power_soc.label;
                        ttt.fsample=1000;
                        ttt.time=1;
                        ttt.dimord= 'subj_chan_time';
                        % specify parameters for plotting
                        cfg = [];   
                        cfg.parameter           = 'avg';                     
                        cfg.layout              = 'elec1005.lay';     
                        cfg.colorbar            = 'EastOutside';
                        cfg.comment             = 'no';
                        cfg.zlim                = [zLimLB zLimUB];
                        ft_topoplotER(cfg,ttt)  
                        title(sprintf('Theta [4-6Hz] Soc'))
                subplot(4,5,4) 
                        ttt.avg=Power_toy.logpow_theta; % each row is 1 channel
                        ttt.label=Power_toy.label;
                        ttt.fsample=1000;
                        ttt.time=1;
                        ttt.dimord= 'subj_chan_time';
                        % specify parameters for plotting
                        cfg = [];   
                        cfg.parameter           = 'avg';                     
                        cfg.layout              = 'elec1005.lay';     
                        cfg.colorbar            = 'EastOutside';
                        cfg.comment             = 'no';
                        cfg.zlim                = [zLimLB zLimUB];
                        ft_topoplotER(cfg,ttt)  
                        title(sprintf('Theta [4-6Hz] Toy') )
                        clear zLimLB zLimUB
                subplot(4,5,5) 
                        ttt.avg=Power_soc.logpow_theta - Power_toy.logpow_theta; % each row is 1 channel
                        ttt.label= Power_soc.label;
                        ttt.fsample=1000;
                        ttt.time=1;
                        ttt.dimord= 'subj_chan_time';
                        % specify parameters for plotting
                        cfg = [];   
                        cfg.parameter           = 'avg';                     
                        cfg.layout              = 'elec1005.lay';     
                        cfg.colorbar            = 'EastOutside';
                        cfg.comment             = 'no';
                        ft_topoplotER(cfg,ttt)  
                        title(sprintf('Theta [4-6Hz] Diff (S - T)') )
                        clear zLimLB zLimUB

                elseif Neps_soc == 0 && Neps_toy > 0 % no soc trials but toy trials

                        Th_ind_low =  find(Power_toy.freq == 4);
                        Th_ind_high = find(Power_toy.freq == 6);
                        Power_toy.logpow_theta = mean(Power_toy.logpowspctrm_avg(:,Th_ind_low:Th_ind_high),2);
                        clear Th_ind_low Th_ind_high
                    minPow = min(Power_toy.logpow_theta);
                    maxPow = max(Power_toy.logpow_theta);
                    zLimUB = ceil(maxPow*10)/10; zLimLB = floor(minPow*10)/10;
                    figure(Spectra_Fig)
                    subplot(4,5,3) 
                            plot(1)
                            txt = 'No social trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                    subplot(4,5,4) 
                            ttt.avg=Power_toy.logpow_theta; % each row is 1 channel
                            ttt.label=Power_toy.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Theta [4-6Hz] Toy') )
                            clear zLimLB zLimUB
                    subplot(4,5,5) 
                            plot(1)
                            txt = 'No social trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);


                elseif Neps_soc > 0 && Neps_toy == 0 % soc trials but no toy trials
                        Th_ind_low =  find(Power_soc.freq == 4);
                        Th_ind_high = find(Power_soc.freq == 6);
                        Power_soc.logpow_theta = mean(Power_soc.logpowspctrm_avg(:,Th_ind_low:Th_ind_high),2);
                        clear Th_ind_low Th_ind_high
                    minPow = min(Power_soc.logpow_theta);
                    maxPow = max(Power_soc.logpow_theta);
                    zLimUB = ceil(maxPow*10)/10; zLimLB = floor(minPow*10)/10;
                    figure(Spectra_Fig)
                    subplot(4,5,3) 
                            ttt.avg=Power_soc.logpow_theta; % each row is 1 channel
                            ttt.label=Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Theta [4-6Hz] Soc'))
                    subplot(4,5,4) 
                            plot(1)
                            txt = 'No toy trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                    subplot(4,5,5) 
                            plot(1)
                            txt = 'No toy trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);

                else 
                    warning('Unexpected situation for N trials')

                end

                % alpha topoplots
                if Neps_soc > 0 && Neps_toy > 0
                        Al_ind_low =  find(Power_soc.freq == 7);
                        Al_ind_high = find(Power_soc.freq == 12);
                        Power_soc.logpow_alpha = mean(Power_soc.logpowspctrm_avg(:,Al_ind_low:Al_ind_high),2);
                        clear Th_ind_low Th_ind_high
                        Al_ind_low =  find(Power_toy.freq == 7);
                        Al_ind_high = find(Power_toy.freq == 12);
                        Power_toy.logpow_alpha = mean(Power_toy.logpowspctrm_avg(:,Al_ind_low:Al_ind_high),2);
                        clear Th_ind_low Th_ind_high           
                    minPow = min([Power_soc.logpow_alpha; Power_toy.logpow_alpha]);
                    maxPow = max([Power_soc.logpow_alpha; Power_toy.logpow_alpha]);
                    zLimUB = ceil(maxPow*10)/10; zLimLB = floor(minPow*10)/10;
                    subplot(4,5,8) 
                            ttt.avg=Power_soc.logpow_alpha; % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify paramters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Soc') )
                    subplot(4,5,9) 
                            ttt.avg=Power_toy.logpow_alpha; % each row is 1 channel
                            ttt.label= Power_toy.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify paramters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Toy') )       
                    subplot(4,5,10) 
                            ttt.avg=Power_soc.logpow_alpha - Power_toy.logpow_alpha; % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Diff (S - T)') )
                            clear zLimLB zLimUB
                elseif Neps_soc > 0 && Neps_toy == 0 % no toy trials
                            Al_ind_low =  find(Power_soc.freq == 7);
                            Al_ind_high = find(Power_soc.freq == 12);
                            Power_soc.logpow_alpha = mean(Power_soc.logpowspctrm_avg(:,Al_ind_low:Al_ind_high),2);
                            clear Th_ind_low Th_ind_high   
                        minPow = min(Power_soc.logpow_alpha);
                        maxPow = max(Power_soc.logpow_alpha);
                        zLimUB = ceil(maxPow*10)/10; zLimLB = floor(minPow*10)/10;
                        subplot(4,5,8) 
                                ttt.avg=Power_soc.logpow_alpha; % each row is 1 channel
                                ttt.label= Power_soc.label;
                                ttt.fsample=1000;
                                ttt.time=1;
                                ttt.dimord= 'subj_chan_time';
                                % specify paramters for plotting
                                cfg = [];   
                                cfg.parameter           = 'avg';                     
                                cfg.layout              = 'elec1005.lay';     
                                cfg.colorbar            = 'EastOutside';
                                cfg.comment             = 'no';
                                cfg.zlim                = [zLimLB zLimUB];
                                ft_topoplotER(cfg,ttt)  
                                title(sprintf('Alpha [7-12Hz] Soc') )
                        subplot(4,5,9) 
                                plot(1)
                                txt = 'No toy trials';
                                text(0.8,1,txt, 'FontSize', 16);
                                xlim([0 4]); ylim([0 2])
                                set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                                set(gca,'ytick',[]); set(gca,'xtick',[]);
                        subplot(4,5,10) 
                                plot(1)
                                txt = 'No toy trials';
                                text(0.8,1,txt, 'FontSize', 16);
                                xlim([0 4]); ylim([0 2])
                                set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                                set(gca,'ytick',[]); set(gca,'xtick',[]);

                elseif Neps_soc == 0 && Neps_toy > 0 % no social trials
                            Al_ind_low =  find(Power_toy.freq == 7);
                            Al_ind_high = find(Power_toy.freq == 12);
                            Power_toy.logpow_alpha = mean(Power_toy.logpowspctrm_avg(:,Al_ind_low:Al_ind_high),2);
                            clear Th_ind_low Th_ind_high           
                        minPow = min(Power_toy.logpow_alpha);
                        maxPow = max(Power_toy.logpow_alpha);
                        zLimUB = ceil(maxPow*10)/10; zLimLB = floor(minPow*10)/10;
                        subplot(4,5,8) 
                                plot(1)
                                txt = 'No social trials';
                                text(0.8,1,txt, 'FontSize', 16);
                                xlim([0 4]); ylim([0 2])
                                set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                                set(gca,'ytick',[]); set(gca,'xtick',[]);
                        subplot(4,5,9) 
                                ttt.avg=Power_toy.logpow_alpha; % each row is 1 channel
                                ttt.label= Power_toy.label;
                                ttt.fsample=1000;
                                ttt.time=1;
                                ttt.dimord= 'subj_chan_time';
                                % specify paramters for plotting
                                cfg = [];   
                                cfg.parameter           = 'avg';                     
                                cfg.layout              = 'elec1005.lay';     
                                cfg.colorbar            = 'EastOutside';
                                cfg.comment             = 'no';
                                cfg.zlim                = [zLimLB zLimUB];
                                ft_topoplotER(cfg,ttt)  
                                title(sprintf('Alpha [7-12Hz] Toy') )       
                        subplot(4,5,10) 
                                plot(1)
                                txt = 'No trials';
                                text(0.8,1,txt, 'FontSize', 16);
                                xlim([0 4]); ylim([0 2])
                                set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                                set(gca,'ytick',[]); set(gca,'xtick',[]);
                else 
                    warning('Unexpected situation for N trials')  
                end

                % plot connectivity findings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % all
                GlobFC_freqs_all = zeros(size(FC_all.dbWPLI,1),1);
                for ff = 1:size(FC_all.dbWPLI,1)
                    FCmat = squeeze(FC_all.dbWPLI(ff,:,:));
                    Mat1 = tril(FCmat,-1);
                    Mat1(Mat1 == 0) = NaN;
                    GlobFC_freqs_all(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                    clear Mat1
                end
                FC_all.GlobFC_freqs_all = GlobFC_freqs_all;
                % soc
                if Neps_soc > 0
                    GlobFC_freqs_soc = zeros(size(FC_soc.dbWPLI,1),1);
                    for ff = 1:size(FC_soc.dbWPLI,1)
                        FCmat = squeeze(FC_soc.dbWPLI(ff,:,:));
                        Mat1 = tril(FCmat,-1);
                        Mat1(Mat1 == 0) = NaN;
                        GlobFC_freqs_soc(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                        clear Mat1
                    end
                    FC_soc.GlobFC_freqs_all = GlobFC_freqs_soc;
                end
                % toy
                if Neps_toy > 0
                    GlobFC_freqs_toy = zeros(size(FC_toy.dbWPLI,1),1);
                    for ff = 1:size(FC_toy.dbWPLI,1)
                        FCmat = squeeze(FC_toy.dbWPLI(ff,:,:));
                        Mat1 = tril(FCmat,-1);
                        Mat1(Mat1 == 0) = NaN;
                        GlobFC_freqs_toy(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                        clear Mat1
                    end
                    FC_toy.GlobFC_freqs_all = GlobFC_freqs_toy;
                end

                subplot(4,5,[11 12]) % FC spectrum across all channels
                    h1 = plot(FC_all.Freqs_Hz', GlobFC_freqs_all,'Color',[0 0 .5],'LineWidth',2);
                    hold on 
                    if Neps_soc > 0
                        h2 = plot(FC_soc.Freqs_Hz', GlobFC_freqs_soc,'Color',[.4 0 .6],'LineWidth',2,'LineStyle',':');
                    end
                    if Neps_toy > 0
                        h3 = plot(FC_toy.Freqs_Hz', GlobFC_freqs_toy,'Color',[.9 .5 .1],'LineWidth',2,'LineStyle','-.');
                    end
                    xlim([1 32]);
                    props = gcf; 
                    if Neps_soc == 0
                        h2 = plot(FC_all.Freqs_Hz',ones(1,size(FC_all.Freqs_Hz,2))* props.CurrentAxes.YLim(1,1),'Color',[.4 0 .6],'LineWidth',2,'LineStyle',':'); %soc
                    end
                    if Neps_toy == 0
                        h3 = plot(FC_all.Freqs_Hz',ones(1,size(FC_all.Freqs_Hz,2))* props.CurrentAxes.YLim(1,1),'Color',[.9 .5 .1],'LineWidth',2,'LineStyle','-.'); %toy
                    end
                    % add lines for fois
                    h4 = plot([3.8 3.8],props.CurrentAxes.YLim,'k-');
                    h5 = plot([6.2 6.2],props.CurrentAxes.YLim,'k-');
                    h6 = plot([6.8 6.8],props.CurrentAxes.YLim,'k:');
                    h7 = plot([12.2 12.2],props.CurrentAxes.YLim,'k:');
                    title('FC dbWPLI spectrumm (all ChsxChs)');
                    legend([h1 h2 h3],{'All','Soc', 'Toy'});
                    xlabel('Frequency (Hz)'); ylabel('dbWPLI');


                % average for theta and alpha FC
                if Neps_soc > 0 && Neps_toy > 0
                    % theta
                    Th_ind_low =  find(FC_all.Freqs_Hz == 4);
                    Th_ind_high = find(FC_all.Freqs_Hz == 6);
                    FC_all.Theta_FCmat_avg = squeeze(mean(FC_all.dbWPLI([(Th_ind_low):(Th_ind_high)],:,:),1));
                    FC_soc.Theta_FCmat_avg = squeeze(mean(FC_soc.dbWPLI([(Th_ind_low):(Th_ind_high)],:,:),1));
                    FC_toy.Theta_FCmat_avg = squeeze(mean(FC_toy.dbWPLI([(Th_ind_low):(Th_ind_high)],:,:),1));
                    clear Th_ind_low Th_ind_high
                    % alpha
                    Al_ind_low =  find(FC_all.Freqs_Hz == 7);
                    Al_ind_high = find(FC_all.Freqs_Hz == 12);
                    FC_all.Alpha_FCmat_avg = squeeze(mean(FC_all.dbWPLI([(Al_ind_low):(Al_ind_high)],:,:),1));
                    FC_soc.Alpha_FCmat_avg = squeeze(mean(FC_soc.dbWPLI([(Al_ind_low):(Al_ind_high)],:,:),1));
                    FC_toy.Alpha_FCmat_avg = squeeze(mean(FC_toy.dbWPLI([(Al_ind_low):(Al_ind_high)],:,:),1));
                    clear Al_ind_low Al_ind_high
                elseif Neps_soc > 0 && Neps_toy == 0 % no toy trials
                    % theta
                    Th_ind_low =  find(FC_all.Freqs_Hz == 4);
                    Th_ind_high = find(FC_all.Freqs_Hz == 6);
                    FC_all.Theta_FCmat_avg = squeeze(mean(FC_all.dbWPLI([(Th_ind_low):(Th_ind_high)],:,:),1));
                    FC_soc.Theta_FCmat_avg = squeeze(mean(FC_soc.dbWPLI([(Th_ind_low):(Th_ind_high)],:,:),1));
                    clear Th_ind_low Th_ind_high
                    % alpha
                    Al_ind_low =  find(FC_all.Freqs_Hz == 7);
                    Al_ind_high = find(FC_all.Freqs_Hz == 12);
                    FC_all.Alpha_FCmat_avg = squeeze(mean(FC_all.dbWPLI([(Al_ind_low):(Al_ind_high)],:,:),1));
                    FC_soc.Alpha_FCmat_avg = squeeze(mean(FC_soc.dbWPLI([(Al_ind_low):(Al_ind_high)],:,:),1));
                    clear Al_ind_low Al_ind_high
                elseif Neps_soc == 0 && Neps_toy > 0 % no soc trials
                    % theta
                    Th_ind_low =  find(FC_all.Freqs_Hz == 4);
                    Th_ind_high = find(FC_all.Freqs_Hz == 6);
                    FC_all.Theta_FCmat_avg = squeeze(mean(FC_all.dbWPLI([(Th_ind_low):(Th_ind_high)],:,:),1));
                    FC_toy.Theta_FCmat_avg = squeeze(mean(FC_toy.dbWPLI([(Th_ind_low):(Th_ind_high)],:,:),1));
                    clear Th_ind_low Th_ind_high
                    % alpha
                    Al_ind_low =  find(FC_all.Freqs_Hz == 7);
                    Al_ind_high = find(FC_all.Freqs_Hz == 12);
                    FC_all.Alpha_FCmat_avg = squeeze(mean(FC_all.dbWPLI([(Al_ind_low):(Al_ind_high)],:,:),1));
                    FC_toy.Alpha_FCmat_avg = squeeze(mean(FC_toy.dbWPLI([(Al_ind_low):(Al_ind_high)],:,:),1));
                    clear Al_ind_low Al_ind_high
                end

                % theta topoplots
                if Neps_soc > 0 && Neps_toy > 0 
                    ThMat_S = FC_soc.Theta_FCmat_avg; ThMat_S(ThMat_S == 0) = NaN;
                    ThMat_T = FC_toy.Theta_FCmat_avg; ThMat_T(ThMat_T == 0) = NaN;
                    minPow = min([mean(ThMat_S, 2, 'omitnan'); mean(ThMat_T, 2, 'omitnan')],[],'all','omitnan');
                    maxPow = max([mean(ThMat_S, 2, 'omitnan'); mean(ThMat_T, 2, 'omitnan')],[],'all','omitnan');
                    zLimUB = ceil(maxPow*100)/100; zLimLB = floor(minPow*100)/100;
                    subplot(4,5,13) 
                            ttt.avg=mean(ThMat_S, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Theta [4-6Hz] Soc') )
                    subplot(4,5,14) 
                            ttt.avg=mean(ThMat_T, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Theta [4-6Hz] Toy') )
                            clear zLimLB zLimUB
                    subplot(4,5,15) 
                            ttt.avg=mean(ThMat_S, 2, 'omitnan') - mean(ThMat_T, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Theta [4-6Hz] Diff (S - T)') )
                            clear zLimLB zLimUB

                            clear ThMat_S ThMat_T
                elseif Neps_soc > 0 && Neps_toy == 0 % no toy trials
                    ThMat_S = FC_soc.Theta_FCmat_avg; ThMat_S(ThMat_S == 0) = NaN;
                    minPow = min(mean(ThMat_S, 2, 'omitnan'),[],'all','omitnan');
                    maxPow = max(mean(ThMat_S, 2, 'omitnan'),[],'all','omitnan');
                    zLimUB = ceil(maxPow*100)/100; zLimLB = floor(minPow*100)/100;
                    subplot(4,5,13) 
                            ttt.avg=mean(ThMat_S, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Theta [4-6Hz] Soc') )
                    subplot(4,5,14) 
                            plot(1)
                            txt = 'No toy trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                    subplot(4,5,15) 
                           plot(1)
                            txt = 'No trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);

                            clear ThMat_S 
                elseif Neps_soc == 0 && Neps_toy > 0 % no soc trials
                    ThMat_T = FC_toy.Theta_FCmat_avg; ThMat_T(ThMat_T == 0) = NaN;
                    minPow = min(mean(ThMat_T, 2, 'omitnan'),[],'all','omitnan');
                    maxPow = max(mean(ThMat_T, 2, 'omitnan'),[],'all','omitnan');
                    zLimUB = ceil(maxPow*100)/100; zLimLB = floor(minPow*100)/100;
                    subplot(4,5,13) 
                            plot(1)
                            txt = 'No social trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                    subplot(4,5,14) 
                            ttt.avg=mean(ThMat_T, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Theta [4-6Hz] Toy') )
                            clear zLimLB zLimUB
                    subplot(4,5,15) 
                            plot(1)
                            txt = 'No trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);

                            clear ThMat_S ThMat_T
                else
                    warning('Unexpected situation for N trials')  
                end

                % alpha topoplots
                if Neps_soc > 0 && Neps_toy > 0
                    AlMat_S = FC_soc.Alpha_FCmat_avg; AlMat_S(AlMat_S == 0) = NaN;
                    AlMat_T = FC_toy.Alpha_FCmat_avg; AlMat_T(AlMat_T == 0) = NaN;
                    minPow = min([mean(AlMat_S, 2, 'omitnan'); mean(AlMat_T, 2, 'omitnan')],[],'all','omitnan');
                    maxPow = max([mean(AlMat_S, 2, 'omitnan'); mean(AlMat_T, 2, 'omitnan')],[],'all','omitnan');
                    zLimUB = ceil(maxPow*100)/100; zLimLB = floor(minPow*100)/100;
                    subplot(4,5,18) 
                            ttt.avg=mean(AlMat_S, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify paramters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Soc') )
                    subplot(4,5,19) 
                            ttt.avg=mean(AlMat_T, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_toy.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify paramters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Toy') )       
                    subplot(4,5,20) 
                            ttt.avg=mean(AlMat_S, 2, 'omitnan') - mean(AlMat_T, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify parameters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Diff (S - T)') )
                            clear zLimLB zLimUB
                elseif Neps_soc > 0 && Neps_toy == 0 % no toy trials
                    AlMat_S = FC_soc.Alpha_FCmat_avg; AlMat_S(AlMat_S == 0) = NaN;
                    minPow = min(mean(AlMat_S, 2, 'omitnan') ,[],'all','omitnan');
                    maxPow = max(mean(AlMat_S, 2, 'omitnan'),[],'all','omitnan');
                    zLimUB = ceil(maxPow*100)/100; zLimLB = floor(minPow*100)/100;
                    subplot(4,5,18) 
                            ttt.avg=mean(AlMat_S, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_soc.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify paramters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Soc') )
                    subplot(4,5,19) 
                            plot(1)
                            txt = 'No toy trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                    subplot(4,5,20) 
                            plot(1)
                            txt = 'No trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                            clear zLimLB zLimUB
                elseif Neps_soc == 0 && Neps_toy > 0 % no social trials
                    AlMat_T = FC_toy.Alpha_FCmat_avg; AlMat_T(AlMat_T == 0) = NaN;
                    minPow = min(mean(AlMat_T, 2, 'omitnan'),[],'all','omitnan');
                    maxPow = max(mean(AlMat_T, 2, 'omitnan'),[],'all','omitnan');
                    zLimUB = ceil(maxPow*100)/100; zLimLB = floor(minPow*100)/100;
                    subplot(4,5,18) 
                            plot(1)
                            txt = 'No social trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                    subplot(4,5,19) 
                            ttt.avg=mean(AlMat_T, 2, 'omitnan'); % each row is 1 channel
                            ttt.label= Power_toy.label;
                            ttt.fsample=1000;
                            ttt.time=1;
                            ttt.dimord= 'subj_chan_time';
                            % specify paramters for plotting
                            cfg = [];   
                            cfg.parameter           = 'avg';                     
                            cfg.layout              = 'elec1005.lay';     
                            cfg.colorbar            = 'EastOutside';
                            cfg.comment             = 'no';
                            cfg.zlim                = [zLimLB zLimUB];
                            ft_topoplotER(cfg,ttt)  
                            title(sprintf('Alpha [7-12Hz] Toy') )       
                    subplot(4,5,20) 
                            plot(1)
                            txt = 'No trials';
                            text(0.8,1,txt, 'FontSize', 16);
                            xlim([0 4]); ylim([0 2])
                            set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                            set(gca,'ytick',[]); set(gca,'xtick',[]);
                            clear zLimLB zLimUB
                else
                    warning('Unexpected situation for N trials') 
                end
                clear AlMat_S AlMat_T


                % add info on N trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                subplot(4,5,[6 7])
                txt = 'N trials:';
                text(0.2,1,txt, 'FontSize', 16);
                hold on
                txt2 = sprintf('Social: %i',Neps_soc);
                text(0.5,.8,txt2, 'FontSize', 16);
                txt3 = sprintf('Toy: %i', Neps_toy);
                text(0.5,.6,txt3, 'FontSize', 16);
                txt4 = sprintf('Total: %i', Neps_tot);
                text(0.2,.4,txt4, 'FontSize', 16);
                set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
                set(gca,'ytick',[]); set(gca,'xtick',[]);
                xlim([0 3])
                ylim([0 1.5])


                sgtitle({'Subject: ', convertCharsToStrings(ID)})

                % clear up
                clear h1 h2 h3 h4 h5 h6 h7 
                clear txt txt2 txt3 txt4 txt5 ttt
                clear GlobFC_freqs_all GlobFC_freqs_soc GlobFC_freqs_toy
                clear minPow maxPow

            end % end if there are no trials
        end % end plot data


        %% save data and add to current subject to report table 
        % keep track of saved data location and subject preprocessed
        % Info
        if exist('EEG','var') && isfield(EEG,'euaims')
            Site{subject} = EEG.euaims.site;
        else
            Site{subject} = 'NA';
        end
        DateCur{subject} = datetime('now');
        Ft_file{subject} = strcat(output_location, filesep, 'processed_data', filesep, ID, '_fieldtrip_data.mat');
        save(Ft_file{subject},'ft_EEG','ft_EEG_soc','ft_EEG_toy','ft_EEG_orig')
        % Power
        Power_file{subject} = strcat(Pow_location, filesep, ID, '_Power_data.mat');
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
        FC_file{subject} = strcat(FC_location, filesep, ID,  '_FC_data.mat');
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

        % figure
        if exist('Spectra_Fig','var')
                % Save the figure
                Fig_name = strcat(Fig_location, filesep, ID,  '_Report_figure.png');
                saveas(Spectra_Fig, Fig_name)
                clear Spectra_Fig
                close
        end


        % add to tracking table
        cd(output_location)
        if exist(strcat(InfoPath, '/HAPPILEE_MeasuresOfInterest_report.mat'),'file')
            load HAPPILEE_MeasuresOfInterest_report.mat
        end
        data_table_newrow = table({ID}, Site(1,subject), DateCur(1,subject), Ft_file(1,subject), {Neps_tot}, {Neps_soc}, {Neps_toy}, ...
            {Freqs}, Power_file(1,subject), {Power_global_all},{Power_global_soc},{Power_global_toy},...
            FC_file(1,subject), {dbWPLI_global_all},{dbWPLI_global_soc},{dbWPLI_global_toy});
        data_table_newrow.Properties.VariableNames={'ID','Site', 'Date','Fieldtrip_file','Neps_tot','Neps_soc','Neps_toy', ...
            'Freqs','Power_path','GlobPow_all','GlobPow_soc','GlobPow_toy',...
            'FC_path','GlobFC_all','GlobFC_soc','GlobFC_toy'};
        HAPPILEE_MOI_table(subject,:) = data_table_newrow;
        save('HAPPILEE_MeasuresOfInterest_report.mat','HAPPILEE_MOI_table');
        clear HAPPILEE_MOI_table
    else
        % Info
        Site{subject} = 'NA';
        DateCur{subject} = datetime('now');
        Ft_file{subject} = 'NA';
        Neps_tot = 0;
        Neps_soc = 0; 
        Neps_toy = 0;
        % Power
        Power_file{subject} = 'NA';
            Power_global_all = NaN(1,32);
            Freqs = NaN(1,32);
            Power_global_soc = NaN(1,32);
            Power_global_toy = NaN(1,32);
        % FFT and FC
        FC_file{subject} = 'NA';
            dbWPLI_global_all = NaN(1,32);
            dbWPLI_global_soc = NaN(1,32);
            dbWPLI_global_toy = NaN(1,32);

        % add to tracking table
        cd(output_location)
        if exist(strcat(InfoPath, '/HAPPILEE_MeasuresOfInterest_report.mat'),'file')
            load HAPPILEE_MeasuresOfInterest_report.mat
        end
        data_table_newrow = table({ID}, Site(1,subject), DateCur(1,subject), Ft_file(1,subject), {Neps_tot}, {Neps_soc}, {Neps_toy}, ...
            {Freqs}, Power_file(1,subject), {Power_global_all},{Power_global_soc},{Power_global_toy},...
            FC_file(1,subject), {dbWPLI_global_all},{dbWPLI_global_soc},{dbWPLI_global_toy});
        data_table_newrow.Properties.VariableNames={'ID','Site', 'Date','Fieldtrip_file','Neps_tot','Neps_soc','Neps_toy', ...
            'Freqs','Power_path','GlobPow_all','GlobPow_soc','GlobPow_toy',...
            'FC_path','GlobFC_all','GlobFC_soc','GlobFC_toy'};
        HAPPILEE_MOI_table(subject,:) = data_table_newrow;
        save('HAPPILEE_MeasuresOfInterest_report.mat','HAPPILEE_MOI_table');
        clear HAPPILEE_MOI_table
        
    end 
    
    % clear up variables
    clear ft_EEG ft_EEG_soc ft_EEG_toy ft_EEG_orig ft_EEG_soc_new ft_EEG_toy_new
    clear Power_all Power_soc Power_toy 
    clear FFT_all FFT_orig FFT_soc FFT_toy FC_all FC_soc FC_toy 
    clear Neps_tot Neps_soc Neps_toy
    clear OVERLAP_all OVERLAP_soc OVERLAP_toy Freqs ff F FCmat 
    clear GlobdbWPLI_all GlobdbWPLI_soc GlobdbWPLI_toy
    clear dbWPLI_global_all dbWPLI_global_soc dbWPLI_global_toy
    clear Power_global_all Power_global_soc Power_global_toy
    clear Site FC_file Pow_file data_table_newrow
    clear cfg cfgfft cfgpow CurrID DateCur
    clear EEG Ft_file
    clear filename_path_preproc CurrFileName ID
    
    close all
    

end

