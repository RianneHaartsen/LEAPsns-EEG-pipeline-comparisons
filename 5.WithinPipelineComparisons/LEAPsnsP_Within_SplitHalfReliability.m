%% Analysis of social and non-social videos in LEAP: Comparison within pipelines

% This script calculates split-half reliability for the EEG metrics of 
% interest for each of the 8 pipelines: Manual, MADE, BOND-MADE, HAPPEv1, 
% HAPPEv4, MADE-BONDld, HAPPILEE, and miniMADE.

% This script reads in the datasets with power from the different
% pipelines. Then it takes the power across all trials and creates dataset
% A and B with alternating trials. This method is repeated for the social
% and toy trials. Next, power is averaged across the trials in the datasets
% A and B and the power spectra are saved.

% In the second step, ICCs for the split half are calculated for the
% different frequency bands (delta, theta, alpha, beta) and pipelines. 

% The following metrics are extracted:

% 1) Whole brain power values: for delta, theta, alpha, and beta frequencies
% 2) Whole brain connectivity values: for delta, theta, alpha, and beta frequencies

% Data is exported to csv files for further analysis and visualisation in Python

% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')

% Created by Rianne Haartsen, PhD.; 04-2023 
% Birkbeck College, University of London

% This script is released under the GNU General Public License version 3.



%% Power for subset A and subset B
% 1) Roche 1 manual pipeline
cd xxx/Preproc_Manual
load('Manual_MeasuresOfInterest_report.mat')
Power_folder = 'xxx/Preproc_Manual/A_Power_data';
for ii = 1:height(Manual_MOI_table)
    ID = Manual_MOI_table.ID{ii};
    load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
    Freqs = 1:1:32;
    % across all trials
    if ~isempty(Power_all)
        Ntrls_tot = size(Power_all.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_all = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_all = nan(1,32);
        PowB_all = nan(1,32);
        Nsplit_all = 0;  
    end
    
    % social trials
    if ~isempty(Power_soc)
        Ntrls_tot = size(Power_soc.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_soc = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_soc = nan(1,32);
        PowB_soc = nan(1,32);
        Nsplit_soc = 0;  
    end
    
    % toy trials
    if ~isempty(Power_toy)
        Ntrls_tot = size(Power_toy.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_toy = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_toy = nan(1,32);
        PowB_toy = nan(1,32);
        Nsplit_toy = 0;  
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            Manual_SplitHalf = data_table_newrow;
            cd xxx/Preproc_Manual
            save('Manual_SplitHalf_reliability.mat','Manual_SplitHalf')
            clear Manual_SplitHalf data_table_newrow 
        else
            cd xxx/Preproc_Manual
            load Manual_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            Manual_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_Manual
            save('Manual_SplitHalf_reliability.mat','Manual_SplitHalf')
            clear Manual_SplitHalf data_table_newrow 
        end

        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear Manual_MOI_table ii
disp('Manual split half reliability done')

% 2) MADE pipeline
cd xxx/Preproc_MADE
load('MADE_MeasuresOfInterest_report.mat')
Power_folder = 'xxx/Preproc_MADE/A_Power_data';
for ii = 1:height(MADE_MOI_table)
    ID = MADE_MOI_table.ID{ii};
    load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
    Freqs = 1:1:32;
    % across all trials
    if ~isempty(Power_all)
        Ntrls_tot = size(Power_all.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_all = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_all = nan(1,32);
        PowB_all = nan(1,32);
        Nsplit_all = 0;  
    end
    
    % social trials
    if ~isempty(Power_soc)
        Ntrls_tot = size(Power_soc.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_soc = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_soc = nan(1,32);
        PowB_soc = nan(1,32);
        Nsplit_soc = 0;  
    end
    
    % toy trials
    if ~isempty(Power_toy)
        Ntrls_tot = size(Power_toy.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_toy = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_toy = nan(1,32);
        PowB_toy = nan(1,32);
        Nsplit_toy = 0;  
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADE_SplitHalf = data_table_newrow;
            cd xxx/Preproc_MADE
            save('MADE_SplitHalf_reliability.mat','MADE_SplitHalf')
            clear MADE_SplitHalf data_table_newrow
        else
            cd xxx/Preproc_MADE
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADE_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_MADE
            save('MADE_SplitHalf_reliability.mat','MADE_SplitHalf')
            clear MADE_SplitHalf
        end
        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear MADE_MOI_table ii
disp('MADE split half reliability done')

% 3) MADE-BOND pipeline
cd xxx/Preproc_MADEBOND
load('MADEBOND_MeasuresOfInterest_report.mat')
Power_folder = 'xxx/Preproc_MADEBOND/A_Power_data';
for ii = 1:height(MADE_MOI_table)
    ID = MADE_MOI_table.ID{ii};
    load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
    Freqs = 1:1:32;
    % across all trials
    if ~isempty(Power_all)
        Ntrls_tot = size(Power_all.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_all = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_all = nan(1,32);
        PowB_all = nan(1,32);
        Nsplit_all = 0;  
    end
    
    % social trials
    if ~isempty(Power_soc)
        Ntrls_tot = size(Power_soc.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_soc = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_soc = nan(1,32);
        PowB_soc = nan(1,32);
        Nsplit_soc = 0;  
    end
    
    % toy trials
    if ~isempty(Power_toy)
        Ntrls_tot = size(Power_toy.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_toy = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_toy = nan(1,32);
        PowB_toy = nan(1,32);
        Nsplit_toy = 0;  
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADEBOND_SplitHalf = data_table_newrow;
            cd xxx/Preproc_MADEBOND
            save('MADEBOND_SplitHalf_reliability.mat','MADEBOND_SplitHalf')
            clear MADEBOND_SplitHalf data_table_new
        else
            cd xxx/Preproc_MADEBOND
            load MADEBOND_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADEBOND_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_MADEBOND
            save('MADEBOND_SplitHalf_reliability.mat','MADEBOND_SplitHalf')
            clear MADEBOND_SplitHalf data_table_new
        end
        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear MADE_MOI_table ii
disp('BOND-MADE split half reliability done')

% 4) HAPPE v1 pipeline
cd xxx/Preproc_HAPPEv1
load('HAPPEv1_MeasuresOfInterest_report.mat')
Power_folder = 'xxx/Preproc_HAPPEv1/A_Power_data';
for ii = 1:height(HAPPEv1_MOI_table)
    ID = HAPPEv1_MOI_table.ID{ii};
    if exist(strcat(Power_folder, '/',ID,'_Power_data.mat'),'file')
        load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
        Freqs = 1:1:32;
        % across all trials
        if ~isempty(Power_all)
            Ntrls_tot = size(Power_all.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_all = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_all = nan(1,32);
            PowB_all = nan(1,32);
            Nsplit_all = 0;  
        end

        % social trials
        if ~isempty(Power_soc)
            Ntrls_tot = size(Power_soc.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_soc = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_soc = nan(1,32);
            PowB_soc = nan(1,32);
            Nsplit_soc = 0;  
        end

        % toy trials
        if ~isempty(Power_toy)
            Ntrls_tot = size(Power_toy.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_toy = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_toy = nan(1,32);
            PowB_toy = nan(1,32);
            Nsplit_toy = 0;  
        end
        
    else
        Freqs = nan(1,32);
        Nsplit_all = 0; PowA_all = nan(1,32); PowB_all = nan(1,32);
        Nsplit_soc = 0; PowA_soc = nan(1,32); PowB_soc = nan(1,32);
        Nsplit_toy = 0; PowA_toy = nan(1,32); PowB_toy = nan(1,32);
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPEv1_SplitHalf = data_table_newrow;
            cd xxx/Preproc_HAPPEv1
            save('HAPPEv1_SplitHalf_reliability.mat','HAPPEv1_SplitHalf')
            clear HAPPEv1_SplitHalf data_table_newrow
        else
            cd xxx/Preproc_HAPPEv1
            load HAPPEv1_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPEv1_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_HAPPEv1
            save('HAPPEv1_SplitHalf_reliability.mat','HAPPEv1_SplitHalf')
            clear HAPPEv1_SplitHalf data_table_newrow
        end
        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear HAPPEv1_MOI_table ii
disp('HAPPEv1 split half reliability done')

% 5) HAPPE v4 pipeline
cd xxx/Preproc_HAPPEv4
load HAPPEv4_MeasuresOfInterest_report.mat
Power_folder = 'xxx/Preproc_HAPPEv4/A_Power_data';
for ii = 1:height(HAPPEv4_MOI_table)
    ID = HAPPEv4_MOI_table.ID{ii};
    if exist(strcat(Power_folder, '/',ID,'_Power_data.mat'),'file')
        load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
        Freqs = 1:1:32;
        % across all trials
        if ~isempty(Power_all)
            Ntrls_tot = size(Power_all.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_all = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_all = nan(1,32);
            PowB_all = nan(1,32);
            Nsplit_all = 0;  
        end

        % social trials
        if ~isempty(Power_soc)
            Ntrls_tot = size(Power_soc.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_soc = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_soc = nan(1,32);
            PowB_soc = nan(1,32);
            Nsplit_soc = 0;  
        end

        % toy trials
        if ~isempty(Power_toy)
            Ntrls_tot = size(Power_toy.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_toy = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_toy = nan(1,32);
            PowB_toy = nan(1,32);
            Nsplit_toy = 0;  
        end
        
    else
        Freqs = nan(1,32);
        Nsplit_all = 0; PowA_all = nan(1,32); PowB_all = nan(1,32);
        Nsplit_soc = 0; PowA_soc = nan(1,32); PowB_soc = nan(1,32);
        Nsplit_toy = 0; PowA_toy = nan(1,32); PowB_toy = nan(1,32);
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPEv4_SplitHalf = data_table_newrow;
            cd xxx/Preproc_HAPPEv4
            save('HAPPEv4_SplitHalf_reliability.mat','HAPPEv4_SplitHalf')
            clear HAPPEv4_SplitHalf data_table_newrow 
        else
            cd xxx/Preproc_HAPPEv4
            load HAPPEv4_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPEv4_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_HAPPEv4
            save('HAPPEv4_SplitHalf_reliability.mat','HAPPEv4_SplitHalf')
            clear HAPPEv4_SplitHalf data_table_newrow 
        end
    
        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear HAPPEv4_MOI_table ii
disp('HAPPEv4 split half reliability done')

% 6) BOND-MADE ld pipeline
cd xxx/Preproc_MADEBONDld
load MADEBOND_ld_MeasuresOfInterest_report.mat
Power_folder = 'xxx/Preproc_MADEBONDld/A_Power_data';
for ii = 1:height(MADEBOND_ld_MOI_table)
    ID = MADEBOND_ld_MOI_table.ID{ii};
    if exist(strcat(Power_folder, '/',ID,'_Power_data.mat'),'file')
        load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
        Freqs = 1:1:32;
        % across all trials
        if ~isempty(Power_all)
            Ntrls_tot = size(Power_all.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_all = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_all = nan(1,32);
            PowB_all = nan(1,32);
            Nsplit_all = 0;  
        end

        % social trials
        if ~isempty(Power_soc)
            Ntrls_tot = size(Power_soc.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_soc = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_soc = nan(1,32);
            PowB_soc = nan(1,32);
            Nsplit_soc = 0;  
        end

        % toy trials
        if ~isempty(Power_toy)
            Ntrls_tot = size(Power_toy.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_toy = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_toy = nan(1,32);
            PowB_toy = nan(1,32);
            Nsplit_toy = 0;  
        end
        
    else
        Freqs = nan(1,32);
        Nsplit_all = 0; PowA_all = nan(1,32); PowB_all = nan(1,32);
        Nsplit_soc = 0; PowA_soc = nan(1,32); PowB_soc = nan(1,32);
        Nsplit_toy = 0; PowA_toy = nan(1,32); PowB_toy = nan(1,32);
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADEBOND_ld_SplitHalf = data_table_newrow;
            cd xxx/Preproc_MADEBONDld
            save('MADEBOND_ld_SplitHalf_reliability.mat','MADEBOND_ld_SplitHalf')
            clear MADEBOND_ld_SplitHalf data_table_newrow
        else
            cd xxx/Preproc_MADEBONDld
            load MADEBOND_ld_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADEBOND_ld_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_MADEBONDld
            save('MADEBOND_ld_SplitHalf_reliability.mat','MADEBOND_ld_SplitHalf')
            clear MADEBOND_ld_SplitHalf data_table_newrow
        end
        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear MADEBOND_ld_MOI_table ii
disp('MADEBONDld split half reliability done')

% 7) HAPPILEE pipeline
cd xxx/Preproc_HAPPILEE
load HAPPILEE_MeasuresOfInterest_report.mat
Power_folder = 'xxx/Preproc_HAPPILEE/A_Power_data';
for ii = 1:height(HAPPILEE_MOI_table)
    ID = HAPPILEE_MOI_table.ID{ii};
    if exist(strcat(Power_folder, '/',ID,'_Power_data.mat'),'file')
        load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
        Freqs = 1:1:32;
        % across all trials
        if ~isempty(Power_all)
            Ntrls_tot = size(Power_all.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_all = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_all = nan(1,32);
            PowB_all = nan(1,32);
            Nsplit_all = 0;  
        end

        % social trials
        if ~isempty(Power_soc)
            Ntrls_tot = size(Power_soc.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_soc = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_soc = nan(1,32);
            PowB_soc = nan(1,32);
            Nsplit_soc = 0;  
        end

        % toy trials
        if ~isempty(Power_toy)
            Ntrls_tot = size(Power_toy.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_toy = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_toy = nan(1,32);
            PowB_toy = nan(1,32);
            Nsplit_toy = 0;  
        end
        
    else
        Freqs = nan(1,32);
        Nsplit_all = 0; PowA_all = nan(1,32); PowB_all = nan(1,32);
        Nsplit_soc = 0; PowA_soc = nan(1,32); PowB_soc = nan(1,32);
        Nsplit_toy = 0; PowA_toy = nan(1,32); PowB_toy = nan(1,32);
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPILEE_SplitHalf = data_table_newrow;
            cd xxx/Preproc_HAPPILEE
            save('HAPPILEE_SplitHalf_reliability.mat','HAPPILEE_SplitHalf')
            clear HAPPILEE_SplitHalf data_table_newrow 
        else
            cd xxx/Preproc_HAPPILEE
            load HAPPILEE_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPILEE_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_HAPPILEE
            save('HAPPILEE_SplitHalf_reliability.mat','HAPPILEE_SplitHalf')
            clear HAPPILEE_SplitHalf data_table_newrow 
        end
    
        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear HAPPILEE_MOI_table ii 
disp('HAPPILEE split half reliability done')

% 8) miniMADE pipeline
cd xxx/Preproc_miniMADE
load miniMADE_MeasuresOfInterest_report.mat
Power_folder = 'xxx/Preproc_miniMADE/A_Power_data';
for ii = 1:height(miniMADE_MOI_table)
    ID = miniMADE_MOI_table.ID{ii};
    if exist(strcat(Power_folder, '/',ID,'_Power_data.mat'),'file')
        load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
        Freqs = 1:1:32;
        % across all trials
        if ~isempty(Power_all)
            Ntrls_tot = size(Power_all.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_all = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_all = nan(1,32);
            PowB_all = nan(1,32);
            Nsplit_all = 0;  
        end

        % social trials
        if ~isempty(Power_soc)
            Ntrls_tot = size(Power_soc.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_soc = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_soc = nan(1,32);
            PowB_soc = nan(1,32);
            Nsplit_soc = 0;  
        end

        % toy trials
        if ~isempty(Power_toy)
            Ntrls_tot = size(Power_toy.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_toy = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_toy = nan(1,32);
            PowB_toy = nan(1,32);
            Nsplit_toy = 0;  
        end
        
    else
        Freqs = nan(1,32);
        Nsplit_all = 0; PowA_all = nan(1,32); PowB_all = nan(1,32);
        Nsplit_soc = 0; PowA_soc = nan(1,32); PowB_soc = nan(1,32);
        Nsplit_toy = 0; PowA_toy = nan(1,32); PowB_toy = nan(1,32);
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            miniMADE_SplitHalf = data_table_newrow;
            cd xxx/Preproc_miniMADE
            save('miniMADE_SplitHalf_reliability.mat','miniMADE_SplitHalf')
            clear miniMADE_SplitHalf data_table_newrow
        else
            cd xxx/Preproc_miniMADE
            load miniMADE_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_table_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            miniMADE_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_miniMADE
            save('miniMADE_SplitHalf_reliability.mat','miniMADE_SplitHalf')
            clear miniMADE_SplitHalf data_table_newrow  
        end

        clear Nsplit_all Nsplit_soc Nsplit_toy PowA_all PowA_soc PowA_toy PowB_all PowB_soc PowB_toy 
        clear Power_all Power_soc Power_toy ID Freqs

end
clear miniMADE_MOI_table ii
disp('miniMADE split half reliability done')

clear Power_folder


%% Correlations for split half %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check Ntrls for all and condition differences
cd xxx/Preproc_Manual
load Manual_SplitHalf_reliability.mat
cd xxx/Preproc_MADE
load MADE_SplitHalf_reliability.mat
cd xxx/Preproc_MADEBOND
load MADEBOND_SplitHalf_reliability.mat
cd xxx/Preproc_HAPPEv1
load HAPPEv1_SplitHalf_reliability.mat
cd xxx/Preproc_HAPPEv4
load HAPPEv4_SplitHalf_reliability.mat
cd xxx/Preproc_MADEBONDld
load MADEBOND_ld_SplitHalf_reliability.mat
cd xxx/Preproc_HAPPILEE
load HAPPILEE_SplitHalf_reliability.mat
cd xxx/Preproc_miniMADE
load miniMADE_SplitHalf_reliability.mat

% set those below threshold to NaNs
Ntrls = [cell2mat(Manual_SplitHalf.Nsplit_all) cell2mat(Manual_SplitHalf.Nsplit_soc) cell2mat(Manual_SplitHalf.Nsplit_toy) ...
    cell2mat(MADE_SplitHalf.Nsplit_all) cell2mat(MADE_SplitHalf.Nsplit_soc) cell2mat(MADE_SplitHalf.Nsplit_toy) ...
    cell2mat(MADEBOND_SplitHalf.Nsplit_all) cell2mat(MADEBOND_SplitHalf.Nsplit_soc) cell2mat(MADEBOND_SplitHalf.Nsplit_toy) ...
    cell2mat(HAPPEv1_SplitHalf.Nsplit_all) cell2mat(HAPPEv1_SplitHalf.Nsplit_soc) cell2mat(HAPPEv1_SplitHalf.Nsplit_toy) ...
    cell2mat(HAPPEv4_SplitHalf.Nsplit_all) cell2mat(HAPPEv4_SplitHalf.Nsplit_soc) cell2mat(HAPPEv4_SplitHalf.Nsplit_toy) ...
    cell2mat(MADEBOND_ld_SplitHalf.Nsplit_all) cell2mat(MADEBOND_ld_SplitHalf.Nsplit_soc) cell2mat(MADEBOND_ld_SplitHalf.Nsplit_toy) ...
    cell2mat(HAPPILEE_SplitHalf.Nsplit_all) cell2mat(HAPPILEE_SplitHalf.Nsplit_soc) cell2mat(HAPPILEE_SplitHalf.Nsplit_toy) ...
    cell2mat(miniMADE_SplitHalf.Nsplit_all) cell2mat(miniMADE_SplitHalf.Nsplit_soc) cell2mat(miniMADE_SplitHalf.Nsplit_toy)];
Ind_incl_SH = find(sum(Ntrls >= 20, 2) == 24);
Ind_incl_SH_fc = find(sum(Ntrls >= 90, 2) == 24);

% select and organise data
DATA_8pipelines_SH = struct();
DATA_8pipelines_SH.pipelines_names = {'Manual','MADE','MADE-BOND','HAPPEv1', ...
    'HAPPEv4','MADE-BONDld','HAPPILEE','miniMADE'};
DATA_8pipelines_SH.pipelines_SH_tables = {Manual_SplitHalf, MADE_SplitHalf, MADEBOND_SplitHalf, HAPPEv1_SplitHalf, ...
    HAPPEv4_SplitHalf, MADEBOND_ld_SplitHalf, HAPPILEE_SplitHalf, miniMADE_SplitHalf};


    % Get data - all trials
    for ii = 1:length(DATA_8pipelines_SH.pipelines_names)
        Sel_powerspectra = cell2mat(table2array(DATA_8pipelines_SH.pipelines_SH_tables{1,ii}(Ind_incl_SH,4)));
        DATA_8pipelines_SH.PowerSpectraAll_datasetA{ii} = Sel_powerspectra; 
        clear Sel_powerspectra
        Sel_powspectra = cell2mat(table2array(DATA_8pipelines_SH.pipelines_SH_tables{1,ii}(Ind_incl_SH,5)));
        DATA_8pipelines_SH.PowerSpectraAll_datasetB{ii} = Sel_powspectra; 
        clear Sel_powspectra
    end
    % Get data - soc trials
    for ii = 1:length(DATA_8pipelines_SH.pipelines_names)
        Sel_powerspectra = cell2mat(table2array(DATA_8pipelines_SH.pipelines_SH_tables{1,ii}(Ind_incl_SH,7)));
        DATA_8pipelines_SH.PowerSpectraSoc_datasetA{ii} = Sel_powerspectra; 
        clear Sel_powerspectra
        Sel_powspectra = cell2mat(table2array(DATA_8pipelines_SH.pipelines_SH_tables{1,ii}(Ind_incl_SH,8)));
        DATA_8pipelines_SH.PowerSpectraSoc_datasetB{ii} = Sel_powspectra; 
        clear Sel_powspectra
    end
    % Get data - toy trials
    for ii = 1:length(DATA_8pipelines_SH.pipelines_names)
        Sel_powerspectra = cell2mat(table2array(DATA_8pipelines_SH.pipelines_SH_tables{1,ii}(Ind_incl_SH,10)));
        DATA_8pipelines_SH.PowerSpectraToy_datasetA{ii} = Sel_powerspectra; 
        clear Sel_powerspectra
        Sel_powspectra = cell2mat(table2array(DATA_8pipelines_SH.pipelines_SH_tables{1,ii}(Ind_incl_SH,11)));
        DATA_8pipelines_SH.PowerSpectraToy_datasetB{ii} = Sel_powspectra; 
        clear Sel_powspectra
    end
    % Get IDs
    DATA_8pipelines_SH.PowerSpectra_IDxFreq = {table2array(DATA_8pipelines_SH.pipelines_SH_tables{1,1}(Ind_incl_SH,1)), 1:1:32};
    
% save data
cd xxx/DataForComparisons
save('SH_Incl_indices.mat','Ind_incl_SH','Ind_incl_SH_fc')
save('SH_DATA_8pipelines.mat','DATA_8pipelines_SH')


% Calculate correlations per frequency band for each pipeline
Freqs = 1:1:32;
Delta_ind = [find(Freqs == 2),find(Freqs == 3)];
Theta_ind = [find(Freqs == 4),find(Freqs == 6)];
Alpha_ind = [find(Freqs == 7),find(Freqs == 12)];
Beta_ind = [find(Freqs == 13),find(Freqs == 30)];

Freqs_inds = [Delta_ind; Theta_ind; Alpha_ind; Beta_ind];

% Across all trials
SH_all_r_vals = zeros(8,4);
SH_all_p_vals = zeros(8,4);
for pp = 1:length(DATA_8pipelines_SH.pipelines_names)
        for ff = 1:size(SH_all_r_vals,2)
            PowA = mean(DATA_8pipelines_SH.PowerSpectraAll_datasetA{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
            PowB = mean(DATA_8pipelines_SH.PowerSpectraAll_datasetB{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
            [rvalcurr, pvalcurr] = corr(PowA,PowB, 'type', 'Spearman');
            SH_all_r_vals(pp,ff) = rvalcurr;
            SH_all_p_vals(pp,ff) = pvalcurr;
            clear PowA PowB rvalcurr rpvalcurr
        end
end

% Condition differences: soc - toy
SH_cdiffs_r_vals = zeros(8,4);
SH_cdiffs_p_vals = zeros(8,4);
for pp = 1:length(DATA_8pipelines_SH.pipelines_names)
        for ff = 1:size(SH_cdiffs_r_vals,2)
            % difference for datasets A
            Pows_A = mean(DATA_8pipelines_SH.PowerSpectraSoc_datasetA{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
            Powt_A = mean(DATA_8pipelines_SH.PowerSpectraToy_datasetA{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
            Powd_A = Pows_A - Powt_A;
            % difference for datasets B
            Pows_B = mean(DATA_8pipelines_SH.PowerSpectraSoc_datasetB{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
            Powt_B = mean(DATA_8pipelines_SH.PowerSpectraToy_datasetB{1,pp}(:,[Freqs_inds(ff,1):Freqs_inds(ff,2)]),2);
            Powd_B = Pows_B - Powt_B;

            [rvalcurr, pvalcurr] = corr(Powd_A,Powd_B, 'type', 'Spearman');
            SH_cdiffs_r_vals(pp,ff) = rvalcurr;
            SH_cdiffs_p_vals(pp,ff) = pvalcurr;
            clear Pows_A Powt_A Pows_B Powt_B
            clear Powd_A Powd_B rvalcurr rpvalcurr
        end
end

% add correlation values to data structure
DATA_8pipelines_SH.SHcorrsAll_rvals = SH_all_r_vals;
DATA_8pipelines_SH.SHcorrsAll_pvals = SH_all_p_vals;
DATA_8pipelines_SH.SHcorrsCdiffs_rvals = SH_cdiffs_r_vals;
DATA_8pipelines_SH.SHcorrsCdiffs_pvals = SH_cdiffs_p_vals;
DATA_8pipelines_SH.SHcorrs_PsxFreq = {DATA_8pipelines_SH.pipelines_names, Freqs_inds};

cd xxx/DataForComparisons
save('SH_DATA_8pipelines.mat','DATA_8pipelines_SH')


%% save data in csv for python visualisation
cd xxx/DataForComparisons/data_csv/
% all trials
writematrix(SH_all_r_vals, 'SplitHalf_8P_all_rvals.csv')
writematrix(SH_all_p_vals, 'SplitHalf_8P_all_pvals.csv')
% condition differences
writematrix(SH_cdiffs_r_vals, 'SplitHalf_8P_diff_rvals.csv')
writematrix(SH_cdiffs_p_vals, 'SplitHalf_8P_diff_pvals.csv')
