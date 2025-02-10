%% Analysis of social and non-social videos in LEAP: Figures

% This script creates the figures published in the manuscripts:
% Figure 2) Global power and connectivity output metrics across pipelines
% Figure 4) Comparisons of ICCs

% In addition, this scripts extracts the ICC values for between pipeline
% comparisons and saves those into csv files for further generation of
% figures in Python. 

% Figure 2 uses raincloud plots to show the
% distributions of values for each pipeline. This section uses the
% Raincloud plot MATLAB code from Micah Allan and colleagues:
% Allen M, Poggiali D, Whitaker K et al. Raincloud plots: a multi-platform 
% tool for robust data visualization [version 2; peer review: 2 approved]. 
% Wellcome Open Res 2021, 4:63 (https://doi.org/10.12688/wellcomeopenres.15191.2)

% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')

% Created by Rianne Haartsen, PhD.; 06-2024 
% Birkbeck College, University of London

% This script is released under the GNU General Public License version 3.



%% Figure 2) Spectra and values within frequency bands

cd xxx/DataForComparisons
load('Incl_indices.mat')
load('DATA_8pipelines.mat')

% Power and connectivity spectrum for each pipeline %%%%%%%%%%%%%%%%%%%%%%%
spectraFig = figure;
Colours = colormap(parula(8));
% power in Figure 2, panel a
subplot(2,1,1)
freqs = DATA_8pipelines.PowerSpectra_IDxFreq{1,2};
l1 = plot([0,-10],[32,-10]);
hold on
for ff = 1:8
    Power_values = DATA_8pipelines.PowerSpectraA_data{1,ff};
    mnPow = mean(Power_values,1);
    plot(freqs, mnPow, 'Color',Colours(ff,:), 'LineWidth', 2);
    clear mnPow sdPow curve1 curve2 inBetween freqs2
end
xlim([1, 32])
xticks([1 2 3 4 6 7 12 13 30])
l2 = xline(1.5); l3 = xline(3.5); 
l4 = xline(6.5); l5 = xline(12.5); 
l6 = xline(30.5);
xlabel('Frequency (Hz)'); ylabel ('Power (log)')
% add legend
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
legend({'P1','P2', 'P3', 'P4', 'P5', 'P6', 'P7','P8'},'Location','southoutside', 'Orientation','horizontal') 


% FC in Figure 2, panel c
subplot(2,1,2)
freqs = DATA_8pipelines.FunConSpectra_IDxFreq{1,2};
l1 = plot([0,-10],[32,-10]);
hold on
for ff = 1:8
    FC_values = DATA_8pipelines.FunConSpectraA_data{1,ff};
    mnFC = mean(FC_values,1);
    plot(freqs, mnFC, 'Color',Colours(ff,:), 'LineWidth', 2, 'DisplayName',strcat('P',num2str(ff)));
    clear mnFC sdFC curve1 curve2 inBetween freqs2
end
xlim([1, 32])
xticks([1 2 3 4 6 7 12 13 30])
l2 = xline(1.5); l3 = xline(3.5); 
l4 = xline(6.5); l5 = xline(12.5); 
l6 = xline(30.5);
ylim([0 .06])
xlabel('Frequency (Hz)'); ylabel ('Connectivity (dbWPLI)')
% add legend
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
set(get(get(l6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
legend({'P1','P2', 'P3', 'P4', 'P5', 'P6', 'P7','P8'},'Location','eastoutside') 


% Individual values and distribution of power within canonical bands %%%%%%
addpath xxx/RainCloudPlots-master/tutorial_matlab

Colours = colormap(parula(8)); 
Freqs = DATA_8pipelines.PowerSpectra_IDxFreq{1,2};
Delta_ind = [find(Freqs == 2),find(Freqs == 3)];
Theta_ind = [find(Freqs == 4),find(Freqs == 6)];
Alpha_ind = [find(Freqs == 7),find(Freqs == 12)];
Beta_ind = [find(Freqs == 13),find(Freqs == 30)];
Freqs_inds = [Delta_ind; Theta_ind; Alpha_ind; Beta_ind];

RaincloudPlots = figure;
% power in Figure 2, panel b
P_del = subplot(2,4,1); % delta
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.PowerSpectraA_data{1,ff}(:,Delta_ind(1,1):Delta_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-1 .8]);
    box off 
    view([-90 90]);
    clear data_cur
    title('Delta')
    xlabel('Global power (log)')
P_the = subplot(2,4,2); % theta
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.PowerSpectraA_data{1,ff}(:,Theta_ind(1,1):Theta_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-1 .8]);
    box off 
    view([-90 90]);
    title('Theta')
    xlabel('Global power (log)')
P_alp = subplot(2,4,3); % alpha
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.PowerSpectraA_data{1,ff}(:,Alpha_ind(1,1):Alpha_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-1 .8]);
    box off 
    view([-90 90]);
    title('Alpha')
    xlabel('Global power (log)')
P_bet = subplot(2,4,4); % beta
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.PowerSpectraA_data{1,ff}(:,Beta_ind(1,1):Beta_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-1 .8]);
    box off 
    view([-90 90]);
    title('Beta')
    xlabel('Global power (log)')

% connectivity in Figure 2, panel d
F_del = subplot(2,4,5); % delta
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.FunConSpectraA_data{1,ff}(:,Delta_ind(1,1):Delta_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-180 130]);
    box off 
    view([-90 90]);
    clear data_cur
    xlabel('Global connectivity (dbWPLI)')
F_the = subplot(2,4,6); % theta
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.FunConSpectraA_data{1,ff}(:,Theta_ind(1,1):Theta_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-140 110]);
    box off 
    view([-90 90]);
F_alp = subplot(2,4,7); % alpha
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.FunConSpectraA_data{1,ff}(:,Alpha_ind(1,1):Alpha_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-90 75]);
    box off 
    view([-90 90]);
F_bet = subplot(2,4,8); % beta
    for ff = 1:8
        Pow_vals = mean(DATA_8pipelines.FunConSpectraA_data{1,ff}(:,Beta_ind(1,1):Beta_ind(1,2)),2);
        raincloud_plot(Pow_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear Pow_vals
    end
    set(gca, 'YLim', [-200 160]);
    box off 
    view([-90 90]);




    

%% Figure 4) Comparisons of ICCs

cd xxx/DataForComparisons
load Stats_BetweenPipelines.mat

% parameters
pval_thresh = .05/9; % Bonferroni correction
pval_ind = 2; % 2 = uncorrected, 7 = FDR bh adjusted
minmaxYAxis = [-22 22]; 
symbolsize = 8;
pipelinecombinations = [1:1:9]; % number of ICC comparisons made = 9 for on the x-axis


BP_stats_figure = figure;
SP11 = subplot(2,2,1); % power all
    d1 = plot(pipelinecombinations, BP_stats.Power_all.Delta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor','none','MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    hold on
    d12 = plot(pipelinecombinations((BP_stats.Power_all.Delta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.Power_all.Delta.All_minimalvals((BP_stats.Power_all.Delta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    
    d2 = plot(pipelinecombinations, BP_stats.Power_all.Theta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor','none','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    d22 = plot(pipelinecombinations((BP_stats.Power_all.Theta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.Power_all.Theta.All_minimalvals((BP_stats.Power_all.Theta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    
    d3 = plot(pipelinecombinations, BP_stats.Power_all.Alpha.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor','none','MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    d32 = plot(pipelinecombinations((BP_stats.Power_all.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.Power_all.Alpha.All_minimalvals((BP_stats.Power_all.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    
    d4 = plot(pipelinecombinations, BP_stats.Power_all.Beta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor','none','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    d42 = plot(pipelinecombinations((BP_stats.Power_all.Beta.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.Power_all.Beta.All_minimalvals((BP_stats.Power_all.Beta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    
    xticks(1:1:9)
    xticklabels({'12vs13', '12vs14','12vs15', '13vs14', '13vs15', '14vs15', ...
        '28vs57', '38vs57', '36vs38'});
    xtickangle(45)
    xline(6.5)
    ylabel({'Test statistic'; 'ICC1 < ICC2                            ICC1 > ICC2'})
    title('Pow All'); xlim([0 (size(pipelinecombinations,2)+1)]); yline(0); ylim(minmaxYAxis)
    legend([d12 d22 d32 d42],{'Delta','Theta','Alpha','Beta'}, 'Location','southoutside', 'Orientation', 'horizontal')
    clear d1 d2 d3 d4 d12 d22 d32 d42


SP12 = subplot(2,2,2); % power cdiffs
    d1 = plot(pipelinecombinations, BP_stats.Power_cdiffs.Delta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor','none','MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    hold on
    d12 = plot(pipelinecombinations((BP_stats.Power_cdiffs.Delta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.Power_cdiffs.Delta.All_minimalvals((BP_stats.Power_cdiffs.Delta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    
    d2 = plot(pipelinecombinations, BP_stats.Power_cdiffs.Theta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor','none','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    d22 = plot(pipelinecombinations((BP_stats.Power_cdiffs.Theta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.Power_cdiffs.Theta.All_minimalvals((BP_stats.Power_cdiffs.Theta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    
    d3 = plot(pipelinecombinations, BP_stats.Power_cdiffs.Alpha.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor','none','MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    d32 = plot(pipelinecombinations((BP_stats.Power_cdiffs.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.Power_cdiffs.Alpha.All_minimalvals((BP_stats.Power_cdiffs.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    
    d4 = plot(pipelinecombinations, BP_stats.Power_cdiffs.Beta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor','none','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    d42 = plot(pipelinecombinations((BP_stats.Power_cdiffs.Beta.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.Power_cdiffs.Beta.All_minimalvals((BP_stats.Power_cdiffs.Beta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    
    xticks(1:1:9)
    xticklabels({'12vs13', '12vs14','12vs15', '13vs14', '13vs15', '14vs15', ...
        '28vs57', '38vs57', '36vs38'});
    xtickangle(45)
    xline(6.5)
    ylabel({'Test statistic'; 'ICC1 < ICC2                            ICC1 > ICC2'})
    title('Pow Cond diffs'); xlim([0 (size(pipelinecombinations,2)+1)]); yline(0); ylim(minmaxYAxis)
    legend([d12 d22 d32 d42],{'Delta','Theta','Alpha','Beta'}, 'Location','southoutside', 'Orientation', 'horizontal')
    clear d1 d2 d3 d4 d12 d22 d32 d42


SP21 = subplot(2,2,3); % fc all
    d1 = plot(pipelinecombinations, BP_stats.FunCon_all.Delta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor','none','MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    hold on
    d12 = plot(pipelinecombinations((BP_stats.FunCon_all.Delta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.FunCon_all.Delta.All_minimalvals((BP_stats.FunCon_all.Delta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    
    d2 = plot(pipelinecombinations, BP_stats.FunCon_all.Theta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor','none','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    d22 = plot(pipelinecombinations((BP_stats.FunCon_all.Theta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.FunCon_all.Theta.All_minimalvals((BP_stats.FunCon_all.Theta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    
    d3 = plot(pipelinecombinations, BP_stats.FunCon_all.Alpha.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor','none','MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    d32 = plot(pipelinecombinations((BP_stats.FunCon_all.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.FunCon_all.Alpha.All_minimalvals((BP_stats.FunCon_all.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    
    d4 = plot(pipelinecombinations, BP_stats.FunCon_all.Beta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor','none','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    d42 = plot(pipelinecombinations((BP_stats.FunCon_all.Beta.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.FunCon_all.Beta.All_minimalvals((BP_stats.FunCon_all.Beta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    
    xticks(1:1:9)
    xticklabels({'12vs13', '12vs14','12vs15', '13vs14', '13vs15', '14vs15', ...
        '28vs57', '38vs57', '36vs38'});
    xtickangle(45)
    xline(6.5)
    ylabel({'Test statistic'; 'ICC1 < ICC2                            ICC1 > ICC2'})
    title('FunCon All'); xlim([0 (size(pipelinecombinations,2)+1)]); yline(0); ylim(minmaxYAxis)
    legend([d12 d22 d32 d42],{'Delta','Theta','Alpha','Beta'}, 'Location','southoutside', 'Orientation', 'horizontal')
    clear d1 d2 d3 d4 d12 d22 d32 d42


SP22 = subplot(2,2,4); % fc cdiffs
    d1 = plot(pipelinecombinations, BP_stats.FunCon_cdiffs.Delta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor','none','MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    hold on
    d12 = plot(pipelinecombinations((BP_stats.FunCon_cdiffs.Delta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.FunCon_cdiffs.Delta.All_minimalvals((BP_stats.FunCon_cdiffs.Delta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','o','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerSize',symbolsize);
    
    d2 = plot(pipelinecombinations, BP_stats.FunCon_cdiffs.Theta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor','none','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    d22 = plot(pipelinecombinations((BP_stats.FunCon_cdiffs.Theta.All_minimalvals(:,pval_ind) <= pval_thresh)), ...
        BP_stats.FunCon_cdiffs.Theta.All_minimalvals((BP_stats.FunCon_cdiffs.Theta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','square', 'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',symbolsize);
    
    d3 = plot(pipelinecombinations, BP_stats.FunCon_cdiffs.Alpha.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor','none','MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    d32 = plot(pipelinecombinations((BP_stats.FunCon_cdiffs.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.FunCon_cdiffs.Alpha.All_minimalvals((BP_stats.FunCon_cdiffs.Alpha.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','pentagram','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',symbolsize);
    
    d4 = plot(pipelinecombinations, BP_stats.FunCon_cdiffs.Beta.All_minimalvals(:,1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor','none','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    d42 = plot(pipelinecombinations((BP_stats.FunCon_cdiffs.Beta.All_minimalvals(:,pval_ind) <= pval_thresh)) , ...
        BP_stats.FunCon_cdiffs.Beta.All_minimalvals((BP_stats.FunCon_cdiffs.Beta.All_minimalvals(:,pval_ind) <= pval_thresh),1),...
        'LineStyle','none', 'Marker','diamond','MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerSize',symbolsize);
    
    xticks(1:1:9)
    xticklabels({'12vs13', '12vs14','12vs15', '13vs14', '13vs15', '14vs15', ...
        '28vs57', '38vs57', '36vs38'});
    xtickangle(45)
    xline(6.5)
    ylabel({'Test statistic'; 'ICC1 < ICC2                            ICC1 > ICC2'})
    title('FunCon Cond diffs'); xlim([0 (size(pipelinecombinations,2)+1)]); yline(0); ylim(minmaxYAxis)
    legend([d12 d22 d32 d42],{'Delta','Theta','Alpha','Beta'}, 'Location','southoutside', 'Orientation', 'horizontal')
    clear d1 d2 d3 d4 d12 d22 d32 d42




%% Preparation for figure 3

cd xxx/DataForComparisons
load('DATA_8pipelines.mat','DATA_8pipelines')

ICC_rvals_all_pow = DATA_8pipelines.Powerband_all.ICCs_r_LB_UB_pvals{1,1};
ICC_rvals_cdiff_pow = DATA_8pipelines.Powerband_cdiffs.ICCs_r_LB_UB_pvals{1,1};
% power - all
    aP_del = ICC_rvals_all_pow(1:8,1:8);
    aP_the = ICC_rvals_all_pow(9:16,9:16);
    aP_alp = ICC_rvals_all_pow(17:24,17:24);
    aP_bet = ICC_rvals_all_pow(25:32,25:32);
% power - condition differences
    cP_del = ICC_rvals_cdiff_pow(1:8,1:8);
    cP_the = ICC_rvals_cdiff_pow(9:16,9:16);
    cP_alp = ICC_rvals_cdiff_pow(17:24,17:24);
    cP_bet = ICC_rvals_cdiff_pow(25:32,25:32);

ICC_rvals_all_fc = DATA_8pipelines.FunConband_all.ICCs_r_LB_UB_pvals{1,1};
ICC_rvals_cdiff_fc = DATA_8pipelines.FunConband_cdiffs.ICCs_r_LB_UB_pvals{1,1};
% FC - all
    aF_del = ICC_rvals_all_fc(1:8,1:8);
    aF_the = ICC_rvals_all_fc(9:16,9:16);
    aF_alp = ICC_rvals_all_fc(17:24,17:24);
    aF_bet = ICC_rvals_all_fc(25:32,25:32);
% FC - condition differences
    cF_del = ICC_rvals_cdiff_fc(1:8,1:8);
    cF_the = ICC_rvals_cdiff_fc(9:16,9:16);
    cF_alp = ICC_rvals_cdiff_fc(17:24,17:24);
    cF_bet = ICC_rvals_cdiff_fc(25:32,25:32);

% save values
    cd xxx/DataForComparisons/data_csv
    % power
    writematrix(aP_del, 'Pow_ICCs_alltrls_rvals_delta.csv')
    writematrix(aP_the, 'Pow_ICCs_alltrls_rvals_theta.csv')
    writematrix(aP_alp, 'Pow_ICCs_alltrls_rvals_alpha.csv')
    writematrix(aP_bet, 'Pow_ICCs_alltrls_rvals_beta.csv')
    writematrix(cP_del, 'Pow_ICCs_cdiffs_rvals_delta.csv')
    writematrix(cP_the, 'Pow_ICCs_cdiffs_rvals_theta.csv')
    writematrix(cP_alp, 'Pow_ICCs_cdiffs_rvals_alpha.csv')
    writematrix(cP_bet, 'Pow_ICCs_cdiffs_rvals_beta.csv')
    % FC
    writematrix(aF_del, 'FC_ICCs_alltrls_rvals_delta.csv')
    writematrix(aF_the, 'FC_ICCs_alltrls_rvals_theta.csv')
    writematrix(aF_alp, 'FC_ICCs_alltrls_rvals_alpha.csv')
    writematrix(aF_bet, 'FC_ICCs_alltrls_rvals_beta.csv')
    writematrix(cF_del, 'FC_ICCs_cdiffs_rvals_delta.csv')
    writematrix(cF_the, 'FC_ICCs_cdiffs_rvals_theta.csv')
    writematrix(cF_alp, 'FC_ICCs_cdiffs_rvals_alpha.csv')
    writematrix(cF_bet, 'FC_ICCs_cdiffs_rvals_beta.csv')
    
