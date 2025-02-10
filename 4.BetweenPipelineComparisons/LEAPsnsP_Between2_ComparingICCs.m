%% Analysis of social and non-social videos in LEAP: Comparison between pipelines

% This script tests for differences between ICCs between different pipeline.
% First, ICCs between the manual and automated pipelines for high-density
% layouts are compared. Next, the ICC between the MADE pipelines for the
% low and high-density layouts are compared with ICC between the HAPPE
% pipelines for the low and high-density layouts.  Finally, the consistency
% between the MADE-BOND pipeline and the MADE-BONDld pipline is compared
% with the consistency between the MADE-BOND pipeline and the miniMADE
% pipeline. 

% For comparisons of ICCs including common pipelines (e.g. 12 vs 13, where 
% P1 is common) we used Hotelling’s T2 tests, whereas we used Z-difference 
% scores between ICCs for comparisons without common pipelines (e.g. 28 vs 57). 
% All tests were 2-tailed. Results were considered significant if the p-value 
% exceeded .0056 (Bonferroni correction: .05/9)

% The comparisons are made for power and connectivity, calculated across
% all epochs and for condition differences, for each of the frequency
% bands. 

% The current script calls to the ICC function which can be downloaded
% from:
% https://uk.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc

% The bbcorrdiff function and folder can be downloaded from:
% https://uk.mathworks.com/matlabcentral/fileexchange/23783-bbcorrdiff-bootstrap-statistics-for-the-difference-of-correlation-coefficients

% The r_test_paired function and folder can be downloaded from:
% https://uk.mathworks.com/matlabcentral/fileexchange/25984-r_test_paired


% Note; folder paths commented out where appropriate for sharing on github
% (substituted by 'xxx')

% Created by Rianne Haartsen, PhD.; 09-2024 
% Birkbeck College, University of London

% This script is released under the GNU General Public License version 3.

%% Set paths and read in data

clear
cd xxx/DataForComparisons
addpath('xxx/ICC/')
addpath('xxx/r_test_paired/')
addpath('xxx/bbcorrdiff/')

load("DATA_8pipelines.mat")


%% Power across trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_all.ForICC(:,1:8);
    % for p1 vs other high density
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,2);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,3);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,4);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = Power_values_mat(:,3); % p3
        k = Power_values_mat(:,6); % p6
        h = Power_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_delta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted 
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_delta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_delta(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_all.Delta.All_minimalvals = MinimalsVals_all;
BP_stats.Power_all.Delta.Seperate.HotellingsT2 = HotellingsT2_delta;
BP_stats.Power_all.Delta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_delta ICCdifs_other MinimalsVals_all Zvals1 Zvals2

% Theta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_all.ForICC(:,9:16);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,2);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,3);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,4);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = Power_values_mat(:,3); % p3
        k = Power_values_mat(:,6); % p6
        h = Power_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_theta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear Power_values_mat HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_theta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_theta(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_all.Theta.All_minimalvals = MinimalsVals_all;
BP_stats.Power_all.Theta.Seperate.HotellingsT2 = HotellingsT2_theta;
BP_stats.Power_all.Theta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_theta ICCdifs_other MinimalsVals_all Zvals1 Zvals2



% Alpha power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_all.ForICC(:,17:24);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,2);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,3);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,4);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = Power_values_mat(:,3); % p3
        k = Power_values_mat(:,6); % p6
        h = Power_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_alpha = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear Power_values_mat HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_alpha(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_alpha(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_all.Alpha.All_minimalvals = MinimalsVals_all;
BP_stats.Power_all.Alpha.Seperate.HotellingsT2 = HotellingsT2_alpha;
BP_stats.Power_all.Alpha.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_alpha ICCdifs_other MinimalsVals_all Zvals1 Zvals2


% Beta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_all.ForICC(:,25:32);
% for p1 vs others
otherpipelines = [3; 4; 5];
HotellingsT2_ICC = zeros(length(otherpipelines),6);
for pp = 1:3 % Roche - MADE vs Roche - others
    j = Power_values_mat(:,1);
    k = Power_values_mat(:,2);
    h = Power_values_mat(:,otherpipelines(pp));
    [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
end
HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
for pp = 2:3 % Roche - MADE-BOND vs Roche - others
    j = Power_values_mat(:,1);
    k = Power_values_mat(:,3);
    h = Power_values_mat(:,otherpipelines(pp));
    [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
end
HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
for pp = 3 % Roche - HAPPEv1 vs Roche - others
    j = Power_values_mat(:,1);
    k = Power_values_mat(:,4);
    h = Power_values_mat(:,otherpipelines(pp));
    [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
end

% other comparisons
% 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
X = Power_values_mat(:,[2,8,5,7]);
    [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
    pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
Zvals1 = [Zdiff, Zpc, pval_2tail];
clear Zdiff Zpc pval_2tail X
% 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
X = Power_values_mat(:,[3,8,5,7]);
    [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
    pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
Zvals2 = [Zdiff, Zpc, pval_2tail];
clear Zdiff Zpc pval_2tail X
% 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
    j = Power_values_mat(:,3); % p3
    k = Power_values_mat(:,6); % p6
    h = Power_values_mat(:,8); % p8
    [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
clear p T2 df r_jk r_jh r_kh j k h

HotellingsT2_beta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
clear Power_values_mat HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted

ICCdifs_other = [Zvals1; Zvals2];

MinimalsVals_all = [HotellingsT2_beta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_beta(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_all.Beta.All_minimalvals = MinimalsVals_all;
BP_stats.Power_all.Beta.Seperate.HotellingsT2 = HotellingsT2_beta;
BP_stats.Power_all.Beta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_beta ICCdifs_other MinimalsVals_all Zvals1 Zvals2


%% Power differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_cdiffs.ForICC(:,1:8);
    % for p1 vs other high density
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,2);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,3);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,4);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = Power_values_mat(:,3); % p3
        k = Power_values_mat(:,6); % p6
        h = Power_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_delta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted 
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_delta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_delta(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_cdiffs.Delta.All_minimalvals = MinimalsVals_all;
BP_stats.Power_cdiffs.Delta.Seperate.HotellingsT2 = HotellingsT2_delta;
BP_stats.Power_cdiffs.Delta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_delta ICCdifs_other MinimalsVals_all Zvals1 Zvals2


% Theta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_cdiffs.ForICC(:,9:16);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,2);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,3);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,4);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = Power_values_mat(:,3); % p3
        k = Power_values_mat(:,6); % p6
        h = Power_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_theta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear Power_values_mat HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_theta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_theta(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_cdiffs.Theta.All_minimalvals = MinimalsVals_all;
BP_stats.Power_cdiffs.Theta.Seperate.HotellingsT2 = HotellingsT2_theta;
BP_stats.Power_cdiffs.Theta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_theta ICCdifs_other MinimalsVals_all Zvals1 Zvals2



% Alpha power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_cdiffs.ForICC(:,17:24);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,2);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,3);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,4);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = Power_values_mat(:,3); % p3
        k = Power_values_mat(:,6); % p6
        h = Power_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_alpha = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear Power_values_mat HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_alpha(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_alpha(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_cdiffs.Alpha.All_minimalvals = MinimalsVals_all;
BP_stats.Power_cdiffs.Alpha.Seperate.HotellingsT2 = HotellingsT2_alpha;
BP_stats.Power_cdiffs.Alpha.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_alpha ICCdifs_other MinimalsVals_all Zvals1 Zvals2


% Beta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Power_values_mat = DATA_8pipelines.Powerband_cdiffs.ForICC(:,25:32);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,2);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,3);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = Power_values_mat(:,1);
        k = Power_values_mat(:,4);
        h = Power_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = Power_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = Power_values_mat(:,3); % p3
        k = Power_values_mat(:,6); % p6
        h = Power_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_beta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear Power_values_mat HotellingsT2_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_beta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_beta(7,[2, 1])];

% collate into 1 structure
BP_stats.Power_cdiffs.Beta.All_minimalvals = MinimalsVals_all;
BP_stats.Power_cdiffs.Beta.Seperate.HotellingsT2 = HotellingsT2_beta;
BP_stats.Power_cdiffs.Beta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_beta ICCdifs_other MinimalsVals_all Zvals1 Zvals2




%% Functional connectivity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta FC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_all.ForICC(:,1:8);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_delta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_delta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_delta(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_all.Delta.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_all.Delta.Seperate.HotellingsT2 = HotellingsT2_delta;
BP_stats.FunCon_all.Delta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_delta ICCdifs_other MinimalsVals_all Zvals1 Zvals2



% Theta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_all.ForICC(:,9:16);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_theta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_theta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_theta(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_all.Theta.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_all.Theta.Seperate.HotellingsT2 = HotellingsT2_theta;
BP_stats.FunCon_all.Theta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_theta ICCdifs_other MinimalsVals_all Zvals1 Zvals2



% Alpha power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_all.ForICC(:,17:24);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_alpha = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_alpha(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_alpha(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_all.Alpha.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_all.Alpha.Seperate.HotellingsT2 = HotellingsT2_alpha;
BP_stats.FunCon_all.Alpha.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_alpha ICCdifs_other MinimalsVals_all Zvals1 Zvals2


% Beta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_all.ForICC(:,25:32);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end

    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_beta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_beta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_beta(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_all.Beta.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_all.Beta.Seperate.HotellingsT2 = HotellingsT2_beta;
BP_stats.FunCon_all.Beta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_beta ICCdifs_other MinimalsVals_all Zvals1 Zvals2



%% FC differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta FC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_cdiffs.ForICC(:,1:8);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_delta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_delta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_delta(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_cdiffs.Delta.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_cdiffs.Delta.Seperate.HotellingsT2 = HotellingsT2_delta;
BP_stats.FunCon_cdiffs.Delta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_delta ICCdifs_other MinimalsVals_all Zvals1 Zvals2

% Theta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_cdiffs.ForICC(:,9:16);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_theta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_theta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_theta(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_cdiffs.Theta.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_cdiffs.Theta.Seperate.HotellingsT2 = HotellingsT2_theta;
BP_stats.FunCon_cdiffs.Theta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_theta ICCdifs_other MinimalsVals_all Zvals1 Zvals2


% Alpha power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_cdiffs.ForICC(:,17:24);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_alpha = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_alpha(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_alpha(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_cdiffs.Alpha.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_cdiffs.Alpha.Seperate.HotellingsT2 = HotellingsT2_alpha;
BP_stats.FunCon_cdiffs.Alpha.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_alpha ICCdifs_other MinimalsVals_all Zvals1 Zvals2

% Beta power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FunCon_values_mat = DATA_8pipelines.FunConband_cdiffs.ForICC(:,25:32);
    % for p1 vs others
    otherpipelines = [3; 4; 5];
    HotellingsT2_ICC = zeros(length(otherpipelines),6);
    for pp = 1:3 % Roche - MADE vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,2);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC(pp,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC2 = zeros(length(otherpipelines)-1,6);
    for pp = 2:3 % Roche - MADE-BOND vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,3);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC2(pp-1,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    HotellingsT2_ICC3 = zeros(length(otherpipelines)-2,6);
    for pp = 3 % Roche - HAPPEv1 vs Roche - others
        j = FunCon_values_mat(:,1);
        k = FunCon_values_mat(:,4);
        h = FunCon_values_mat(:,otherpipelines(pp));
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
        HotellingsT2_ICC3(pp-2,:) = [p, T2, df, r_jk, r_jh, r_kh];
        clear p T2 df r_jk r_jh r_kh j k h
    end
    % other comparisons
    % 28 vs 57: MADE-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[2,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals1 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 38 vs 57: MADE-BOND-miniMADE vs HAPPEv4-HAPPILEE
    X = FunCon_values_mat(:,[3,8,5,7]);
        [~,~,Zdiff,Zpc,~,~] = bbICCdiff(X,1,5000,.95,1,1);
        pval_2tail = (1 - normcdf(abs(Zdiff)))*2;
    Zvals2 = [Zdiff, Zpc, pval_2tail];
    clear Zdiff Zpc pval_2tail X
    % 36 vs 38: MADE-BOND - MADE-BONDld vs MADE-BOND - miniMADE
        j = FunCon_values_mat(:,3); % p3
        k = FunCon_values_mat(:,6); % p6
        h = FunCon_values_mat(:,8); % p8
        [p, T2, df, r_jk, r_jh, r_kh] = rICC_test_paired(j,k,h,0);
    HotellingsT2_ICCadapted = [p, T2, df, r_jk, r_jh, r_kh];
    clear p T2 df r_jk r_jh r_kh j k h
    
    HotellingsT2_beta = [HotellingsT2_ICC; HotellingsT2_ICC2; HotellingsT2_ICC3; HotellingsT2_ICCadapted];
    clear FunCon_values_mat HotellingsT2_ICC_ICC HotellingsT2_ICC2 HotellingsT2_ICC3 HotellingsT2_ICCadapted
    ICCdifs_other = [Zvals1; Zvals2];
    MinimalsVals_all = [HotellingsT2_beta(1:6,[2, 1]); ICCdifs_other(:,[1, 4]); HotellingsT2_beta(7,[2, 1])];

% collate into 1 structure
BP_stats.FunCon_cdiffs.Beta.All_minimalvals = MinimalsVals_all;
BP_stats.FunCon_cdiffs.Beta.Seperate.HotellingsT2 = HotellingsT2_beta;
BP_stats.FunCon_cdiffs.Beta.Seperate.Zdiffs = ICCdifs_other;

clear HotellingsT2_beta ICCdifs_other MinimalsVals_all Zvals1 Zvals2

% Save data
cd xxx/DataForComparisons
save('Stats_BetweenPipelines.mat','BP_stats')



%% For reporting in text
cd xxx/DataForComparisons
load Stats_BetweenPipelines.mat

% Power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Power all
All_Ts = [BP_stats.Power_all.Delta.All_minimalvals(:,1), ...
    BP_stats.Power_all.Theta.All_minimalvals(:,1), ...
    BP_stats.Power_all.Alpha.All_minimalvals(:,1), ...
    BP_stats.Power_all.Beta.All_minimalvals(:,1)];
All_ps = [BP_stats.Power_all.Delta.All_minimalvals(:,2), ...
    BP_stats.Power_all.Theta.All_minimalvals(:,2), ...
    BP_stats.Power_all.Alpha.All_minimalvals(:,2), ...
    BP_stats.Power_all.Beta.All_minimalvals(:,2)];

All_Ts_corr = All_Ts;
All_Ts_corr(All_ps > .05/9) = NaN;

All_ps_corr = All_ps;
All_ps_corr(All_ps > .05/9) = NaN;

% Power cdiffs
Diffs_Ts = [BP_stats.Power_cdiffs.Delta.All_minimalvals(:,1), ...
    BP_stats.Power_cdiffs.Theta.All_minimalvals(:,1), ...
    BP_stats.Power_cdiffs.Alpha.All_minimalvals(:,1), ...
    BP_stats.Power_cdiffs.Beta.All_minimalvals(:,1)];
Diffs_ps = [BP_stats.Power_cdiffs.Delta.All_minimalvals(:,2), ...
    BP_stats.Power_cdiffs.Theta.All_minimalvals(:,2), ...
    BP_stats.Power_cdiffs.Alpha.All_minimalvals(:,2), ...
    BP_stats.Power_cdiffs.Beta.All_minimalvals(:,2)];

Diffs_Ts_corr = Diffs_Ts;
Diffs_Ts_corr(Diffs_ps > .05/9) = NaN;

Diffs_ps_corr = Diffs_ps;
Diffs_ps_corr(Diffs_ps > .05/9) = NaN;

Inds1 = [2, 4, 6]; % HAPPEv4
max(abs([All_Ts_corr(Inds1,:) Diffs_Ts_corr(Inds1,:)]),[],'all')
min(abs([All_Ts_corr(Inds1,:) Diffs_Ts_corr(Inds1,:)]),[],'all')
max([All_ps_corr(Inds1,:) Diffs_ps_corr(Inds1,:)], [], 'all')

Inds2 = [1, 3, 5]; Inds2f = [3, 4]; % 12-13 vs 12/13-15
max(abs([All_Ts(Inds2, Inds2f)]),[],'all')
min(abs([All_Ts(Inds2, Inds2f)]),[],'all')
min([All_ps(Inds2, Inds2f)], [], 'all')

Inds1 = [1]; % 12 vs 13 all
max(abs([All_Ts(Inds1,:)]),[],'all')
min(abs([All_Ts(Inds1,:)]),[],'all')
min([All_ps(Inds1,:) ], [], 'all')

Inds1 = [1]; Inds1f = [2, 3, 4];% 12 vs 13 condition differences
max(abs([Diffs_Ts_corr(Inds1,Inds1f)]),[],'all')
min(abs([Diffs_Ts_corr(Inds1,Inds1f)]),[],'all')
max([Diffs_ps_corr(Inds1,Inds1f)], [], 'all')

Inds2 = [7, 8]; % low vs high density MADE vs HAPPE family
max(abs([All_Ts(Inds2, :) Diffs_Ts(Inds2,:)]),[],'all')
min(abs([All_Ts(Inds2, :) Diffs_Ts(Inds2,:)]),[],'all')
min([All_ps(Inds2, :) Diffs_ps(Inds2, :)], [], 'all')

Inds2 = [9]; Inds2f = [1, 2, 3]; % 36 vs 38 all
max(abs([All_Ts(Inds2, Inds2f)]),[],'all')
min(abs([All_Ts(Inds2, Inds2f)]),[],'all')
max([All_ps(Inds2, Inds2f)], [], 'all')

Inds2 = [9]; Inds2f = [4]; % 36 vs 38 all
max(abs([Diffs_Ts(Inds2, Inds2f)]),[],'all')
min(abs([Diffs_Ts(Inds2, Inds2f)]),[],'all')
max([Diffs_ps(Inds2, Inds2f)], [], 'all')


% FunCon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FunCon all
All_Ts = [BP_stats.FunCon_all.Delta.All_minimalvals(:,1), ...
    BP_stats.FunCon_all.Theta.All_minimalvals(:,1), ...
    BP_stats.FunCon_all.Alpha.All_minimalvals(:,1), ...
    BP_stats.FunCon_all.Beta.All_minimalvals(:,1)];
All_ps = [BP_stats.FunCon_all.Delta.All_minimalvals(:,2), ...
    BP_stats.FunCon_all.Theta.All_minimalvals(:,2), ...
    BP_stats.FunCon_all.Alpha.All_minimalvals(:,2), ...
    BP_stats.FunCon_all.Beta.All_minimalvals(:,2)];

All_Ts_corr = All_Ts;
All_Ts_corr(All_ps > .05/9) = NaN;

All_ps_corr = All_ps;
All_ps_corr(All_ps > .05/9) = NaN;

% FunCon cdiffs
Diffs_Ts = [BP_stats.FunCon_cdiffs.Delta.All_minimalvals(:,1), ...
    BP_stats.FunCon_cdiffs.Theta.All_minimalvals(:,1), ...
    BP_stats.FunCon_cdiffs.Alpha.All_minimalvals(:,1), ...
    BP_stats.FunCon_cdiffs.Beta.All_minimalvals(:,1)];
Diffs_ps = [BP_stats.FunCon_cdiffs.Delta.All_minimalvals(:,2), ...
    BP_stats.FunCon_cdiffs.Theta.All_minimalvals(:,2), ...
    BP_stats.FunCon_cdiffs.Alpha.All_minimalvals(:,2), ...
    BP_stats.FunCon_cdiffs.Beta.All_minimalvals(:,2)];

Diffs_Ts_corr = Diffs_Ts;
Diffs_Ts_corr(Diffs_ps > .05/9) = NaN;

Diffs_ps_corr = Diffs_ps;
Diffs_ps_corr(Diffs_ps > .05/9) = NaN;

Inds1 = [2, 4, 6]; % HAPPEv4
max(abs([All_Ts_corr(Inds1,:) Diffs_Ts_corr(Inds1,[3, 4])]),[],'all')
min(abs([All_Ts_corr(Inds1,:) Diffs_Ts_corr(Inds1,[3, 4])]),[],'all')
max([All_ps_corr(Inds1,:) Diffs_ps_corr(Inds1,:)], [], 'all')

Inds2 = [1]; Inds2f = [1, 2]; % 12-13  delta theta all
max(abs([All_Ts(Inds2, Inds2f)]),[],'all')
min(abs([All_Ts(Inds2, Inds2f)]),[],'all')
max([All_ps(Inds2, Inds2f)], [], 'all')

Inds2 = [1]; Inds2f = [1, 3]; % 12-13  delta alpha diffs
max(abs([Diffs_Ts(Inds2, Inds2f)]),[],'all')
min(abs([Diffs_Ts(Inds2, Inds2f)]),[],'all')
max([Diffs_ps(Inds2, Inds2f)], [], 'all')

Inds1 = [6]; % 14 vs 15 newer HAPPE version
max(abs([All_Ts_corr(Inds1,:) Diffs_Ts_corr(Inds1,:)]),[],'all')
min(abs([All_Ts_corr(Inds1,:) Diffs_Ts_corr(Inds1,:)]),[],'all')
max([All_ps_corr(Inds1,:) Diffs_ps_corr(Inds1,:)], [], 'all')

Diffs_Ts(Inds1,2)
Diffs_ps(Inds1,2)


Inds2 = [7, 8]; % low vs high density MADE vs HAPPE family
max(abs([All_Ts(Inds2, :) Diffs_Ts(Inds2,:)]),[],'all')
min(abs([All_Ts(Inds2, :) Diffs_Ts(Inds2,:)]),[],'all')
min([All_ps(Inds2, :) Diffs_ps(Inds2, :)], [], 'all')

