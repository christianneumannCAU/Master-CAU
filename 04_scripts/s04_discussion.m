%% Startup

clear all;  % remove all variables from current workspace
close all;  % close all plots
clc;        % clear all text from command window 

% add subfolders and initiate fieldtrip (addpath(genpath(MAIN)) is not
% possible, because fieldtrip needs to be added seperately
MAIN = [fileparts(pwd) '\'];
addpath(genpath([MAIN '101_software\matlab functions']));
addpath(genpath([MAIN '02_data\']));
addpath(genpath([MAIN '04_scripts\']));
addpath([MAIN '101_software\fieldtrip-20210411\']);
ft_defaults;

% Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

% Go to Folder with fooof_results 
%(needs to be in the same folder as the script-folder)
%(see git for the structure of the folders)
PATHIN_conv = [MAIN '02_data' filesep '04_final' filesep];
cd([PATHIN_conv]);
load('00_fooof_results.mat')

%% create new table for correlation
fooof = [];
x = 1;
for p = 1:size(fooof_results,1)
    for d = 1:size(fooof_results,2)
        if isempty(fooof_results{p,d})                                      % ignore empty depths
            continue
        else
            for c = 1:length(fooof_results{p,d})
                if sum(fooof_results{p,d}(c).ap_fit) ~= 0                   % ignore empty channels
                    fooof{x,1} = p;                                         % id
                    fooof{x,2} = SIDE{p,d};  
                    fooof{x,3} = DEPTH{p,d}; 
                    fooof{x,4} = label{p,d}(c);                             % name of channels                 
                    fooof{x,5} = samples{p,d}(2);                           % number of samplepoints 
                    
                    % theta Power
                    fooof{x,6} = mean(fooof_results{p,d}(c).power_spectrum(fooof_results{p,d}(c).freqs>=4&fooof_results{p,d}(c).freqs<8));
                    % alpha Power
                    fooof{x,7} = mean(fooof_results{p,d}(c).power_spectrum(fooof_results{p,d}(c).freqs>=8&fooof_results{p,d}(c).freqs<13));
                    % beta Power
                    fooof{x,8} = mean(fooof_results{p,d}(c).power_spectrum(fooof_results{p,d}(c).freqs>=13&fooof_results{p,d}(c).freqs<30));
                    % low beta
                    fooof{x,9} = mean(fooof_results{p,d}(c).power_spectrum(fooof_results{p,d}(c).freqs>=13&fooof_results{p,d}(c).freqs<20));
                    % high beta
                    fooof{x,10} = mean(fooof_results{p,d}(c).power_spectrum(fooof_results{p,d}(c).freqs>=20&fooof_results{p,d}(c).freqs<30));

                    x = x+1;
                    
                end
            end
        end
    end
end

% convert to table
T = cell2table(fooof,'VariableNames',{'ID','SIDE','DEPTH','CHANNEL','SAMPLES','THETA_POWER','ALPHA_POWER','BETA_POWER','low_beta','high_beta'});
clear 'fooof';

% delete channels at depth 10 (starting point of operation, still
% calibrating)
errorcount_4 = 0;
c = 1;
while c ~= height(T)
    if str2double(T.DEPTH(c)) >= 10 
        T(c,:) = [];
        errorcount_4 = errorcount_4 + 1;
    else
        c = c+1;
    end
end

% delete channels at depth < -3 (end point)
errorcount_5 = 0;
c = 1;
while c ~= height(T)
    if str2double(T.DEPTH(c)) < -3 
        T(c,:) = [];
        errorcount_5 = errorcount_5 + 1;
    else
        c = c+1;
    end
end

%% add z-transformed data + build 2 groups for near/far from target

nms_ids = unique(T.ID); % find all participants
for s = 1:numel(nms_ids)
    idx_sub = T.ID == nms_ids(s);
    
    % z-transformations
    T.z_theta(idx_sub)      = zscore(T.THETA_POWER(idx_sub));
    T.z_alpha(idx_sub)      = zscore(T.ALPHA_POWER(idx_sub));
    T.z_beta(idx_sub)       = zscore(T.BETA_POWER(idx_sub));
    T.z_lbeta(idx_sub)      = zscore(T.low_beta(idx_sub));
    T.z_hbeta(idx_sub)      = zscore(T.high_beta(idx_sub));
    
    theta_id{s}     = T.z_theta(idx_sub);
    alpha_id{s}     = T.z_alpha(idx_sub);
    beta_id{s}      = T.z_beta(idx_sub);
    lbeta_id{s}     = T.z_lbeta(idx_sub);
    hbeta_id{s}     = T.z_hbeta(idx_sub);
    depth_id{s}     = str2double(T.DEPTH(idx_sub));
    
    
    % create new table with values for every variable near 0 and far from 0
    target{s,1}     = theta_id{s}(dsearchn(depth_id{s},0));
    target{s,2}     = theta_id{s}(dsearchn(depth_id{s},10));
    target{s,3}     = alpha_id{s}(dsearchn(depth_id{s},0));
    target{s,4}     = alpha_id{s}(dsearchn(depth_id{s},10));
    
    nf_beta{s,1}    = depth_id{s}(dsearchn(depth_id{s},0));
    target{s,7}     = beta_id{s}(dsearchn(depth_id{s},0));
    nf_beta{s,2}    = depth_id{s}(dsearchn(depth_id{s},10));
    target{s,8}     = beta_id{s}(dsearchn(depth_id{s},10));
    
    target{s,7}    = lbeta_id{s}(dsearchn(depth_id{s},0));
    target{s,8}    = lbeta_id{s}(dsearchn(depth_id{s},10));
    target{s,9}    = hbeta_id{s}(dsearchn(depth_id{s},0));
    target{s,10}   = hbeta_id{s}(dsearchn(depth_id{s},10));
    
end
target_dist_idx = cell2table(target,'VariableNames',{'near_theta','far_theta','near_alpha','far_alpha','near_beta','far_beta','near_lbeta','far_lbeta','near_hbeta','far_hbeta'}); 
nf_beta         = cell2table(nf_beta,'VariableNames',{'near_beta','far_beta'});
clear 'target';

%% safe for R 
writetable(T,'regression_table_discussion.csv');
writetable(target_dist_idx,'ttest_table_discussion.csv');
writetable(nf_beta,'or_beta_depth_nf.csv');

%% cleaning up 
clear 'c' 'd' 'dif_neg' 'l' 'label' 'MAIN' 'p' 'PATHIN_conv' 'SIDE' 'str' 'x' 's' 'exp_id' 'nms_ids' 'beta_id' 'depth_id' 'idx_sub' 'rms_id' 'theta_id' 'alpha_id' 'hbeta_id' 'lbeta_id' 'neg_m' 'dif_pow'