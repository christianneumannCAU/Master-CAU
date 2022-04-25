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

%% counting errors from s01_pipeline

errorcount_1 = 0;
errorcount_2 = 0; % for later; counts channel that got deleted because of a bad fit
for p = 1:size(error,1)
    for d = 1:size(error,2)
        for c = 1:5
            if isempty(error{p,d})
                continue
            else
                if isempty(error{p,d}{1,c})
                    continue
                else
                    errorcount_1 = errorcount_1 +1;
                end
            end
        end
    end
end

%% delete channels with bad fit

for p = 1:size(fooof_results,1)                                 % p = patient
    l = 1;
    for d = 1:size(fooof_results,2)                             % d = depth
        for c = 1:length(fooof_results{p,d})                    % c = channel
            if isempty(fooof_results{p,d}(c).power_spectrum)    %ignore empty cells
                continue
            else
                % calculate difference between first power value of original
                % spectrum and aperiodic component
                dif_pow{p,d}(c) = abs(fooof_results{p,d}(c).power_spectrum(1) - fooof_results{p,d}(c).ap_fit(1));

                % calculate mean of negative space when the ap_fit is
                % substracted from the powerspectrum
                dif_neg{p,d}(c,:) = fooof_results{p,d}(c).power_spectrum - fooof_results{p,d}(c).ap_fit;
                neg_m{p,d}(c) = mean(dif_neg{p,d}(dif_neg{p,d}<0));

                % delete suspicious data
                if dif_pow{p,d}(c) > 0.3 & neg_m{p,d}(c) < -0.15
                    subplot(5,10,l);
                    figure(p);
                    plot(fooof_results{p,d}(c,:).freqs,fooof_results{p,d}(c).power_spectrum);
                    %plot(fooof_results{p,d}(c,:).freqs,fooof_results{p,d}(c).fooofed_spectrum);
                    hold on;
                    plot(fooof_results{p,d}(c,:).freqs,fooof_results{p,d}(c).ap_fit);
%                     str = append(num2str(p),' ',DEPTH{p,d},' ',label{p,d}(c));
%                     title(str);
%                     xlabel 'Frequency [Hz]';
%                     ylabel 'power';
                    hold off;
                    l = l+1;
                    fooof_results{p,d}(c).aperiodic_params = [];
                    fooof_results{p,d}(c).peak_params = [];
                    fooof_results{p,d}(c).gaussian_params = [];
                    fooof_results{p,d}(c).error = [];
                    fooof_results{p,d}(c).r_squared = [];
                    fooof_results{p,d}(c).freqs = [];
                    fooof_results{p,d}(c).power_spectrum = [];
                    fooof_results{p,d}(c).fooofed_spectrum = [];
                    fooof_results{p,d}(c).ap_fit = [];
                    errorcount_2 = errorcount_2 + 1;
                    bad_fit{p,d}(c) = label{p,d}(c);
                end    
            end
        end
    end
    clear 'l';
end

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
                    fooof{x,6} = fooof_results{p,d}(c).aperiodic_params(2); % Exponent von der aperiodischen Komponente

                    % create powerspectrum without aperiodic component 
                    spectrum_wo_ap{x,1} = fooof_results{p,d}(c).power_spectrum - fooof_results{p,d}(c).ap_fit;

                    % theta Power
                    fooof{x,7} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>=4&fooof_results{p,d}(c).freqs<8));
                    % alpha Power
                    fooof{x,8} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>=8&fooof_results{p,d}(c).freqs<13));
                    % beta Power
                    fooof{x,9} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>=13&fooof_results{p,d}(c).freqs<30));
                    % root_mean_square
                    fooof{x,10} = rmsd{p,d}(c);
                    % low beta
                    fooof{x,11} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>=13&fooof_results{p,d}(c).freqs<20));
                    % high beta
                    fooof{x,12} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>=20&fooof_results{p,d}(c).freqs<30));

                    x = x+1;
                    
                end
            end
        end
    end
end

% convert to table
T = cell2table(fooof,'VariableNames',{'ID','SIDE','DEPTH','CHANNEL','SAMPLES','AP_EXPONENT','THETA_POWER','ALPHA_POWER','BETA_POWER','root_mean_square','low_beta','high_beta'});
clear 'fooof';

% delete channels with negative powers
errorcount_3 = 0;
c = 1;
while c ~= height(T)
    if (T.THETA_POWER(c) < 0) | (T.ALPHA_POWER(c) < 0) | (T.BETA_POWER(c) < 0)
        T(c,:) = [];
        errorcount_3 = errorcount_3 + 1;
    else
        c = c+1;
    end
end

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

% delete channels with too big root mean square
errorcount_6 = 0;
c = 1;
while c ~= height(T)
    if T.root_mean_square(c) > 32 
        T(c,:) = [];
        errorcount_6 = errorcount_6 + 1;
    else
        c = c+1;
    end
end

%% add z-transformed data + build 2 groups for near/far from target

nms_ids = unique(T.ID); % find all participants
for s = 1:numel(nms_ids)
    idx_sub = T.ID == nms_ids(s);
    
    % z-transformations
    T.z_exp(idx_sub)        = zscore(T.AP_EXPONENT(idx_sub));
    T.z_theta(idx_sub)      = zscore(T.THETA_POWER(idx_sub));
    T.z_alpha(idx_sub)      = zscore(T.ALPHA_POWER(idx_sub));
    T.z_beta(idx_sub)       = zscore(T.BETA_POWER(idx_sub));
    T.z_rms(idx_sub)        = zscore(T.root_mean_square(idx_sub));
    T.z_lbeta(idx_sub)      = zscore(T.low_beta(idx_sub));
    T.z_hbeta(idx_sub)      = zscore(T.high_beta(idx_sub));
    
    exp_id{s}       = T.z_exp(idx_sub);
    theta_id{s}     = T.z_theta(idx_sub);
    alpha_id{s}     = T.z_alpha(idx_sub);
    beta_id{s}      = T.z_beta(idx_sub);
    rms_id{s}       = T.z_rms(idx_sub);
    lbeta_id{s}     = T.z_lbeta(idx_sub);
    hbeta_id{s}     = T.z_hbeta(idx_sub);
    depth_id{s}     = str2double(T.DEPTH(idx_sub));
    
    
    % create new table with values for every variable near 0 and far from 0
    target{s,1}     = exp_id{s}(dsearchn(depth_id{s},0));
    target{s,2}     = exp_id{s}(dsearchn(depth_id{s},10));
    target{s,3}     = theta_id{s}(dsearchn(depth_id{s},0));
    target{s,4}     = theta_id{s}(dsearchn(depth_id{s},10));
    target{s,5}     = alpha_id{s}(dsearchn(depth_id{s},0));
    target{s,6}     = alpha_id{s}(dsearchn(depth_id{s},10));
    
    nf_beta{s,1}    = depth_id{s}(dsearchn(depth_id{s},0));
    target{s,7}     = beta_id{s}(dsearchn(depth_id{s},0));
    nf_beta{s,2}    = depth_id{s}(dsearchn(depth_id{s},10));
    target{s,8}     = beta_id{s}(dsearchn(depth_id{s},10));
    
    target{s,9}     = rms_id{s}(dsearchn(depth_id{s},0));
    target{s,10}    = rms_id{s}(dsearchn(depth_id{s},10));
    target{s,11}    = lbeta_id{s}(dsearchn(depth_id{s},0));
    target{s,12}    = lbeta_id{s}(dsearchn(depth_id{s},10));
    target{s,13}    = hbeta_id{s}(dsearchn(depth_id{s},0));
    target{s,14}    = hbeta_id{s}(dsearchn(depth_id{s},10));
    
end
target_dist_idx = cell2table(target,'VariableNames',{'near_exp','far_exp','near_theta','far_theta','near_alpha','far_alpha','near_beta','far_beta','near_rms','far_rms','near_lbeta','far_lbeta','near_hbeta','far_hbeta'});
nf_beta         = cell2table(nf_beta,'VariableNames',{'near_beta','far_beta'});
clear 'target';

%% safe for R 
writetable(T,'regression_table.csv');
writetable(target_dist_idx,'ttest_table.csv');
writetable(nf_beta,'beta_depth_nf.csv');

%% safe for visual inspection
save([MAIN '02_data' filesep '04_final' filesep 'T' '.mat'], 'T');

%% cleaning up 
clear 'c' 'd' 'dif_neg' 'l' 'label' 'MAIN' 'p' 'PATHIN_conv' 'SIDE' 'str' 'x' 's' 'nms_ids' 'idx_sub' 'rms_id' 'theta_id' 'alpha_id' 'hbeta_id' 'lbeta_id' 'neg_m' 'dif_pow'