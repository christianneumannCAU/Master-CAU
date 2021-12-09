%% Startup

clear all;  %remove all variables from current workspace
close all;  %close all plots
clc;        %clear all text from command window 

%add subfolders and initiate fieldtrip (addpath(genpath(MAIN)) is not
%possible, because fieldtrip needs to be added seperately
MAIN = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
PATHIN_stats = fullfile(MAIN,'02_data','04_final');
cd(MAIN)

%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

%% load data
tab_reg = readtable(fullfile(PATHIN_stats,'regression_table.csv'));

% normalise data per participant

nms_ids = unique(tab_reg.ID); % find all participants

for s = 1:numel(nms_ids)
    idx_sub = tab_reg.ID == nms_ids(s);
    
    % z-transform data per participant
    tab_reg.z_exp(idx_sub)      = zscore(tab_reg.AP_EXPONENT(idx_sub));
    tab_reg.z_alpha(idx_sub)    = zscore(tab_reg.ALPHA_POWER(idx_sub));
    tab_reg.z_beta(idx_sub)     = zscore(tab_reg.BETA_POWER(idx_sub));

end

tab_reg.CHANNEL     = lower(strrep(extractAfter(tab_reg.CHANNEL,'CRAW_'),'__','')); 
tab_reg.abs_depth   = abs(tab_reg.DEPTH);

writetable(tab_reg,fullfile(PATHIN_stats,'norm_data.csv'));
