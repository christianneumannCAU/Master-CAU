%% Startup
clear all;  %remove all variables from current workspace
close all;  %close all plots
clc;        %clear all text from command window 

%add subfolders and initiate fieldtrip (addpath(genpath(MAIN)) is not
%possible, because fieldtrip needs to be added seperately
MAIN = [fileparts(pwd) '\'];
addpath(genpath([MAIN '101_software\matlab functions']));
addpath(genpath([MAIN '02_data\']));
addpath(genpath([MAIN '04_scripts\']));
addpath([MAIN '101_software\fieldtrip-20210411\']);
ft_defaults;

%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

% Go to Folder with data 
%(needs to be in the same folder as the script-folder)
%(see git for the structure of the folders)
PATHIN_conv = [MAIN '02_data' filesep '01_converted' filesep];
cd([PATHIN_conv]);
patient     = dir;

% clean and sort patient structure
patient([1,2],:) = [];
for s = 1:length(patient)
    patient(s).name = str2num(patient(s).name);
end
T       = struct2table(patient);
patient = sortrows(T);
patient = table2struct(patient);
for s = 1:length(patient)
    patient(s).name = num2str(patient(s).name);
end
clear T s;

%% define variables
% define lower boundary for variance
vlim_l          = 0.001;

% empty structures
data            = [];
data_FFT        = [];
TFR             = [];
m               = [];
fooof_results   = [];
rootmeansquare  = [];

%% loop through every patient
for p = 1:1
    
    cd([PATHIN_conv patient(p).name filesep]); % switch to a patient
    indat = dir('*.mat');     
    DEPTH(p,1:length(indat)) = extractBetween({indat.name},'D','F'); % extract Depth from filename
    SIDE(p,1:length(indat)) = extract({indat.name},1); % extract Side from filename (does not work with Matlab R2018b, use extractBefore)
    TRAJECTORY(1:length(indat)) = extractBetween({indat.name},2,3); % extract Trajectory from filename
    
    %% loop every file (Depth) in folder for one patient
    for d = 1:length(DEPTH)
        
        % preparing search for errors
        error{p,d}      = cell(1,5);
        
        %% read data
        cfg         = [];
        cfg.dataset = [PATHIN_conv patient(p).name filesep indat(d).name];
        data{d}     = ft_preprocessing(cfg);        % read data unfiltered
        
        %% raw data and rms for spike-activity
        cfg             = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [300 3000];           %bandpass-filter
        cfg.bpfilttype  = 'firws';
        
        try
            data_spikes{d}   = ft_preprocessing(cfg,data{d});
        catch ME
            continue
        end
    end
end