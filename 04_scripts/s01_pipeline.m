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
PATHIN_conv = [MAIN '02_data' filesep '02_test' filesep];
cd([PATHIN_conv]);
patient     = dir;

%% define variables
% define lower boundary for variance
vlim_l          = 0.001;
% empty structures
data            = [];
data_FFT        = [];
TFR             = [];
m               = [];
fooof_results   = [];

%% loop through every patient
for p = 3:length(patient)
    
    cd([PATHIN_conv patient(p).name filesep]); % switch to a patient
    indat = dir('*.mat'); 
    DEPTH(p-2,1:length(indat)) = extractBetween({indat.name},'D','F'); % extract Depth from filename
    SIDE(p-2,1:length(indat)) = extract({indat.name},1); % extract Side from filename (does not work with Matlab R2018b, use extractBefore)
    TRAJECTORY(1:length(indat)) = extractBetween({indat.name},2,3); % extract Trajectory from filename
    
    %% loop every file (Depth) in folder for one patient
    for d = 1:length(indat)
        
        error{p-2,d} = []; % prepare empty cell for error structure
        
        %% read data + reject trials without data or with artefacts
        % read data
        cfg         = [];
        cfg.dataset = [PATHIN_conv patient(p).name filesep indat(d).name];
        data{d}     = ft_preprocessing(cfg);        % read data unfiltered        
        
        % downsampling
        cfg             = [];
        cfg.resamplefs  = 512;
        data{d}         = ft_resampledata(cfg,data{d});
        
        
        for c = 1:length(data{d}.label)                 % c = channels
            mnm = min(data{d}.trial{1,1}(c,:));         % lowest point 
            mxm = max(data{d}.trial{1,1}(c,:));         % hightest point

            % normalize data
            if mxm - mnm ~= 0
                norm_raw{d}(c,:) = (data{d}.trial{1,1}(c,:) - mnm) / (mxm - mnm);
            else
                norm_raw{d}(c,:) = 1;
            end

            % split data and calculate variance

            vrc{p-2,d}(c,1) = var(norm_raw{d}(c,1:ceil(end/4)));               % first quarter
            vrc{p-2,d}(c,2) = var(norm_raw{d}(c,ceil(end/4):ceil(end/2)));     % second quarter
            vrc{p-2,d}(c,3) = var(norm_raw{d}(c,ceil(end/2):ceil(end*0.75)));  % third quarter
            vrc{p-2,d}(c,4) = var(norm_raw{d}(c,ceil(end*0.75):end));          % fourth quarter 

            % replace data with nan if variance is smaller than 0.001
            if (vrc{p-2,d}(c,1) < vlim_l) || (vrc{p-2,d}(c,2) < vlim_l) || (vrc{p-2,d}(c,3) < vlim_l) || (vrc{p-2,d}(c,4) < vlim_l)             
                data{d}.trial{1}(c,:) = nan(size(data{d}.trial{1}(c,:)));
                error{p-2,d}          = 'artefacts found';
            end
            
            % replace data with nan if there are less than 1280 samplepoints in
            % Trial
            if length(data{d}.trial{1}) < 1280
                if isempty(error) 
                    data{d}.trial{1}(c,:)   = nan(size(data{d}.trial{1}(c,:)));
                    error{p-2,d}            = 'not enough samplepoints';
                end
            end
        end

        %% FFT
        % preprocessing
        cfg             = [];
        cfg.demean      = 'yes';    %remove DC offset
        cfg.hpfilter    = 'yes';
        cfg.hpfreq      = .5;       %high-pass filter, cutting everything under .5 Hz
        cfg.hpfilttype  = 'firws';
        cfg.lpfilter    = 'yes';
        cfg.lpfreq      = 45;       %low-pass filter, cutting everything over 45 Hz
        cfg.lpfilttype  = 'firws';
        cfg.pad         = 'nextpow2'

        try
            data_FFT{d}             = ft_preprocessing(cfg,data{d}); 
        catch ME
            if isempty(error)
                error{p-2,d}        = 'filtering didnt work';
            end
            continue
        end

        % FFT, Hanning taper
        cfg             = [];
        cfg.method      = 'mtmconvol'; 
        cfg.output      = 'pow';        % Output parameter
        cfg.foi         = [2:.05:30];   % Frequency resolution
        cfg.toi         = [0:.01: 5];   % Temporal resolution
        cfg.t_ftimwin   = 5./cfg.foi;
        cfg.taper       = 'hanning';    % Frequency-Adaptive Smoothing

        TFR{d}          = ft_freqanalysis(cfg,data_FFT{d});
        
        % save frequency bins from original spectrum
        or_freq{p-2,d}  = TFR{d}.freq;
        
        %% calculate powerspectrum
        m{p-2,d}            = mean(TFR{d}.powspctrm,[3],'omitnan'); %avarege powerspectrum over time
        
        %% FOOOF
        %make sure, that the version of python in cmd and in matlab match each
        %other (python 3.8 works optimally, python 3.6 from psychopy or 
        %python from anaconda does not work)

        settings                    = []; 
        settings.peak_width_limits  = [1.5 12]; %minimum and maximum widths of etracted peaks
        settings.peak_threshold     = 2; %standard deviation of the aperiodic-removed powerspectrum, above which a data point must pass to be considered a candidate peak
        f_range                     = [2 30]; %fitting range
        return_model                = 1; 
        freqs{d}                    = TFR{d}.freq;

        for c = 1:length(data_FFT{d}.label)
            power_spectrum{d}(c,:) = m{p-2,d}(c,:);
            try
                fooof_results{p-2,d}(c,:) = fooof(freqs{d}, power_spectrum{d}(c,:) , f_range ,settings , return_model);
            catch ME
                if isempty(error)
                    error{p-2,d}         = 'fooof didnt work';
                end
                continue
            end
        end
        
        %% create structure for label
        label{p-2,d} = data{d}.label;
    end
    %% Save in a patient-file
    save([MAIN '02_data' filesep '03_processed' filesep int2str(p-2) '_' patient(p).name '.mat'],'data','data_FFT','DEPTH','SIDE','TRAJECTORY','TFR','error','fooof_results');
    %% clear for next loop
    clearvars -except MAIN PATHIN_conv patient vlim_l fooof_results DEPTH label SIDE vrc error m or_freq TFR
end

%% Save fooof results
save([MAIN '02_data' filesep '03_processed' filesep '00_fooof_results.mat'],'fooof_results','DEPTH','label','SIDE','vrc','error','m','or_freq');