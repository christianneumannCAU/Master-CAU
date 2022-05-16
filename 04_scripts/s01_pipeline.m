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

% Go to Folder with data 
%(needs to be in the same folder as the script-folder)
%(see git for the structure of the folders)
PATHIN_conv = [MAIN '02_data' filesep '01_converted' filesep];
cd([PATHIN_conv]);
patient     = dir;

% sort patient structure (patient = dir read the titles of the folders as
% characters and therefore did not sort them by 1-30 
patient([1,2],:)    = [];
for s               = 1:length(patient)
    patient(s).name = str2num(patient(s).name);
end
T       = struct2table(patient);
patient = sortrows(T);
patient = table2struct(patient);
for s   = 1:length(patient)
    patient(s).name = num2str(patient(s).name);
end
clear T s;

%% define variables

% define boundaries for variance
vlim_l          = 0.003;
vlim_m          = 0.075;

% empty structures
data            = [];
data_FFT        = [];
TFR             = [];
m               = [];
fooof_results   = [];
    
%% loop through every patient

for p = 1:length(patient)
    
    cd([PATHIN_conv patient(p).name filesep]);                          % switch to current patientfolder
    indat                       = dir('*.mat');     
    DEPTH(p,1:length(indat))    = extractBetween({indat.name},'D','F'); % extract Depth from filename
    SIDE(p,1:length(indat))     = extract({indat.name},1);              % extract Side from filename (does not work with Matlab R2018b, use extractBefore)
    TRAJECTORY(1:length(indat)) = extractBetween({indat.name},2,3);     % extract Trajectory from filename
    
    %% loop every file in the current patientfolder
    
    for d = 1:length(indat)
        
        % preparing search for errors
        error{p,d}  = cell(1,5);

        %% read LFP-data
        
        cfg         = [];
        cfg.dataset = [PATHIN_conv patient(p).name filesep indat(d).name];
        data{d}     = ft_preprocessing(cfg);        % read data unfiltered
        
        %% read raw data and calculate rms for spike-activity
        
        load(indat(d).name);
        raw{d}      = [];
        for c       = 1:length(data{d}.label)
            chans{d}(c)     = append('CSPK',extractAfter(data{d}.label(c),4));
            name            = string(chans{d}(c));
            raw{d}(c,:)     = eval(name);
            rmsd{p,d}(c)    = rms(raw{d}(c,:));
        end

        %% downsampling to 512 Hz
        
        cfg             = [];
        cfg.resamplefs  = 512;
        data{d}         = ft_resampledata(cfg,data{d});
        
        %% reject trials without data, with artefacts or with enough samplepoints
        
        for c   = 1:length(data{d}.label)               % c = channels
            mnm = min(data{d}.trial{1,1}(c,:));         % lowest point 
            mxm = max(data{d}.trial{1,1}(c,:));         % hightest point

            % normalize data
            if mxm - mnm ~= 0
                norm_raw{d}(c,:) = (data{d}.trial{1,1}(c,:) - mnm) / (mxm - mnm);
            else
                norm_raw{d}(c,:) = 1;
            end

            % split data and calculate variance
            vrc{p,d}(c,1) = var(norm_raw{d}(c,1:ceil(end/4)));               % first quarter
            vrc{p,d}(c,2) = var(norm_raw{d}(c,ceil(end/4):ceil(end/2)));     % second quarter
            vrc{p,d}(c,3) = var(norm_raw{d}(c,ceil(end/2):ceil(end*0.75)));  % third quarter
            vrc{p,d}(c,4) = var(norm_raw{d}(c,ceil(end*0.75):end));          % fourth quarter 
            
            % replace data with nan if variance is smaller or equal 0.003
            if (vrc{p,d}(c,1) <= vlim_l) || (vrc{p,d}(c,2) <= vlim_l) || (vrc{p,d}(c,3) <= vlim_l) || (vrc{p,d}(c,4) <= vlim_l)             
                data{d}.trial{1}(c,:)   = nan(size(data{d}.trial{1}(c,:)));
                error{p,d}(c)           = cellstr('artefacts found');
            end
            
            % replace data with nan if variance is greater than 0.055
            if (vrc{p,d}(c,1) > vlim_m) || (vrc{p,d}(c,2) > vlim_m) || (vrc{p,d}(c,3) > vlim_m) || (vrc{p,d}(c,4) > vlim_m)             
                if isempty(error{p,d}{1,c})
                    data{d}.trial{1}(c,:)   = nan(size(data{d}.trial{1}(c,:)));
                    error{p,d}(c)           = cellstr('artefacts found');
                end
            end
            
            % replace data with nan if there are less than 1280 samplepoints in
            % Trial
            if length(data{d}.trial{1}) < 1280
                if isempty(error{p,d}{1,c})
                    data{d}.trial{1}(c,:)   = nan(size(data{d}.trial{1}(c,:)));   
                    error{p,d}(c)           = cellstr('not enough samplepoints');
                end
            end
        end
        %% FFT
        % preprocessing
        cfg             = [];
        cfg.demean      = 'yes';    % remove DC offset
        cfg.hpfilter    = 'yes';
        cfg.hpfreq      = .5;       % high-pass filter, cutting everything under .5 Hz
        cfg.hpfilttype  = 'firws';
        cfg.lpfilter    = 'yes';
        cfg.lpfreq      = 45;       % low-pass filter, cutting everything over 45 Hz
        cfg.lpfilttype  = 'firws';
        cfg.pad         = 'nextpow2';

        try
            data_FFT{d} = ft_preprocessing(cfg,data{d}); 
        catch ME
            continue
        end

        % FFT, Hanning taper
        cfg             = [];
        cfg.method      = 'mtmconvol'; 
        cfg.output      = 'pow';        % Output parameter
        cfg.foi         = [2:.05:35];   % Frequency resolution
        cfg.toi         = [0:.01: 5];   % Temporal resolution
        cfg.t_ftimwin   = 5./cfg.foi;
        cfg.taper       = 'hanning';    % Frequency-Adaptive Smoothing

        TFR{d}          = ft_freqanalysis(cfg,data_FFT{d});
        
        % save frequency bins from original spectrum
        or_freq{p,d}    = TFR{d}.freq;
        
        %% calculate powerspectrum
        m{p,d}          = mean(TFR{d}.powspctrm,[3],'omitnan'); % avarege powerspectrum over time
        
        %% FOOOF
        % make sure, that the version of python in cmd and in matlab match each
        % other (python 3.8 works optimally, python 3.6 from psychopy or 
        % python from anaconda does not work)

        settings                    = []; 
        settings.peak_width_limits  = [0.5 12]; % minimum and maximum widths of etracted peaks
        settings.peak_threshold     = 2;        % standard deviation of the aperiodic-removed powerspectrum, above which a data point must pass to be considered a candidate peak
        f_range                     = [4 30];   % fitting range
        return_model                = 1; 
        freqs{d}                    = TFR{d}.freq;

        for c = 1:length(data_FFT{d}.label)
            power_spectrum{d}(c,:) = m{p,d}(c,:);
            try
                fooof_results{p,d}(c,:) = fooof(freqs{d}, power_spectrum{d}(c,:) , f_range ,settings , return_model);
            catch ME
                continue
            end
        end
        
        %% create structures for label and samplesize
        label{p,d}      = data{d}.label;
        samples{p,d}    = data_FFT{d}.sampleinfo;
    end
    %% Save in a patient-file
    save([MAIN '02_data' filesep '03_processed' filesep patient(p).name '.mat'],'data','data_FFT','DEPTH','SIDE','TRAJECTORY','TFR','error','fooof_results');
    %% clear for next loop
    clearvars -except MAIN PATHIN_conv patient vlim_l vlim_m fooof_results DEPTH label SIDE vrc error m or_freq e samples rmsd 
end

%% Save fooof results
save([MAIN '02_data' filesep '04_final' filesep '00_fooof_results.mat'],'fooof_results','DEPTH','label','SIDE','vrc','error','or_freq','samples','rmsd');