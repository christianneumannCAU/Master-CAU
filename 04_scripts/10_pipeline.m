%% Startup
clear all;  %remove all variables from current workspace
close all;  %close all plots
clc;        %clear all text from command window 

%add toolboxes and initiate fieldtrip
MAIN = [fileparts(pwd) '\'];
addpath(genpath(MAIN));
addpath([userpath '\toolboxes\eeglab_current\']);
addpath(genpath([userpath '\toolboxes\fieldtrip-20210411\']));
ft_defaults;

%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

% Go to Folder with data 
%(needs to be in the same folder as the script-folder)
%(see git for the structure of the folders)
PATHIN_conv = [MAIN '02_Data' filesep '02_test' filesep];
cd([PATHIN_conv])
indat = dir('*.mat');
DEPTH = extractBetween({indat.name},'D','F'); % extract Depth from filename

%% define variables 
vlim_l = 0.001; % define lower boundary for variance 

%% loop every file in folder for one patient

for v = 1:length(indat)     
    %% read data + reject trials without data or with artefacts
    % read data
    cfg         = [];
    cfg.dataset = [PATHIN_conv indat(v).name];
    data{v}     = ft_preprocessing(cfg);        % read data unfiltered
    
    for c = 1:length(data{v}.label)                 % c = channels
        mnm = min(data{v}.trial{1,1}(c,:));         % lowest point 
        mxm = max(data{v}.trial{1,1}(c,:));         % hightest point

        % normalize data
        if mxm - mnm ~= 0
            norm_raw{v}(c,:) = (data{v}.trial{1,1}(c,:) - mnm) / (mxm - mnm);
        else
            norm_raw{v}(c,:) = 1;
        end

        % split data and calculate variance

        vrc1 = var(norm_raw{v}(c,1:ceil(end/4)));               % first quarter
        vrc2 = var(norm_raw{v}(c,ceil(end/4):ceil(end/2)));     % second quarter
        vrc3 = var(norm_raw{v}(c,ceil(end/2):ceil(end*0.75)));  % third quarter
        vrc4 = var(norm_raw{v}(c,ceil(end*0.75):end));          % fourth quarter 

        % replace data with nan if variance is smaller than 0.001
        if (vrc1 < vlim_l) || (vrc2 < vlim_l) || (vrc3 < vlim_l) || (vrc4 < vlim_l)         
            data{v}.trial{1,1}(c,:) = nan(size(data{v}.trial{1}(c,:)));
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

    try
        data_FFT{v}         = ft_preprocessing(cfg,data{v}); 
    catch ME
        continue
    end

    % FFT, Hanning taper
    cfg             = [];
    cfg.method      = 'mtmconvol'; 
    cfg.output      = 'pow';        % Output parameter
    cfg.foi         = [1:.05:30];   % Frequency resolution
    cfg.toi         = [0:.01: 5];   % Temporal resolution
    cfg.t_ftimwin   = 5./cfg.foi;
    cfg.taper       = 'hanning';    % Frequency-Adaptive Smoothing
    
    TFR{v}          = ft_freqanalysis(cfg,data_FFT{v});

    % Plotting Option 1

    m{v} = mean(TFR{v}.powspctrm,[3],'omitnan'); %avarege over time

    %normalize data
    for c = 1:length(data_FFT{v}.label)
        mnm = min(m{v}(c,:),[],'omitnan')
        mxm = max(m{v}(c,:),[],'omitnan')
        norm_FFT{v}(c,:) = (m{v}(c,:)- mnm)/(mxm - mnm);
    end

    figure; hold;
    plot(TFR{v}.freq,norm_FFT{v});
    str = ['powerspectrum FFT (Hanning)' ' ' DEPTH(v)];
    title(str) ;
    xlabel 'Frequency [Hz]';
    ylabel 'Power';
    lgd = legend(extractAfter(data_FFT{v}.label,10));
    lgd.NumColumns = length(data_FFT{v}.label);
    
%% Extracting Spikes (Rey, Pedreira & Quiroga, 2015)
    cfg             = [];
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [300 3000];           %bandpass-filter
    cfg.bpfilttype  = 'firws';
    cfg.dataset     = [PATHIN_conv indat(v).name];
    
    try
        data_spikes{v}   = ft_preprocessing(cfg,data{v});
    catch ME
        continue
    end

    for c = 1:length(data_spikes{v}.label)   % c = channel

        %extract time of spikes
        threshold = (median(abs(data_spikes{v}.trial{1,1}(c,:))))/0.6745;
        spike_time{v} = spike_detection(data_spikes{v}.trial{1,1}(c,:),threshold);

        %plot spikes and mark every spike with a star
        try
            figure; hold;
            plot(data_spikes{v}.time{1,1},data_spikes{v}.trial{1,1}(c,:));
            plot(data_spikes{v}.time{1,1}(spike_time{v}),0,'*');
            str = [data_spikes{v}.label(c) ' ' DEPTH(v)];
            title(str) ;
        catch ME
            continue
        end
    end
end