%% Startup
clear all;  %remove all variables from current workspace
close all;  %close all plots
clc;        %clear all text from command window 

%add and initiate toolboxes
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
PATHIN_conv = [MAIN '02_Data' filesep '01_converteddata' filesep 'Andreas_Arndt' filesep];
cd ([PATHIN_conv])
indat = dir('*.mat');
DEPTH = extractBetween({indat.name},'D','F'); % extract Depth from filename

%% read in the data + preprocessing
cfg = [];
cfg.demean = 'yes';             %remove DC offset
cfg.hpfilter = 'yes';
cfg.hpfreq = .5;                %high-pass filter, cutting everything under .5 Hz
cfg.hpfilttype = 'firws';
cfg.lpfilter = 'yes';
cfg.lpfreq = 45;                %low-pass filter, cutting everything over 45 Hz
cfg.lpfilttype = 'firws';
cfg.dataset = [PATHIN_conv indat(1).name];
data = ft_preprocessing(cfg);

%% Option 1: FFT, Hanning taper
    cfg=[];
    cfg.method='mtmconvol'; 
    cfg.output='pow'; % Output parameter
    cfg.foi=[1:2:30]; % Frequency resolution
    cfg.toi=[0:.01: 5]; % Temporal resolution
    cfg.t_ftimwin = 3./cfg.foi;
    cfg.taper = 'hanning'; % Frequency-Adaptive Smoothing

    TFR1=ft_freqanalysis(cfg,data);

%% Option 2: Wavelet
    cfg=[];
    cfg.method='wavelet'; % Method: Wavelet Transformation
    cfg.output='pow'; % Output parameter
    cfg.foi=[1:2:30]; % Frequency resolution
    cfg.toi=[0:.01: 5]; % Temporal resolution
    cfg.width = 3;
    TFR2 = ft_freqanalysis(cfg,data);
    
%% plots
%raw data
plot(data.time{1,1},data.trial{1,1});

%% for Option 1
% avarage over all timesamples;
m1 = mean(TFR1.powspctrm,[3],'omitnan'); 
plot(TFR1.freq,m1);

%% for Option 2
% avarage over all timesamples;
m2 = mean(TFR2.powspctrm,[3],'omitnan'); 
plot(TFR2.freq,m2);

%% TF plots
subplot(2,1,1)

% FFT
imagesc(TFR1.time,TFR1.freq,squeeze(TFR1.powspctrm(1,:,:)));
    % legend
    title 'FFT (Hanning) TFR'
    xlabel 'Time [s]';
    ylabel 'Frequency [Hz]';
    cb = colorbar;
    cb.Label.String = 'Power [?] :)';
    
% Wavelets
subplot(2,1,2)
imagesc(TFR2.time,TFR2.freq,squeeze(TFR2.powspctrm(1,:,:)));
    % legend
    title 'Wavelet TFR'
    xlabel 'Time [s]';
    ylabel 'Frequency [Hz]';
    cb = colorbar;
    cb.Label.String = 'Power [?] :)';

