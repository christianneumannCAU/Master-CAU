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

for v = 1:length(indat)
    try
        cfg = [];
        cfg.demean = 'yes';             %remove DC offset
        cfg.hpfilter = 'yes';
        cfg.hpfreq = .5;                %high-pass filter, cutting everything under .5 Hz
        cfg.hpfilttype = 'firws';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 45;                %low-pass filter, cutting everything over 45 Hz
        cfg.lpfilttype = 'firws';
        cfg.dataset = [PATHIN_conv indat(v).name];
        data{v} = ft_preprocessing(cfg);
    catch ME
        display(['ERROR IN' ' ' PATHIN_conv indat(v).name])  
        continue
    end
end
%% Plotting raw data
% plot(data.time{1,1},data.trial{1,1});

%% Extracting Spikes (Rey, Pedreira & Quiroga, 2015)
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [300 3000];
cfg.bpfilttype = 'firws';
cfg.dataset = [PATHIN_conv indat(1).name];
spikes_raw = ft_preprocessing(cfg);

for v = 1:2
    threshold = (median(abs(spikes_raw.trial{1,1}(v,:))))/0.6745;
    spike_time = spike_detection(spikes_raw.trial{1,1}(v,:),threshold);
    
    figure; hold;
    plot(spikes_raw.time{1,1},spikes_raw.trial{1,1}(v,:));
    plot(spikes_raw.time{1,1}(spike_time),0,'*');
end
    
%% Option 1: FFT, Hanning taper
% Calculating multiple FFTs for a different number of Cycles (3 to 7)

x = 1;
for v = 3:7 
    cfg=[];
    cfg.method='mtmconvol'; 
    cfg.output='pow'; % Output parameter
    cfg.foi=[1:.05:30]; % Frequency resolution
    cfg.toi=[0:.01: 5]; % Temporal resolution
    cfg.t_ftimwin = v./cfg.foi;
    cfg.taper = 'hanning'; % Frequency-Adaptive Smoothing
    TFR1(x)=ft_freqanalysis(cfg,data);
    x = x+1;
end

%% Plotting Option 1
n = 3; %number of Cycles
for v = 1:5
    subplot(5,1,v);
    m = mean(TFR1(v).powspctrm,[3],'omitnan'); %avarege over time
    
    %normalize data
    for u = 1:2
        norm_data(u,:) = (m(u,:)-min(m(u,:),[],'omitnan'))/(max(m(u,:),[],'omitnan')-min(m(u,:),[],'omitnan'));
    end
    
    plot(TFR1(v).freq,norm_data);
    str = ['powerspectrum FFT (Hanning)' ' ' num2str(n) ' ' 'Cycles'];
    title(str) ;
    xlabel 'Frequency [Hz]';
    ylabel 'Power';
    lgd = legend('Central','Anterior');
    lgd.NumColumns = 2;
    n = n+1;
end

%% Option 2: Wavelet
% Calculating multiple FFTs for a different number of Cycles (3 to 7)

x = 1;
for v = 3:7
    cfg=[];
    cfg.method='wavelet'; % Method: Wavelet Transformation
    cfg.output='pow'; % Output parameter
    cfg.foilim=[1 30]; % Frequency resolution
    cfg.toi=[0:.01: 5]; % Temporal resolution
    cfg.width = v;
    TFR2(x) = ft_freqanalysis(cfg,data);
    x = x+1;
end

%% Plotting Option 2
n = 3; %number of Cycles
for v = 1:5
    subplot(5,1,v);
    m = mean(TFR2(v).powspctrm,[3],'omitnan'); %avarege over time
    
    %normalize data
    for u = 1:2
        norm_data(u,:) = (m(u,:)-min(m(u,:),[],'omitnan'))/(max(m(u,:),[],'omitnan')-min(m(u,:),[],'omitnan'));
    end
    
    plot(TFR2(v).freq,norm_data);
    str = ['powerspectrum FFT (Wavelet)' ' ' num2str(n) ' ' 'Cycles'];
    title(str) ;
    xlabel 'Frequency [Hz]';
    ylabel 'Power';
    lgd = legend('Central','Anterior');
    lgd.NumColumns = 2;
    n = n+1;
end

%% TF plots
subplot(2,1,1)

% FFT
imagesc(TFR1(1).time,TFR1(1).freq,squeeze(TFR1(1).powspctrm(1,:,:)));
    % legend
    title 'FFT (Hanning) TFR'
    xlabel 'Time [s]';
    ylabel 'Frequency [Hz]';
    cb = colorbar;
    cb.Label.String = 'Power [?] :)';
    
% Wavelets
subplot(2,1,2)
imagesc(TFR2(1).time,TFR2(1).freq,squeeze(TFR2(1).powspctrm(1,:,:)));
    % legend
    title 'Wavelet TFR'
    xlabel 'Time [s]';
    ylabel 'Frequency [Hz]';
    cb = colorbar;
    cb.Label.String = 'Power [?] :)';