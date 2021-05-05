%% Startup
clear all;  %remove all variables from current workspace
close all;  %close all plots
clc;        %clear all text from command window 

%add toolboxes and initiate fieltrip
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
cd ([PATHIN_conv])
indat = dir('*.mat');
DEPTH = extractBetween({indat.name},'D','F'); % extract Depth from filename

%% loop every file in folder for one patient

for v = 1:length(indat) 
    %% read data + preprocessing
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
        display(['ERROR IN' ' ' PATHIN_conv indat(v).name]);  
        continue
    end 

%% excluding files with no data
    if length(data{v}) & sum(data{v}.trial{1,1},'all') ~= 0
%% plot data
%         figure; hold;
%         plot(data{1,v}.time{1,1},data{1,v}.trial{1,1});
%% Extracting Spikes (Rey, Pedreira & Quiroga, 2015)
%         cfg = [];
%         cfg.bpfilter = 'yes';
%         cfg.bpfreq = [300 3000]; %bandpass-filter
%         cfg.bpfilttype = 'firws';
%         cfg.dataset = [PATHIN_conv indat(v).name];
%         spikes_raw{v} = ft_preprocessing(cfg);
% 
%         for c = 1:length(spikes_raw{v}.label) %loop for channels
%             %extract time of spikes
%             threshold = (median(abs(spikes_raw{v}.trial{1,1}(c,:))))/0.6745;
%             spike_time{v} = spike_detection(spikes_raw{v}.trial{1,1}(c,:),threshold);
% 
%             %mark every spike with a star
%             figure; hold;
%             plot(spikes_raw{v}.time{1,1},spikes_raw{v}.trial{1,1}(c,:));
%             plot(spikes_raw{v}.time{1,1}(spike_time{v}),0,'*');
%         end
       
%% FFT, Hanning taper

        cfg=[];
        cfg.method='mtmconvol'; 
        cfg.output='pow'; % Output parameter
        cfg.foi=[1:.05:30]; % Frequency resolution
        cfg.toi=[0:.01: 5]; % Temporal resolution
        cfg.t_ftimwin = 5./cfg.foi;
        cfg.taper = 'hanning'; % Frequency-Adaptive Smoothing
        TFR{v}=ft_freqanalysis(cfg,data{v});
    
%% Plotting Option 1
% 
%         m{v} = mean(TFR{v}.powspctrm,[3],'omitnan'); %avarege over time
% 
%         %normalize data
%         for u = 1:length(data{v}.label)
%             norm_data{v}(u,:) = (m{v}(u,:)-min(m{v}(u,:),[],'omitnan'))/(max(m{v}(u,:),[],'omitnan')-min(m{v}(u,:),[],'omitnan'));
%         end
% 
%         figure; hold;
%         plot(TFR{v}.freq,norm_data{v});
%         str = ['powerspectrum FFT (Hanning)' ' ' DEPTH(v)];
%         title(str) ;
%         xlabel 'Frequency [Hz]';
%         ylabel 'Power';
%         lgd = legend(extractAfter(data{v}.label,10));
%         lgd.NumColumns = length(data{v}.label);
    else
        continue
    end
end