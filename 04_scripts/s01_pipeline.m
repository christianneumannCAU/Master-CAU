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
patient = dir;

%% define variables 
% counter for errors
e               = 1;
r               = 1;
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
    error           = [];
    
    %% loop every file in folder for one patient
    for v = 1:length(indat)     
        %% read data + reject trials without data or with artefacts
        % read data
        cfg         = [];
        cfg.dataset = [PATHIN_conv patient(p).name filesep indat(v).name];
        data{v}     = ft_preprocessing(cfg);        % read data unfiltered        
        
        % downsampling
        cfg             = []
        cfg.resamplefs  = 512
        data{v}         = ft_resampledata(cfg,data{v})
        
        
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
                data{v}.trial{1}(c,:) = nan(size(data{v}.trial{1}(c,:)));
                error{e,r}              = DEPTH{v};
                error{e,r+1}            = data{v}.label(c);
                error{e,r+2}            = 'artefacts found';
                e                       = e+1;
            end
            
            % replace data with nan if there are less than 1280 samplepoints in
            % Trial
            if length(data{v}.trial{1}) < 1280
                data{v}.trial{1}(c,:) = nan(size(data{v}.trial{1}(c,:)));
                error{e,r}              = DEPTH{v};
                error{e,r+1}            = data{v}.label(c);
                error{e,r+2}            = 'not enough samplepoints';
                e                       = e+1;
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
            data_FFT{v}         = ft_preprocessing(cfg,data{v}); 
        catch ME
            error{e,r}          = DEPTH{v};
            error{e,r+1}        = data{v}.label(c);
            error{e,r+2}        = 'filtering didnt work';
            e                   = e+1;
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

        TFR{v}          = ft_freqanalysis(cfg,data_FFT{v});

        m{v} = mean(TFR{v}.powspctrm,[3],'omitnan'); %avarege powerspectrum over time

    %     % Plotting FFT
    % 
    %     %normalize data
    %     for c = 1:length(data_FFT{v}.label)
    %         mnm = min(m{v}(c,:),[],'omitnan');
    %         mxm = max(m{v}(c,:),[],'omitnan');
    %         norm_FFT{v}(c,:) = (m{v}(c,:)- mnm)/(mxm - mnm);
    %     end
    % 
    %     figure; hold;
    %     plot(TFR{v}.freq,norm_FFT{v});
    %     str = ['powerspectrum FFT (Hanning)' ' ' DEPTH(v)];
    %     title(str) ;
    %     xlabel 'Frequency [Hz]';
    %     ylabel 'Power';
    %     lgd = legend(extractAfter(data_FFT{v}.label,10));
    %     lgd.NumColumns = length(data_FFT{v}.label);

    %% FOOOF
    %make sure, that the version of python in cmd and in matlab match each
    %other (python 3.8 works optimally, python 3.6 from psychopy or 
    %python from anaconda does not work)

        settings = []; 
        settings.peak_width_limits = [1.5 12]; %minimum and maximum widths of etracted peaks
        settings.peak_threshold = 2; %standard deviation of the aperiodic-removed powerspectrum, above which a data point must pass to be considered a candidate peak
        f_range = [1 35]; %fitting range
        return_model = 1; 
        freqs{v} = TFR{v}.freq;

        for c = 1:length(data_FFT{v}.label)
            power_spectrum{v}(c,:) = m{v}(c,:);
            try
                fooof_results{p-2,v}(c,:) = fooof(freqs{v}, power_spectrum{v}(c,:) , f_range ,settings , return_model);
            catch ME
                error{e,r}              = DEPTH{v};
                error{e,r+1}            = data{v}.label(c);
                error{e,r+2}            = 'fooof didnt work';
                e                       = e+1;
                continue
            end
        end
    %% create structure for label
        label{p-2,v} = data{v}.label;    
    end
    %% Save in a patient-file
    save([MAIN '02_data' filesep '03_processed' filesep int2str(p-2) '_' patient(p).name '.mat'],'data','data_FFT','DEPTH','SIDE','TRAJECTORY','TFR','error','fooof_results');
    %% clear for next loop
    clearvars -except MAIN PATHIN_conv patient vlim_l e r fooof_results DEPTH label SIDE 
end

%% Save fooof results
save([MAIN '02_data' filesep '03_processed' filesep '00_fooof_results.mat'],'fooof_results','DEPTH','label','SIDE');