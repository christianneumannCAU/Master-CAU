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
PATHIN_conv = [MAIN '02_data' filesep '02_test' filesep];
cd([PATHIN_conv]);
patient = dir;
%% define variables 
vlim_l = 0.001; % define lower boundary for variance 

%% loop through every patient

for p = 3:length(patient)
    cd([PATHIN_conv patient(p).name filesep]); % switch to a patient
    indat = dir('*.mat'); 
    DEPTH(1:length(indat),p-2) = extractBetween({indat.name},'D','F'); % extract Depth from filename
    SIDE(1:length(indat),p-2) = extract({indat.name},1); % extract Side from filename 
    TRAJECTORY(1:length(indat),p-2) = extractBetween({indat.name},2,3); % extract Trajectory from filename

    %% loop every file in folder for one patient
    for v = 1:length(indat)     
        %% read data + reject trials without data or with artefacts
        % read data
        cfg         = [];
        cfg.dataset = [PATHIN_conv patient(p).name filesep indat(v).name];
        data{v,p-2}     = ft_preprocessing(cfg);        % read data unfiltered

        for c = 1:length(data{v,p-2}.label)                 % c = channels
            mnm = min(data{v,p-2}.trial{1,1}(c,:));         % lowest point 
            mxm = max(data{v,p-2}.trial{1,1}(c,:));         % hightest point

            % normalize data
            if mxm - mnm ~= 0
                norm_raw{v}(c,:) = (data{v,p-2}.trial{1,1}(c,:) - mnm) / (mxm - mnm);
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
                data{v,p-2}.trial{1,1}(c,:) = nan(size(data{v,p-2}.trial{1}(c,:)));
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
            data_FFT{v,p-2}         = ft_preprocessing(cfg,data{v,p-2}); 
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

        TFR{v,p-2}          = ft_freqanalysis(cfg,data_FFT{v,p-2});

        m{v,p-2} = mean(TFR{v,p-2}.powspctrm,[3],'omitnan'); %avarege powerspectrum over time

    %     % Plotting FFT
    % 
    %     %normalize data
    %     for c = 1:length(data_FFT{v,p-2}.label)
    %         mnm = min(m{v,p-2}(c,:),[],'omitnan');
    %         mxm = max(m{v,p-2}(c,:),[],'omitnan');
    %         norm_FFT{v}(c,:) = (m{v,p-2}(c,:)- mnm)/(mxm - mnm);
    %     end
    % 
    %     figure; hold;
    %     plot(TFR{v,p-2}.freq,norm_FFT{v});
    %     str = ['powerspectrum FFT (Hanning)' ' ' DEPTH(v)];
    %     title(str) ;
    %     xlabel 'Frequency [Hz]';
    %     ylabel 'Power';
    %     lgd = legend(extractAfter(data_FFT{v,p-2}.label,10));
    %     lgd.NumColumns = length(data_FFT{v,p-2}.label);

    %% FOOOF
    %make sure, that the version of python in cmd and in python match each
    %other

        settings = []; 
        settings.peak_width_limits = [1.5 12]; %minimum and maximum widths of etracted peaks
        settings.peak_threshold = 2; %standard deviation of the aperiodic-removed powerspectrum, above which a data point must pass to be considered a candidate peak
        f_range = [1 35]; %fitting range
        return_model = 1; 
        freqs{v,p-2} = TFR{v,p-2}.freq;

        for c = 1:length(data_FFT{v,p-2}.label)
            power_spectrum{v,p-2}(c,:) = m{v,p-2}(c,:);
            try
                fooof_results{v,p-2}(c,:) = fooof(freqs{v,p-2}, power_spectrum{v,p-2}(c,:) , f_range ,settings , return_model);
            catch ME
                continue
            end
        end

    %% Extracting Spikes (Rey, Pedreira & Quiroga, 2015)
    %     cfg             = [];
    %     cfg.bpfilter    = 'yes';
    %     cfg.bpfreq      = [300 3000];           %bandpass-filter
    %     cfg.bpfilttype  = 'firws';
    %     cfg.dataset     = [PATHIN_conv indat(v).name];
    %     
    %     try
    %         data_spikes{v}   = ft_preprocessing(cfg,data{v});
    %     catch ME
    %         continue
    %     end
    % 
    %     for c = 1:length(data_spikes{v}.label)   % c = channel
    % 
    %         %extract time of spikes
    %         threshold = (median(abs(data_spikes{v}.trial{1,1}(c,:))))/0.6745;
    %         spike_time{v} = spike_detection(data_spikes{v}.trial{1,1}(c,:),threshold);
    % 
    %         %plot spikes and mark every spike with a star
    %         try
    %             figure; hold;
    %             plot(data_spikes{v}.time{1,1},data_spikes{v}.trial{1,1}(c,:));
    %             plot(data_spikes{v}.time{1,1}(spike_time{v}),0,'*');
    %             str = [data_spikes{v}.label(c) ' ' DEPTH(v)];
    %             title(str) ;
    %         catch ME
    %             continue
    %         end
    %     end
    end
    %% clear for next loop
    clear norm_raw;
    clear power_spectrum
end

%% Save for later
% convert Matlab variables, name and save them
DEPTH = array2table(DEPTH,'VariableNames',{patient(3:end).name});
SIDE = array2table(SIDE,'VariableNames',{patient(3:end).name});
TRAJECTORY = array2table(TRAJECTORY,'VariableNames',{patient(3:end).name});
data = array2table(data,'VariableNames',{patient(3:end).name});
data_FFT = array2table(data_FFT,'VariableNames',{patient(3:end).name});
fooof_results = array2table(fooof_results,'VariableNames',{patient(3:end).name});
TFR = array2table(TFR,'VariableNames',{patient(3:end).name});
save([MAIN '02_data' filesep '03_processed' filesep,'allData_preproc.mat'],'data','data_FFT','fooof_results','DEPTH','SIDE','TRAJECTORY','TFR');