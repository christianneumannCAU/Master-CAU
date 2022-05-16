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

% read remaining channels
PATHIN_conv = [MAIN '02_data' filesep '04_final' filesep];
cd([PATHIN_conv]);
load('T.mat');
load('00_fooof_results.mat')

% Go to Folder with LFP-data 
%(needs to be in the same folder as the script-folder)
%(see git for the structure of the folders)
PATHIN_conv = [MAIN '02_data' filesep '01_converted' filesep];
cd([PATHIN_conv]);
patient     = dir;

% sort patient structure (patient = dir read the titles of the folders as
% characters and therefore did not sort them by 1-30) 
patient([1,2],:)    = [];
for s               = 1:length(patient)
    patient(s).name = str2num(patient(s).name);
end
O       = struct2table(patient);
patient = sortrows(O);
patient = table2struct(patient);
for s   = 1:length(patient)
    patient(s).name = num2str(patient(s).name);
end
clear O s;

%% switch to patientfolder
% give p a value according to the patientfolder you want to look at

p                           = 1;   
cd([PATHIN_conv patient(p).name filesep]);                          % switch to current patientfolder
indat                       = dir('*.mat');     
DEPTH(p,1:length(indat))    = extractBetween({indat.name},'D','F'); % extract Depth from filename
SIDE(p,1:length(indat))     = extract({indat.name},1);              % extract Side from filename (does not work with Matlab R2018b, use extractBefore)
TRAJECTORY(1:length(indat)) = extractBetween({indat.name},2,3);     % extract Trajectory from filename

%% read LFP-data

for d = 1:length(indat)
    cfg         = [];
    cfg.dataset = [PATHIN_conv patient(p).name filesep indat(d).name];
    data{p,d}     = ft_preprocessing(cfg);        % read data unfiltered
end

%% file information

n_chans = 0;
for p = 1:length(patient)
    mt = 0;
    for d = 1:size(rmsd,2)
        if isempty(rmsd{p,d})
            mt = mt + 1;
        else 
            c_count{p,d} = length(rmsd{p,d});
            n_chans = n_chans + c_count{p,d}; % how many channels?
        end
    end
    f_count{p,1} = size(rmsd,2) - mt; 
end

f_count = cell2mat(f_count);

% how many files?
n_files = sum(f_count);
% avarage? 
avg_files = mean(f_count);
% standard deviation?
sd_files = std(f_count);
% least files?
min_files = min(f_count);
%most files?
max_files = max(f_count);

%% normalize data
p = 1;
for d   = 1:length(indat)
    for c   = 1:length(data{p,d}.label)               % c = channels
                mnm = min(data{p,d}.trial{1,1}(c,:));         % lowest point 
                mxm = max(data{p,d}.trial{1,1}(c,:));         % hightest point

                % normalize data
                if mxm - mnm ~= 0
                    norm_raw{p,d}(c,:) = (data{p,d}.trial{1,1}(c,:) - mnm) / (mxm - mnm);
                else
                    norm_raw{p,d}(c,:) = 1;
                end
    end
end

%% plot raw data to inspect artefacts

figure(1);
for d = 1:size(data,2)
    if length(norm_raw{p,d}) > 6
        subplot(10,16,d)
        plot(data{p,d}.time{1},norm_raw{p,d});
    else
        continue
    end
end

%% plot remaining channels

figure(2);
for d = 1:length(indat)
    subplot(10,16,d)
    if sum(DEPTH{p,d}==string(T.DEPTH(T.ID == p))) == 0
        continue
    else
        for c = 1:length(data{p,d}.label)
            if sum(data{p,d}.label(c) == string(T.CHANNEL(DEPTH{p,d}==string(T.DEPTH)))) == 1 & string(T.SIDE(DEPTH{p,d}==string(T.DEPTH))) == SIDE{p,d}
                plot(data{p,d}.time{1},norm_raw{p,d}(c,:));
                hold on
            else
                continue
            end
        end
        hold off
    end
end

%% plot histrogram for all variances
figure(3);
% variance 
vrc = cat(1,vrc{:});
histogram(vrc,0.00001:0.00098:0.15);
xlabel('Varianz');

%% plot aperiodic component and powerspectrum for channels near target and far from target
nms_ids = unique(T.ID); % find all participants

for s = 1:numel(nms_ids)
    idx_sub = T.ID == nms_ids(s);
    
    depth_id{s}     = str2double(T.DEPTH(idx_sub));
     
    nf_beta{s,1}    = depth_id{s}(dsearchn(depth_id{s},0));
    nf_beta{s,2}    = depth_id{s}(dsearchn(depth_id{s},10));
    
end

nf_beta         = cell2table(nf_beta,'VariableNames',{'near_beta_ID','far_beta_ID'});

figure(4)
for p = 1:length(patient)
    for d = 1:size(fooof_results,2)
        for c = 1:length(fooof_results{p,d})
            if str2double(DEPTH{p,d}) == depth_id{p}(dsearchn(depth_id{p},0))
                if isempty(fooof_results{p,d}(c).freqs)
                    continue
                elseif isempty(get(subplot(6,5,p), 'children'))
                    subplot(6,5,p)
                    plot(fooof_results{p,d}(c).freqs, fooof_results{p,d}(c).power_spectrum)
                    hold on
                    plot(fooof_results{p,d}(c).freqs, fooof_results{p,d}(c).ap_fit)
                    str = append('Patient',' ',num2str(p),';',' ','Tiefe',' ',DEPTH{p,d});
                    title(str);
                    xlabel 'Frequenz (Hz)';
                    ylabel 'Power (µV)';
                    hold off
                else 
                    continue
                end
            end
        end
    end
end

figure(5)
for p = 1:length(patient)
    for d = 1:size(fooof_results,2)
        for c = 1:length(fooof_results{p,d})
            if str2double(DEPTH{p,d}) == depth_id{p}(dsearchn(depth_id{p},10))
                if isempty(fooof_results{p,d}(c).freqs)
                    continue
                elseif isempty(get(subplot(6,5,p), 'children'))
                    subplot(6,5,p)
                    plot(fooof_results{p,d}(c).freqs, fooof_results{p,d}(c).power_spectrum)
                    hold on
                    plot(fooof_results{p,d}(c).freqs, fooof_results{p,d}(c).ap_fit)
                    str = append('Patient',' ',num2str(p),';',' ','Tiefe',' ',DEPTH{p,d});
                    title(str);
                    xlabel 'Frequenz (Hz)';
                    ylabel 'Power (µV)';
                    hold off
                else 
                    continue
                end
            end
        end
    end
end