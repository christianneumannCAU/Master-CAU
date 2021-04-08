%Startup
clear all;
close all;
clc;

MAIN = [fileparts(pwd) '\'];
addpath(genpath(MAIN));
addpath([userpath '\toolboxes\eeglab_current\']);
addpath([userpath '\toolboxes\fieldtrip-20190611\']);
ft_defaults;

%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

% Go to Current Folder
PATHIN_conv = [MAIN '02_Data' filesep '01_converteddata' filesep 'Andreas_Arndt' filesep];
cd ([PATHIN_conv])
indat = dir('*.mat');
DEPTH = extractBetween({indat.name},'D','F'); % Depth is a cell array containing all depth in one participant folder


%% read in the data
for v = 1:length(indat)
    data(v).name = load(indat(v).name);
    nd(v).central.data = data(v).name.CLFP_01___Central;
    nd(v).anterior.data = data(v).name.CLFP_02___Anterior;
    nd(v).posterior = data(v).name.CLFP_01___Posterior;
end
%%
v = 1
data = importdata(indat(v).name);
data.depth = extractBetween(indat(v).name,'D','F');

