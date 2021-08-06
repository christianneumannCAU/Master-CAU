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
cd([PATHIN_conv])
indat = dir('*.mat');
DEPTH = extractBetween({indat.name},'D','F'); % extract Depth from filename