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

% Go to Folder with fooof_results 
%(needs to be in the same folder as the script-folder)
%(see git for the structure of the folders)
PATHIN_conv = [MAIN '02_data' filesep '03_processed' filesep];
cd([PATHIN_conv]);
load('00_fooof_results.mat')

%% fooof plot
% for p = 1:size(fooof_results,1) % p = patient
%     subplot(size(fooof_results,1),1,p)
%     for d = 1:size(fooof_results,2) % d = depth
%         for c = 1:length(fooof_results{p,d}) % c = channel
%             fooof_plot = plot(fooof_results{p,d}(c,:).freqs,fooof_results{p,d}(c,:).fooofed_spectrum);
%             hold on;
%         end
%     end
%     str = p;
%     title(str);
%     xlabel 'Frequency [Hz]';
%     ylabel 'fooofed spectrum';
%     hold off;
% end

%% descriptive statistic
% %variance 
% vrc_hist        = cell2mat(vrc);
% histogram(vrc_hist);

%% counting errors
errorcount = 0;
for p = 1:size(error,1)
    for d = 1:size(error,2)
        if isempty(error{p,d})
            continue
        else
            errorcount = errorcount +1;
        end
    end
end

%% create new table for regression
fooof = [];
x = 1;
for p = 1:size(fooof_results,1)
    for d = 1:size(DEPTH,2)
        for c = 1:length(fooof_results{p,d})
            if sum(fooof_results{p,d}(c).ap_fit) ~= 0
                fooof{x,1} = p; %id
                fooof{x,2} = SIDE{p,d};  
                fooof{x,3} = DEPTH{p,d}; 
                fooof{x,4} = label{p,d}(c); 
                fooof{x,5} = fooof_results{p,d}(c).aperiodic_params(2); %Exponent von der aperiodischen Komponente
                 
                %save first frequency bin from original spectrum
                first_freq_o{p,d} = or_freq{p,d}(1);
                % save first frequency bin from fooofed spectrum
                first_freq_f{p,d} = fooof_results{p,d}(c).freqs(1);
                dif_freq{p,d} = first_freq_f{p,d} - first_freq_o{p,d};
                
                % create powerspectrum without aperiodic component with
                % fooof 
                fooofed_peaks{x,1} = fooof_results{p,d}(c).fooofed_spectrum - fooof_results{p,d}(c).ap_fit;
                
                % theta Power
                fooof{x,6} = mean(fooofed_peaks{x,1}(fooof_results{p,d}(c).freqs>5&fooof_results{p,d}(c).freqs<7));
                % alpha Power
                fooof{x,7} = mean(fooofed_peaks{x,1}(fooof_results{p,d}(c).freqs>7&fooof_results{p,d}(c).freqs<12));
                % beta Power
                fooof{x,8} = mean(fooofed_peaks{x,1}(fooof_results{p,d}(c).freqs>12&fooof_results{p,d}(c).freqs<30));
                x = x+1;
            end
        end
    end
end
T = cell2table(fooof,'VariableNames',{'ID','SIDE','DEPTH','CHANNEL','AP_EXPONENT','THETA_POWER','ALPHA_POWER','BETA_POWER'}); 
%% safe for R 
writetable(T,'regression_table.csv');