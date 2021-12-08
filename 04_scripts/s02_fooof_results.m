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

%% counting errors from s01_pipeline

errorcount_1 = 0;
errorcount_2 = 0; % for later; counts channel that got deleted because of a bad fit
for p = 1:size(error,1)
    for d = 1:size(error,2)
        for c = 1:3
            if isempty(error{p,d})
                continue
            else
                if isempty(error{p,d}{1,c})
                    continue
                else
                    errorcount_1 = errorcount_1 +1;
                end
            end
        end
    end
end

%% delete channels with bad fit

l = 1;
for p = 1:size(fooof_results,1) % p = patient
    for d = 1:size(fooof_results,2) % d = depth
        for c = 1:length(fooof_results{p,d}) % c = channel
            if isempty(fooof_results{p,d}(c).power_spectrum)
                continue
            else
                % calculate difference between first power value of original
                % spectrum and aperiodic component
                dif_pow{p,d}(c) = abs(fooof_results{p,d}(c).power_spectrum(1) - fooof_results{p,d}(c).ap_fit(1));

                % calculate mean of negative space when the ap_fit is
                % substracted from the powerspectrum
                dif_neg{p,d}(c,:) = fooof_results{p,d}(c).power_spectrum - fooof_results{p,d}(c).ap_fit;
                neg_m{p,d}(c) = mean(dif_neg{p,d}(dif_neg{p,d}<0));

                % delete suspicious data
                if dif_pow{p,d}(c) > 0.3|neg_m{p,d}(c) < -0.1
%                     subplot(size(DEPTH,2),3,l)
%                     plot(fooof_results{p,d}(c,:).freqs,fooof_results{p,d}(c).power_spectrum);
%                     %plot(fooof_results{p,d}(c,:).freqs,fooof_results{p,d}(c).fooofed_spectrum);
%                     hold on;
%                     plot(fooof_results{p,d}(c,:).freqs,fooof_results{p,d}(c).ap_fit);
%                     str = append(num2str(p),' ',DEPTH{p,d},' ',label{p,d}(c));
%                     title(str);
%                     xlabel 'Frequency [Hz]';
%                     ylabel 'power';
%                     hold off;
%                     l = l+1;
                    fooof_results{p,d}(c).aperiodic_params = [];
                    fooof_results{p,d}(c).peak_params = [];
                    fooof_results{p,d}(c).gaussian_params = [];
                    fooof_results{p,d}(c).error = [];
                    fooof_results{p,d}(c).r_squared = [];
                    fooof_results{p,d}(c).freqs = [];
                    fooof_results{p,d}(c).power_spectrum = [];
                    fooof_results{p,d}(c).fooofed_spectrum = [];
                    fooof_results{p,d}(c).ap_fit = [];
                    errorcount_2 = errorcount_2 + 1;
                    bad_fit{p,d}(c) = label{p,d}(c);
                end    
            end
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
                
                % create powerspectrum without aperiodic component with
                % fooof 
                spectrum_wo_ap{x,1} = fooof_results{p,d}(c).power_spectrum - fooof_results{p,d}(c).ap_fit;
                
                % theta Power
                fooof{x,6} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>5&fooof_results{p,d}(c).freqs<7));
                % alpha Power
                fooof{x,7} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>7&fooof_results{p,d}(c).freqs<12));
                % beta Power
                fooof{x,8} = mean(spectrum_wo_ap{x,1}(fooof_results{p,d}(c).freqs>12&fooof_results{p,d}(c).freqs<30));
                x = x+1;
            end
        end
    end
end

T = cell2table(fooof,'VariableNames',{'ID','SIDE','DEPTH','CHANNEL','AP_EXPONENT','THETA_POWER','ALPHA_POWER','BETA_POWER'}); 
clear 'fooof';

% delete Channels with negative powers
errorcount_3 = 1;
c = 1;
while c ~= height(T)
    if (T.THETA_POWER(c) < 0) | (T.ALPHA_POWER(c) < 0) | (T.BETA_POWER(c) < 0)
        T(c,:) = [];
        errorcount_3 = errorcount_3 + 1;
    else
        c = c+1;
    end
end

%% safe for R 
writetable(T,'regression_table.csv');

%% descriptive statistic
% variance 
vrc = cat(1,vrc{:});
histogram(vrc,0.00001:0.00098:0.15);

%% cleaning up 
clear 'c' 'd' 'DEPTH' 'dif_neg' 'l' 'label' 'MAIN' 'p' 'PATHIN_conv' 'SIDE' 'str' 'x'


