%make sure, that the version of python in cmd and in python match each
%other
for v = 1:length(TFR)
    settings = [];
    settings.peak_width_limits = [1.5 12]; %minimum and maximum widths of etracted peaks
    settings.peak_threshold = 2; %standard deviation of the aperiodic-removed powerspectrum, above which a data point must pass to be considered a candidate peak
    f_range = [1 35]; %fitting range
    return_model = 1;
    freqs{v} = TFR{v}.freq;
    
    for c = 1:length(data_FFT{v}.label)
        power_spectrum{v}(c,:) = m{v}(c,:);
        try
            fooof_results{v}(c,:) = fooof(freqs{v}, power_spectrum{v}(c,:) , f_range ,settings , return_model);
        catch ME
            continue
        end
    end

end