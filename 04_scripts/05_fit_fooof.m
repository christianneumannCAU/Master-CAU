settings = [];
settings.peak_width_limits = [1 12]
freqs = TFR{1, 1}.freq;
power_spectrum = m{1,1}(1,:);
missing = find(isnan(power_spectrum)==1);
power_spectrum(isnan(power_spectrum))= power_spectrum(missing+1);
f_range = [1 30];
return_model = 1;

fooof_results = fooof(freqs, power_spectrum , f_range ,settings , return_model);