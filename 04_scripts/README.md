# File structure

## 01_old scripts *(multiple scripts that were made to fulfill specific tasks)*

* explore_spectra_FOOOF.py *(python-script by Julius to compare python and matlab)*
* 01_startup.m *(prepare following loops)* 
* 02_read_data.m *(read and clean up data)*
* 03_FFT_Hanning_Taper.m *(Filter data, compute and plot FFT)*
* 04_extracting_spikes.m *(extract and plot spikes)*
* 05_fit_fooof.m *(compute fooof results)*
* 06_save_R.m *(script by Julius to save data for R)*
* 07_stats.m *(script by Julius to prepare analysis)*
* test_lfp_rms_jw.m *(script by Julius to look at root mean square and lfp for depth near 0 and depth near 10)*

## s00_compare_FTT.m 

* computing and plotting FFTs with different number of cycles for Hinning Tapers vs Wavelet

## s01_pipeline

* read data
* clean up data
* filter data
* FFT with Hanning taper
* calculate powerspectrum
* FOOOF-Algorithm to seperate aperiodic and periodic components
* save data

## s02_fooof_results

* counting errors from s01_pipeline
* delete channels with bad fit
* create new table for correlation/regression
* delete channels with negative powers, depth >= 10, depth < -3 and root mean sqaure > 32
* add z-transformed data 
* build new table with values near 0 and far from 0 for t-tests
* for betapower: also extract at which depth the values near 0 and far from 0 were extracted
* save data

## s03_exploration

* calculating electrophysiological variables again, but with the original powerspectrum (without using FOOOF)

## s04_plots

* plot raw data to look for artifacts
* plot remaining channels to see if criteria for cleaning data are chosen well
* plot histogram for all variances to see if cut-off was chosen well
* plot aperiodic component and powerspectrum for channels near target and far from target to evaluate the estimate of the aperiodic component

## s05_statistics.R

* set libraries and read data
* visualize data to check normal distribution (aperiodic exponent, root mean square, theta-, alpha- and betapower)
* visualize z-transformed data to compare
* visualize log-transformed data to compare (only theta-, alpha- and betapower)
* H1: Paired t-tests (compare values near 0 with values near 10 for betapower and root mean square)
* H2: Paired t-test between Depth and aperiodic exponent and regression with linear mixed model to look for correlation between aperiodic exponent and depth of electrode
* Exploration: Regression for every predictor possible (full model), paired t-tests (Theta-, Alpha-, low-beta- and high-betapower) and compare beta near target with far from target again but with original powerspectrum (without using FOOOF)
* Discussion: Compare depth of channels near target for original powerspectrum vs powerspectrum after fooof, boxplot for means when using fooof but without deleting channels with bad fit

## s06_discussion

* calculating electrophysiological variables again, but without deleting channels with a bad fit between aperiodic component and original powerspectrum
 

