# Folder structure

## 01_old scripts *(multiple scripts that were made to fulfill specific tasks)*

* explore_spectra_FOOOF.py *(python-script by Julius to compare python and matlab)*
* 00_compare_FTT.m *(computing and plotting FFTs with different number of cycles for Hinning Tapers vs Wavelet)*
* 01_startup.m *(prepare following loops)* 
* 02_read_data.m *(read and clean up data)*
* 03_FFT_Hanning_Taper.m *(Filter data, compute and plot FFT)*
* 04_extracting_spikes.m *(extract and plot spikes)*
* 05_fit_fooof.m *(compute fooof results)*
* 06_save_R.m *(script by Julius to save data for R)*
* 07_stats.m *(script by Julius to prepare analysis)*
* test_lfp_rms_jw.m *(script by Julius to look at root mean square and lfp for depth near 0 and depth near 10)*


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
* create new table for correlation
* delete channels with negative powers, depth >= 10, depth < -3 and root mean sqaure > 32
* add z-transformed data 
* build new table with values near 0 and far from 0 for ttests
* save data

## s03_statistics.R

* visualize data to check normal distribution (aperiodic exponent, root mean square, theta-, alpha- and betapower)
* visualize z-transformed data to compare
* visualize log-transformed data (only theta-, alpha- and betapower)
* paired t-tests for H1 (compare values near 0 with values near 10 for betapower and root mean square)
* paired t-test and correlation test between Depth and aperiodic exponent for H2
* paired t-tests for exploration (Theta-, Alpha-, low-beta- and high-betapower)
* correlation matrix for exploration (Depth, aperiodic exponent, theta-, alpha-, betapower and root mean square
* visualize found correlations


 

