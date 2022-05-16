# Legend

## 00_raw

* format: .mpx
* no raw data included, because they are too big

## 01_converted

* converted raw data to .mat

## 02_test

* few folders with few files to try out the scripts

## Output from s01_pipeline in 03_processed (patients are coded with numbers 1-25)

* data: Cleaned up but unfiltered data
* data_FFT: Filtered data (demean, high pass, low pass)
* error: Struct of Every rejected channel and their errors in pipeline
* TFR: results of the computed FFT using Hanning taper with 5 cycles
* fooof_results: parameters for the peaks and aperiodic component, that got extracted (for more explanation see https://fooof-tools.github.io/fooof/auto_tutorials/plot_04-MoreFOOOF.html#sphx-glr-auto-tutorials-plot-04-morefooof-py)

* DEPTH: position of the electrode (usually starts with 10, 0 is the targetpoint of the STN)
* SIDE: information in which hemisphere the electrode is (left(L) or right(R))
* TRAJECTORY: information in which trajectory the electrode is (T1: central; T2: anterior; T3: left side: medial, right side: lateral; T4: posterior; T5: left side: lateral, right side: medial)

## Output from s01_pipeline in 04_final (00_fooof_results.mat)

* error, fooof_results, DEPTH, SIDE for every file and every patient (rows: patient; columns: files)
* label: names of the channels for every file and every patient (rows: patient; columns: files)
* or_freq: frequency bins from original spectrum for every file and every patient (rows: patient; columns: files)
* rmsd: root mean square of raw signal for every file and every patient (rows: patient; columns: files)
* samples: number of samplepoints for every file and every patient (rows: patient; columns: files)
* vrc: variances of the splitted data for every file and every patient (rows: patient; columns: files)

## Output from s02_fooof_results in 04_final

* bad_fit: name of channels that have bad accordance between estimated aperiodic component and original powerspectrum 
* spectrum_wo_ap: original powerspectrum minus the estimated aperiodic power
* errorcount_1: number of deleted channels in s01_pipeline
* errorcount_2: number of deleted channels because of bad fit
* errorcount_3: number of deleted channels because of negative power (spectrum_wo_ap)
* errorcount_4: number of deleted channels with Depth = 10
* errorcount_5: number of deleted channels with Dept < -3
* errorcount_6: number of deleted channels with root mean square > 32

* T (saved in regression_table.csv and T.mat): table with ID, Side, Depth, Channel, samples, aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared and z-transformation of electrophysiolical variables (aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared) 
* target_dist_idx (saved in ttest_table.csv): values near target (closest to 0) and far away from target (closest to 10) for every electrophysiolical variable (aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared)
* nf_beta (saved in beta_depth_nf.csv): depth of the channels that were chosen for betapower near target and far away from target 

## Output from s03_exploration in 04_final

* T (saved in regression_table_exploring.csv): table with ID, Side, Depth, Channel, samples, aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared and z-transformation of electrophysiolical variables (aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared) but without using FOOOF
* target_dist_idx (saved in ttest_table_exploring.csv): values near target (closest to 0) and far away from target (closest to 10) for every electrophysiolical variable (aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared) but without using FOOOF
* nf_beta (saved in or_beta_depth_nf.csv): depth of the channels that were chosen for betapower near target and far away from target but without using FOOOF

## Output from s06_discussion in 04_final

* T (saved in regression_table_wobadfit.csv): table with ID, Side, Depth, Channel, samples, aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared and z-transformation of electrophysiolical variables (aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared) but without deleting channels with bad fit between aperiodic component and original powerspectrum
* target_dist_idx (saved in ttest_table_wobadfit.csv): values near target (closest to 0) and far away from target (closest to 10) for every electrophysiolical variable (aperiodic exponent, thetapower, alphapower, betapower, low-betapower, high-betapower, root mean squared) but without deleting channels with bad fit between aperiodic component and original powerspectrum