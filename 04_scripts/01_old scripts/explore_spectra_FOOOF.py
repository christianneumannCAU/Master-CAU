# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 19:41:19 2021

@author: User
"""

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


import numpy as np
data = np.loadtxt('data.txt')


import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=1.2)

# Define sampling frequency and time vector
sf = 100.
time = np.arange(data.size) / sf

# Plot the signal
fig, ax = plt.subplots(1, 1, figsize=(12, 4))
plt.plot(time, data, lw=1.5, color='k')
plt.xlabel('Time (seconds)')
plt.ylabel('Voltage')
plt.xlim([time.min(), time.max()])
plt.title('N3 sleep EEG data (F3)')
sns.despine()


from scipy import signal

# Define window length (4 seconds)
win = 4 * sf
freqs, psd = signal.welch(data, sf, nperseg=win)

# Plot the power spectrum
sns.set(font_scale=1.2, style='white')
plt.figure(figsize=(8, 4))
plt.plot(freqs, psd, color='k', lw=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectral density (V^2 / Hz)')
plt.ylim([0, psd.max() * 1.1])
plt.title("Welch's periodogram")
plt.xlim([0, freqs.max()])
sns.despine()

from fooof import FOOOF
from fooof.plts.spectra import plot_spectrum

fm = FOOOF(peak_width_limits = [1.5,12],peak_threshold = 2.0)
fm.fit(freqs,psd,[1,35])

fm.plot(plot_peaks='shade', peak_kwargs={'color' : 'green'},plt_log=False)

idx_f_start = find_nearest_idx(freqs,1)
idx_f_end = find_nearest_idx(freqs,35)

plot_spectrum(freqs[idx_f_start:idx_f_end], np.log10(psd[idx_f_start:idx_f_end]), log_freqs=False, log_powers=False,
              color='black', label='Original Spectrum')