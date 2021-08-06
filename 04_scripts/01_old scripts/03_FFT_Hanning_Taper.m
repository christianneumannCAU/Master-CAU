for v = 1:length(data) % read and clean data first 
    %% FFT
    % preprocessing
    cfg             = [];
    cfg.demean      = 'yes';    %remove DC offset
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = .5;       %high-pass filter, cutting everything under .5 Hz
    cfg.hpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 45;       %low-pass filter, cutting everything over 45 Hz
    cfg.lpfilttype  = 'firws';
    
    try
        data_FFT{v}     = ft_preprocessing(cfg,data{v});
    catch ME
        continue
    end
    
    % FFT, Hanning taper
    cfg             = [];
    cfg.method      = 'mtmconvol'; 
    cfg.output      = 'pow';        % Output parameter
    cfg.foi         = [2:.05:35];   % Frequency resolution
    cfg.toi         = [0:.01: 5];   % Temporal resolution
    cfg.t_ftimwin   = 5./cfg.foi;
    cfg.taper       = 'hanning';    % Frequency-Adaptive Smoothing
    
    TFR{v}          = ft_freqanalysis(cfg,data_FFT{v});
    
    m{v} = mean(TFR{v}.powspctrm,[3],'omitnan'); %avarege powerspectrum over time
    
    % Plotting Option 1
 

    %normalize data
    for c = 1:length(data_FFT{v}.label)
        mnm = min(m{v}(c,:),[],'omitnan');
        mxm = max(m{v}(c,:),[],'omitnan');
        norm_FFT{v}(c,:) = (m{v}(c,:)- mnm)/(mxm - mnm);
    end

    figure; hold;
    plot(TFR{v}.freq,norm_FFT{v});
    str = ['powerspectrum FFT (Hanning)' ' ' DEPTH(v)];
    title(str) ;
    xlabel 'Frequency [Hz]';
    ylabel 'Power';
    lgd = legend(extractAfter(data_FFT{v}.label,10));
    lgd.NumColumns = length(data_FFT{v}.label);
    
end