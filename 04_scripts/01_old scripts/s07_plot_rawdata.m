    %     % Plotting FFT
    % 
    %     %normalize data
    %     for c = 1:length(data_FFT{v}.label)
    %         mnm = min(m{v}(c,:),[],'omitnan');
    %         mxm = max(m{v}(c,:),[],'omitnan');
    %         norm_FFT{v}(c,:) = (m{v}(c,:)- mnm)/(mxm - mnm);
    %     end
    % 
    %     figure; hold;
    %     plot(TFR{v}.freq,norm_FFT{v});
    %     str = ['powerspectrum FFT (Hanning)' ' ' DEPTH(v)];
    %     title(str) ;
    %     xlabel 'Frequency [Hz]';
    %     ylabel 'Power';
    %     lgd = legend(extractAfter(data_FFT{v}.label,10));
    %     lgd.NumColumns = length(data_FFT{v}.label);