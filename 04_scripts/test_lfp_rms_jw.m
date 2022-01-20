close all
nms_data = {'LT1D0.008F0001','LT1D8.791F0001'}
cp = 1;

for f = 1:numel(nms_data)
    
    load([nms_data{f} '.mat']);
    
    lfp = CMacro_LFP_01___Central;
    raw = CSPK_01___Central;

    srate_lfp = CMacro_LFP_01___Central_KHz * 1000;
    srate_raw = CSPK_01___Central_KHz * 1000;

    tvec_lfp    = linspace(0,numel(lfp)/srate_lfp,numel(lfp));
    tvec_raw    = linspace(0,numel(raw)/srate_raw,numel(raw));
    subplot(2,2,cp)
    rms = rmsd_lfp(raw);
    plot(tvec_raw,raw)
        xlabel 'Time [s]'
        ylabel '\muV'
        title (strjoin(['Depth: ' extractBetween(nms_data(f),'D','F') ', RMS ' num2str(rms)],''))
        axis tight
        cp = cp+1;
    
    subplot(2,2,cp)
    plot(tvec_lfp,lfp)
        xlabel 'Time [s]'
        ylabel '\muV'
        title 'LFP'
        axis tight
        cp = cp+1;
        
end