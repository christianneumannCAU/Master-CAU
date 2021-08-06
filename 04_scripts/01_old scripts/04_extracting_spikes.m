for v = 1:length(data) % read and clean data first 
%% Extracting Spikes (Rey, Pedreira & Quiroga, 2015)
        cfg             = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [300 3000];           %bandpass-filter
        cfg.bpfilttype  = 'firws';
        
        try
            data_spikes{v}   = ft_preprocessing(cfg,data{v});
        catch ME
            continue
        end

        for c = 1:length(data_spikes{v}.label)   % c = channel
            %extract time of spikes
            threshold = (median(abs(data_spikes{v}.trial{1,1}(c,:))))/0.6745;
            spike_time{v} = spike_detection(data_spikes{v}.trial{1,1}(c,:),threshold);

            %plot spikes and mark every spike with a star
            try
                figure; hold;
                plot(data_spikes{v}.time{1,1},data_spikes{v}.trial{1,1}(c,:));
                plot(data_spikes{v}.time{1,1}(spike_time{v}),0,'*');
                str = [data_spikes{v}.label(c) ' ' DEPTH(v)];
                title(str) ;
            catch ME
                continue
            end
        end
end