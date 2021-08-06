vlim_l = 0.001; % define lower boundary for variance 

for v = 1:length(indat) % startup first 
%% read data + reject trials without data or with artefacts
    % read data
    cfg         = [];
    cfg.dataset = [PATHIN_conv indat(v).name];
    data{v}     = ft_preprocessing(cfg);        % read data unfiltered
    
    for c = 1:length(data{v}.label)                 % c = channels
        mnm = min(data{v}.trial{1,1}(c,:));         % lowest point 
        mxm = max(data{v}.trial{1,1}(c,:));         % hightest point

        % normalize data
        if mxm - mnm ~= 0
            norm_raw{v}(c,:) = (data{v}.trial{1,1}(c,:) - mnm) / (mxm - mnm);
        else
            norm_raw{v}(c,:) = 1;
        end

        % split data and calculate variance

        vrc1 = var(norm_raw{v}(c,1:ceil(end/4)));               % first quarter
        vrc2 = var(norm_raw{v}(c,ceil(end/4):ceil(end/2)));     % second quarter
        vrc3 = var(norm_raw{v}(c,ceil(end/2):ceil(end*0.75)));  % third quarter
        vrc4 = var(norm_raw{v}(c,ceil(end*0.75):end));          % fourth quarter 

        % replace data with nan if variance is smaller than 0.001
        if (vrc1 < vlim_l) || (vrc2 < vlim_l) || (vrc3 < vlim_l) || (vrc4 < vlim_l)        
            data{v}.trial{1,1}(c,:) = nan(size(data{v}.trial{1,1}(c,:)));
        end
    end
end