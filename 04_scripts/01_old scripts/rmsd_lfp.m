function rmsd = rmsd_lfp(dat)
% rmse_lfp: Function to calculate the root mean square deviation from the
% data mean
%     
% INPUT: 
%     lfp: any 1 dim vector
% OUTPUT:
%     rmse: average absolute deviation from data mean
% 
% 
% Author: (Julius Welzel, University of Kiel, 2022)

sz_dat = size(dat);
if sum(sz_dat > 1) > 1;
    error 'Data not 1 dimensional'
end

dat_mean    = mean(dat)
rmsd        = sqrt(mean((dat - dat_mean).^2));

end

