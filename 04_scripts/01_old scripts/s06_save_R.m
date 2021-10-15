% save for R
        try
            T       = struct2table(fooof_results{v},'RowNames',data{v}.label,'AsArray',true); 
        catch ME
            continue
        end
        T       = removevars(T,{'gaussian_params','freqs','power_spectrum','fooofed_spectrum','ap_fit'});
        title   = [strcat(SIDE(v),TRAJECTORY(v),DEPTH(v),'.txt')];
        writetable(T,string(title));