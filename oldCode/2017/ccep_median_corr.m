function data = ccep_median_corr(data,times,ccep,bs_period)
% Baseline correction of data
% 
% function subtracts the average difference from the median response
% data = baseline(data,times,ccep,bs_period);
% data - data matrix
% times - time-points vector
% ccep - structure for which epochs belong to stimulated trials
%       should have:
%           ccep.els (stimulated electrodes)
%           ccep.epochs (epochs when els were stimulated)
% bs_period - period to get median from e.g [10 600]
%
%
% Hermes and Vassileva, 2016, UMC Utrecht


% set the baseline period
% bs_period = [10 600];
time_bs = times>bs_period(1) & times<bs_period(2);

for el=1:size(data,1)
    if numel(ccep(el).epochs)==0
        disp(['el ' int2str(el) ' not stimulated - skip'])
    else
        for elm = 1:size(data)
            median_ccep = nanmedian(data(elm,time_bs,ccep(el).epochs),3);
            for k = 1:length(ccep(el).epochs)
                trial_ccep = data(elm,time_bs,ccep(el).epochs(k));
                
                % what is the average difference between the trial and the
                % median? (average across time-points)
                av_difference = nanmean(trial_ccep - median_ccep);
                
                % correct the trial with this offset
                % only use this line for debuggin / plotting:
                corr_ccep = trial_ccep - av_difference;
                
                % correct the data with this offset
                data(elm,:,ccep(el).epochs(k)) = data(elm,:,ccep(el).epochs(k)) - av_difference;
            end
        end
        
%     disp(['el ' int2str(el)])
    end
end
disp('corrected for the difference with the median response; done')

     