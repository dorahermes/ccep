function data = ccep_baseline_corr(data,times)
% Baseline correction of data
% subtracts the average of the baseline period (bs_period) from the data
% data = baseline(data,times);
% data - data matrix
% times - time-points vector
%
% Hermes and Vassileva, 2016, UMC Utrecht

% set the baseline period
bs_period = [-200 -100];

for el=1:size(data,1) % do baseline subtraction for each electrode
    base_epochs = squeeze(mean(data(el,times>bs_period(1) & times<bs_period(2),:),2)); %overall mean voltage between bs_period
    for l=1:size(data,3) % baseline subtracted from each epoch
        data(el,:,l) = data(el,:,l) -  base_epochs(l);
    end
end
disp('baseline corrected; done')
