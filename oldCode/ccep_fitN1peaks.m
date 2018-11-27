    
function [n1_parms,test_model_error,error_meas,corr_test,fit_curve, signal_train,signal_test,t_use] = fit_n1peaks(data,t,el_m,epochs,srate)
%
% function fit_n1peaks
%
% use:
% [n1_parms,rel_rmse,corr_test] = fit_n1peaks(data,t,el_m,epochs)
%
% inputs:
% data      = electrodes X time X epochs
% t         = time
% el_m      = electrode to get N1 from
% epochs    = epochs that another electrode was stimulated
%
% outputs:
% n1_parms  = amplitude, latency (ind of t_use), width, latency (ms)
% rel_rmse  = relative root mean squared error (mode error / test-train error)
% corr_tes  = correlation between test-set and prediction
%
% DH, April 2016, UMC Utrecht

% get time to use
t_ind = t>=10 & t<100;
t_use = t(t_ind);

% get signal to fit:
signal_use = double(squeeze(data(el_m,t_ind,epochs))); % time X epochs

% z-score:
signal_use = (signal_use-mean(signal_use(:)))./std(signal_use(:));

% % detrend signal
% for k=1:size(signal_use,2)
%     signal_use(:,k) = detrend(signal_use(:,k));
% end

% find the max and set to zero
for k=1:size(signal_use,2)
    [m_val,m_ind,m_w]=findpeaks(signal_use(:,k),...
        'MinPeakDistance',round(.08*srate),'SortStr','descend');
    signal_use(:,k) = bsxfun(@minus,signal_use(:,k),m_val(1));
end
% % set first value to zero (subtract the first time-point):
% signal_use = bsxfun(@minus,signal_use,signal_use(1,:));

% training/testing data
signal_train = signal_use(:,1:2:size(signal_use,2));
signal_test = signal_use(:,2:2:size(signal_use,2));

% now start the 'fitting'
% find 1 negative peak in the training data:
x = -median(signal_train,2);
x = smooth(x,round(srate*.010),'lowess');
% we need to both smooth the signal and set a minimum peak distance:
% Smoothing: makes the peak width detected ok
% Distance: only selects the maximum peak in a range

[m_val,m_ind,m_w]=findpeaks(x,t_use,...
    'Annotate','extents','WidthReference','halfheight',...
    'MinPeakDistance',50);

if ~isempty(m_val)
    n1_parms = [m_val(1),m_ind(1),m_w(1)];
    % FITTED CURVE SHAPE:
    n = 1;
    peak_width = m_w(n)/2;%m_w(n)/2;
%     fit_curve = m_val(n)*exp(-(([1:length(t_use)] - m_ind(n))/peak_width).^2);
    fit_curve = m_val(n)*exp(-(t_use - m_ind(n)).^2/(2*peak_width.^2));
    fit_curve = -fit_curve';

    % for the error calculations, select the N1 duration, otherwise we are
    % including the N2
%     start_point = round(max([1 n1_parms(2)-2*n1_parms(3)])); % first time point or latency - 2sd
%     end_point = round(min([length(signal_train) n1_parms(2)+2*n1_parms(3)])); % last time point or latency + 2sd
    start_point = 1; % this is the first time-point in t_use 
    end_point = find(t_use>=60,1); % use t = 60
    signal_test_4error = median(signal_test(start_point:end_point,:),2);
    signal_train_4error = median(signal_train(start_point:end_point,:),2);
    fit_curve_4error = fit_curve(start_point:end_point);
    
    % error
    test_train_error = sqrt(sum((signal_test_4error - signal_train_4error).^2,1))/length(signal_test_4error);
    % model error
    test_model_error = sqrt(sum((signal_test_4error - fit_curve_4error).^2,1))/length(signal_test_4error);
    
    % relative RMS error FOR N1 PART:
    rel_rmse = test_model_error./test_train_error;

    error_meas(1) = rel_rmse;
    
    % total area under the curve:
    total_undercurve = sum(abs(signal_test_4error));
    % total absolute error: fitted curve - test_signal
    total_abserror = sum(abs(signal_test_4error-fit_curve_4error));
    % how much error wrt data (normalized absolute error
    error_meas(2) = total_abserror./total_undercurve;

    % dot product (can translate to angle between curves)
    error_meas(3) = signal_test_4error'*fit_curve_4error./...
        (sqrt(sum((signal_test_4error.^2)))*sqrt(sum((fit_curve_4error.^2))));

    % correlation with test set:
    corr_test = corr(signal_test_4error,fit_curve_4error);

else
    n1_parms = [NaN NaN NaN NaN];
    error_meas = [NaN NaN NaN];
    corr_test = NaN;
end




