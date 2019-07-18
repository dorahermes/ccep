
function [n1_mat,n2_mat,rmseN1,rmseN2,fitted_line_n1,fitted_line_n2] = ccep_fitN1N2(paramsCCEPfit,ccep,data,times,el,elm)



%% N1 and N2 component detection. Fits a line plus gaussian with preset parameters.
% Each component is detected separately. Best fit is estimated by
% non-linear least squares estimate. Outputs to structs that contain N1 and
% N2 data:
% 1st dimension - stimulated electrode
% 2nd dimension - measured electrode
% 3rd dimension - (1st value - peak amplitude (microV), 2nd value - width (ms), 3rd value - latency (ms))
%
% Additionally, for each component a root mean square error is
% created, to evaluate the goodness of fit for each component:
% 1st dimension - stimulated electrode
% 2nd dimension - measured electrode
%
% Vassileva and Hermes, 2016, UMC Utrecht


%% Define data and times to fit for each component:

% N1:
data2fitN1 = squeeze(nanmean(data(elm,times>paramsCCEPfit.t_low_N1 & times<paramsCCEPfit.t_up_N1,ccep(el).epochs),3));
data2fitN1 = double(data2fitN1);

% % N2:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data2fitN2 = squeeze(nanmean(data(elm,times>paramsCCEPfit.t_low_N2 & times<paramsCCEPfit.t_up_N2,ccep(el).epochs),3));
data2fitN2 = double(data2fitN2);

if find(isnan(data2fitN1),1)>0
    disp(['WARNING: input data N1 has NANs, fitting function is going to crash at electrode ' int2str(elm)])
end
if find(isnan(data2fitN2),1)>0
    disp(['WARNING: input data N2 has NANs, fitting function is going to crash at electrode ' int2str(elm)])
end

%% Detect N1

% calculate least squares
[a, resNorm] =...
  lsqnonlin(@(x) LinePlusGauss(x,data2fitN1,paramsCCEPfit.times2fitN1),paramsCCEPfit.X01,paramsCCEPfit.LB1,paramsCCEPfit.UB1,paramsCCEPfit.my_options);

% fit a line
fitted_line_n1 = a(1)*sqrt(2*pi)*normpdf(paramsCCEPfit.times2fitN1,a(2),a(3));

% calculate outputs
amp = a(1);                  % amplitude
funWidth = a(3); % width
peakLat = a(2);                   % peak latency
% amp = a(1)/a(3);                  % amplitude
% funWidth = 2*sqrt(2*log(2))*a(3); % width
% peakLat = a(2);                   % peak latency

% write stuff down:
n1_mat = [amp funWidth peakLat];
rmseN1 = sqrt(resNorm/paramsCCEPfit.alltN1); % root mean square error
% 
clear a amp funWidth peakLat resNorm

%% Plot N1 data + fitted line
% figure,
% hold on
% plot(times2fitN1,data2fitN1,'k')
% plot(times2fitN1,fitted_line_n1,'r')
% xlabel('ms')
% ylabel('mV')
% title(['N1 ' int2str(el) ' to ' int2str(elm) ' subj ' s])

%% Detect N2

% calculate least squares
[b, resNorm] =...
    lsqnonlin(@(x) LinePlusGauss(x,data2fitN2,paramsCCEPfit.times2fitN2),paramsCCEPfit.X02,paramsCCEPfit.LB2,paramsCCEPfit.UB2,paramsCCEPfit.my_options);

% fit a line
fitted_line_n2 = b(1)*sqrt(2*pi)*normpdf(paramsCCEPfit.times2fitN2,b(2),b(3));

% calculate parameters
amp = b(1)/b(3);                  % amplitude
funWidth = 2*sqrt(2*log(2))*b(3); % width
peakLat = b(2);                   % peak latency

% write stuff down again
n2_mat = [amp funWidth peakLat];
rmseN2 = sqrt(resNorm/paramsCCEPfit.alltN2);  % root mean square

clear b amp funWidth peakLat resNorm X0 LB UB

%% Plot N2 data + fitted line
% figure,
% hold on
% plot(times2fitN2,data2fitN2,'k')
% plot(times2fitN2,fitted_line_n2,'g')
% xlabel('ms')
% ylabel('mV')
% title(['N2 ' int2str(el) ' to ' int2str(elm)])
