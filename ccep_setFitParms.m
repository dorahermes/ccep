function ccep_parmsN1fit = ccep_setFitParms(t)

%% Function fitting parameters
% Pick t for each component
% t_low_N1 - lower bounds for N1
% t_up_N1 - upper bound for N1
% t_low_N2 - lower bounds for N2
% t_up_N2 - upper bounds for N2

% N1:
ccep_parmsN1fit.t_low_N1 = 10;
ccep_parmsN1fit.t_up_N1 = 50;
% N2:
ccep_parmsN1fit.t_low_N2 = 50;
ccep_parmsN1fit.t_up_N2 = 550;

%%
ccep_parmsN1fit.t2fitN1 = t(1,t>ccep_parmsN1fit.t_low_N1 & t<ccep_parmsN1fit.t_up_N1);
ccep_parmsN1fit.alltN1 = size(ccep_parmsN1fit.t2fitN1,2);

ccep_parmsN1fit.t2fitN2 = t(1,t>ccep_parmsN1fit.t_low_N2 & t<ccep_parmsN1fit.t_up_N2);
ccep_parmsN1fit.alltN2 = size(ccep_parmsN1fit.t2fitN2,2);

% options for lsqnonlin:
ccep_parmsN1fit.my_options=optimset('Display','off','Algorithm','trust-region-reflective'); 

%% N1 Parameters for the fit:
%  amplitude, time, width, offset
ccep_parmsN1fit.X01 = [-1 30 1 0];    % starting points
ccep_parmsN1fit.LB1 = [-Inf 10 0 -Inf]; % lower bounds
ccep_parmsN1fit.UB1 = [ 0 40 150 Inf]; % upper bounds

%% N2 parameters for the fit
%  amplitude, time, width
% ccep_parmsN1fit.X02 = [-1 40 1];   % starting points
% ccep_parmsN1fit.LB2 = [-Inf 0 50]; % lower bounds
% ccep_parmsN1fit.UB2 = [0 Inf 700]; % upper bounds