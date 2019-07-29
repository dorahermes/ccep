function ccep_parms = ccep_setFitParms(times)

%% Function fitting parameters
% Pick Times for each component
% t_low_N1 - lower bounds for N1
% t_up_N1 - upper bound for N1

ccep_parms.t_low_N1 = 10;
ccep_parms.t_up_N1 = 80;

%%
ccep_parms.times2fitN1 = times(1,times>ccep_parms.t_low_N1 & times<ccep_parms.t_up_N1);
ccep_parms.alltN1 = size(ccep_parms.times2fitN1,2);

% options for lsqnonlin:
ccep_parms.my_options=optimset('Display','off','Algorithm','trust-region-reflective'); 

%% N1 Parameters for the fit:
%  amplitude, time, width, offset
ccep_parms.X01 = [-1 30 1 0];    % starting points
ccep_parms.LB1 = [-Inf 10 0 -Inf]; % lower bounds
ccep_parms.UB1 = [ 0 40 150 Inf]; % upper bounds

%% N2 parameters for the fit
%  amplitude, time, width
ccep_parms.X02 = [-1 40 1];   % starting points
ccep_parms.LB2 = [-Inf 0 50]; % lower bounds
ccep_parms.UB2 = [0 Inf 700]; % upper bounds


