
% first run the preprocessing in master02

%% %%%%%%%
%% %%%%%%% Now, we are going to fit some data
%% %%%%%%%

%% N1 and N2 detection
close all

% Define electrodes
el = 71; % stimulated electrode

% Preallocate matrices to contain components' amplitude, width and peak latency.
% 1st dim - stimulated electrode
% 2nd dim - measured electrode
% 3rd dim - 4 fitted line parms: a(1)*sqrt(2*pi)*normpdf(t,a(2),a(3)) + a(4); 
% a(1) - uncorrected Amplitude, a(2) - Latency, a(3) uncorrected Width,...
% a(4) - Offset, 4 - corrAmplitude, 5 - corrWidth)
n1_mat_all = NaN(size(data,1),size(data,1),6); % parameters matrix


% Preallocate matrices for error measurements:
rmseN1_all = NaN(size(data,1),size(data,1),2);

%% And now run fitting procedure

ccep_parms = ccep_setFitParms(t); % load fit parameters (times to fit, starting points, lower/upper bounds)

tic
for elec_ms = 1:size(data,1)
    if mod(elec_ms,10) == 0, disp(['el ' int2str(elec_ms) ' of ' int2str(size(data,1))]),end
    data2fit = squeeze(nanmean(data(elec_ms,:,ccep(el).epochs),3));
    if sum(isnan(data2fit))==0 % check whether this epoch was good 
        [n1_mat,rmseN1,fitted_line_n1] = ...
            ccep_fitN1(ccep_parms,ccep,data,t,el,elec_ms);
        n1_mat_all(el,elec_ms,:) = n1_mat;
    end
end
toc

clear data2fit fitted_line_n1 

%% plot fitting result for one measured electrode

elm = 13;
t_min = -100;

t_n1 = t>ccep_parms.t_low_N1 & t<ccep_parms.t_up_N1;

% fitted line
n1_elm = squeeze(n1_mat_all(el,elm,:));
fitted_line_n1 = n1_elm(1)*sqrt(2*pi)*normpdf(ccep_parms.times2fitN1,n1_elm(2),n1_elm(3)) + n1_elm(4);

figure,hold on
t_plot = t(t>t_min);
plot(t_plot,zeros(size(t_plot)),'k:')

data_plot = squeeze(data(elm,t>t_min,ccep(el).epochs));
plot(t_plot,data_plot,'Color',[.9 .9 .9])

data_plot = squeeze(mean(data(elm,t>t_min,ccep(el).epochs),3));
plot(t_plot,data_plot,'k')

plot(t(t_n1),fitted_line_n1,'r')

%% render figure
figure
ctmr_gauss_plot(cortex,[0 0 0],0)

n1_amp = n1_mat_all(el,:,5);
n1_amp(isnan(n1_amp)) = 0;
n1_lat = n1_mat_all(el,:,2);
n1_lat(isnan(n1_lat)) = 0;

el_add_sizecolor(elecmatrix,n1_amp,n1_lat,[200],40)
el_add(elecmatrix(ccep(el).els,:),'k',50) % stimulated electrodes


%% plot fitting result for many measured electrodes
figure,hold on
t_min = 5;
t_plot = t(t>t_min);
t_n1 = t>ccep_parms.t_low_N1 & t<ccep_parms.t_up_N1;

for l = 1:length(chans2inc)
    elm = chans2inc(l);
    subplot(8,ceil(length(chans2inc)/8),l),hold on
    % fitted line
    n1_elm = squeeze(n1_mat_all(el,elm,:));
    fitted_line_n1 = n1_elm(1)*sqrt(2*pi)*normpdf(ccep_parms.times2fitN1,n1_elm(2),n1_elm(3)) + n1_elm(4);

%     plot(t_plot,zeros(size(t_plot)),'k:')

%     data_plot = squeeze(data(elm,t>t_min,ccep(el).epochs));
%     plot(t_plot,data_plot,'Color',[.9 .9 .9])

    data_plot = squeeze(mean(data(elm,t>t_min,ccep(el).epochs),3));
    plot(t_plot,data_plot,'k')

    plot(t(t_n1),fitted_line_n1,'r')
end
%%
for l = 1:length(chans2inc)
    subplot(8,ceil(length(chans2inc)/8),l),hold on
    xlim([0 100])
    ylim([-1000 1000])
end