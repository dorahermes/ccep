
% first run the preprocessing in master02

%% %%%%%%%
%% %%%%%%% Now, we are going to fit some data
%% %%%%%%%

el = 70;

%% show all CCEPs in 1 plot
% get time to use
t_ind = t>=10 & t<1000;
t_use = t(t_ind);

figure,hold on
for el_m = 1:size(data,1)
    signal_use = double(squeeze(data(el_m,t_ind,ccep(el).epochs))); % time X epochs
    signal_plot = mean(signal_use,2);
    signal_plot = zscore(signal_plot);    
    plot(t_use,signal_plot)
end

%% fit 
% define output matrices
% n1_mat_all: parameters matrix: els X 3
n1_mat_all = NaN(size(data,1),3); % amplitude, latency (ms), width
% rmseN1_all: fit quality matrix: els X 2
rmseN1_all = NaN(size(data,1),5); 
 
[~,~,~,~,fit_curve] = ccep_fitN1peaks(data,t,1,ccep(el).epochs,srate);    
fitted_curves = NaN(size(data,1),length(fit_curve));

for el_m = 1:size(data,1)
    [n1_parms,test_model_error,error_meas,corr_test,fit_curve] =...
        ccep_fitN1peaks(data,t,el_m,ccep(el).epochs,srate);
    n1_mat_all(el_m,1:3)    = n1_parms;
    rmseN1_all(el_m,1)      = test_model_error;
    rmseN1_all(el_m,2)      = corr_test;
    rmseN1_all(el_m,3:5)    = error_meas;
    fitted_curves(el_m,:)   = fit_curve;
end
clear el_m

%% plot fitted curve for one electrode set:

el_m = 74;

[n1_parms,test_model_error,error_meas,corr_test,fit_curve, signal_train,signal_test,t_use] = ...
     ccep_fitN1peaks(data,t,el_m,ccep(el).epochs,srate);

figure
subplot(2,1,1),hold on %%%% SIGNAL TO FIT
plot(t_use,signal_train,'k')
plot(t_use,signal_test,'k')
xlabel('time (ms)')
ylabel('amplitude (z-score)')

title(['stim el ' int2str(el) ' measure el ' int2str(el_m) ' all trials'])
% title(['rel RMSE ' num2str(rmseN1_all(el_m,3)) ...
%     ' Norm error ' num2str(rmseN1_all(el_m,4)) ...
%     ' corr ' num2str(rmseN1_all(el_m,2))])

subplot(2,1,2),hold on %%%% SIGNAL TO FIT

% FITTED CURVE SHAPE:
plot(t_use,median(signal_train,2),'k')
plot(t_use,median(signal_test,2),'k--')
% plot(t_use,fitted_curves(el_m,:),'r')
plot(t_use,fit_curve,'m')
xlabel('time (ms)')
ylabel('amplitude (z-score)')
legend({'training data','testing data','fitted curve'})

title(['R^2 = ' num2str(rmseN1_all(el_m,2).^2) ...
    ' Norm error = ' num2str(rmseN1_all(el_m,4)) ...
    ' Amp = ' num2str(n1_mat_all(el_m,1)) ', lat = ' num2str(n1_mat_all(el_m,2))])

set(gcf, 'PaperPositionMode', 'auto');
% print('-painters','-r300','-dpng',strcat(['./figures/fast_fit']));
% print('-painters','-r300','-depsc',strcat(['./figures/fast_fit']));


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

for l = 1:length(chans2inc)
    el_m = chans2inc(l);
    
    [n1_parms,test_model_error,error_meas,corr_test,fit_curve, signal_train,signal_test,t_use] = ...
        ccep_fitN1peaks(data,t,el_m,ccep(el).epochs,srate);

    subplot(8,ceil(length(chans2inc)/8),l),hold on
    
    % FITTED CURVE SHAPE:
    plot(t_use,mean(signal_train,2),'k')
    plot(t_use,mean(signal_test,2),'k--')
    plot(t_use,fitted_curves(el_m,:),'r')
    ylim([-5 5])
end

%% plot CCEPs accroding to visual hierarchy

figure('Position',[0 0 800 300]),hold on

cmap = hsv(28);

for el_m = 1:size(n1_mat_all,1) % electrodes
    if ~isnan(labels_mode(el_m)) % if this electrode has a label
        plot(labels_mode(el_m),n1_mat_all(el_m,2),'.','MarkerSize',20,'Color',cmap(round(labels_mode(el_m)),:))
    end
end
xlim([0 26])
ylabel('latency')
set(gca,'XTick',[1:25],'XTickLabel',Wang_ROI_Names)

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',strcat(['./figures/VisualHierarchy/firstLatency']));
% print('-painters','-r300','-depsc',strcat(['./figures/VisualHierarchy/firstLatency']));


figure('Position',[0 0 800 600])
for el_m = 1:size(n1_mat_all,1) % electrodes
    if ~isnan(labels_mode(el_m)) % if this electrode has a label
        subplot(5,5,round(labels_mode(el_m))),hold on
        [n1_parms,test_model_error,error_meas,corr_test,fit_curve, signal_train,signal_test,t_use] = ...
            ccep_fitN1peaks(data,t,el_m,ccep(el).epochs,srate);
        plot(t_use,median(signal_train,2),'k')
%         plot(t_use,mean(signal_test,2),'k--')
        plot(t_use,fitted_curves(el_m,:),'r')
    end
end
for k = 1:25
    subplot(5,5,k)
    title(Wang_ROI_Names{k})
end

%     find(labels_mode == k)
% end