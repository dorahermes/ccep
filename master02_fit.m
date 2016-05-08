%% Master script to run fitting procedure
% First run the preprocessing in master01.

%% Run fitting procedure for each trial separately

% Load fit parameters (t to fit, starting points, lower/upper bounds):
ccep_parmsN1fit = ccep_setFitParms(t);

% This script creates vectors named after the anatomical locations. Each
% vector contains the electrode numbers that correspond to that location.
% To change electrode numbers, edit script.
[IFG, SMG, STG, els2check] = ccep_anat_locs(s);

% Preallocate matrices:
n1_mat_epochs = NaN(size(data,1),10,6);
n1_ttest = NaN(size(data,1),2);
rmse_ep = NaN(size(data,1),10);

% This data structure contains the number of each electrode, its anatomical
% location, the fitting data, the t-test data and the fitting error :
n1_data = struct('el', [], 'loc',[], 'fit', [], 'ttest', [], 'rmse', []);
for k = 1:size(data,1)
    n1_data(k).el = k;
end
% Write down the location corresponding to each electrode in the data structure:
for loc = IFG
    n1_data(loc).loc = 'IFG';
end
for loc = STG
    n1_data(loc).loc = 'STG';
end
% for loc = MTG
%     n1_data(loc).loc = 'MTG';
% end
% for loc = AG
%     n1_data(loc).loc = 'AG';
% end
for loc = SMG
    n1_data(loc).loc = 'SMG';
end

clear k loc

% set alpha for t-test
alphaVal = 0.05/length(els2check);

tic

%%% loop through stimulated electrodes
for el = els2check
    
    if numel(ccep(el).epochs) == 0 % skip this electrode, if it wasn't stimulated.
        disp(['el ' int2str(el) ' not stimulated'])
    else
        
%%%%%%%%%%% loop through measured electrodes
        for elm = els2check 
            
            % Calculate a threshold to check against in the t-test.
            % 2 standard deviations of the pre-stimulation baseline between
            % -500 and -100 ms
            data2check = squeeze(nanmean(data(elm,:,ccep(el).epochs),3));
            signThr = -2*std(data2check(1,t > -500 & t < -100),[],2);
            
%%%%%%%%%%%%%%% loop through epochs
            for ep = 1:10 
                
                % check if epoch has NaNs and if not run fit and test procedures
                epochCheck = squeeze(nanmean(data(elm,:,ccep(el).epochs(ep)),3));
                if sum(isnan(epochCheck))==0
                    
                %%% fit each epoch individually %%%
                    [n1_vec,rmseN1,~] = ...
                        ccep_fitN1_ind_epochs(ccep_parmsN1fit,ccep,data,t,el,elm,ep);
                    
                    n1_mat_epochs(elm,ep,:) = n1_vec; % write down ccep data
                    rmse_ep(elm,ep) = rmseN1; % write down error
                    
                end
            end
            % Check the fitting errors:
            errorCheck = std(rmse_ep(elm,:)); % get standard deviation of fitting errors
            a = find(rmse_ep(elm,:)>2*errorCheck); % check if any epoch fit exceeds fit error SD
            % If any epochs do not pass the error check and they're less
            % than half the total number of epochs, fill them with NaNs
            if ~isempty(a)==1 && numel(a)<(size(n1_mat_epochs,2))/2
                n1_mat_epochs(elm,a,:) = NaN;
            end
            clear errorCheck a
            if nanmean(n1_mat_epochs(elm,:,2),2) > 15 % Check if CCEP peaks after 15 ms.
                % t-test for the amplitudes of the epochs against the calculated significance threshold
                % and store the p-value and t-statistic:
                [~,p,~,stats] = ttest(n1_mat_epochs(elm,:,5), signThr, 'Dim', 2, 'Tail', 'left', 'Alpha', alphaVal);
                n1_ttest(elm,:) = [stats.tstat p];
            end
            % Write data in structure. If you don't want to fit and test all electrodes,
            % it skips them, but keeps structure with the same size as data size
            % (i.e. number of electrodes). Structure fields for untested electrodes
            % remain empty
            n1_data(el).el = el;
            n1_data(el).fit = n1_mat_epochs;
            n1_data(el).ttest = n1_ttest;
            n1_data(el).rmse = rmse_ep;
        end
    end
    disp(['fitting and testing el ' int2str(el)])
end
clear n1_mat_epochs n1_ttest rmse_ep el elm ep stats p n1_vec rmseN1 epochCheck data2check signThr
disp('data fitted')

%% Save data with subject number in file

switch s
    case '99'
        n1_99 = n1_data;
        save('m-files/detectN1/data/n1_99', 'n1_99');
        clear n1_99
    case '31'
        n1_31 = n1_data;
        save('m-files/detectN1/data/n1_31', 'n1_31');
        clear n1_31
    case '78'
        n1_78 = n1_data;
        save('m-files/detectN1/data/n1_78', 'n1_78');
        clear n1_78
    case '81'
        n1_81 = n1_data;
        save('m-files/detectN1/data/n1_81', 'n1_81');
        clear n1_81
    case '88'
        n1_88 = n1_data;
        save('m-files/detectN1/data/n1_88', 'n1_88');
        clear n1_88
    case '50'
        n1_50 = n1_data;
        save('m-files/detectN1/data/n1_50', 'n1_50');
        clear n1_50
end
disp(strcat(['Fitting data for subject ' s ' was saved.']))

%% Check goodness of fit
ccep_goodnessOfFit;

%% Fitting and goodness of fit are done, now plot some stuff.
%% Plot fitting result for one measured electrode.
%Plots all trials, the mean and the fitted line (calculated as the mean of all trial fits)
%Choose electrodes:
% 
% el = 61; % stimulated
% elm = 64; % measured
% 
% % time-window to plot:
% t_min = 5;
% t_max = 300;
% 
% % N1 time-window
% t_n1 = t>ccep_parmsN1fit.t_low_N1 & t<ccep_parmsN1fit.t_up_N1;
% 
% % % calculate fitted line with the mean parameters across epochs
% n1_elm = squeeze(mean(n1_data(el).fit(elm,:,:),2));
% fitted_line_n1 = n1_elm(1)*sqrt(2*pi)*normpdf(ccep_parmsN1fit.t2fitN1,n1_elm(2),n1_elm(3)) + n1_elm(4);
% 
% % calculate the line for every epoch and then take the mean
% n1_elm = squeeze(n1_data(el).fit(elm,:,:));
% fitted_line_n1_epoch = zeros(size(fitted_line_n1,2),size(n1_elm,1)); % time X epochs
% for k=1:size(n1_elm,1) % epochs
%     n1_elm_curr = n1_elm(k,:);
%     fitted_line_n1_epoch(:,k) = n1_elm_curr(1)*sqrt(2*pi)*...
%         normpdf(ccep_parmsN1fit.t2fitN1,n1_elm_curr(2),n1_elm_curr(3)) + n1_elm_curr(4);
% end
% clear n1_elm_curr
% fitted_line_n1_mean = mean(fitted_line_n1_epoch,2);
% 
% % plot
% figure,hold on
% 
% t_plot = t(t>t_min & t<t_max); % time to plot
% 
% plot(t_plot,zeros(size(t_plot)),'k:') % plot zero line
% 
% data_plot = squeeze(data(elm,t>t_min & t<t_max,ccep(el).epochs)); % plot all epochs
% plot(t_plot,data_plot,'Color',[.9 .9 .9])
% 
% data_plot = squeeze(mean(data(elm,t>t_min & t<t_max,ccep(el).epochs),3)); % plot mean of all epochs
% plot(t_plot,data_plot,'k')
% 
% % plot(t(t_n1),fitted_line_n1,'r') % plot fitted line % mean parameters
% % plot(t(t_n1),fitted_line_n1_mean,'m') % plot fitted line % mean line
% plot(t(t_n1),fitted_line_n1_epoch,'m')
% 
% 
% title({['el ' int2str(el) ' to ' int2str(elm)]...
%     'mean (of all trials) fit'})

%% Render significant CCEPs on the cortex for all electrodes of interest
% % Run goodness of fit script first. It classifies good fits and bad fits
% % based on relative error (should be smaller than 2.5), spearman's and its
% % p-value
% 
% % Get alpha value as a string (as .05):
% alphaVal = 0.05;
% dispAlpha = num2str(alphaVal);
% dispAlpha = strsplit(dispAlpha, '.');
% dispAlpha = dispAlpha{2};
% dispAlpha = strcat(['.' dispAlpha]);
% 
% figure('Position',[0 0 600 600]);
% set(gcf,'PaperPositionMode', 'auto')
% 
% for el = 51%els2check
%     if numel(ccep(el).epochs) == 0 % skip this electrode, if it wasn't stimulated.
%         disp(['el ' int2str(el) ' not stimulated'])
%     else
% 
%         % t-test p-value
%         n1_p = n1_data(el).ttest(:,2);
%         % amplitude
%         n1_amp = zeros(size(data,1),1);
%         n1_lat = zeros(size(data,1),1);
%         n1_amp(els2check,1) = median(n1_data(el).fit(els2check,:,5),2);
%         % latency
%         n1_lat(els2check,1) = mean(n1_data(el).fit(els2check,:,2),2);
%         % alpha Bonferroni corrected
%         n1_amp(n1_p>=alphaVal) = 0; % if its not significant, we're setting it to zero
%         % set as 0 if we're not fitting a real n1
%         n1_fit = goodnessOfFitClassified(el,:);
%         n1_amp((n1_fit~=1),1) = NaN;
%         % render
%         ctmr_gauss_plot(cortex,[0 0 0],0)
%         el_add(elecmatrix((1:size(data,1)),:),'b',5)
%         el_add(elecmatrix(ccep(el).els,:),[.9 .9 .9],20)
%         el_add_sizecolor(elecmatrix(1:size(data,1),:),n1_amp,n1_lat,200,40)
%         
%         title({['significant CCEPs when ' int2str(el) ' was stimulated']...
%             ['alpha = ' dispAlpha ' uncorrected']})
% %                 print(gcf,'-r300', '-dpng', strcat(['figures/' s '/cortexCCEPs/cortex_cceps_el' int2str(el) '_s' s]))
% %         clf(gcf,'reset');
% 
%     end
% end
% %
% % close all
% 
% % break %#ok Stop the script before it runs last section (plots a lot of figures)

%% Save fitting plots for CCEPs
% % This plots and saves A LOT of figures! All combinations of
% % stimulated-measured electrode (of interest) pairs!
% 
% t_min = 5;
% t_max = 300;
% 
% t_n1 = t>ccep_parmsN1fit.t_low_N1 & t<ccep_parmsN1fit.t_up_N1;
% 
% figure('Position',[0 0 500 500])
% set(gcf,'PaperPositionMode', 'auto')
% % cd(['figures/' s '/CCEPs']);
% for el = els2check
%     
%     if numel(ccep(el).epochs) == 0 % skip this electrode, if it wasn't stimulated.
%         
%         disp(['el ' int2str(el) ' not stimulated'])
%     else
%         %%%% This bit gets the saved files into the correct directories %%%
%         % Make sure we're in the right directory
%         
%         % check if there already is a directory for this electrode
%         dirCheck = int2str(el);
%         a = exist('dirCheck', 'dir');
%         if a == 1 % if a directory already exist, go to that directory
%             cd(int2str(el));
%         else
%             mkdir(int2str(el)); % if there's no directory for this el, make one and go there
%             cd(int2str(el));
%         end        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         for elm = els2check
%             if find(ccep(el).els==elm)
%                 continue
%             else
%                 hold on
%                 n1_elm = median(n1_data(el).fit(elm,:,:),2);
%                 fitted_line_n1 = n1_elm(1)*sqrt(2*pi)*normpdf(ccep_parmsN1fit.t2fitN1,n1_elm(2),n1_elm(3)) + n1_elm(4);
%                 
%                 t_plot = t(t>t_min & t<t_max);
%                 plot(t_plot,zeros(size(t_plot)),'k:')
%                 
%                 data_plot = squeeze(data(elm,t>t_min & t<t_max,ccep(el).epochs));
%                 plot(t_plot,data_plot,'Color',[.9 .9 .9])
%                 
%                 data_plot = squeeze(mean(data(elm,t>t_min & t<t_max,ccep(el).epochs),3));
%                 plot(t_plot,data_plot,'k')
%                 
%                 plot(t(t_n1),fitted_line_n1,'r')
%                 
%                 title(['el ' int2str(el) ' to ' int2str(elm)])
%                 
%                 print(gcf,'-r300', '-dpng', strcat(['ccep_fit_el' int2str(el) '_to_' int2str(elm) '_s' s]))
%                 hold off
%                 clf(gcf,'reset');
%             end
%         end
%         cd('../')
%     end
%     
% end

close all
disp('all done!')
