% This script is used to find the best threshold and range for the detectER
% script. It does this by running through multiple values of the parameters
% and compares their sensitivity, specificity, ROC-curve and total score.

%% parameters to investigate 

% thresh default is 2.5, which gives default treshhold of 125 microV
% now set to between 50 and 250 with steps of 10 microV
thresh = [1:.2:5];   
% default is 50, kept the same because already taken into account by thresh
minSD = 50;   

% N1 peak default is 10 and 50 ms. Also try 40, 50, 60, 70 and 80ms as end

n1_samples_start = find(tt>0.01,1);
n1_samples_end = [find(tt>0.04,1), find(tt>0.05,1), find(tt>0.06,1), find(tt>0.07,1), find(tt>0.08,1), find(tt>0.09,1), find(tt>0.1,1), find(tt>0.11,1)]; % can be higher for RESP0768


%%
% initialize parameters matrix to save results of different parameters
parameters_optimalize_mat = NaN(length(thresh),length(n1_samples_end),3);

%% Loop with the given parameters through the detectER loop 
tic;

% for all different thresholds
for th = 1:length(thresh)
    % for all different N1 peak ranges
    for ss = 1:length(n1_samples_end)
        
        % run through the detect ER script, but only part that finds N1
        
        % initialize temporarily matrix to store sample of peak
        temp_mat = NaN(size(cc_epoch_sorted_avg,1),size(cc_epoch_sorted_avg,2));
        
        % for every averaged stimulation
        for jj = 1:size(cc_epoch_sorted_avg,2)
            % for every channel
            for ii = 1:size(cc_epoch_sorted_avg,1)

                % baseline subtraction: take median of part of the averaged signal for
                % this stimulation pair before stimulation, which is the half of the
                % epoch
                baseline_tt = tt>-2 & tt<-.1;
                signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt),3);

                % subtract median baseline from signal
                new_signal = squeeze(cc_epoch_sorted_avg(ii,jj,:)) - signal_median;
                % testplot new signal: plot(tt,squeeze(new_signal))

                % take area before the stimulation of the new signal and calculate its SD
                pre_stim_sd = std(new_signal(baseline_tt));

                % if the pre_stim_sd is smaller that the minimally needed SD, use this
                % the minSD as pre_stim_sd
                if pre_stim_sd < minSD
                    pre_stim_sd = minSD;
                end

                % when the electrode is stimulated
                if ii == cc_stimsets(jj,1) || ii == cc_stimsets(jj,2)
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;

                    % in other electrode
                else
                    % look in stimsets which electrodes are
                    % Next: only do analysis on the channels that are not stimulated
                    % Should we e.g. channel F01 and F02 in stimulation of
                    % F01-F02 load the data in as NaN?
                    % so if member, do nothing or assign NaN
                    % elseif: run rest of the loop and find peak

                    % use peakfinder to find all positive and negative peaks and their
                    % amplitude.
                    % tt are the samples of the epoch based on the Fs and -2.5 to 2.5
                    % seconds sample of the total epoch
                    % As tt use first sample after timepoint 0  (+ extra samples to not take artifact)
                    % till first sample after 0,5 seconds (rougly 1000 samples)
                    [all_samppos, all_amplpos] = ccep_peakfinder(new_signal(find(tt>0,1)+extrasamps:find(tt>0.5,1)),sel,[],1);
                    [all_sampneg, all_amplneg] = ccep_peakfinder(new_signal(find(tt>0,1)+extrasamps:find(tt>0.5,1)),sel,[],-1);

                    % If the first selected sample is a peak, this is not a real peak,
                    % so delete
                    all_amplpos(all_samppos==1) = [];
                    all_samppos(all_samppos==1) = [];
                    all_amplneg(all_sampneg==1) = [];
                    all_sampneg(all_sampneg==1) = [];

                    % convert back timepoints based on tt, substract 1 because before
                    % the first sample after stimulation is taken
                    all_samppos = all_samppos + find(tt>0,1) + extrasamps -1;
                    all_sampneg = all_sampneg + find(tt>0,1) + extrasamps -1;


                    % first do analysis on the N1, only if this is found, continue with
                    % other properties

                    % for N1, first select the range in which the N1 could appear, and
                    % select the peaks found in this range
                    temp_n1_peaks_samp = all_sampneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end(ss)));
                    temp_n1_peaks_ampl = all_amplneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end(ss)));

                    % if peak(s) found, select biggest peak
                    if ~isempty(temp_n1_peaks_samp)
                        max_n1_ampl = find(abs(temp_n1_peaks_ampl) == max(abs(temp_n1_peaks_ampl)));
                        n1_peak_sample = temp_n1_peaks_samp(max_n1_ampl(1));
                        n1_peak_amplitude = temp_n1_peaks_ampl(max_n1_ampl(1));
                        % otherwise give the amplitude the value 0
                    elseif isempty(temp_n1_peaks_samp)
                        n1_peak_amplitude = NaN;
                    end

                    % when peak amplitude is saturated, it is deleted
                    if abs(n1_peak_amplitude) > 3000
                        n1_peak_sample = NaN;
                        n1_peak_amplitude = NaN;
                    end

                    % if the peak is not big enough to consider as a peak, assign NaN
                    if abs(n1_peak_amplitude) < thresh(th)* abs(pre_stim_sd)
                        n1_peak_sample = NaN;
                        n1_peak_amplitude = NaN;
                    end
                end
                   
                % write n1_peak_sample to temporarily matrix
                temp_mat(ii,jj) = n1_peak_sample;
                
                if ~isnan(temp_mat(ii,jj))
                temp_mat(ii,jj) = 1;
                elseif isnan(temp_mat(ii,jj))
                temp_mat(ii,jj) = 0;
                end 
                
            end
        end
    
        class_1_visual = reshape(validation_matrix(:,:,1),1,size(cc_epoch_sorted_avg,1)*size(cc_epoch_sorted_avg,2));
        class_2_code = reshape(temp_mat(:,:),1,size(cc_epoch_sorted_avg,1)*size(cc_epoch_sorted_avg,2));
        
        
        TP = length(find(class_1_visual == 1 & class_2_code == 1));
        TN = length(find(class_1_visual == 0 & class_2_code == 0));
        FP = length(find(class_1_visual == 0 & class_2_code == 1));
        FN = length(find(class_1_visual == 1 & class_2_code == 0));

        sensitivity = TP/ (TP + FN);
        specificity = TN / (TN + FP);
        
        [c,cm,ind,per] = confusion(class_1_visual,class_2_code);

        parameters_optimalize_mat(th,ss,1) = sensitivity;
        parameters_optimalize_mat(th,ss,2) = specificity;
        parameters_optimalize_mat(th,ss,3) = 1-c;

    end
end

toc;

save([fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
['sub-' sub_label '_ses-' ses_label '_run-' run_label '_parameters_optimalization.mat'])],...
'parameters_optimalize_mat')

%% ROC plots with different parameters - different ROCs for time

figure
subplot(1,4,1:3)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(parameters_optimalize_mat,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(parameters_optimalize_mat,2)

    sens_plot = parameters_optimalize_mat(:,time_th,1);
    spes_plot = parameters_optimalize_mat(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1-specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])

% plot legenda for all ROCs
time_end = tt(n1_samples_end);

subplot(1,4,4),hold on
for time_th = 1:size(parameters_optimalize_mat,2)
    plot(0,time_end(time_th),'.','MarkerSize',20,'Color',my_colors(time_th,:))
end
ylabel('time end (s)')

%% ROC plots with different parameters - different ROCs for amplitude

figure
subplot(1,4,1:3)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(parameters_optimalize_mat,1));

% for all different ranges plot a ROC- curve
for ampl_th = 1:size(parameters_optimalize_mat,1)

    sens_plot = parameters_optimalize_mat(ampl_th,:,1);
    spes_plot = parameters_optimalize_mat(ampl_th,:,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(ampl_th,:))
end
xlabel('1-specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])

% plot legenda for all ROCs
thresh = [1:.2:5]*50; 

subplot(1,4,4),hold on
for ampl_th = 1:size(parameters_optimalize_mat,1)
    plot(0,thresh(ampl_th),'.','MarkerSize',20,'Color',my_colors(ampl_th,:))
end
ylabel('significant amplitude threshold (uV)')