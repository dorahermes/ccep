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
                        % otherwise give the amplitude the value NaN
                    elseif isempty(temp_n1_peaks_samp)
                        n1_peak_amplitude = NaN;
                    end

                    % if N1 exceeds positive threshold, it is deleted
                    if temp_n1_peaks_ampl > 0
                       n1_peak_sample = NaN;
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
    
        class_1_visual = reshape(validation_matrix(:,:,1),1,size(validation_matrix(:,:,1),1)*size(validation_matrix(:,:,1),2));
        class_2_code = reshape(temp_mat(1:81,1:37),1, ...
           size(validation_matrix(:,:,1),1)*size(validation_matrix(:,:,1),2));
        
        
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

% write parameters_optimalization.mat to folder
working_dir = fullfile('/Fridge','users','jaap','ccep','dataBIDS');
if ~exist(fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_run-' run_label '_parameters_optimalization.mat']))
    disp('writing output parameters_optimalize_mat.mat')
    save([fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_run-' run_label '_parameters_optimalization.mat'])],...
    'parameters_optimalize_mat_0458')
else
    disp(['ERROR: can not overwrite, output file already exists '])
end

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


%% Combine the validation scores of multiple datasets + ROCcurve

% The validation_matrix and parameters_optimalize_mat of all validated data
% are loaded (the variables in the mat-files are renamed so they do not overwrite)
load(fullfile(working_dir,'validation', 'sub-RESP0458_ses-1_run-011714_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0458_ses-1_run-011714_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0468_ses-1_run-031729_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0468_ses-1_run-031729_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0621_ses-1_run-021147_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0621_ses-1_run-021147_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0706_ses-1_run-041501_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0706_ses-1_run-041501_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0733_ses-1b_run-050941_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0733_ses-1b_run-050941_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0768_ses-1_run-021704_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0768_ses-1_run-021704_validation_matrix.mat'))



% The scores of the different datasets are based on a different amount of averaged epochs
% To calculate the average scores, it first needs recalculation.

% Calculate total amount of epochs over all datasets
total_validated_epochs = (size(validation_matrix_0458,1)*size(validation_matrix_0458,2)) ...
    + (size(validation_matrix_0468,1)*size(validation_matrix_0468,2)) ...
    + (size(validation_matrix_0621,1)*size(validation_matrix_0621,2)) ...
    + (size(validation_matrix_0706,1)*size(validation_matrix_0706,2)) ...
    + (size(validation_matrix_0733,1)*size(validation_matrix_0733,2)) ... 
    + (size(validation_matrix_0768,1)*size(validation_matrix_0768,2));

% multiply each parameters_optimalize_mat by the fraction on total epochs
parameters_optimalize_mat_0458 = parameters_optimalize_mat_0458 * ...
    (size(validation_matrix_0458,1)*size(validation_matrix_0458,2)/total_validated_epochs);
parameters_optimalize_mat_0468 = parameters_optimalize_mat_0468 * ...
    (size(validation_matrix_0468,1)*size(validation_matrix_0468,2)/total_validated_epochs);
parameters_optimalize_mat_0621 = parameters_optimalize_mat_0621 * ...
    (size(validation_matrix_0621,1)*size(validation_matrix_0621,2)/total_validated_epochs);
parameters_optimalize_mat_0706 = parameters_optimalize_mat_0706 * ...
    (size(validation_matrix_0706,1)*size(validation_matrix_0706,2)/total_validated_epochs);
parameters_optimalize_mat_0733 = parameters_optimalize_mat_0733 * ...
    (size(validation_matrix_0733,1)*size(validation_matrix_0733,2)/total_validated_epochs); 
parameters_optimalize_mat_0768 = parameters_optimalize_mat_0768 * ...
    (size(validation_matrix_0768,1)*size(validation_matrix_0768,2)/total_validated_epochs); 


% NOTE TO SELF: add extra option in which there is no correction for the
% amount of epochs per patient

% add scores to get 1 averaged matrix with parameter scores
averaged_parameter_scores = parameters_optimalize_mat_0458 + parameters_optimalize_mat_0468 + parameters_optimalize_mat_0621 ...
    + parameters_optimalize_mat_0706 + parameters_optimalize_mat_0733 + parameters_optimalize_mat_0768;

figure
subplot(1,4,1:3)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(averaged_parameter_scores,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(averaged_parameter_scores,2)

    sens_plot = averaged_parameter_scores(:,time_th,1);
    spes_plot = averaged_parameter_scores(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1-specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
set(gca,'XTick',[0:.05:1],'YTick',[0:.05:1],'FontName','Arial','FontSize',10) % get(gca) to see properties to change

% plot legenda for all ROCs
time_end = tt(n1_samples_end);

subplot(1,4,4),hold on
for time_th = 1:size(averaged_parameter_scores,2)
    plot(0,time_end(time_th),'.','MarkerSize',20,'Color',my_colors(time_th,:))
end
ylabel('time end (s)')

working_dir = fullfile('/Fridge','users','jaap','ccep','dataBIDS');
if ~exist(fullfile(working_dir,'validation', 'averaged_parameters_optimalization.mat'))
    disp('writing output averaged_parameters_optimalization.mat')
    save(fullfile(working_dir,'validation', 'averaged_parameters_optimalization.mat'),...
    'averaged_parameter_scores')
else
    disp(['ERROR: can not overwrite, output file already exists '])
end


%% ROC plots with different parameters - different ROCs for amplitude

figure('Position',[0 0 450 450])
subplot(1,4,1:3)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(averaged_parameter_scores,1));

% for all different ranges plot a ROC- curve
for ampl_th = 1:size(averaged_parameter_scores,1)

    sens_plot = averaged_parameter_scores(ampl_th,:,1);
    spes_plot = averaged_parameter_scores(ampl_th,:,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(ampl_th,:))
end
axis square
xlabel('1-specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
set(gca,'XTick',[0:.05:1],'YTick',[0:.05:1],'FontName','Arial','FontSize',10) % get(gca) to see properties to change

% plot legenda for all ROCs
thresh = [1:.2:5]*50; 

subplot(1,4,4),hold on
for ampl_th = 1:size(averaged_parameter_scores,1)
    plot(0,thresh(ampl_th),'.','MarkerSize',20,'Color',my_colors(ampl_th,:))
end
ylabel('significant amplitude threshold (uV)')

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/Fridge/users/jaap/temp2')


%% looking at total ROC without RESP0768

% The validation_matrix and parameters_optimalize_mat of all validated data
% are loaded (the variables in the mat-files are renamed so they do not overwrite)
load(fullfile(working_dir,'validation', 'sub-RESP0621_ses-1_run-021147_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0621_ses-1_run-021147_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0706_ses-1_run-041501_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0706_ses-1_run-041501_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0733_ses-1b_run-050941_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0733_ses-1b_run-050941_validation_matrix.mat'))

% The scores of the different datasets are based on a different amount of averaged epochs
% To calculate the average scores, it first needs recalculation.

% Calculate total amount of epochs over all datasets
total_validated_epochs_corrected = (size(validation_matrix_0621,1)*size(validation_matrix_0621,2)) ...
    + (size(validation_matrix_0706,1)*size(validation_matrix_0706,2)) ...
    + (size(validation_matrix_0733,1)*size(validation_matrix_0733,2));

% multiply each parameters_optimalize_mat by the fraction on total epochs
parameters_optimalize_mat_0621 = parameters_optimalize_mat_0621 * ...
    (size(validation_matrix_0621,1)*size(validation_matrix_0621,2)/total_validated_epochs_corrected);
parameters_optimalize_mat_0706 = parameters_optimalize_mat_0706 * ...
    (size(validation_matrix_0706,1)*size(validation_matrix_0706,2)/total_validated_epochs_corrected);
parameters_optimalize_mat_0733 = parameters_optimalize_mat_0733 * ...
    (size(validation_matrix_0733,1)*size(validation_matrix_0733,2)/total_validated_epochs_corrected); 

% add scores to get 1 averaged matrix with parameter scores
averaged_parameter_scores_corrected = parameters_optimalize_mat_0621 + parameters_optimalize_mat_0706 ...
    + parameters_optimalize_mat_0733;

figure
subplot(1,4,1:3)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(averaged_parameter_scores_corrected,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(averaged_parameter_scores_corrected,2)

    sens_plot = averaged_parameter_scores_corrected(:,time_th,1);
    spes_plot = averaged_parameter_scores_corrected(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1-specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
set(gca,'XTick',[0:.05:1],'YTick',[0:.05:1],'FontName','Arial','FontSize',10) % get(gca) to see properties to change

% plot legenda for all ROCs
time_end = tt(n1_samples_end);

subplot(1,4,4),hold on
for time_th = 1:size(averaged_parameter_scores_corrected,2)
    plot(0,time_end(time_th),'.','MarkerSize',20,'Color',my_colors(time_th,:))
end
ylabel('time end (s)')

working_dir = fullfile('/Fridge','users','jaap','ccep','dataBIDS');
if ~exist(fullfile(working_dir,'validation', 'averaged_parameters_optimalization_corrected.mat'))
    disp('writing output averaged_parameters_optimalization_corrected.mat')
    save(fullfile(working_dir,'validation', 'averaged_parameters_optimalization_corrected.mat'),...
    'averaged_parameter_scores_corrected')
else
    disp(['ERROR: can not overwrite, output file already exists '])
end

%% Finding optimal parameters


% change all not-optimal endpoints to NaN (all but 70 and 80ms)
averaged_parameter_scores_corrected(:,1:3,:) = NaN;
averaged_parameter_scores_corrected(:,6:8,:) = NaN;


% Set necessary specificity on 95%
test = averaged_parameter_scores_corrected(:,:,2) >= 0.95;

% loop through results and change 
for opt_len = 1:size(averaged_parameter_scores_corrected(:,:,2),1)
    for opt_wid = 1:size(averaged_parameter_scores_corrected(:,:,2),2)
        if test(opt_len,opt_wid) == 0 
            averaged_parameter_scores_corrected(opt_len,opt_wid,:) = NaN;
        end
    end
end

% find optimal sensitivity under condition of specificitity >= 95%
optimal_sensi = max(max(averaged_parameter_scores_corrected(:,:,1)));

% find corresponding parameters location in matrix
optimum_mat = find(averaged_parameter_scores_corrected(:,:,1)==optimal_sensi);


% convert this position to the parameters
optimum_mat_x = rem(optimum_mat,size(averaged_parameter_scores_corrected(:,:,1),1));
optimum_mat_y = floor(optimum_mat/size(averaged_parameter_scores_corrected(:,:,1),1))+1;
optimum_threshold = thresh(optimum_mat_x);
optimum_endpoint = tt(n1_samples_end(optimum_mat_y));
