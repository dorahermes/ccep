% this script is used to validate the ccep_detectER.m 
% it plots the epochs and standard deviations.  By visual assessment is decided
% whether there is an significant N1 peak.
% This assesment is then compared to what the code found to validate the code

% load master_preprocessCCEP script for the participant first,
% to get all necessary variables


% By J. van der Aar & D. Hermes, UMC Utrech 05-2019
% inspired by code and paper of D. van Blooijs (2018)
%% load in peakfinder
addpath(genpath('/Fridge/users/jaap/github/ccep/functions')) 



%% set defined parameters
% pre-allocation: variables determined in Dorien's paper
thresh = 2.5;   %
minSD = 50;     % in microV: minimum standard deviation of prestimulus baseline in case actual std of prestim baseline is too small
sel = 20;       % how many samples around peak not considered as another peak
extrasamps = 5; % set to 5 now, to read properties of N1 onset
SDfactor=[];
SD = [];

%% set timeframes

% set timeframes in which you expect the parameter to show up

% P1/N1 onset: between 2 and 20 ms
p1_samples_start = find(tt>0.002,1);
p1_samples_end = find(tt>0.02,1);

% N1 peak: 10 and 50 ms
n1_samples_start = find(tt>0.01,1);
n1_samples_end = find(tt>0.05,1);

% P2/N2 onset: 50 and 150 ms
p2_samples_start = find(tt>0.05,1);
p2_samples_end = find(tt>0.15,1);

% N2 peak: between 50 and 300 ms
n2_samples_start =  find(tt>0.05,1);
n2_samples_end = find(tt>0.3,1);

% End N2: ... % possible to add this as well

%% Loop and plot

% output in channels X stimulations X [p1 sample, p1 amplitude , n1 sample
% n1 amplitude, p2 sample, p2 amplitude, n2 sample, n2 amplitude]
output_ER_all = NaN(size(cc_epoch_sorted_avg,1),size(cc_epoch_sorted_avg,2),8);

validation_matrix = NaN(size(cc_epoch_sorted_avg,1),size(cc_epoch_sorted_avg,2),2);

% for every channel
for ii = 1:size(cc_epoch_sorted_avg,1)
    % for every averaged stimulation
    for jj = 1:size(cc_epoch_sorted_avg,2)
        
        % baseline subtraction: take median of part of the averaged signal for
        % this stimulation pair before stimulation, which is the half of the
        % epoch
        baseline_tt = tt>-1 & tt<-.01;
        signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt),3);
        
        % subtract median baseline from signal
        new_signal = squeeze(cc_epoch_sorted_avg(ii,jj,:)) - signal_median;
        % testplot new signal: plot(tt,squeeze(new_signal))
        
        % take area before the stimulation of the new signal and calculate its SD
        pre_stim_sd_orig = std(new_signal(baseline_tt));
        
        % if the pre_stim_sd is smaller that the minimally needed SD, use this
        % the minSD as pre_stim_sd
        if pre_stim_sd_orig < minSD
            pre_stim_sd = minSD;
        else 
            pre_stim_sd =  pre_stim_sd_orig;
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
            temp_n1_peaks_samp = all_sampneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
            temp_n1_peaks_ampl = all_amplneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
            
            % if peak(s) found, select biggest peak
            if ~isempty(temp_n1_peaks_samp)
                max_n1_ampl = find(abs(temp_n1_peaks_ampl) == max(abs(temp_n1_peaks_ampl)));
                n1_peak_sample = temp_n1_peaks_samp(max_n1_ampl(1));
                n1_peak_amplitude = temp_n1_peaks_ampl(max_n1_ampl(1));
                % otherwise give the amplitude the value 0
            elseif isempty(temp_n1_peaks_samp)
                n1_peak_amplitude = 0;
            end
            
            % when peak amplitude is saturated, it is deleted
            if abs(n1_peak_amplitude) > 3000
                n1_peak_sample = NaN;
                n1_peak_amplitude = NaN;
            end
            
            % if the peak is not big enough to consider as a peak, assign NaN
            if abs(n1_peak_amplitude) < thresh* abs(pre_stim_sd)
                n1_peak_sample = NaN;
                n1_peak_amplitude = NaN;
            end
        end
        % If there is no significant n1 peak, assign NaN to all other
        % properties. If there is a significant peak, find other properties
        if isnan(n1_peak_sample)
            p1_peak_sample = NaN;
            p1_peak_amplitude = NaN;
            n2_peak_sample = NaN;
            n2_peak_amplitude = NaN;
            p2_peak_sample = NaN;
            p2_peak_amplitude = NaN;
            % If there is a significant peak, find other properties
        elseif isnumeric(n1_peak_sample)
            
            % for P1, first select the range in which the P1 could appear, and
            % select the peaks found in this range
            temp_p1_peaks_samp = all_samppos((p1_samples_start <= all_samppos) & (all_samppos <= p1_samples_end));
            temp_p1_peaks_ampl = all_amplpos((p1_samples_start <= all_samppos) & (all_samppos <= p1_samples_end));
            
            % if peak(s) found, select biggest peak
            if ~isempty(temp_p1_peaks_samp)
                max_p1_ampl = find(abs(temp_p1_peaks_ampl) == max(abs(temp_p1_peaks_ampl)));
                p1_peak_sample = temp_p1_peaks_samp(max_p1_ampl(1));
                p1_peak_amplitude = temp_p1_peaks_ampl(max_p1_ampl(1));
                % otherwise give the sample and amplitude NaN
            elseif isempty(temp_p1_peaks_samp)
                p1_peak_sample = NaN;
                p1_peak_amplitude = NaN;
            end
            
            % for N2, first select the range in which the N2 could appear, and
            % select the peaks found in this range
            temp_n2_peaks_samp = all_sampneg((n2_samples_start <= all_sampneg) & (all_sampneg <= n2_samples_end));
            temp_n2_peaks_ampl = all_amplneg((n2_samples_start <= all_sampneg) & (all_sampneg <= n2_samples_end));
            
            % if peak(s) found, select biggest peak
            if ~isempty(temp_n2_peaks_samp)
                max_n2_ampl = find(abs(temp_n2_peaks_ampl) == max(abs(temp_n2_peaks_ampl)));
                n2_peak_sample = temp_n2_peaks_samp(max_n2_ampl(1));
                n2_peak_amplitude = temp_n2_peaks_ampl(max_n2_ampl(1));
                % otherwise give the sample and amplitude NaN
            elseif isempty(temp_n2_peaks_samp)
                n2_peak_sample = NaN;
                n2_peak_amplitude = NaN;
            end
            
            % for P2, first select the range in which the P2 could appear, and
            % select the peaks found in this range
            temp_p2_peaks_samp = all_samppos((p2_samples_start <= all_samppos) & (all_samppos <= p2_samples_end));
            temp_p2_peaks_ampl = all_amplpos((p2_samples_start <= all_samppos) & (all_samppos <= p2_samples_end));
            
            % if peak(s) found, select biggest peak
            if ~isempty(temp_p2_peaks_samp)
                max_p2_ampl = find(abs(temp_p2_peaks_ampl) == max(abs(temp_p2_peaks_ampl)));
                p2_peak_sample = temp_p2_peaks_samp(max_p2_ampl(1));
                p2_peak_amplitude = temp_p2_peaks_ampl(max_p2_ampl(1));
                % otherwise give the sample and amplitude NaN
            elseif isempty(temp_p2_peaks_samp)
                p2_peak_sample = NaN;
                p2_peak_amplitude = NaN;
            end
            
        end
        
        % add properties to output frame - careful -
        % [p1 sample, p1 amplitude , n1 sample n1 amplitude, p2 sample, p2 amplitude, n2 sample, n2 amplitude]
        output_ER_all(ii,jj,1) = p1_peak_sample;
        output_ER_all(ii,jj,2) = p1_peak_amplitude;
        output_ER_all(ii,jj,3) = n1_peak_sample;
        output_ER_all(ii,jj,4) = n1_peak_amplitude;
        output_ER_all(ii,jj,5) = p2_peak_sample;
        output_ER_all(ii,jj,6) = p2_peak_amplitude;
        output_ER_all(ii,jj,7) = n2_peak_sample;
        output_ER_all(ii,jj,8) = n2_peak_amplitude;
        
        
        
        % plotting and validation part 
        
        % display figure in specific size
        figh = figure(1);
        pos = get(figh,'position');
        set(figh,'position',[pos(1:2)/4 pos(3:4)*3])
        % plot relevant part of the signal (1s before till 500ms after stimulation)
        subplot(1,2,1)
        plot(tt(tt>-1 & tt<.5),squeeze(new_signal(tt>-1 & tt<.5)))
        
        xlabel('time(s)')
        ylabel('amplitude(uV)')
        title(['elec ' data_hdr.label{ii} ' for stimulation of ' data_hdr.label{cc_stimsets(jj,1)} ' and ' data_hdr.label{cc_stimsets(jj,2)} ])
        ylim([-1000 1000])
        
        hold on
        % plot(tt(output_ER_all(ii,jj,1)),output_ER_all(ii,jj,2),'b*')
        % plot(tt(output_ER_all(ii,jj,3)),output_ER_all(ii,jj,4),'r*')
        % plot(tt(output_ER_all(ii,jj,5)),output_ER_all(ii,jj,6),'g*')
        % plot(tt(output_ER_all(ii,jj,7)),output_ER_all(ii,jj,8),'y*')
        
        % plot calculated baseline standard deviation
        plot(tt(baseline_tt), pre_stim_sd_orig+zeros(size(tt(baseline_tt))), 'r-')
        plot(tt(baseline_tt), -pre_stim_sd_orig+zeros(size(tt(baseline_tt))), 'r-')
        
        % plot adjusted baseline (when calculated < minSD)
        plot(tt(baseline_tt), pre_stim_sd+zeros(size(tt(baseline_tt))), 'g-')
        plot(tt(baseline_tt), -pre_stim_sd+zeros(size(tt(baseline_tt))), 'g-')
        hold off
        
        % plot part with n1 peak to inspect
        subplot(1,2,2)
        plot(tt(tt>0.01 & tt<.1),squeeze(new_signal(tt>0.01 & tt<.1)))
        title(['elec ' data_hdr.label{ii} ' for stimulation of ' data_hdr.label{cc_stimsets(jj,1)} ' and ' data_hdr.label{cc_stimsets(jj,2)} ])
        xlabel('time(s)')
        ylabel('amplitude(uV)')
        ylim([-1000 1000])
        hold on
        % plot SD * threshhold to see if it exceeds and is significant
        plot(tt(tt>0 & tt<.1), thresh*pre_stim_sd+zeros(size(tt(tt>0 & tt<.1))), 'g-')
        plot(tt(tt>0 & tt<.1), thresh*-pre_stim_sd+zeros(size(tt(tt>0 & tt<.1))), 'g-')
        hold off
        
        % ask for input, do you see a significant n1 peak?
        visually_found_n1_peak = input('n1 peak? [y/n] ','s');
        close(figh)
        
        % save in matrix
        % !!!!!!!!!!         110 = 'N', 121 = 'Y'          !!!!!!!!!!
        validation_matrix(ii,jj,1) = visually_found_n1_peak;
        % change NaNs in output_ER_all to 'n' and # to 'y'
        % and add these to matrix
        if ~isnan(output_ER_all(ii,jj,3))
            validation_matrix(ii,jj,2) = 'y';
        elseif isnan(output_ER_all(ii,jj,3))
            validation_matrix(ii,jj,2) = 'n';
        end
        
    end 
    
    % all stimulations of a channel, save data, so continuing later is
    % possible
    working_dir = fullfile('/Fridge','users','jaap','ccep','dataBIDS');
    save([fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_run-' run_label '_validation_matrix.mat'])],...
    'validation_matrix')
    
    
end

%% statistics & ROC
% write code what finds:
% true positives
% true negatives
% false positives
% false negatives