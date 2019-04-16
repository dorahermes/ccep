%% load in peakfinder
addpath(genpath('/Fridge/users/jaap/github/ccep/functions')) 

%% set defined parameters
% pre-allocation: variables determined in Dorien's thesis
thresh = 2.5;   %
minSD = 50;     % in microV: minimum standard deviation of prestimulus baseline in case actual std of prestim baseline is too small
sel = 20;       % how many samples around peak not considered as another peak
extrasamps = 5; % set to 5 now, to read properties of N1 onset
SDfactor=[];
SD = [];

%% Run through all avg epochs of all channels to find peaks

% this channel + stimulation gives nice N1 to test with
% ii = 3; jj = 1; 

% output in channels X stimulations X [onset(seconds) & amplitude]
output_ER = NaN(size(cc_epoch_sorted_avg,1),size(cc_epoch_sorted_avg,2),2);

% for every channel
for ii = 3%1:size(cc_epoch_sorted_avg,1)
    % for every averaged stimulation
    for jj = 1%:size(cc_epoch_sorted_avg,2)
    
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
        [temp_samppos, temp_amplpos] = ccep_peakfinder(new_signal(find(tt>0,1)+extrasamps:find(tt>0.5,1)),sel,[],1);
        [temp_sampneg, temp_amplneg] = ccep_peakfinder(new_signal(find(tt>0,1)+extrasamps:find(tt>0.5,1)),sel,[],-1);

        % If the first selected sample is a peak, this is not a real peak,
        % so delete
        temp_amplpos(temp_samppos==1) = [];
        temp_samppos(temp_samppos==1) = [];
        temp_amplneg(temp_sampneg==1) = [];
        temp_sampneg(temp_sampneg==1) = [];
        
        % Delete all data after the 0.1s after stimulation, because these cannot be 
        % the N1 peak. Set this to bigger if you also want to find the
        % delayed response
        temp_amplpos(temp_samppos >= find(tt>0.1,1) - (find(tt>0,1) + extrasamps)) = []; 
        temp_samppos(temp_samppos >= find(tt>0.1,1) - (find(tt>0,1) + extrasamps)) = [];
        temp_amplneg(temp_sampneg >= find(tt>0.1,1) - (find(tt>0,1) + extrasamps)) = [];
        temp_sampneg(temp_sampneg >= find(tt>0.1,1) - (find(tt>0,1) + extrasamps)) = [];
        
        %%%% NEXT LINES ARE TO CALCULATE THE N1 PEAK
        %%%% FIRST CALCULATE BOTH NEGATIVE AND POSITIVE PEAKS
        %%%% CONVERT BACK TO THE SAMPLES BASED ON TT
        %%%% THEN LOOK IF IT IS ACTUALLY A POSITIVE OR NEGATIVE PEAK % WHY?
        
        % positive peaks find max and get latency 
        temp_n1peak_samplepos = [];
        temp_n1peak_amplipos = [];
        if ~isempty(temp_samppos)
            % find peak with maximum amplitude
            maxamplpos = find(abs(temp_amplpos) == max(abs(temp_amplpos)));
            % recalculate sample of maxamplpos based on the sampling with
            % tt, necessary because peakfinder adjusted this
            temp_n1peak_samplepos = temp_samppos(maxamplpos(1)) + find(tt>0,1) + extrasamps-1;
            temp_n1peak_amplipos = temp_amplpos(maxamplpos(1));
        end  
        
        % negative peaks find max and get latency
        temp_n1peak_sampleneg = [];
        temp_n1peak_amplineg = [];
        if ~isempty(temp_sampneg)
            % find the maximum amplitude
            maxamplneg = find(abs(temp_amplneg) == max(abs(temp_amplneg)));
            % recalculate sample of maxamplneg based on the sampling with
            % tt, necessary because peakfinder adjusted this
            temp_n1peak_sampleneg = temp_sampneg(maxamplneg(1)) + find(tt>0,1) + extrasamps-1;
            temp_n1peak_amplineg = temp_amplneg(maxamplneg(1));
        end  

        
        % now it is known which datapoints are the peaks, negative and
        % positive peaks, now check which of those is the real n1 peak
        % ------ but why? ------
        
        % if only a positive peak is found, this is the n1
        if ~isempty(temp_n1peak_samplepos) && isempty(temp_n1peak_sampleneg) 
            sample = temp_n1peak_samplepos;
            ampl =  temp_n1peak_amplipos;
        % if only a negative peak is found, this is the n1
        elseif isempty(temp_n1peak_samplepos) && ~isempty(temp_n1peak_sampleneg) 
            sample = temp_n1peak_sampleneg;
            ampl = temp_n1peak_amplineg;
        % if both are found, the one with biggest amplitude is the n1
        elseif ~isempty(temp_n1peak_samplepos) && ~isempty(temp_n1peak_sampleneg) 
            if abs( temp_n1peak_amplipos) > abs(temp_n1peak_amplineg) 
                sample = temp_n1peak_samplepos;
                ampl =  temp_n1peak_amplipos;
            elseif abs(temp_n1peak_amplineg) > abs( temp_n1peak_amplipos) 
                sample = temp_n1peak_sampleneg;
                ampl = temp_n1peak_amplineg;
            % if both amplitudes are as big, take first one in time as 
            elseif abs(temp_n1peak_amplineg) == abs(temp_n1peak_amplipos) 
                if temp_n1peak_sampleneg < temp_n1peak_samplepos 
                    sample = temp_n1peak_sampleneg;
                    ampl = temp_n1peak_amplineg;
                elseif temp_n1peak_sampleneg > temp_n1peak_samplepos
                    sample = temp_n1peak_samplepos;
                    ampl =  temp_n1peak_amplipos;
                end
            end
        % if there are no peaks found 
        elseif isempty(temp_n1peak_samplepos) && isempty(temp_n1peak_sampleneg)
            sample = 0;
            ampl = 0;
        end      
        
        output_ER(ii,jj,1) = sample;
        output_ER(ii,jj,2) = ampl;
        
        %%%% AND THEN HERE ANALYZE WHETHER THIS IS A SIGNIFICANT PEAK 
        %%%% SO IF IT CAN BE SEEN AS A REAL EARLY RESPONS (USE PRE_STIM_SD)
        
    end
end


%% plot to test

%%%% TESTING FOUND OUT ITS ONE SAMPLE OFF/ SHOULD BE SAMPLE 5515

plot((find(tt>0,1)+extrasamps:find(tt>0.1,1)),new_signal(find(tt>0,1)+extrasamps:find(tt>0.1,1)))
hold on

plot(sample,ampl,'r*')

timepoint_n1_peak = tt(1,sample);

%% Try 2: selecting P1/N1/P2/N2

% P1/N1 onset: between 2 and 20 ms
p1_samples_start = find(tt>0.002,1);
p1_samples_end = find(tt>0.02,1);

% N1 peak: 10 and 50 ms
n1_samples_start = find(tt>0.01,1);
n1_samples_finish = find(tt>0.05,1);

% P2/N2 onset: 50 and 150 ms
p2_samples_start = find(tt>0.05,1);
p2_samples_finish = find(tt>0.15,1);

% N2 peak: between 50 and 300 ms
n2_samples_start =  find(tt>0.05,1);
n2_samples_finish = find(tt>0.3,1);

% End N2: ... 

for ii = 1:size(cc_epoch_sorted_avg,1)
    % for every averaged stimulation
    for jj = 1:size(cc_epoch_sorted_avg,2)
    
    % baseline subtraction: take median of part of the averaged signal for
    % this stimulation pair before stimulation, which is the half of the
    % epoch
    baseline_tt = tt>-2 & tt<-.1;
    signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt));
    
    % subtract median baseline from signal
    new_signal = cc_epoch_sorted_avg(ii,jj,:) - signal_median; 
    % testplot new signal: plot(tt,squeeze(new_signal))
    
    % take area before the stimulation of the new signal and calculate its SD
    pre_stim_sd = std(new_signal(baseline_tt));
    
    % if the pre_stim_sd is smaller that the minimally needed SD, use this
    % the minSD as pre_stim_sd
    if pre_stim_sd < minSD
       pre_stim_sd = minSD;
    end  
    
    % Next: only do analysis on the channels that are not stimulated
    % Should we e.g. channel F01 and F02 in stimulation of
    % F01-F02 load the data in as NaN?
    % so if member, do nothing or assign NaN
    % elif: run rest of the loop and find peak
    
        % use peakfinder to find all positive and negative peaks and their
        % amplitude. 
        % tt are the samples of the epoch based on the Fs and -2.5 to 2.5
        % seconds sample of the total epoch 
        % As tt use first sample after timepoint 0  (+ extra samples to not take artifact)
        % till first sample after 0.5 seconds (1024 samples)
        [all_samppos, all_amplpos] = peakfinder(new_signal(1,find(tt>0,1)+extrasamps:find(tt>0.5,1)),sel,[],1);
        [all_sampneg, all_amplneg] = peakfinder(new_signal(1,find(tt>0,1)+extrasamps:find(tt>0.5,1)),sel,[],-1);

        % If the first selected sample is a peak, this is not a real peak,
        % so delete
        all_amplpos(all_samppos==1) = [];
        all_samppos(all_samppos==1) = [];
        all_amplneg(all_sampneg==1) = [];
        all_sampneg(all_sampneg==1) = [];
        
        % convert back timepoints based on tt
        all_samppos = all_samppos + find(tt>0,1) + extrasamps;
        all_sampneg = all_sampneg + find(tt>0,1) + extrasamps;
        
        % find possible p1 peaks by selecting sample range
        % IN PROGRESS
        temp_p1_peaks = [];
        (p1_samples_start <= all_samppos) && (all_samppos <= p1_samples_end); 
        
        % same for n1/p2/n2
        
        
        

        
    end
end