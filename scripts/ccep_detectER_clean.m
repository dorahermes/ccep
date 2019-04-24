%% load in peakfinder
addpath(genpath('/Fridge/users/jaap/github/ccep/functions')) 


%% shortcut for loading in data without using preprocessCCEP script

% loads cc_epoch_sorted_avg and tt
cd /Fridge/users/jaap/github/ccep/scripts/
load('/home/jaap/data_avg_epoch_and_timevector')
%% set defined parameters
% pre-allocation: variables determined in Dorien's thesis
thresh = 2.5;   %
minSD = 50;     % in microV: minimum standard deviation of prestimulus baseline in case actual std of prestim baseline is too small
sel = 20;       % how many samples around peak not considered as another peak
extrasamps = 5; % set to 5 now, to read properties of N1 onset
SDfactor=[];
SD = [];

%% Option 1: Run through all avg epochs of all channels to find peaks

%%%% this script is fully based on Dorien's function

% this channel + stimulation gives nice N1 to test with
% ii = 3; jj = 1; 

% output in channels X stimulations X [onset(seconds) & amplitude]
output_ER = NaN(size(cc_epoch_sorted_avg,1),size(cc_epoch_sorted_avg,2),2);

% for every channel
for ii = 1:size(cc_epoch_sorted_avg,1)
    % for every averaged stimulation
    for jj = 1:size(cc_epoch_sorted_avg,2)
    
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
            % tt, necessary because peakfinder adjusted this. Because find
            % function takes first sample after stimulation, subtract 1
            % sample
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
            % tt, necessary because peakfinder adjusted this. Because find
            % function takes first sample after stimulation, subtract 1
            % sample
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
        
        % when peak amplitude is saturated, it is deleted
        if abs(ampl) > 3000
            sample = NaN;
            ampl = NaN;
        end
    
        % is the peak big enough to consider as peak
        if abs(ampl) < thresh* abs(pre_stim_sd)
            sample = NaN;
            ampl = NaN;
        end
    
        output_ER(ii,jj,1) = sample;
        output_ER(ii,jj,2) = ampl;

        
    end
end


%% plot to test

% makes it possible to plot multiple CCEPs to check
for ii = [3,11]
    for jj = 1:5
        
        % recalculate new_signal for this ii and jj to plot
        signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt),3);
        new_signal = squeeze(cc_epoch_sorted_avg(ii,jj,:)) - signal_median; 
    
        plot((find(tt>0,1)+extrasamps:find(tt>0.1,1)),new_signal(find(tt>0,1)+extrasamps:find(tt>0.1,1)))
        hold on

        plot(output_ER(ii,jj,1),output_ER(ii,jj,2),'r*')
    end
end


%% Option 2: selecting P1/N1/P2/N2 - variables 

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

%% Script to first find N1 and if significant, find other properties

% output in channels X stimulations X [p1 sample, p1 amplitude , n1 sample
% n1 amplitude, p2 sample, p2 amplitude, n2 sample, n2 amplitude]
output_ER_all = NaN(size(cc_epoch_sorted_avg,1),size(cc_epoch_sorted_avg,2),8);

for ii = 1:size(cc_epoch_sorted_avg,1)
    % for every averaged stimulation
    for jj = 1:size(cc_epoch_sorted_avg,2)
    
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
    end
        
end

%% plot to test

% makes it possible to plot multiple CCEPs to check
for ii = 3
    for jj = [1,8,10]
        
        if ~isnan(output_ER_all(ii,jj,3)) % exclude this if, to see how the ones without N1 look 
            
        % recalculate new_signal for this ii and jj to plot
        signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt),3);
        new_signal = squeeze(cc_epoch_sorted_avg(ii,jj,:)) - signal_median; 
    
        plot((find(tt>-0.1,1)+extrasamps:find(tt>0.3,1)),new_signal(find(tt>-0.1,1)+extrasamps:find(tt>0.3,1)))
        hold on

        plot(output_ER_all(ii,jj,1),output_ER_all(ii,jj,2),'b*')
        plot(output_ER_all(ii,jj,3),output_ER_all(ii,jj,4),'r*')
        plot(output_ER_all(ii,jj,5),output_ER_all(ii,jj,6),'g*')
        plot(output_ER_all(ii,jj,7),output_ER_all(ii,jj,8),'y*')
        
        end 
    end
end

%% render brain + effects



%% Render plain with used electrodes
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';
subjects = {'RESP0706'};
hemi_cap = {'R'};

% pick a viewing angle:
v_dirs = [90 0];%;90 0;90 -60;270 -60;0 0];

stim_pair = 1;
% which electrodes are stimulated? can we render these also?)
n1_plot = squeeze(output_ER_all(:,stim_pair,3:4)); % measured electrodes X latency/ampl

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
    % electrode locations name:
    dataLocName = dir(fullfile(dataRootPath,['sub-' subj],'ses-1','ieeg',...
        ['sub-' subj '_ses-1_task-SPESclin_run-*_electrodes.tsv']));
    dataLocName = fullfile(dataLocName(1).folder,dataLocName(1).name);
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
%         ecog_Label(els,30,12) % add electrode positions
        el_add(els,[0 0 0],20)
        ccep_eladd_sizecolor(els,n1_plot(:,1),n1_plot(:,2),1,100)

        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
%             ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end
