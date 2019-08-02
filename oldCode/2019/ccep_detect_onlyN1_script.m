%% load in peakfinder
addpath(genpath('/Fridge/users/giulio/github/ccep/functions')) 

%% set defined parameters
% pre-allocation: variables determined in Dorien's thesis
thresh = 2.5;   %
minSD = 50;     % in microV: minimum standard deviation of prestimulus baseline in case actual std of prestim baseline is too small
sel = 20;       % how many samples around peak not considered as another peak
extrasamps = 5; % set to 5 now, to read properties of N1 onset
SDfactor=[];
SD = [];

%% Run through all avg epochs of all channels to find peaks

%%%% this script is fully based on Dorien's function

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
    % when the electrode is stimulated
    if ii == cc_stimsets(jj,1) || ii == cc_stimsets(jj,2)
            n1_peak_sample = NaN;
            n1_peak_amplitude = NaN;
    else 
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
            sample = NaN;
            ampl = NaN;
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
end 


%% plot to test in case E1-E2 is thought as different from E2-E1 
figure

% pick one stimulation pair
%odd numbers will give a stimulation pair and the following even number
%will show the stimulation in the opposite direction. 
jj = 39;

suptitle(['Electrode by electrode analysis of the stimulation of' {cc_stimsets(jj,1)} ' to ' {cc_stimsets(jj,2)}])

% plot all measured channels
nr_column = 10;
nr_rows = ceil(length(channel_table.name)/nr_column);

% makes it possible to plot multiple CCEPs to check
for ii = 1:length(channel_table.name)
    if isequal(channel_table.type{ii},'ECOG')
        subplot(nr_rows,nr_column,ii),hold on
        % recalculate new_signal for this ii and jj to plot
        signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt),3);
        new_signal = squeeze(cc_epoch_sorted_avg(ii,jj,:)) - signal_median; 
        plot(tt(find(tt>0,1)+extrasamps:find(tt>0.1,1)),new_signal(find(tt>0,1)+extrasamps:find(tt>0.1,1)))

        hold on
        if ~isnan(output_ER(ii,jj,1)) % if a N1 is detected
            plot(tt(output_ER(ii,jj,1)),output_ER(ii,jj,2),'r*')
        end
    end
end


%elecpair (jj,:) %to see exactly which electrodes are stimulated 

%% render brain + effects

%% Render plain with used electrodes (setting the parameters)
addpath('/Fridge/users/jaap/github/ecogBasicCode/render/')
dataRootPath = '/Fridge/users/giulio/ccep/dataBIDS/';
subjects = {'RESP0751'};
hemi_cap = {'R'};

% set stimulated pair you want to render
%stim_pair = 2;

s = 1; %1:length(subjects)
% subject code
subj = subjects{s};
hemi_small = {'r'};

% gifti file name:
dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
% load gifti:
g = gifti(dataGiiName);

% electrode locations name:
dataLocName = dir(fullfile(dataRootPath,['sub-' subj],'ses-1','ieeg',...
    ['sub-' subj '_ses-1_electrodes.tsv']));
dataLocName = fullfile(dataLocName(1).folder,dataLocName(1).name);
% load electrode locations
loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% channels locations name:
channelsName = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_channels.tsv']);
% load channels info
channel_info = readtable(channelsName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% surface labels 
surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'surf',...
    [hemi_small{s} 'h.wang15_mplbl.mgz']);
surface_labels = MRIread(surface_labels_name);
vert_label = surface_labels.vol(:);


% put electrode positions in correct order of data
elecmatrix = NaN(size(data_epoch,1),3); % number of channels in data X xyz
for kk = 1:size(data_epoch,1) % number of channels in data
    % find matching electrode for this channel
    for mm = 1:length(loc_info.name) % loop through electrodes untill we find a match
        if isequal(channel_info.name{kk},loc_info.name{mm}) %channel name matches electrode name
            elecmatrix(kk,:) = [loc_info.x(mm) loc_info.y(mm) loc_info.z(mm)];
        end
    end
end

%% Rendering 

% pick a viewing angle:
v_dirs = [45 0]; %;90 0;90 -60;270 -60;0 0];
% make a colormap for the labelss
cmap = lines(max(vert_label));
Wang_ROI_Names = {...
        'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
        'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
        'IPS5' 'SPL1' 'FEF'};    
for stim_pair = 39
    
    % select significant peaks in the other channels
    n1_plot = squeeze(output_ER(:,stim_pair,:)); % measured electrodes X latency/ampl of N1 peak
    % remove positive N1
    n1_plot(n1_plot(:,2)>0,2) = NaN;
    
    
   % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
  
        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle     
        title(['stimulation of ' data_hdr.label{cc_stimsets(stim_pair,1)} ' and ' data_hdr.label{cc_stimsets(stim_pair,2)}])

        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);


        ecog_Label(els, 30,0.0001) % add electrode positions
        % add electrodes
        ccep_el_add(els,[1 1 1 ])
        % give stimulated electrodes other color
        ccep_el_add(els(cc_stimsets(stim_pair,1:2),:),[0 0 0],50)
        % set size and color of channels with significant peak 
        % based on sample (from stimulation on, so -5120) and the amplitude
        % color = latency, size = amplitude
        % minus N1 --> the more negative the bigger the circle
        ccep_el_add_size_and_color(els,-n1_plot(:,2),(n1_plot(:,1)-5120),500,35)

        set(gcf,'PaperPositionMode','auto')
        % print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
        % ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end


%%
figure('Position',[0 0 400 300]),hold on
plot([0.002 0.002],[-550 350],'Color',[.5 .5 .5],'LineWidth',7)
% plot([0.010 0.010],[-550 350],'Color',[.5 .5 .5])
% pick one stimulation pair
jj = 39;
% makes it possible to plot multiple CCEPs to check
for ii = 55
    
    if isequal(channel_info.type{ii},'ECOG')
        % recalculate new_signal for this ii and jj to plot
        signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt),3);
        new_signal = squeeze(cc_epoch_sorted_avg(ii,jj,:)) - signal_median; 

        plot(tt,new_signal,'k')
        hold on
        if ~isnan(output_ER(ii,jj,1)) % if a N1 is detected
            plot(tt(output_ER(ii,jj,1)),output_ER(ii,jj,2),'r.','MarkerSize',30)
        end
    end
end
xlim([-.05 .25]),ylim([-550 350])
xlabel('Time (s)'),ylabel('Amplitude (uV)')
box off
set(gca,'FontName','Lato','FontSize',14)

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/home/giulio/figures/Ccep_electrode')

%%
%% parula map
%
cm = parula(99);
cm = [0 0 0; cm];

elsize=[15:(45-15)/(100-1):60];

figure('Color',[1 1 1],'Position',[30 50 150 300]),hold on
for k=1:100
    plot(1,k/2,'.','MarkerSize',elsize(k),'Color',cm(k,:))
end
ylim([-5 55])

box off
set(gca,'FontName','Lato','FontSize',14,'XTick',[])

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/home/giulio/figures/colormap_renderLatency')

%% grey scale for amplitude

figure('Color',[1 1 1],'Position',[30 10 150 300]),hold on
for k=1:100
    plot(1,5*k,'.','MarkerSize',elsize(k),'Color',[.5 .5 .5])
end
ylim([-20 520])

box off
set(gca,'FontName','Lato','FontSize',14,'XTick',[])

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/home/giulio/figures/colormap_renderAmplitude')
