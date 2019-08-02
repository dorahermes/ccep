% Function to load BIDS data
%
% Depends on Fieldtrip functions to load data
%
% D. Hermes & J. van der Aar & Giulio Castegnaro 2019

clear all

addpath(genpath('/Fridge/users/giulio/github/ccep/'))
addpath('/Fridge/users/giulio/github/fieldtrip/')
ft_defaults


%% load raw data (brainvision format):

%
%      %%%%% load data as BrainVision BIDS %%%%%
%
% Sample script that calls Fieldtrip functions to write a Brainvision dataset
% Added fields are examples, read these from the raw data
%
% Fieldtrip has to be in the path!

% Setting data root path
dataRootPath = '/Fridge/users/giulio/ccep/dataBIDS/';

% Subject information
subjects = {'RESP0751'}; % subject number
sessions = {'1'}; % session number
task_labels = {'SPESclin'};
run_labels = {'021513'};
data_nr = 1; % dataset number in the above cell array

sub_label = subjects{data_nr};
ses_label = sessions{data_nr};
task_label= task_labels{data_nr};
run_label = run_labels{data_nr};

ieeg_name = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_ieeg.eeg']);
        
% load the data
data = ft_read_data(ieeg_name,'dataformat','brainvision_eeg');
data_hdr = ft_read_header(ieeg_name,'dataformat','brainvision_eeg');

%% look at the signal from one channel
% % make a time vector t
% t = [1:size(data,2)]./data_hdr.Fs;
% for ii = 1
%     figure
%     plot(t,data(ii,:))
%     axis([0, 1000, -4000, 4000])
%     xlabel('time(s)')
%     ylabel('amplitude(uV)')
%     
% end
%% re-reference (if we want)

%% epoch into channel X epoch X time

% load the _events.tsv
events_name = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_events.tsv']);

cc_events = readtable(events_name,'FileType','text','Delimiter','\t');


electrodes_file = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],...
    'ieeg',['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);

electrodes_table = readtable(electrodes_file,'Filetype','text','Delimiter','\t');

% might want to add a check here if the stimulations are not closer than 5
% seconds  to eachother (might be the case in older patients)

% set up index that iterates over artefact rows 
% recalculate to samples and change these samples into NaN
for qq = find(strcmp(cc_events.trial_type,'artefact'))' 
    % if artefact is in all channels
    if strcmp(cc_events.electrodes_involved_onset(qq), 'all')
        data(:,(cc_events.sample_start(qq):str2double(cc_events.sample_end(qq)))) = NaN; 
    % if artefact is only in specific channels. 
    % WARNING: it might be possible that there are channel specific
    % artifacts with the same onset but another offset. Therefore, there is
    % an extra 'else' included that will capture this problem if a
    % participant has this. Then this can be fixed in this loop. Same for
    % multiple channel specific artifacts at the same time
    elseif isequal(cc_events.electrodes_involved_onset(qq),cc_events.electrodes_involved_offset(qq))
        arti_channel = find(strcmp(electrodes_table.name,cc_events.electrodes_involved_onset(qq)));
        data(arti_channel,(cc_events.sample_start(qq):str2double(cc_events.sample_end(qq)))) = NaN;
    else
        disp('Warning, there might be multiple channel specific artifacts')
        disp('These could have another onset, or offset')
        disp('Add this option to the loop to process these')
    end
end

% Exclude bad channels and channels that are not included
% Load channels.tsv
channel_table = readtable(fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_channels.tsv']),...
    'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
% iterate over all channels
for channel = 1:height(channel_table)
    % if a channel is marked as bad, change data to NaN
    if strcmp(channel_table.status(channel),'bad') 
       data(channel,:) = NaN; 
    % if a channel is not included, change data to NaN
    elseif strcmp(channel_table.status_description(channel), 'not included')
       data(channel,:) = NaN;
    elseif strcmp(channel_table.status_description(channel), 'included')
       % do nothing
    else
        disp('check channels.tsv, because there is another annotation')  
    end  
end

% set epoch parameters
epoch_length = 5; % in seconds, -2.5:2.5
epoch_prestim = 2.5;

% count how much stimulations there are 
total_stim_count = sum(strcmp(cc_events.trial_type,'electrical_stimulation'));

% define the output structure
data_epoch = zeros(size(data,1),total_stim_count,round(epoch_length*data_hdr.Fs));

% define time vector for all epochs (for plotting and knowing when tt==0 for onset stimulation)
tt = [1:epoch_length*data_hdr.Fs]/data_hdr.Fs - epoch_prestim; 

% set size of table with only stimlation events:
% - question - 
% how to initiate the size of the table, because then it needs the headers already
% also, is there an easier way to do this? e.g. Select only part of the
% original table instead of iterating over the rows and creating new one

% create cc_events_onlystims, which makes table of only the
% stimulations and not e.g. artifact data
ll_counter = 1;
for ll = 1:length(cc_events.onset)
    if strcmp(cc_events.trial_type(ll),'electrical_stimulation')
        cc_events_onlystims(ll_counter,:) = cc_events(ll,:);
        ll_counter = ll_counter + 1;
    end
end

% loop through all stimulations and add to the output structure
for elec = 1:size(data,1) % for all channels
    for ll = 1:length(cc_events_onlystims.onset) % for all epochs
        % if data contain artifacts, which are adjusted to NaN's
        % consider that epoch as unreliable and change epoch to NaN's
        if sum(isnan(data_epoch(elec,ll,(tt>-1 & tt<.5)))) > 0
            data_epoch(elec,ll,:) = NaN;
        else
            data_epoch(elec,ll,:) = data(elec,cc_events_onlystims.sample_start(ll)-round((epoch_prestim*data_hdr.Fs))+1:cc_events_onlystims.sample_start(ll)+round(((epoch_length-epoch_prestim)*data_hdr.Fs)));
        end
    end
end

%% Make figures for specific epoch
figure
plot(tt,squeeze(data_epoch(3,1,:)))
xlabel('time(s)')
ylabel('amplitude(uV)')

% %% Averaging the epochs - assuming F01-F02 is different from F02-F01
% %%%% First assume that F01-F02 is different from F02-F01
% % get the unique number of stimulated pairs:
% [cc_stimsets,IA,IC] =  unique(cc_events_onlystims.electrical_stimulation_site, 'stable');
% % run through these and see how many there are
% cc_nroftimes = zeros(length(cc_stimsets),1); % number of times each pair is stimulated
% for kk = 1:length(cc_stimsets)
%     cc_nroftimes(kk) = sum(ismember(cc_events_onlystims.electrical_stimulation_site,cc_stimsets{kk}));
% end
% 
% 
% % check whether there is consistent stimulations, if all stimulations pairs
% % are done 5 times. 
% if cc_nroftimes == 5*(ones(round(height(cc_events_onlystims)/5),1))
%     disp('All stimulations are done 5 times')
%     % cc_epoch_sorted_if_all5 = squeeze(mean(reshape(data_epoch,size(data_epoch,1),
%     % max_amount_same_stimulation,(size(cc_events,1)/amount_same_stimulation),size(data_epoch,3)),2));
% else
%     disp('Caution: not all stimulations are done 5 times')
%     for ll = 1:length(cc_nroftimes)
%         if cc_nroftimes(ll) ~= 5
%             xx = [cc_stimsets(ll),' is stimulated ',cc_nroftimes(ll), ' times.'];
%             disp(xx)
%         end
%     end
% end
% %% Averaging the epochs - assuming F01-F02 is different as F02-F01 - different method - now with a double of electrodes number as output 
% elecpair =  unique(cc_events_onlystims.electrical_stimulation_site, 'stable'); 
% 
% 
% % first extract stimulation sites from the cc_events_onlystims
% % for all rows in the table
% for ee = 1:size(cc_events_onlystims,1)
%     % extract the first stimulation site and add as column
%     cc_events_onlystims.electrical_stimulation_site_num_1(ee) = str2double(extractBefore(cc_events_onlystims.electrical_stimulation_site_num(ee),'  '));
%     % extract the second stimulation site and add as column
%     cc_events_onlystims.electrical_stimulation_site_num_2(ee) = str2double(extractAfter(cc_events_onlystims.electrical_stimulation_site_num(ee),'  '));
% end
% 
% %%%% Assume that F01-F02 is DIFFERENT than F02-F01
% 
% % get the unique number of stimulated pairs, but now it is a cell array 
% [C,IA,IC] =  unique(cc_events_onlystims.electrical_stimulation_site_num, 'stable'); %C pairs of stimulated electrodes, represented by number 
% 
% CC = split(C);
% cc_stimsets = str2double(CC); % now we have a stimsets in which E1-E2 is different than E2-E1
% 
% %count how many times a stimulation is performed 
% cc_nroftimes = accumarray (IC,1); % how many times the stim is performed 
% %value_counts = [C, stim_counts] %trying to display which stimulation is performed for how many times
% 
% % check whether there is consistent stimulations, if all stimulations pairs
% % are done 5 times. 
% if cc_nroftimes == 5*(ones(round(height(cc_events_onlystims)/5),1))
%     disp('All stimulations are done 5 times')
%     % cc_epoch_sorted_if_all5 = squeeze(mean(reshape(data_epoch,size(data_epoch,1),
%     % max_amount_same_stimulation,(size(cc_events,1)/amount_same_stimulation),size(data_epoch,3)),2));
% else
%     disp('Caution: not all stimulations are done 5 times')
%     for ll = 1:length(cc_nroftimes)
%         if cc_nroftimes(ll) ~= 5
%             xx = ['The pair ' cc_stimsets(ll,:),' is stimulated ',cc_nroftimes(ll), ' times.'];
%             disp(xx)
%         end
%     end
% end

%% Averaging the epochs - assuming F01-F02 is the same as F02-F01

% first extract stimulation sites from the cc_events_onlystims
% for all rows in the table
for ee = 1:size(cc_events_onlystims,1)
    % extract the first stimulation site and add as column
    cc_events_onlystims.electrical_stimulation_site_num_1(ee) = str2double(extractBefore(cc_events_onlystims.electrical_stimulation_site_num(ee),'  '));
    % extract the second stimulation site and add as column
    cc_events_onlystims.electrical_stimulation_site_num_2(ee) = str2double(extractAfter(cc_events_onlystims.electrical_stimulation_site_num(ee),'  '));
end

%%%% Assume that F01-F02 is the same as F02-F01
unique_pairs = sort([cc_events_onlystims.electrical_stimulation_site_num_1 cc_events_onlystims.electrical_stimulation_site_num_2],2);
% IC has length all stimulations 
[cc_stimsets,~,IC] = unique(unique_pairs,'rows'); % make sure that we only have a vector of nX1 not nX2
% cc_stimsets has length of unique stimulation pairs (assuming F02-F01=F01-F2)
cc_nroftimes = zeros(size(cc_stimsets,1),1); % number of times each pair is stimulated
for kk = 1:size(cc_stimsets,1)
    cc_nroftimes(kk) = sum(ismember(unique_pairs,cc_stimsets(kk,:),'rows'));
end

%% Taking average of epochs 

% maximum stimulations of 1 unique pair 
max_amount_same_stimulation = max(cc_nroftimes);

% first create NaN matrix based on the max amount of stimulations in the
% data set. For every unique stimulation 
% Size: elec X nr of stimulations X stimtype X  time
cc_epoch_sorted = NaN(size(data_epoch,1),max_amount_same_stimulation,size(cc_nroftimes,1),size(data_epoch,3));

% for loop that runs through all unique stimulation pairs and for every pair
% it runs through the amount of stimulations to fill in the NaN matrix
% if stimulations of that pair < max_amount_stimulations, this should leave NaN's 
for yy = 1:length(cc_nroftimes) % loop through epoch types (stimulated pairs)
    % find the epochs of the current type
    events_for_this_set = find(yy==IC);
    % see how many there are
    nr_events = length(events_for_this_set);
    
    % put things back in the matrix
    cc_epoch_sorted(:,1:nr_events,yy,:) = data_epoch(:,events_for_this_set,:);
    % if 'All stimulations are done 5 times' is displayed earlier, isequal(cc_epoch_sorted,data_epoch) == True 
end

% Take mean of these stimulations and squeeze into 3D
cc_epoch_sorted_avg = squeeze(nanmean(cc_epoch_sorted,2));

%% plot averaged epochs
stim_pair_nr = 39
ccep_elec = 55;
figure
plot(tt,squeeze(cc_epoch_sorted_avg(ccep_elec,stim_pair_nr,:)))
xlabel('time(s)')
% set(gca,'Ydir','reverse')
ylabel('amplitude(uV)')
title(['elec ' data_hdr.label{ccep_elec} ' for stimulation of ' data_hdr.label{cc_stimsets(stim_pair_nr,1)} ' and ' data_hdr.label{cc_stimsets(stim_pair_nr,2)} ])

