% Function to load BIDS data
%
% Depends on Fieldtrip functions to load data
%
% D. Hermes & J. van der Aar & Giulio Castegnaro 2019

clear all

addpath(genpath('/Fridge/users/jaap/github/ccep/'))
addpath('/Fridge/users/jaap/github/fieldtrip/')
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
dataRootPath = '/Fridge/CCEP/';

% Subject information
subjects = {'RESP0706'}; % subject number
sessions = {'1'}; % session number
task_labels = {'SPESclin'};
run_labels = {'041501'};
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
% make a time vector t
t = [1:size(data,2)]./data_hdr.Fs;
for ii = 1
    figure
    plot(t,data(ii,:))
    axis([0, 1000, -4000, 4000])
    xlabel('time(s)')
    ylabel('amplitude(uV)')
    
end
%% re-reference (if we want)

%% epoch into channel X epoch X time

% load the _events.tsv
events_name = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_events.tsv']);

cc_events = readtable(events_name,'FileType','text','Delimiter','\t');

% the cc_events already has a cc_events.sample
% If I would use round(cc_events.onset*data_hdr.FS) it would give slightly
% different numbers.

% might want to add a check here if the stimulations are not closer than 5
% seconds  to eachother (might be the case in older patients)

% set epoch parameters
epoch_length = 5; % in seconds, -2.5:2.5
epoch_prestim = 2.5;

% define the output structure
data_epoch = zeros(size(data,1),length(cc_events.onset),round(epoch_length*data_hdr.Fs));

% define time vector for all epochs (for plotting and knowing when tt==0 for onset stimulation)
tt = [1:epoch_length*data_hdr.Fs]/data_hdr.Fs - epoch_prestim; 

% loop through all events and add to the output structure
for elec = 1:size(data,1) % for all channels
    for ll = 1:length(cc_events.onset) % for all epochs
        data_epoch(elec,ll,:) = data(elec,cc_events.sample(ll)-round((epoch_prestim*data_hdr.Fs))+1:cc_events.sample(ll)+round(((epoch_length-epoch_prestim)*data_hdr.Fs)));
    end
end

%% Make figures for specific epoch
figure
plot(tt,squeeze(data_epoch(5,1,:)))
xlabel('time(s)')
ylabel('amplitude(uV)')

%% Average the epochs with the same stimulation site

% Now we just assume this is always 5 stimulations, should we build
% in a check, so the 5 is based on if cc_events.electrical_stimulation_site
% is the same?
amount_same_stimulation = 5; 

% get the unique number of stimulated pairs:
[cc_stimsets,IA,IC] =  unique(cc_events.electrical_stimulation_site, 'stable');
% run through these and see how many there are
cc_nroftimes = zeros(length(cc_stimsets),1); % number of times each pair is stimulated
for kk = 1:length(cc_stimsets)
    cc_nroftimes(kk) = sum(ismember(cc_events.electrical_stimulation_site,cc_stimsets{kk}));
end

% check whether there is consistent stimulations, if all stimulations pairs
% are done 5 times. 
if cc_nroftimes == 5*(ones((height(cc_events)/5),1))
    disp('All stimulations are done 5 times')
    % cc_epoch_sorted_if_all5 = squeeze(mean(reshape(data_epoch,size(data_epoch,1),
    % max_amount_same_stimulation,(size(cc_events,1)/amount_same_stimulation),size(data_epoch,3)),2));
else
    disp('Caution: not all stimulations are done 5 times')
    for ll = 1:length(cc_nroftimes)
        if cc_nroftimes(ll) ~= 5
            xx = [cc_stimsets(ll),' is stimulated ',cc_nroftimes(ll), ' times.'];
            disp(xx)
        end
    end
end

%% Taking average of epochs 

% maximum stimulations of 1 unique pair 
max_amount_same_stimulation = max(cc_nroftimes);

% first create NaN matrix based on the max amount of stimulations in the
% data set. For every unique stimulation 
cc_epoch_sorted = NaN(size(data_epoch,1),(max_amount_same_stimulation*size(cc_nroftimes,1)),size(data_epoch,3));


% for loop tat runs through all unique stimulation pairs and for every pair
% it runs through the amount of stimulations to fill in the NaN matrix
% if stimulations of that pair < max_amount_stimulations, this should leave NaN's 
for yy = 1:length(cc_nroftimes)
    for zz = 1:cc_nroftimes(yy)
        cc_epoch_sorted(:,(IA(yy)+(zz-1)),:) = data_epoch(:,(IA(yy)+(zz-1)),:);
        % if 'All stimulations are done 5 times' is displayed earlier, isequal(cc_epoch_sorted,data_epoch) == True 
    end
end

% Reshape to split into groups of max_amount_stimulations (5), then take
% mean of these stimulations and squeeze into 3D
cc_epoch_sorted_avg = squeeze(mean(reshape(cc_epoch_sorted,size(cc_epoch_sorted,1),...
    max_amount_same_stimulation,size(cc_nroftimes,1),size(cc_epoch_sorted,3)),2));

%% plot avg epochs
figure
plot(tt,squeeze(cc_epoch_sorted_avg(3,1,:)))
xlabel('time(s)')
% set(gca,'Ydir','reverse')
ylabel('amplitude(uV)')

%% This part is just notes and test code i'll just leave for now
% cc_epoch_sorted_if_all5 = squeeze(mean(reshape(data_epoch,size(data_epoch,1),max_amount_same_stimulation,(size(cc_events,1)/amount_same_stimulation),size(data_epoch,3)),2));

% max_amount_same_stimulation = max(cc_nroftimes);
% stimsets_avg_epoch = NaN(size(data_epoch,1),max_amount_same_stimulation,size(data_epoch,3));
% 
% for zz = 1:length(cc_nroftimes)
%     cc_epoch_sorted =  data_epoch(:,IA(zz):IA(zz+1)-1,:); 
% end

max_amount_same_stimulation = max(cc_nroftimes);
stimsets_avg_epoch = NaN(size(data_epoch,1),max_amount_same_stimulation,size(data_epoch,3));


cc_epoch_sorted = NaN(size(data_epoch,1),max_amount_same_stimulation,size(data_epoch,3));
cc_epoch_sorted(:,1,:) = data_epoch(:,1,:); 

% squeeze(mean(reshape( ... 

cc_epoch_sorted = NaN(size(data_epoch,1),(max_amount_same_stimulation*size(cc_nroftimes,1)),size(data_epoch,3));

for yy = 1:length(cc_nroftimes)
    for zz = 1:cc_nroftimes(yy)
        cc_epoch_sorted(:,(IA(yy)+(zz-1)),:) = data_epoch(:,(IA(yy)+(zz-1)),:);
    end
end


cc_events_timer = 1;
for yy = 1:length(cc_nroftimes)
    cc_epoch_sorted = NaN(size(data_epoch,1),max_amount_same_stimulation,size(data_epoch,3));
    for zz = 1:cc_nroftimes(yy)
        cc_epoch_sorted(:,zz,:) = data_epoch(:,cc_events_timer,:);
        cc_events_timer = cc_events_timer + 1;
    end
end
% % generate a NaN array for all stimulation pairs, NaN, because they may not
% % all fill up and we want to be able to use nanmean, nanstd etc 
% cc_epoch_sorted = NaN(electrodes X max aantal stimulatiepairs(5) X time);
% % overschrijf de NaN voor aantal matches
% 



%%