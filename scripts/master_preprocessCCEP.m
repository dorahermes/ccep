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
figure
plot(t,data(5,:))
xlabel('time(s)')
ylabel('amplitude(uV)')

%% re-reference (if we want)

%% epoch into channel X epoch X time

% load the _events.tsv
events_name = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_events.tsv']);

cc_events = readtable(events_name,'FileType','text','Delimiter','\t');


% the cc_events already has a cc_events.sample
% If I would use round(cc_events.onset*data_hdr.FS) it would give slightly
% different numbers.

epoch_length = 5; % in seconds, -2.5:2.5

% define the output structure
data_epoch = zeros(size(data,1),length(cc_events.onset),round(epoch_length*data_hdr.Fs));

% define starting time for each epoch/ should I do anything with this??
tt = [1:epoch_length*data_hdr.Fs]/data_hdr.Fs - 2.5; 

% loop through all events and add to the output structure
for elec = 1:size(data,1) % for all channels
    for l = 1:length(cc_events.onset) % for all epochs
        data_epoch(elec,l,:) = data(elec,cc_events.sample(l)-(epoch_length/2*data_hdr.Fs)+1:cc_events.sample(l)+(epoch_length/2*data_hdr.Fs));
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

% Reshape matrix to split into every unique combination of stimulations
% sites. Take mean of these 5 stimulations and squeeze into 3D matrix. 
avg_epoch_stimulation = squeeze(mean(reshape(data_epoch,size(data_epoch,1),...
    amount_same_stimulation,(size(cc_events,1)/amount_same_stimulation),size(data_epoch,3)),2));

% Look at averaged data of per stimulation combination. 
figure
plot(tt,squeeze(avg_epoch_stimulation(5,1,:)))
xlabel('time(s)')
% set(gca,'Ydir','reverse')
ylabel('amplitude(uV)')

%%