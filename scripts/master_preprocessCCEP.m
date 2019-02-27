% Function to load BIDS data
%
% Depends on Fieldtrip functions to load data
%
% D. Hermes & J. van der Aar & Giulio Castegnaro 2019

clear all

addpath(genpath('/Fridge/users/dora/github/ccep/'))
addpath('/Fridge/users/dora/github/fieldtrip/')
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
plot(t,data(1,:))
xlabel('time(s)')
ylabel('amplitude(uV)')

%% re-reference (if we want)

%% epoch into channel X epoch X time

% load the _events.tsv
events_name = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_events.tsv']);

cc_events = readtable(events_name,'FileType','text','Delimiter','\t');

% define the output structure

% loop through all events and add to the output structure



