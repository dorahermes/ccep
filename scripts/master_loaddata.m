% Function to load BIDS data
%
% Depends on Fieldtrip functions to load data
%
% D. Hermes & J. van der Aar 2018

clear all

addpath(genpath('/Users/dora/Documents/git/ccep/'))

%% load data

dataRootPath = '/Users/dora/Documents/data/ccep/';

sub_label = 'RESP0733';
ses_label = '1day3x1149';
task_label = 'SPESclin';

dataName = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_ieeg.eeg']);

% give an error in case the data do not exist:
if ~exist(dataName)==2
    disp('Error: dataName incorrect')
end

% %%%% load the data
% ccep_data      = ft_read_data(dataName,'dataformat','brainvision_eeg');
% ccep_header    = ft_read_header(dataName);
%%%% there is a problem in loading the brainvision data

% let's try the .trc file:
dataName = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_ieeg.TRC']);
%%%%% load the data and events correctly using Erik's function
SETTINGS.loadevents.state = 'yes';
SETTINGS.loadevents.type = 'analog';
SETTINGS.loadevents.dig_ch1='';
SETTINGS.loadevents.dig_ch1_label='';
SETTINGS.loadevents.dig_ch2='';
SETTINGS.loadevents.dig_ch2_label='';
SETTINGS.chan_adjust_status=0;
SETTINGS.chan_exclude_status=0;
SETTINGS.chan_adjust='';
SETTINGS.chan_exclude='';
SETTINGS.chans=[];
SETTINGS.filename = dataName;
SETTINGS.verbose = 'no';
SETTINGS.loaddata.state = 'yes';
SETTINGS.loaddata.type = 'uV';
SETTINGS.movenonEEGchannels = 'yes';
d = jun_readtrc(SETTINGS);

% create a matrix in which all stimulation pairs (CO1-CO2) are grouped:
% you will get a matrix with: 
% channel(133) x time(4*srate, 2 sec before and 2 sec after) x events(550)
% sort by event type (stimulation pair)

%% now, save as brainvision:

%
%      %%%%% save data as BrainVision BIDS %%%%%
%
% Sample script that calls Fieldtrip functions to write a Brainvision dataset
% Added fields are examples, read these from the raw data
%
% Fieldtrip has to be in the path!
% 
% D. Hermes 2018
 
ieeg_name_save = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_ieeg.eeg']);
        
%%%% assign header fields:

% sampling frequency
dataStruct.hdr.Fs = d.srate; 
dataStruct.hdr.nChans = size(d.data,1); 

% number of channels
dataStruct.hdr.label = cell(dataStruct.hdr.nChans,1);
for kk = 1:dataStruct.hdr.nChans
    dataStruct.hdr.label{kk} = d.channel(kk).name;
end

% number of samples
dataStruct.hdr.nSamples = size(d.data,2); 

% ?
dataStruct.hdr.nSamplesPre = 0; 

% 1 trial for continuous data
dataStruct.hdr.nTrials = 1; 

% channels type, see BIDS list of types
dataStruct.hdr.chantype = cell(dataStruct.hdr.nChans,1); 
for kk = 1:dataStruct.hdr.nChans
    dataStruct.hdr.chantype{kk} = ['ECOG'];
end

% I still don't how to indicate the mu in BIDS, now using letter u
dataStruct.hdr.chanunit = repmat({'uV'},size(d.data,2),1);
        
% labels again, same as before
dataStruct.label = dataStruct.hdr.label; 

% time vector
dataStruct.time{1} = [1:dataStruct.hdr.nSamples]/dataStruct.hdr.Fs;

% put the data matrix here: electrodes x samples
dataStruct.trial{1} = d.data; 

% sampling freq again, same as before
dataStruct.fsample = dataStruct.hdr.Fs; 

% ?
dataStruct.sampleinfo = [1 dataStruct.hdr.nSamples];

% now fetch a header
hdr_data = ft_fetch_header(dataStruct);

% save the data
ft_write_data(ieeg_name_save,dataStruct.trial{1},'header',hdr_data,'dataformat','brainvision_eeg')

