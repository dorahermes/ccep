
clear all

%% Add paths for database

addpath(genpath('/Fridge/users/jaap/github/ccep/'))

%% load database structure 

% list subjects that you want to include in analysis
subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768', 'RESP0294', 'RESP0306', ...
    'RESP0315', 'RESP0348', 'RESP0368', 'RESP0369', 'RESP0400', 'RESP0401', 'RESP0405', ...
    'RESP0435', 'RESP0449', 'RESP0450', 'RESP0458', 'RESP0467', 'RESP0468', 'RESP0477', ...
    'RESP0478', 'RESP0501', 'RESP0502', 'RESP0703', 'RESP0724', 'RESP0754', 'RESP0295'}; %RESP0295-1



% create top path in which the subjects can be found (and their data within
% the folders)
top_path = '/Fridge/users/jaap/ccep/dataBIDS/';

%% run create database structure function

% careful: this will clear database struct
database = ccep_load_database(subjects, top_path);


%% Add paths for preprocess 

top_path = '/Fridge/users/jaap/ccep/dataBIDS/';
addpath('/Fridge/users/dora/github/fieldtrip/')
ft_defaults

%% Run preprocess data function

stim_status = 2; % alternating polarity
% takes at least a few hours for all subj, so only do once. 
database = ccep_data_preprocess(database, top_path, stim_status);

%% Set parameters for n1 peak detection

amplitude_thresh = 3.6;   
n1_peak_range = 100;


%% Run detect n1 peaks function

database = ccep_detect_n1peak(database, amplitude_thresh, n1_peak_range);


%% Save

save('dbstruct_jaap_130819','database', '-v7.3')