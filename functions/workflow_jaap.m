
clear all

%% Add paths for database

addpath(genpath('/Fridge/users/jaap/github/ccep/'))

%% load database structure 

subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768'};

top_path = '/Fridge/users/jaap/ccep/dataBIDS/';

database = ccep_load_database(subjects, top_path);


%% Add paths for preprocess 

top_path = '/Fridge/users/jaap/ccep/dataBIDS/';
addpath('/Fridge/users/dora/github/fieldtrip/')
ft_defaults

%% Run preprocess data function

database = ccep_data_preprocess(database, top_path);

%% Set parameters for n1 peak detection

amplitude_thresh = 2.8;   
n1_peak_range = 70;


%% Run detect n1 peaks function

database = ccep_detect_n1peak(database, amplitude_thresh, n1_peak_range);
