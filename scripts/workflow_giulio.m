

%% Add paths for database

addpath(genpath('/Fridge/users/giulio/github/ccep/'))

%% load database structure 

subjects = {'RESP0315','RESP0751','RESP0401','RESP0405','RESP0306'};

top_path = '/Fridge/users/giulio/ccep/dataBIDS/';

database = ccep_load_database(subjects, top_path);

stim_status = 2; % 1 for monophasic, 2 for biphasic
 
%% Add paths for preprocess 

top_path = '/Fridge/users/giulio/ccep/dataBIDS/';
addpath('/Fridge/users/dora/github/fieldtrip/')
ft_defaults

%% Run preprocess function
database = ccep_data_preprocess(database, top_path, stim_status);
amplitude_thresh = 2.8; 
n1_peak_range = 70; 
% Run detect N1 function 
database = ccep_detect_N1peak(database, amplitude_thresh, n1_peak_range); 

