
load('database_visual_31719.mat')

% Possibility to load saved database (it is stored in /home/giulio
%!! Attention, database3172019 is without subject 0405 which is now
%excluded(in database 3172019)

%% Add paths for database
clear all

addpath(genpath('/Fridge/users/giulio/github/ccep/'))

%% load database structure 

subjects = {'RESP0315','RESP0751','RESP0401','RESP0306', 'RESP0703'};

top_path = '/Fridge/users/giulio/ccep/dataBIDS/';

database = ccep_load_database(subjects, top_path);

%% Add paths for preprocess 

top_path = '/Fridge/users/giulio/ccep/dataBIDS/';
addpath('/Fridge/users/dora/github/fieldtrip/')
ft_defaults

%% Run preprocess function
stim_status = 2; % 1 for monophasic, 2 for biphasic

database = ccep_data_preprocess(database, top_path, stim_status);

% Run detect N1 function 
amplitude_thresh = 2.8; 
n1_peak_range = 70; % Detected N1 only range between 10 and n1_peak_range (10-70ms after stimulation)

database = ccep_detect_n1peak(database, amplitude_thresh, n1_peak_range); 

%% Save the variables for later faster use 

save('database_visual_31719', '-v7.3') %this database is without 0405!!!
 
