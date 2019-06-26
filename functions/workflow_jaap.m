

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

%% 

database = ccep_data_preprocess(database, top_path);