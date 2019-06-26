

%% Add paths

addpath(genpath('/Fridge/users/jaap/github/ccep/'))

%% load database structure 

subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768'};

top_path = '/Fridge/users/jaap/ccep/dataBIDS/';

database = ccep_load_database(subjects, top_path);


%% 