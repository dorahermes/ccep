% Setting data root path
dataRootPath = '/Fridge/CCEP/';

% Subject information
db.subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768'}; % subject number
db.sessions = {'1', '1', '1b', '1'}; % session number
db.task_labels = {'SPESclin', 'SPESclin', 'SPESclin', 'SPESclin'};
db.run_labels = {'021147', '041501', '050941', '021704'};

% IN FUNCTION
data_nr = 1; % dataset number in the above cell array

% IN PREPROCESS
for zz = 1:size(data_nr,2)

    sub_label = db.subjects{db.data_nr};
    ses_label = db.sessions{db.data_nr};
    task_label = db.task_labels{db.data_nr};
    run_label = db.run_labels{db.data_nr};
    
end
