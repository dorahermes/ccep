
% list of subjects that you want to include in analysis
subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768'};
% main path
topPath = '/Fridge/users/jaap/ccep/dataBIDS/';

dataBase = [];
% for all subjects, create subject struct withing dataBase
for ss = 1:length(subjects)
    dataBase.subj{ss} = subjects{ss};
    
    % First create the sesPath by adding the the subj to topPath. 
    sesPath = fullfile(topPath,['sub-' subjects{ss}]);
    % look in the sesPath. sesList returns only names of the folders and
    % only in that folder and not any deeper because depth = 0
    sesList = dirPlus(sesPath, 'ReturnDirs', true, 'Depth',0);
    
    % for all sessions in the sesPath folder create struct. 
    % because the sesPath folder only contains ses folders, we can iterate
    % over the length of the sesList
    for tt = 1:length(sesList)
        dataBase.subj{ss}.ses{tt} = sesList{tt}
    end
    
    
    
    
    % then another for loop: for all tasks, create struct
    
    % the another for loop: for all runs, create struct
end

%% other methods

% using dirwalk
[pathNames, dirNames, fileNames] = dirwalk(sesPath);

% using dir()
fileinfo = dir('/Fridge/users/jaap/ccep/dataBIDS/');
fileinfo.name

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


db.resp
%%
% We need to create a list of the subjects we want to analyse: this is the
% list used for all Jaap's analyses, because there can be subjects added
% afterwards that we don't include and there can be subjects or datasets we
% exclude. This is necessary because we are still curating the database. 
subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768'}; % subject number
dataBase = [];
for ss = 1:length(subjects)
    if isequal(subjects{ss},'RESP0621')
        dataBase.subj = subjects{ss};
        % read dir
        % get data
        % potentially remove excluded data with a reason
        % if thisdata.dataname=='whatever', make it empty, end
        
    end
end

%% notes on structs

% list of subjects that you want to include in analysis
subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768'};
% main path
topPath = '/Fridge/users/jaap/ccep/dataBIDS/';

dataBase = [];
% for all subjects, create subject struct withing dataBase
for ss = 1:length(subjects)
    dataBase.subj{ss} = subjects{ss};
    
    % First create the sesPath by adding the the subj to topPath. 
    sesPath = fullfile(topPath,['sub-' subjects{ss}]);
    % look in the sesPath. sesList returns only names of the folders and
    % only in that folder and not any deeper because depth = 0
    sesList = dirPlus(sesPath, 'ReturnDirs', true, 'Depth',0);
    
    % for all sessions in the sesPath folder create struct. 
    % because the sesPath folder only contains ses folders, we can iterate
    % over the length of the sesList
    for tt = 1:length(sesList)
        dataBase.subj{ss}.ses{tt} = sesList{tt}
    end
    
    
    
    
    % then another for loop: for all tasks, create struct
    
    % the another for loop: for all runs, create struct
end

[pathNames, dirNames, fileNames] = dirwalk(sesPath);

fileinfo = dir('/Fridge/users/jaap/ccep/dataBIDS/RESP0621/');
fileinfo.name