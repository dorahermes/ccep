function [database] = ccep_load_database(subjects, top_path)

% This function is used to create a structure named 'dataBase' that can hold
% all metadata and data of the subjects. This 'dataBase' structure can be
% used to put in several ccep functions. 

% Structure is created using BIDS-format. The script walks through the
% folders of the subjects to find the right metadata and reads this into
% the structure. Later on, by using other functions, this structure can be
% used to add e.g. averaged epochs. 
% (database(subj).metadata(runs).test = avg_epoch;) 

% INPUT
% subjects  = a cell array holding the subjects you want analyse.
%               e.g. {'RESP0001', 'RESP0002'}
% topPath   = folder in which the subjects can be found (begin of the
%               BIDS-structure. E.g. 'Fridge/CCEP' or '/Fridge/users/jaap/ccep/dataBIDS/'

% OUTPUT
% database  = structure including the subjects, all their sessions, and runs

% IMPORTANT: function uses dirPlus function. Add this to path. 

% By Jaap van der Aar & Dorien van Blooijs, 05-2019, UMC Utrecht


% create dataBase structure
database = struct('subject',cell(size(subjects)));

% for all subjects, create subject struct withing dataBase
for subj = 1:length(subjects)
    database(subj).subject = subjects{subj};
    
    % First create the sesPath by adding the the subj to topPath.
    ses_path = fullfile(top_path,['sub-' subjects{subj}]);
    % look in the sesPath. sesList returns only names of the folders and
    % only in that folder and not any deeper because depth = 0
    ses_list = dirPlus(ses_path, 'ReturnDirs', true, 'Depth',0);
    
    % set run_counter that keeps track of the amount of runs when there are
    % multiple sessions and therefore multiple events files
    run_counter = 0;
    
    % create metaData structure, to ensure it will be cleared every
    % iteration. Size not defined yet, because we are not sure how
    % much we want to add to the structure
    metadata = {};
    
    % for all sessions in the sesPath folder create struct.
    % because the sesPath folder only contains ses folders, we can iterate
    % over the length of the sesList
    for ses = 1:length(ses_list)
        
        % The number of events.tsv files correspond with the number of
        % different runs, so use dirPlus to find this amount
        run_list = dirPlus(ses_list{ses}, 'FileFilter', '\_events.tsv$');
          
        % for all runs found in the scans.tsv
        for runs = 1:length(run_list)
            
            % write subject to metadata (it is already in dataBase, but
            % useful to have it here as well
            metadata(runs + run_counter).subject = subjects{subj};
            
            % extract name of the session
            ses_idx = strfind(ses_list{ses},'ses');
            session = ses_list{ses}(ses_idx(end)+4:end);
            
            % write session into metadata
            metadata(runs + run_counter).session = session;
            
            % extract name of the task
            task_idx_start = strfind(run_list{runs}, 'task-');
            task_idx_end = strfind(run_list{runs}, '_run');
            task = run_list{runs}(task_idx_start+5:task_idx_end-1);
            
            % write task into second data column
            metadata(runs + run_counter).task = task;
            
            % extract name of the run
            run_idx_start = strfind(run_list{runs}, 'run-');
            run_idx_end = strfind(run_list{runs}, '_events.tsv');
            run = run_list{runs}(run_idx_start+4:run_idx_end-1);
            
            % write run into third data column
            metadata(runs + run_counter).run = run;
            
        end
        
        % add amount of runs to run_counter
        run_counter = run_counter + length(run_list);
        
    end
    
    % write metaData structure to dataBase structure
    database(subj).metadata = metadata;
       
end


