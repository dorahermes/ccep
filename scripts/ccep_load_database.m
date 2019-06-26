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
    
    % if there is just one sessions, add metadata that are found in these 
    % folders to the structure
    if length(ses_list) == 1
      % for all sessions in the sesPath folder create struct. 
      % because the sesPath folder only contains ses folders, we can iterate
      % over the length of the sesList
      for ses = 1:length(ses_list)

        % The number of events.tsv files correspond with the number of
        % different runs, so use dirPlus to find this amount
        run_list = dirPlus(ses_list{ses}, 'FileFilter', '\_events.tsv$');

        % create metaData structure, to ensure it will be cleared every
        % iteration. Size not defined yet, because we are not sure how
        % much we want to add to the structure
        metadata = {};

        % for all runs found in the scans.tsv
        for runs = 1:length(run_list)

            % write subject to metadata (it is already in dataBase, but
            % useful to have it here as well
            metadata(runs).subject = subjects{subj};

            % extract name of the session
            ses_idx = strfind(ses_list{ses},'ses');
            session = ses_list{ses}(ses_idx(end)+4:end);

            % write session into metadata
            metadata(runs).session = session;

            % extract name of the task
            task_idx_start = strfind(run_list{runs}, 'task-');
            task_idx_end = strfind(run_list{runs}, '_run');
            task = run_list{runs}(task_idx_start+5:task_idx_end-1);

            % write task into second data column
            metadata(runs).task = task;

            % extract name of the run
            run_idx_start = strfind(run_list{runs}, 'run-');
            run_idx_end = strfind(run_list{runs}, '_events.tsv');
            run = run_list{runs}(run_idx_start+4:run_idx_end-1);

            % write run into third data column
            metadata(runs).run = run;

        end
      end     
      
      % write metaData structure to dataBase structure
      database(subj).metadata = metadata;
    
    % in some cases there are two sessions. In this case the metadata of the second 
    % sesson have to be added in the dataBase struct on the lines underneath
    % the metadata of the first session.
    elseif length(ses_list) >= 2
        
        % create metaData structure, to ensure it will be cleared every
        % iteration. Size not defined yet, because we are not sure how
        % much we want to add to the structure
        metadata = {};
            
        % for every session
        for ses = 1:length(ses_list)
            
            % The number of events.tsv files correspond with the number of
            % different runs, so use dirPlus to find this amount
            run_list = dirPlus(ses_list{ses}, 'FileFilter', '\_events.tsv$');
            
            % for the first session, no adjustments to code are made
            if ses == 1
               
                % save number of runs in first session
                run_counter = length(run_list);
                               
                for runs = 1:length(run_list)
                    
                    % write subject to metadata (it is already in dataBase, but
                    % useful to have it here as well
                    metadata(runs).subject = subjects{subj};
                    
                    % extract name of the session
                    ses_idx = strfind(ses_list{ses},'ses');
                    session = ses_list{ses}(ses_idx(end)+4:end);
                    
                    % write session into metadata
                    metadata(runs).session = session;
                    
                    % extract name of the task
                    task_idx_start = strfind(run_list{runs}, 'task-');
                    task_idx_end = strfind(run_list{runs}, '_run');
                    task = run_list{runs}(task_idx_start+5:task_idx_end-1);
                    
                    % write task into second data column
                    metadata(runs).task = task;
                    
                    % extract name of the run
                    run_idx_start = strfind(run_list{runs}, 'run-');
                    run_idx_end = strfind(run_list{runs}, '_events.tsv');
                    run = run_list{runs}(run_idx_start+4:run_idx_end-1);
                    
                    % write run into third data column
                    metadata(runs).run = run;
                    
                end  
            end 
            
            % for next sessions, add the runCounter to the runs to ensure
            % the data will be written in the rows under the rows under the
            % first session
            if ses > 1
                
                % add runCounter to length of runList
                for runs = (1:length(run_list))+run_counter
                    
                    % write subject to metadata (it is already in dataBase, but
                    % useful to have it here as well
                    metadata(runs).subject = subjects{subj};
                    
                    % extract name of the session
                    ses_idx = strfind(ses_list{ses},'ses');
                    session = ses_list{ses}(ses_idx(end)+4:end);
                    
                    % write session into metadata
                    metadata(runs).session = session;
                    
                    % extract name of the task (runCouter subtracted)
                    task_idx_start = strfind(run_list{runs-run_counter}, 'task-');
                    task_idx_end = strfind(run_list{runs-run_counter}, '_run');
                    task = run_list{runs-run_counter}(task_idx_start+5:task_idx_end-1);
                    
                    % write task into second data column
                    metadata(runs).task = task;
                    
                    % extract name of the run (runCouter subtracted)
                    run_idx_start = strfind(run_list{runs-run_counter}, 'run-');
                    run_idx_end = strfind(run_list{runs-run_counter}, '_events.tsv');
                    run = run_list{runs-run_counter}(run_idx_start+4:run_idx_end-1);
                    
                    % write run into third data column
                    metadata(runs).run = run;             
                    
                    % add latest list to runCounter in case there is a
                    % third session
                    run_counter = run_counter + length(run_list);
                    
                end   
            end    
        end       
        
        % write metaData structure to dataBase structure
        database(subj).metadata = metadata;

    end
    
end

% then you can put in this struct whatever you want e.g. avg epochs
% database(ss).metadata(runs).test = avg_epoch;