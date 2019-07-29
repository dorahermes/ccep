% This script is used to create a structure named 'dataBase' that can hold
% all metadata and data of the subjects. This 'dataBase' structure can be
% used to put in several ccep functions. 

% Structure is created using BIDS-format. The script walks through the
% folders of the subjects to find the right metadata and reads this into
% the structure. Later on, by using other functions, this structure can be
% used to e.g. add averaged epochs to them. 

% INPUT
% subjects = a cell array holding the subjects you want analyse.
%               e.g. {'RESP0001', 'RESP0002'}
% topPath = folder in which the subjects can be found (begin of the
%               BIDS-structure. E.g. 'Fridge/CCEP' or '/Fridge/users/jaap/ccep/dataBIDS/'

% IMPORTANT: function uses dirPlus function. Add this to path. 

% By Jaap van der Aar & Dorien van Blooijs, 05-2019, UMC Utrecht


% list of subjects that you want to include in analysis
subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768'};
% main path
topPath = '/Fridge/users/jaap/ccep/dataBIDS/';

% create dataBase structure
dataBase = struct('subj',cell(size(subjects)));

% for all subjects, create subject struct withing dataBase
for ss = 1:length(subjects)
    dataBase(ss).subj = subjects{ss};
    
    % First create the sesPath by adding the the subj to topPath. 
    sesPath = fullfile(topPath,['sub-' subjects{ss}]);
    % look in the sesPath. sesList returns only names of the folders and
    % only in that folder and not any deeper because depth = 0
    sesList = dirPlus(sesPath, 'ReturnDirs', true, 'Depth',0);
    
    % if there is just one sessions, add metadata that are found in these 
    % folders to the structure
    if length(sesList) == 1
      % for all sessions in the sesPath folder create struct. 
      % because the sesPath folder only contains ses folders, we can iterate
      % over the length of the sesList
      for ses = 1:length(sesList)

        % The number of events.tsv files correspond with the number of
        % different runs, so use dirPlus to find this amount
        runList = dirPlus(sesList{ses}, 'FileFilter', '\_events.tsv$');

        % create metaData structure, to ensure it will be cleared every
        % iteration. Size not defined yet, because we are not sure how
        % much we want to add to the structure
        metaData = {};

        % for all runs found in the scans.tsv
        for runs = 1:length(runList)

            % write subject to metadata (it is already in dataBase, but
            % useful to have it here as well
            metaData(runs).subject = subjects{ss};

            % extract name of the session
            sesIdx = strfind(sesList{ses},'ses');
            session = sesList{ses}(sesIdx(end)+4:end);

            % write session into metadata
            metaData(runs).session = session;

            % extract name of the task
            taskIdxStart = strfind(runList{runs}, 'task-');
            taskIdxEnd = strfind(runList{runs}, '_run');
            task = runList{runs}(taskIdxStart+5:taskIdxEnd-1);

            % write task into second data column
            metaData(runs).task = task;

            % extract name of the run
            runIdxStart = strfind(runList{runs}, 'run-');
            runIdxEnd = strfind(runList{runs}, '_events.tsv');
            run = runList{runs}(runIdxStart+4:runIdxEnd-1);

            % write run into third data column
            metaData(runs).run = run;

        end
      end     
      
      % write metaData structure to dataBase structure
      dataBase(ss).metaData = metaData;
    
    % in some cases there are two sessions. In this case the metadata of the second 
    % sesson have to be added in the dataBase struct on the lines underneath
    % the metadata of the first session.
    elseif length(sesList) >= 2
        
        % create metaData structure, to ensure it will be cleared every
        % iteration. Size not defined yet, because we are not sure how
        % much we want to add to the structure
        metaData = {};
            
        % for every session
        for ses = 1:length(sesList)
            
            % The number of events.tsv files correspond with the number of
            % different runs, so use dirPlus to find this amount
            runList = dirPlus(sesList{ses}, 'FileFilter', '\_events.tsv$');
            
            % for the first session, no adjustments to code are made
            if ses == 1
               
                % save number of runs in first session
                runCounter = length(runList);
                               
                for runs = 1:length(runList)
                    
                    % write subject to metadata (it is already in dataBase, but
                    % useful to have it here as well
                    metaData(runs).subject = subjects{ss};
                    
                    % extract name of the session
                    sesIdx = strfind(sesList{ses},'ses');
                    session = sesList{ses}(sesIdx(end)+4:end);
                    
                    % write session into metadata
                    metaData(runs).session = session;
                    
                    % extract name of the task
                    taskIdxStart = strfind(runList{runs}, 'task-');
                    taskIdxEnd = strfind(runList{runs}, '_run');
                    task = runList{runs}(taskIdxStart+5:taskIdxEnd-1);
                    
                    % write task into second data column
                    metaData(runs).task = task;
                    
                    % extract name of the run
                    runIdxStart = strfind(runList{runs}, 'run-');
                    runIdxEnd = strfind(runList{runs}, '_events.tsv');
                    run = runList{runs}(runIdxStart+4:runIdxEnd-1);
                    
                    % write run into third data column
                    metaData(runs).run = run;
                    
                end  
            end 
            
            % for next sessions, add the runCounter to the runs to ensure
            % the data will be written in the rows under the rows under the
            % first session
            if ses > 1
                
                % add runCounter to length of runList
                for runs = (1:length(runList))+runCounter
                    
                    % write subject to metadata (it is already in dataBase, but
                    % useful to have it here as well
                    metaData(runs).subject = subjects{ss};
                    
                    % extract name of the session
                    sesIdx = strfind(sesList{ses},'ses');
                    session = sesList{ses}(sesIdx(end)+4:end);
                    
                    % write session into metadata
                    metaData(runs).session = session;
                    
                    % extract name of the task (runCouter subtracted)
                    taskIdxStart = strfind(runList{runs-runCounter}, 'task-');
                    taskIdxEnd = strfind(runList{runs-runCounter}, '_run');
                    task = runList{runs-runCounter}(taskIdxStart+5:taskIdxEnd-1);
                    
                    % write task into second data column
                    metaData(runs).task = task;
                    
                    % extract name of the run (runCouter subtracted)
                    runIdxStart = strfind(runList{runs-runCounter}, 'run-');
                    runIdxEnd = strfind(runList{runs-runCounter}, '_events.tsv');
                    run = runList{runs-runCounter}(runIdxStart+4:runIdxEnd-1);
                    
                    % write run into third data column
                    metaData(runs).run = run;             
                    
                    % add latest list to runCounter in case there is a
                    % third session
                    runCounter = runCounter + length(runList);
                    
                end   
            end    
        end       
        
        % write metaData structure to dataBase structure
        dataBase(ss).metaData = metaData;

    end
    
end

% then you can put in this struct whatever you want e.g. avg epochs
 dataBase(ss).metaData(runs).test = avg_epoch
