function [database] = ccep_data_preprocess(database, top_path)

% This function is used to preprocess the data, resulting in averaged
% epochs, which are saved in the database structure. 

% INPUT:
% database structure
    % (containing subjects, sessions, tasks and runs) (run ccep_load_database first) 
% top_path 
    % (path with data, e.g. '/Fridge/users/jaap/ccep/dataBIDS/' or '/Fridge/CCEP)
    
% OUTPUT:
% database structure, with data and name of datafiles added to them. 

% IMPORTANT: 
% Depends on Fieldtrip functions to load data. Add Fieldtrip to path (ft_defaults)

% Function includes:

% - Loading data, events, electrodes, channels
% - Removing annotated artifacts from data
% - Bad channel rejection
% - Does not include re-referencing
% - Splits epochs into 2.5s before and 2.5s after stimulation
% - Adds data, events, electrodes and channels files and filenames to
%       datastucture
% - Adds epoched data, averaged epoched data and stimulated pairs to
%       datastructure


% D. Hermes & J. van der Aar & Giulio Castegnaro, UMC Utrecht, 2019




% Fieldtrip has to be in the path!



% iterate over all subjects in database
for subj = 1:length(database)
    % iterate over all their runs
    for runs = 1:length(database(subj).metadata)
        
                
        % load the data, as BrainVision BIDS format
        ieeg_name = fullfile(top_path,['sub-' database(subj).metadata(runs).subject], ...
            ['ses-' database(subj).metadata(runs).session],'ieeg',...
            ['sub-' database(subj).metadata(runs).subject '_ses-' database(subj).metadata(runs).session ...
            '_task-' database(subj).metadata(runs).task '_run-' database(subj).metadata(runs).run '_ieeg.eeg']);
        data = ft_read_data(ieeg_name,'dataformat','brainvision_eeg');
        data_hdr = ft_read_header(ieeg_name,'dataformat','brainvision_eeg');

        
        % to do a quick check on how the data look, these lines of code
        % plot some data for small visual assessment
%         t = [1:size(data,2)]./data_hdr.Fs;
%         figure
%         plot(t,data(ii,:))
%         axis([0, 1000, -4000, 4000])
%         xlabel('time(s)')
%         ylabel('amplitude(uV)')
%         title(['sub: ' database(subj).metadata(runs).subject ' ses: ' database(subj).metadata(runs).session ...
%             ' task: ' database(subj).metadata(runs).task ' run: ' database(subj).metadata(runs).run])

        
        % load the events.tsv
        events_name = fullfile(top_path,['sub-' database(subj).metadata(runs).subject], ...
            ['ses-' database(subj).metadata(runs).session],'ieeg',...
            ['sub-' database(subj).metadata(runs).subject '_ses-' database(subj).metadata(runs).session ...
            '_task-' database(subj).metadata(runs).task '_run-' database(subj).metadata(runs).run '_events.tsv']);
        ccep_events = readtable(events_name,'FileType','text','Delimiter','\t');


        % load the electrodes.tsv
        electrodes_name = fullfile(top_path,['sub-' database(subj).metadata(runs).subject], ...
            ['ses-' database(subj).metadata(runs).session],'ieeg',...
            ['sub-' database(subj).metadata(runs).subject '_ses-' database(subj).metadata(runs).session ...
            '_electrodes.tsv']);
        electrodes_table = readtable(electrodes_name,'Filetype','text','Delimiter','\t');

        
        % might want to add a check here if the stimulations are not closer than 5
        % seconds  to eachother (might be the case in older patients)
        
        % set up index that iterates over artefact rows
        % recalculate to samples and change these samples into NaN
        for qq = find(strcmp(ccep_events.trial_type,'artefact'))'
            % if artefact is in all channels
            if strcmp(ccep_events.electrodes_involved_onset(qq), 'all')
                data(:,(ccep_events.sample_start(qq):str2double(ccep_events.sample_end(qq)))) = NaN;
                % if artefact is only in specific channels.
                % WARNING: it might be possible that there are channel specific
                % artifacts with the same onset but another offset. Therefore, there is
                % an extra 'else' included that will capture this problem if a
                % participant has this. Then this can be fixed in this loop. Same for
                % multiple channel specific artifacts at the same time
            elseif isequal(ccep_events.electrodes_involved_onset(qq),ccep_events.electrodes_involved_offset(qq))
                arti_channel = find(strcmp(electrodes_table.name,ccep_events.electrodes_involved_onset(qq)));
                data(arti_channel,(ccep_events.sample_start(qq):str2double(ccep_events.sample_end(qq)))) = NaN;
            else
                disp('Warning, there might be multiple channel specific artifacts')
                disp('These could have another onset, or offset')
                disp('Add this option to the loop to process these')
            end
        end

        % Exclude bad channels and channels that are not included, first
        % Load channels.tsv
        channels_name = fullfile(top_path,['sub-' database(subj).metadata(runs).subject], ...
            ['ses-' database(subj).metadata(runs).session],'ieeg',...
            ['sub-' database(subj).metadata(runs).subject '_ses-' database(subj).metadata(runs).session ...
            '_task-' database(subj).metadata(runs).task '_run-' database(subj).metadata(runs).run '_channels.tsv']);    
        channel_table = readtable(channels_name, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
        
        % iterate over all channels
        for channel = 1:height(channel_table)
            % if a channel is marked as bad, change data to NaN
            if strcmp(channel_table.status(channel),'bad')
                data(channel,:) = NaN;
                % if a channel is not included, change data to NaN
            elseif strcmp(channel_table.status_description(channel), 'not included')
                data(channel,:) = NaN;
            elseif strcmp(channel_table.status_description(channel), 'included')
                % do nothing
            else
                disp('check channels.tsv, because there is another annotation')
            end
        end
        
        % set epoch parameters
        epoch_length = 5; % in seconds, -2.5:2.5
        epoch_prestim = 2.5;
        
        % count how much stimulations there are
        total_stim_count = sum(strcmp(ccep_events.trial_type,'electrical_stimulation'));
        
        % define the output structure
        data_epoch = zeros(size(data,1),total_stim_count,round(epoch_length*data_hdr.Fs));
        
        % define time vector for all epochs (for plotting and knowing when tt==0 for onset stimulation)
        tt = [1:epoch_length*data_hdr.Fs]/data_hdr.Fs - epoch_prestim;
        
        % create cc_events_onlystims, which makes table of only the
        % stimulations and not e.g. artifact data
        ll_counter = 1;
        ccep_events_onlystims = table();
        for ll = 1:length(ccep_events.onset)
            if strcmp(ccep_events.trial_type(ll),'electrical_stimulation')
                ccep_events_onlystims(ll_counter,:) = ccep_events(ll,:);
                ll_counter = ll_counter + 1;
            end
        end
        
        % loop through all stimulations and add to the output structure
        for elec = 1:size(data,1) % for all channels
            for ll = 1:length(ccep_events_onlystims.onset) % for all epochs
                % if data contain artifacts, which are adjusted to NaN's
                % consider that epoch as unreliable and change epoch to NaN's
                if sum(isnan(data_epoch(elec,ll,(tt>-1 & tt<.5)))) > 0
                    data_epoch(elec,ll,:) = NaN;
                else
                    data_epoch(elec,ll,:) = data(elec,ccep_events_onlystims.sample_start(ll)-round((epoch_prestim*data_hdr.Fs))+1 ...
                    :ccep_events_onlystims.sample_start(ll)+round(((epoch_length-epoch_prestim)*data_hdr.Fs)));
                end
            end
        end
        
        % HERE AN EXTRA OPTION TO SELECT SPLITTING DATA INTO MONOPHASIC 
        % INSTEAD OF THE CURRENT METHOD, BIPHASIC
        
        % Averaging epochs assuming e.g. F01-F02 is the same as F02-F01
        
        % first extract stimulation sites from the cc_events_onlystims
        % for all rows in the table
        for ee = 1:size(ccep_events_onlystims,1)
            % extract the first stimulation site and add as column
            ccep_events_onlystims.electrical_stimulation_site_num_1(ee) = str2double(extractBefore(ccep_events_onlystims.electrical_stimulation_site_num(ee),'  '));
            % extract the second stimulation site and add as column
            ccep_events_onlystims.electrical_stimulation_site_num_2(ee) = str2double(extractAfter(ccep_events_onlystims.electrical_stimulation_site_num(ee),'  '));
        end
        
        %%%% Assume that F01-F02 is the same as F02-F01
        unique_pairs = sort([ccep_events_onlystims.electrical_stimulation_site_num_1 ccep_events_onlystims.electrical_stimulation_site_num_2],2);
        % IC has length all stimulations
        [ccep_stimsets,~,IC] = unique(unique_pairs,'rows'); % make sure that we only have a vector of nX1 not nX2
        % cc_stimsets has length of unique stimulation pairs (assuming F02-F01=F01-F2)
        ccep_nroftimes = zeros(size(ccep_stimsets,1),1); % number of times each pair is stimulated
        for kk = 1:size(ccep_stimsets,1)
            ccep_nroftimes(kk) = sum(ismember(unique_pairs,ccep_stimsets(kk,:),'rows'));
        end

        % maximum stimulations of 1 unique pair
        max_amount_same_stimulation = max(ccep_nroftimes);
        
        % first create NaN matrix based on the max amount of stimulations in the
        % data set. For every unique stimulation
        % Size: elec X nr of stimulations X stimtype X  time
        ccep_epoch_sorted = NaN(size(data_epoch,1),max_amount_same_stimulation,size(ccep_nroftimes,1),size(data_epoch,3));
        
        % for loop that runs through all unique stimulation pairs and for every pair
        % it runs through the amount of stimulations to fill in the NaN matrix
        % if stimulations of that pair < max_amount_stimulations, this should leave NaN's
        for yy = 1:length(ccep_nroftimes) % loop through epoch types (stimulated pairs)
            % find the epochs of the current type
            events_for_this_set = find(yy==IC);
            % see how many there are
            nr_events = length(events_for_this_set);
            % put things back in the matrix
            ccep_epoch_sorted(:,1:nr_events,yy,:) = data_epoch(:,events_for_this_set,:);
            % if 'All stimulations are done 5 times' is displayed earlier, isequal(cc_epoch_sorted,data_epoch) == True
        end
        
        % might only want to take the mean if there are at least ... stimulations?
        
        % Take mean of these stimulations and squeeze into 3D
        ccep_epoch_sorted_avg = squeeze(nanmean(ccep_epoch_sorted,2));

        
        % Option to plot and test the data:      
%         stim_pair_nr = 1;
%         ccep_elec = 3;
%         figure()
%         plot(tt(tt>-1 & tt<1),zeros(size(tt(tt>-1 & tt<1))),'Color',[.5 .5 .5]) 
%         plot(tt(tt>-1 & tt<1),squeeze(ccep_epoch_sorted_avg(ccep_elec,stim_pair_nr,(tt>-1 & tt<1))))
%         xlabel('time(s)')
%         ylabel('amplitude(uV)')
%         title(['elec ' database(subj).metadata(runs).data_hdr.label{ccep_elec} ' for stimulation of ' ...
%             database(subj).metadata(runs).data_hdr.label{ccep_stimsets(stim_pair_nr,1)} ...
%             ' and ' database(subj).metadata(runs).data_hdr.label{ccep_stimsets(stim_pair_nr,2)} ])
          
        
        % add data paths and data to the database structure:
        % not all might be necessary for further processing
        
        % add data file name and data to database struct
        database(subj).metadata(runs).ieeg_filename = ieeg_name;
        database(subj).metadata(runs).data = data;
        database(subj).metadata(runs).data_hdr = data_hdr;
        
        % add events file name and events to database struct
        database(subj).metadata(runs).events_filename = events_name;
        database(subj).metadata(runs).events = ccep_events;
        
        % add selection with stimulations only from events.tsv
        database(subj).metadata(runs).events_onlystims = ccep_events_onlystims;
   
        % add electrodes file name and electrodes tabel to database struct
        database(subj).metadata(runs).electrodes_filename = electrodes_name;
        database(subj).metadata(runs).electrodes = electrodes_table;
        
        % add chan        % loop through all stimulations and add to the output structure
        for elec = 1:size(data,1) % for all channels
            for ll = 1:length(ccep_events_onlystims.onset) % for all epochs
                % if data contain artifacts, which are adjusted to NaN's
                % consider that epoch as unreliable and change epoch to NaN's
                if sum(isnan(data_epoch(elec,ll,(tt>-1 & tt<.5)))) > 0
                    data_epoch(elec,ll,:) = NaN;
                else
                    data_epoch(elec,ll,:) = data(elec,ccep_events_onlystims.sample_start(ll)-round((epoch_prestim*data_hdr.Fs))+1 ...
                    :ccep_events_onlystims.sample_start(ll)+round(((epoch_length-epoch_prestim)*data_hdr.Fs)));
                end
            end
        endnels file name and electrodes tabel to database struct
        database(subj).metadata(runs).channels_filename = channels_name;
        database(subj).metadata(runs).channels = channel_table;

        % add epoched data to the database structure
        database(subj).metadata(runs).epoched_data = data_epoch;

        % add stimulation pairs and amount of stimulation of those pairs to
        % database structure
        database(subj).metadata(runs).stimulated_pairs = ccep_stimsets;
        database(subj).metadata(runs).stimulated_nroftimes = ccep_nroftimes;
        
        % add averaged epochs to database structure
        database(subj).metadata(runs).epoched_data_avg = ccep_epoch_sorted_avg;
        
    end
end





