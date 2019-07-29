%% EVERYTHING FOR CONNECTIVITY MATRICES 

%% Add electrode atlases (without the MKR1+! and unneccesary channels) to structure

top_path = '/Fridge/users/jaap/ccep/dataBIDS/';
% iterate over all subjects in database
for subj = 1:length(database)
    % iterate over all their runs
    for runs = 1:length(database(subj).metadata)
        
             
        % load electrodes file labeled with atlases
        electrodes_file_atlases = fullfile(top_path,['sub-' database(subj).metadata(runs).subject], ...
            ['ses-' database(subj).metadata(runs).session],'ieeg',...
            ['sub-' database(subj).metadata(runs).subject '_ses-'  ...
            database(subj).metadata(runs).session '_electrode_positions_fouratlases.tsv']);
        
        electrodes_table_atlases = readtable(electrodes_file_atlases,'ReadVariableNames', true,'Filetype','text','Delimiter','\t');
        
        database(subj).metadata(runs).atlases_filename = electrodes_file_atlases;
        
        if mod(length(database(subj).total_grid),2) == 0 
            
            database(subj).metadata(runs).atlases_electrodes = electrodes_table_atlases(1:(length(database(subj).total_grid)),:);
              
        elseif mod(length(database(subj).total_grid),2) == 1
            
            database(subj).metadata(runs).atlases_electrodes = electrodes_table_atlases([1:64,66:(length(database(subj).total_grid))],:);
            
            
        end
    end
end

%% Make histogram of represented brain regions (subj based only)
destrieux_labels = nan(133,length(database));

% only include the first run
runs = 1; 

% iterate over all subjects in database
for subj = 1:length(database)

   
    
    % extract labels from table and select only the ones that are not NaN
    destrieux_labels(1:length(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),subj) = ...
        str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label);
  
end

% plot histogram of represented electrodes and their amount
figure,hist(destrieux_labels(:),[1:74])
title('Amount of electrodes on brain region using Destrieux Atlas')
xlabel('Brain region number')
ylabel('Amount of electrodes')
colormap(hot)
xlim([1 75])

%% Plot connectivity matrix 

% how many connections to each electrode would be possible 
connectivity_mat = zeros(75,75);


% iterate over all subjects in database
for subj = 1:length(database)
    % iterate over all their runs
    for runs = 1:length(database(subj).metadata)
        % for the SPESclin runs
        if strcmp(database(subj).metadata(runs).task, 'SPESclin')
            % run through all stimulations
            for stims = 1:length(database(subj).metadata(runs).stimulated_pairs)

                % add one to every recorded channel except for the ones that
                % are stimulated (e.g. add +1 to channel 3 till 64 one the
                % places (stim1,3:64) and (stim2,3:64). But then converted
                % to the brain re_jgions

                % find locations of the stimulated pairs
                stimnum_1 = database(subj).metadata(runs).stimulated_pairs(stims,1);
                stimnum_2 = database(subj).metadata(runs).stimulated_pairs(stims,1);
                
                if stimnum_1 <= length(database(subj).total_grid) && stimnum_2 <= length(database(subj).total_grid)
               
                    stimloc_1 = str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label(stimnum_1));            
                    stimloc_2 = str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label(stimnum_2));

                    % for all recorded channels
                    for zz = 1:length(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)
                        % if the row does not correspond with the stimulated
                        % electrodes, add +1 to the matrix
                        if zz ~= stimnum_1 || zz ~= stimnum_2

                            rec_loc = str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label(zz));
                            
                            if ~isnan(rec_loc) && ~isnan(stimloc_1) && ~isnan(stimloc_2) 
                                connectivity_mat(stimloc_1,rec_loc) = connectivity_mat(stimloc_1,rec_loc) + 1; 
                                connectivity_mat(stimloc_2,rec_loc) = connectivity_mat(stimloc_2,rec_loc) + 1; 
                            end

                        end
                    end
    %                 
    %                 database(subj).metadata(runs).atlases_electrodes
    %                 
    %                 find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label(stims)) == isnumeric)
    %                 
    %                 
    %                 ~isnan(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label))
    %                 
                end
                
            end
        end
    end
end


%% Connectivity matrix
figure,
imagesc(connectivity_mat(1:75,1:75),[0 10])

title('Connectivity matrix of possible connections in database (27 subj)')
xlabel('Brain region of Destrieux - recorded region ')
ylabel('Brain region of Destrieux - stimulated region ')

axis square
grid on

set(gca,'XMinorTick','on','GridColor', [0.9, 0.9, 0.9])
% grid on
colormap(hot)

%% Connectivity matrix - White is > 1000 stims
figure('Position',[0 0 800 800])
imagesc(connectivity_mat,[0 max(connectivity_mat(:))])
axis square
title('Connectivity matrix of possible connections in database (27 subj)')
xlabel('Brain region of Destrieux - recorded region ')
ylabel('Brain region of Destrieux - stimulated region ')
% xlabel('Target Visual Area'),ylabel('Source Visual Area')
set(gca,'FontName','Lato','FontSize',12,...
    'XDir','normal'); xtickangle(90)
cm = hot(1000);
cm = [0 0 0; cm];
cm(1000:max(connectivity_mat(:)),1) = [1];
cm(1000:max(connectivity_mat(:)),2) = [1];
cm(1000:max(connectivity_mat(:)),3) = [1];
colormap(cm)
hcb = colorbar;
title(hcb,'stimulations')
set(gcf,'PaperPositionMode','auto')

%% finding ROIs with >1000 epochs

% to know which combinations of stim and rec have > 1000 epochs
[stim_region, rec_region] = find(connectivity_mat > 1000);
rec_stim_regions = [rec_region stim_region] - 1;