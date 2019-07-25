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

%% to remove database fields if you want to try other ROI

% select regions with more than 1000 averaged epochs
% top10 = min(maxk(connectivity_mat(:),20))
[stim_region, rec_region] = find(connectivity_mat >= 1000);
rec_stim_regions = [rec_region stim_region] - 1 % -1 for destrieux

%% find ROIs with more then 3 electrode to measure within ROI

ROI_destrieux = [15, 26, 38];

for destrieux_array = 1:length(ROI_destrieux)
    
    
    for subj = 1:length(database)
        database(subj).ROI_within_coverage = [];
        database(subj).ROI_stimpairs_data = [];
    end

    runs = 1;
    ses = 1;
    % find in which patients those areas are represented
    % iterate over all subjects in database
    % in atlas ROI is 1 less then in the labels, so add 1 in code
    
    full_latency_vector = nan(250,27);
    full_age_vector = nan(250,27);

    for subj = 1:length(database)

       % use counter for writing data in the right row of the matrix
       mat_row_counter = 1;
        
       % +1 for destrieux as label in table
        if sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux(destrieux_array) + 1)) > 2

            % print amount of electrodes
            print_num = sprintf('%.f',  sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux(destrieux_array) + 1)));    
            disp(['subject: ' database(subj).subject ' has at least 3 electrodes coverage of ROI. It has: ' print_num ' electrodes']); 


            % next step: which electrodes are this within the subject
            ROI_elec_loc = find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)== ROI_destrieux(destrieux_array) + 1);
            ROI_elec_name = database(subj).metadata(runs).atlases_electrodes.name(ROI_elec_loc);

            % write names of electrodes to the database struct
            database(subj).ROI_within_coverage = ROI_elec_name;


            % look in electrodes table which numbers correpond to these
            % ROI_elec_names. Create array in which the numbers of the 
            % electrodes can be written
            ROI_array = NaN(1,length(database(subj).ROI_within_coverage));  

            % for every stimulated electrode in ROI
            for ROI_elec_stim = 1:length(database(subj).ROI_within_coverage)

                % find the number that corresponds with the name
                ROI_array(1,ROI_elec_stim) = find(strcmp(database(subj).metadata(runs).atlases_electrodes.name,ROI_elec_name(ROI_elec_stim)));       

            end

            % iterate over the ROI array to find if the electrodes are in a 
            % stimulation pair of which the other electrode is also in the
            % pair. 
            for zz = 1:length(ROI_array)

                % find in which stimulation pair combination these electrodes
                % are
                ROI_stim_pair = find(database(subj).metadata(runs).stimulated_pairs(:,1) == (ROI_array(zz))); 

                % if the other electrode of the stimulated pair is also on the
                % ROI, and they are an actual combination of stimulations (this 
                % is only when the second electrodes is one higher than the first)
                % then this electrode pair is part of analysis
                if sum(ismember(ROI_array,database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,2))) >= 1 ...
                        && abs(database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,2) - ROI_array(zz)) == 1

                    % iterate over all other electrodes in ROI that are not stim_pair
                    % to find the latency and amplitude in those electrodes
                    for qq = 1:length(ROI_array)

                        % if an electrode in the ROI_array in not part of the
                        % stimulation, these data can be used. 
                        if sum(ismember(database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,:),ROI_array(qq))) == 0


                            % to ensure this combination of electrodes
                            % stimulation actually happened 
                            if ROI_array <= size(database(subj).all_spesclin_latency,2)

                                % save stim_pair under ROI_stimparis_data

                                % save stim_pair in first column
                                database(subj).ROI_stimpairs_data(mat_row_counter,1) = ROI_stim_pair;
                                % save measured elec in second column
                                database(subj).ROI_stimpairs_data(mat_row_counter,2) = ROI_array(qq);
                                % save latency in third column
                                database(subj).ROI_stimpairs_data(mat_row_counter,3) = database(subj).all_spesclin_latency(ROI_stim_pair,ROI_array(qq));
                                % save amplitude in fourth column
                                database(subj).ROI_stimpairs_data(mat_row_counter,4) = database(subj).all_spesclin_amplitude(ROI_stim_pair,ROI_array(qq));

                                % add one to counter to ensure no overwriting
                                mat_row_counter =  mat_row_counter + 1;

                            end
                        end
                    end                         
                end          
            end
        end 
        
%         database(subj).ROI_coverage_38 = database(subj).ROI_within_coverage;
%         database(subj).ROI_data_38 = database(subj).ROI_stimpairs_data;

    end
end

%% Specific ROI analysis

ROI_destrieux = 38;

% create empty matrix
ROI_plot_matrix = nan(250,length(subjects));


for subj = 1:length(database)
    database(subj).ROI_within_coverage = [];
    database(subj).ROI_stimpairs_data = [];
end

runs = 1;
ses = 1;

% find in which patients those areas are represented
% iterate over all subjects in database
% in atlas ROI is 1 less then in the labels, so add 1 in code
for subj = 1:length(database)
    
    % use counter for writing data in the right row of the matrix
    mat_row_counter = 1;
    
    % +1 for destrieux as label in table
    if sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux(destrieux_array) + 1)) > 2
        
        % print amount of electrodes
        print_num = sprintf('%.f',  sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux(destrieux_array) + 1)));
        disp(['subject: ' database(subj).subject ' has at least 3 electrodes coverage of ROI. It has: ' print_num ' electrodes']);
        
        
        % next step: which electrodes are this within the subject
        ROI_elec_loc = find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)== ROI_destrieux(destrieux_array) + 1);
        ROI_elec_name = database(subj).metadata(runs).atlases_electrodes.name(ROI_elec_loc);
        
        % write names of electrodes to the database struct
        database(subj).ROI_within_coverage = ROI_elec_name;
        
        
        % look in electrodes table which numbers correpond to these
        % ROI_elec_names. Create array in which the numbers of the
        % electrodes can be written
        ROI_array = NaN(1,length(database(subj).ROI_within_coverage));
        
        % for every stimulated electrode in ROI
        for ROI_elec_stim = 1:length(database(subj).ROI_within_coverage)
            
            % find the number that corresponds with the name
            ROI_array(1,ROI_elec_stim) = find(strcmp(database(subj).metadata(runs).atlases_electrodes.name,ROI_elec_name(ROI_elec_stim)));
            
        end
        
        % iterate over the ROI array to find if the electrodes are in a
        % stimulation pair of which the other electrode is also in the
        % pair.
        for zz = 1:length(ROI_array)
            
            % find in which stimulation pair combination these electrodes
            % are
            ROI_stim_pair = find(database(subj).metadata(runs).stimulated_pairs(:,1) == (ROI_array(zz)));
            
            % if the other electrode of the stimulated pair is also on the
            % ROI, and they are an actual combination of stimulations (this
            % is only when the second electrodes is one higher than the first)
            % then this electrode pair is part of analysis
            if sum(ismember(ROI_array,database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,2))) >= 1 ...
                    && abs(database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,2) - ROI_array(zz)) == 1
                
                % iterate over all other electrodes in ROI that are not stim_pair
                % to find the latency and amplitude in those electrodes
                for qq = 1:length(ROI_array)
                    
                    % if an electrode in the ROI_array in not part of the
                    % stimulation, these data can be used.
                    if sum(ismember(database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,:),ROI_array(qq))) == 0
                        
                        
                        % to ensure this combination of electrodes
                        % stimulation actually happened
                        if ROI_array <= size(database(subj).all_spesclin_latency,2)
                            
                            % save stim_pair under ROI_stimparis_data
                            
                            % save stim_pair in first column
                            database(subj).ROI_stimpairs_data(mat_row_counter,1) = ROI_stim_pair;
                            % save measured elec in second column
                            database(subj).ROI_stimpairs_data(mat_row_counter,2) = ROI_array(qq);
                            % save latency in third column
                            database(subj).ROI_stimpairs_data(mat_row_counter,3) = database(subj).all_spesclin_latency(ROI_stim_pair,ROI_array(qq));
                            % save amplitude in fourth column
                            database(subj).ROI_stimpairs_data(mat_row_counter,4) = database(subj).all_spesclin_amplitude(ROI_stim_pair,ROI_array(qq));
                            
                            % add one to counter to ensure no overwriting
                            mat_row_counter =  mat_row_counter + 1;
                            
                        end
                    end
                end
            end
        end
    end
    
    %%%% to transform data to be used for analysis
    
    % if there are data
    if ~isempty(database(subj).ROI_stimpairs_data)
        
        % this matrix can be used for analysis: latencies * subjects
        ROI_plot_matrix(1:length(database(subj).ROI_stimpairs_data),subj) = database(subj).ROI_stimpairs_data(:,3);
        
        
    end
    
    
end


%% Plot latencies within ROI

figure(ROI_destrieux), hold on
% iterate over all subjects in database
for subj = 1:length(database)
    
    % create age-vector for plotting
    age_vector(1,subj) = database(subj).age_ses1;
    
    
    % scatterplot
    scatter(repelem(age_vector(1,subj),length(ROI_plot_matrix)),ROI_plot_matrix(:,subj))

end

xlim([0 50]),ylim([10 90])
xlabel('age subject')
ylabel('latency in ms')
print_label = sprintf('%.f', ROI_destrieux(destrieux_array));
title(['latency within Destieux ROI: ' print_label])
hold off

%% Plot latencies in two groups: 18- and 18+

figure(100+ROI_destrieux), hold on

for gg = 1:length(age_vector)
    if age_vector(gg) <= 18
        age_group(gg) = 1;
    elseif age_vector(gg) > 18
        age_group(gg) = 2;
        
    end
end
   
boxplot(ROI_plot_matrix,age_group)


xlabel('age subject')
ylabel('latency in ms')
title('age-effect on latency two groups, 18- and 18+')
hold off

summed_cceps = (sum(~isnan(ROI_plot_matrix)));

CCEPs_18min = 0;
for aa = find(age_group == 1)

    CCEPs_18min = CCEPs_18min + summed_cceps(aa);
end

CCEPs_18plus = 0;
for bb = find(age_group == 2)

    CCEPs_18plus = CCEPs_18plus + summed_cceps(bb);
end
CCEPs_18min
CCEPs_18plus

%% Scatterplot percentages of CCEPs within ROI

figure(200+ROI_destrieux(destrieux_array)), hold on


for subj = 1:length(database)
    
    if ~isempty(database(subj).ROI_stimpairs_data)
        
        % scatterplot the age by the relative percentage of CCEPs found
        relative_cceps = (sum(~isnan(ROI_plot_matrix(:,subj))) / length(database(subj).ROI_stimpairs_data));
        scatter(age_vector(1,subj),relative_cceps,50,[0 0 0])
        
        
    end
end

xlim([0 50])%,ylim([0 1])
xlabel('age subject')
ylabel('percentage CCEPs')
title('age-effect on percentage of CCEPs in all electrodes')
hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% everything for between ROIs analysis

% optional: 
% select regions with more than 1000 averaged epochs
[stim_region, rec_region] = find(connectivity_mat >= 1000);
rec_stim_regions = [rec_region stim_region] - 1 % -1 for destrieux


% for now: use the same regions as within ROI











%% Scatter all data (full vectors) and plot correlation etc 

figure(ROI_destrieux(destrieux_array)), hold on

full_latency_vector = reshape(full_latency_vector,size(full_latency_vector,1)*size(full_latency_vector,2),1);
full_age_vector = reshape(full_age_vector,size(full_age_vector,1)*size(full_age_vector,2),1);
scatter(full_age_vector, full_latency_vector)
xlabel('age subject')
ylabel('latency in ms')
print_label = sprintf('%.f', ROI_destrieux(destrieux_array));
title(['latency within Destieux ROI: ' print_label])
refline()

%% correlation coefficient???
R = corrcoef(full_age_vector,full_latency_vector);
R_squared=R(2)^2;
scatter(X,Y)
text(x_text, y_text, ['R^2 = ' num2str(R_squared)])
