%%%% This script is used for the analysis for of the percentage of cceps in
%%%% to the amount of stimulations. A measure of propagation

%%%% Includes: 
%%%% Within ROIs percentages - scatter and group level
%%%% Between ROIs percentages - scatter and group level
%%%% Within and between together 

%% Relative percentage CCEPs - within ROIs - scatter

%%%% RUN THE FOLLOWING LOOP (SAME AS IN ANALYSIS_WITHIN_ROIS SCRIPT) TO CREATE THE 
%%%% database(subj).ROI_stimpairs_data FOR A SPECIFIC ROI TO ANALYZE

%%%% after the loop, create figure with percentage of CCEPs relative to the
%%%% amount of stimulations 

% run for 15, 26, 38 or all three
ROI_list = [15, 26, 38];

runs = 1;
ses = 1;

for subj = 1:length(database)
    database(subj).amount_cceps = 0;
    database(subj).total_stims = 0;
end

figure(), hold on


for ROI = 1:size(ROI_list,2)
    
    % create empty matrix
    ROI_plot_matrix = nan(250,length(subjects));
    
    for subj = 1:length(database)
        database(subj).ROI_within_coverage = [];
        database(subj).ROI_stimpairs_data = [];
    end
    
    % find in which patients those areas are represented
    % iterate over all subjects in database
    % in atlas ROI is 1 less then in the labels, so add 1 in code
    for subj = 1:length(database)
        
        % use counter for writing data in the right row of the matrix
        mat_row_counter = 1;
        
        % +1 for destrieux as label in table
        if sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_list(ROI) + 1)) > 2
            
            % print amount of electrodes
            print_num = sprintf('%.f',  sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_list(ROI) + 1)));
            disp(['subject: ' database(subj).subject ' has at least 3 electrodes coverage of ROI. It has: ' print_num ' electrodes']);
            
            
            % next step: which electrodes are this within the subject
            ROI_elec_loc = find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)== ROI_list(ROI) + 1);
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
    

    
    
    for subj = 1:length(database)
        
        if ~isempty(database(subj).ROI_stimpairs_data)
            
            % scatterplot the age by the relative percentage of CCEPs found
            relative_cceps = (sum(~isnan(ROI_plot_matrix(:,subj))) / length(database(subj).ROI_stimpairs_data));
            scatter(age_vector(1,subj),relative_cceps,50,[0 0 0])
            
            
            % create variables and save to analyze as groups later
            % save subject specific in struct. 
            
            % add to current data in variable (because there could already
            % be data from another ROI)
            % save amount of CCEPs
            database(subj).amount_cceps = database(subj).amount_cceps + sum(~isnan(ROI_plot_matrix(:,subj)));
            % save total amount of stimulations 
            database(subj).total_stims = database(subj).total_stims + length(database(subj).ROI_stimpairs_data);
 
            
        end
    end
  
end

xlim([0 50]), ylim([0 1])
xlabel('age subject')
ylabel('percentage CCEPs')
title('age-effect on percentage of CCEPs: ROIs 15, 26, 38')
hold off


%% Relative percentage CCEPs - within ROIs - boxplot group level

figure(), hold on

for gg = 1:length(age_vector)
    if age_vector(gg) < 18
        age_group(gg) = 1;
    elseif age_vector(gg) >= 18
        age_group(gg) = 2;
        
    end
end

% create array of relative percentage CCEPs per subject
rel_perc_cceps_within = nan(1,length(database));

for subj = 1:length(subjects)
    
    rel_perc_cceps_within(subj) = database(subj).amount_cceps / database(subj).total_stims;
    
end

boxplot(rel_perc_cceps_within,age_group)

xlabel('age subject')
ylabel('latency in ms')
title('relative percentage cceps, 18- and 18+, ROI 15, 26 & 38')
hold off

%% Relative percentage CCEPs - between ROIs - scatter

% first concat the three different between ROI lists created in analysis_between_ROIs script 

 for subj = 1:length(database)
     
     database(subj).ROI_between_all = [database(subj).ROI_between_15_26; ...
         database(subj).ROI_between_15_38; database(subj).ROI_between_26_38];
     
 end
 
 rel_perc_cceps_between = nan(1,length(database));
 
 figure(), hold on
 for subj = 1:length(database)

     if ~isempty(database(subj).ROI_between_all)
         
         % scatterplot the age by the relative percentage of CCEPs found
         relative_cceps = (sum(~isnan(database(subj).ROI_between_all(:,3))) / length(database(subj).ROI_between_all(:,3)));
         scatter(age_vector(1,subj),relative_cceps,50,[0 0 0])
         
         % create array of relative cceps for group analysis
         rel_perc_cceps_between(subj) = relative_cceps;
%          % create variables and save to analyze as groups later
%          % save subject specific in struct.
%          
%          % add to current data in variable (because there could already
%          % be data from another ROI)
%          % save amount of CCEPs
%          database(subj).amount_cceps = database(subj).amount_cceps + sum(~isnan(ROI_plot_matrix(:,subj)));
%          % save total amount of stimulations
%          database(subj).total_stims = database(subj).total_stims + length(database(subj).ROI_stimpairs_data);
         
     elseif isempty(database(subj).ROI_between_all)
         
     % write NaN if there are no between ROI data to the array for group
     % analysis
         rel_perc_cceps_between(subj) = NaN;
     
     
     end
 end
 
xlim([0 50]), ylim([0 1])
xlabel('age subject')
ylabel('percentage CCEPs')
title('age-effect on percentage of CCEPs: between ROIs 15, 26, 38')
hold off

%% Boxplot percentage of CCEPs between ROIs in groups, 18+, 18- 

figure(), hold on

for gg = 1:length(age_vector)
    if age_vector(gg) < 18
        age_group(gg) = 1;
    elseif age_vector(gg) >= 18
        age_group(gg) = 2;
        
    end
end
boxplot(rel_perc_cceps_between,age_group,'Labels',{'18-','18+'})

xlabel('age subject')
ylabel('latency in ms')
title('relative percentage cceps, 18- and 18+, between ROI 15, 26 & 38')
hold off


%% Scatter both within and between ROIs together
figure(),hold on

for subj = 1:length(database)
    
    scatter(age_vector(1,subj),rel_perc_cceps_within(subj),50,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1])
    
    
end


for subj = 1:length(subjects)
    
    scatter(age_vector(1,subj),rel_perc_cceps_between(subj),50,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
    
end
xlim([0 50]), ylim([0 1])
xlabel('age subject')
ylabel('latency in ms')
title('relative percentage cceps, 18- and 18+, within and between all ROIs')
hold off
