%% Between: find ROIs with 2 or more electrode in one ROI and 1 or more in another ROI

% write to matrices with all between ROI options, but also to individual
% tracts 

% these are the ROIs
ROI_destrieux = [15, 26, 38];

% clean cells
for subj = 1:length(database)
    database(subj).ROI_stimpairs_data_between = [];
    database(subj).ROI_between_coverage_name = [];
    database(subj).ROI_between_coverage_num = [];
    database(subj).ROI_between_15_26 = [];
    database(subj).ROI_between_15_38 = [];
    database(subj).ROI_between_26_38 = [];
end

% set ses and runs to 1 - because we will not iterate over them
runs = 1;
ses = 1;

% for all subjects 
for subj = 1:length(database)
    
    % count how many 
    ROI_15_sum = sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux(1) + 1));
    ROI_26_sum = sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux(2) + 1));
    ROI_38_sum = sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux(3) + 1));
    ROI_between_total = ROI_15_sum + ROI_26_sum + ROI_38_sum;
    
    % clean cells/arrays and set counter to 1
    ROI_elec_name_all = {};
    ROI_elec_num_all = [];
    mat_row_counter = 1;
    
    % set counters for individual between ROI matrices to 1
    mat_row_counter_15_26 = 1;
    mat_row_counter_15_38 = 1;
    mat_row_counter_26_38 = 1;
    
    % if at least two of the three ROIs have an electrode on the ROI
    % and they are at least 3 electrodes in total
    if ((ROI_15_sum > 0 && ROI_26_sum > 0) ||  (ROI_15_sum > 0 && ROI_38_sum > 0) || (ROI_26_sum > 0 && ROI_38_sum > 0)) && ...
            (ROI_between_total >= 3)
        
        % print amount of electrodes
        print_num = sprintf('%.f', ROI_between_total);
        disp(['On subject: ' database(subj).subject ' between ROI analysis can be done. It has: ' print_num ' electrodes']);
        
        % next step: which electrodes are this within the subject
        ROI_elec_15_num = find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)== ROI_destrieux(1) + 1);
        ROI_elec_26_num = find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)== ROI_destrieux(2) + 1);
        ROI_elec_38_num = find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)== ROI_destrieux(3) + 1);
        
        % create cell matrix with the electrode numbers, for every region other
        % column
        ROI_elec_num_all(1:length(ROI_elec_15_num),1) = ROI_elec_15_num;
        ROI_elec_num_all(1:length(ROI_elec_26_num),2) = ROI_elec_26_num;
        ROI_elec_num_all(1:length(ROI_elec_38_num),3) = ROI_elec_38_num;
        
        % run over matrix to change 0 into NaN
        for aa = 1:size(ROI_elec_num_all,1)
            for bb = 1:size(ROI_elec_num_all,2)
                if ROI_elec_num_all(aa,bb) == 0
                    ROI_elec_num_all(aa,bb) = NaN;
                end
            end
        end
        
        % write nums of electrodes to the database struct
        database(subj).ROI_between_coverage_num = ROI_elec_num_all;
        
        % repeat process for names - might be useful later
        % find the names corresponding to these electrodes to save as well
        ROI_elec_15_name = database(subj).metadata(runs).atlases_electrodes.name(ROI_elec_15_num);
        ROI_elec_26_name = database(subj).metadata(runs).atlases_electrodes.name(ROI_elec_26_num);
        ROI_elec_38_name = database(subj).metadata(runs).atlases_electrodes.name(ROI_elec_38_num);
        
        % create cell matrix with the electrode names, for every region other
        % column
        ROI_elec_name_all(1:length(ROI_elec_15_name),1) = ROI_elec_15_name;
        ROI_elec_name_all(1:length(ROI_elec_26_name),2) = ROI_elec_26_name;
        ROI_elec_name_all(1:length(ROI_elec_38_name),3) = ROI_elec_38_name;
        
        % write names of electrodes to the database struct
        database(subj).ROI_between_coverage_name = ROI_elec_name_all;
        
        
        % iterate over the ROI columns within the matrix to find if the electrodes
        % are in a stimulation pair of which the other electrode is also in the
        % pair.
        
        % for every column
        for pp = 1:size(ROI_elec_num_all,2)
            
            % for every electrode in the column
            for qq = 1:size(ROI_elec_num_all,1)
                
                % find in which stimulation pair combination these electrodes
                % are
                ROI_stim_pair = find(database(subj).metadata(runs).stimulated_pairs(:,1) == (ROI_elec_num_all(qq,pp)));
                
                % if the other electrode of the stimulated pair is also on the
                % ROI, and they are an actual combination of stimulations (this is
                % only possible when the second electrodes is one higher than the first)
                % then this electrode pair is part of analysis
                if sum(ismember(ROI_elec_num_all(:,pp),database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,2))) >= 1 ...
                        && abs(database(subj).metadata(runs).stimulated_pairs(ROI_stim_pair,2) - ROI_elec_num_all(qq,pp)) == 1
                    
                    % iterate over all other electrodes in other ROI column to find
                    % the latency and amplitude in those electrodes
                    
                    % for every column
                    for zz = 1:size(ROI_elec_num_all,2)
                        
                        % if the column is not the column of the stim_pair
                        if zz ~= pp
                            
                            
                            %for every measured electrode in that column
                            for yy = 1:size(ROI_elec_num_all(:,zz),1)
                            
                                % to ensure this combination of electrodes
                                % stimulation actually happened
                                if ROI_elec_num_all(yy,zz) <= size(database(subj).all_spesclin_latency,2)
                                    % save stim_pair under ROI_stimparis_data
                                    
                                    % save stim_pair in first column
                                    database(subj).ROI_stimpairs_data_between(mat_row_counter,1) = ROI_stim_pair;
                                    % save measured elec in second column
                                    database(subj).ROI_stimpairs_data_between(mat_row_counter,2) = ROI_elec_num_all(yy,zz);
                                    % save latency in third column
                                    database(subj).ROI_stimpairs_data_between(mat_row_counter,3) = database(subj).all_spesclin_latency(ROI_stim_pair,ROI_elec_num_all(yy,zz));
                                    % save amplitude in fourth column
                                    database(subj).ROI_stimpairs_data_between(mat_row_counter,4) = database(subj).all_spesclin_amplitude(ROI_stim_pair,ROI_elec_num_all(yy,zz));
                         
                                    % It the step before, all latencies are
                                    % put into same matrix. Here put the
                                    % data in matrices for specific ROIs
                                    
                                    % ROI 15 & 26
                                    if (pp == 1 && zz == 2) || (pp == 2 && zz == 1)
                                        % save stim_pair in first column
                                        database(subj).ROI_between_15_26(mat_row_counter_15_26,1) = ROI_stim_pair;
                                        % save measured elec in second column
                                        database(subj).ROI_between_15_26(mat_row_counter_15_26,2) = ROI_elec_num_all(yy,zz);
                                        % save latency in third column
                                        database(subj).ROI_between_15_26(mat_row_counter_15_26,3) = database(subj).all_spesclin_latency(ROI_stim_pair,ROI_elec_num_all(yy,zz));
                                        % save amplitude in fourth column
                                        database(subj).ROI_between_15_26(mat_row_counter_15_26,4) = database(subj).all_spesclin_amplitude(ROI_stim_pair,ROI_elec_num_all(yy,zz));
                                        
                                        % add 1 to individual counter
                                        mat_row_counter_15_26 = mat_row_counter_15_26 + 1;
                                        
                                    % ROI 15 & 38    
                                    elseif (pp == 1 && zz == 3) || (pp == 3 && zz == 1)
                                        % save stim_pair in first column
                                        database(subj).ROI_between_15_38(mat_row_counter_15_38,1) = ROI_stim_pair;
                                        % save measured elec in second column
                                        database(subj).ROI_between_15_38(mat_row_counter_15_38,2) = ROI_elec_num_all(yy,zz);
                                        % save latency in third column
                                        database(subj).ROI_between_15_38(mat_row_counter_15_38,3) = database(subj).all_spesclin_latency(ROI_stim_pair,ROI_elec_num_all(yy,zz));
                                        % save amplitude in fourth column
                                        database(subj).ROI_between_15_38(mat_row_counter_15_38,4) = database(subj).all_spesclin_amplitude(ROI_stim_pair,ROI_elec_num_all(yy,zz));    
                                       
                                        % add 1 to individual counter
                                        mat_row_counter_15_38 = mat_row_counter_15_38 + 1;
                                        
                                    % ROI 26 & 38    
                                    elseif (pp == 2 && zz == 3) || (pp == 3 && zz == 2)
                                        % save stim_pair in first column
                                        database(subj).ROI_between_26_38(mat_row_counter_26_38,1) = ROI_stim_pair;
                                        % save measured elec in second column
                                        database(subj).ROI_between_26_38(mat_row_counter_26_38,2) = ROI_elec_num_all(yy,zz);
                                        % save latency in third column
                                        database(subj).ROI_between_26_38(mat_row_counter_26_38,3) = database(subj).all_spesclin_latency(ROI_stim_pair,ROI_elec_num_all(yy,zz));
                                        % save amplitude in fourth column
                                        database(subj).ROI_between_26_38(mat_row_counter_26_38,4) = database(subj).all_spesclin_amplitude(ROI_stim_pair,ROI_elec_num_all(yy,zz));
                                        
                                        % add 1 to individual counter 
                                        mat_row_counter_26_38 = mat_row_counter_26_38 + 1;
                                    end   
                                    
                                    % add one to counter to ensure no overwriting
                                    mat_row_counter =  mat_row_counter + 1;
                               end 
                            end                  
                        end
                    end
                end
            end
        end
    end
end



%% Create matrice of data of individuals for analysis

ROI_between_plot_matrix = nan(300,length(subjects));
ROI_between_plot_15_26 = nan(300,length(subjects));
ROI_between_plot_15_38 = nan(300,length(subjects));
ROI_between_plot_26_38 = nan(300,length(subjects));

for subj = 1:length(database)
    
    % if there are data
    if ~isempty(database(subj).ROI_stimpairs_data_between)
        % these matrices can be used for analysis: latencies * subjects
        ROI_between_plot_matrix(1:length(database(subj).ROI_stimpairs_data_between),subj) = database(subj).ROI_stimpairs_data_between(:,3);
    end
    
    % if there are data
    if ~isempty(database(subj).ROI_between_15_26)
        
        ROI_between_plot_15_26(1:length(database(subj).ROI_between_15_26),subj) = database(subj).ROI_between_15_26(:,3);
    end
    
    % if there are data
    if ~isempty(database(subj).ROI_between_15_38)
        
        ROI_between_plot_15_38(1:length(database(subj).ROI_between_15_38),subj) = database(subj).ROI_between_15_38(:,3);       
    end
    
    % if there are data
    if ~isempty(database(subj).ROI_between_26_38)
        
        ROI_between_plot_26_38(1:length(database(subj).ROI_between_26_38),subj) = database(subj).ROI_between_26_38(:,3);
    end
    
    
end

%% Scatterplot of the latencies between ROIs

% list of ROI_between matrices that can be used for plotting:
% - ROI_between_plot_15_26
% - ROI_between_plot_15_38
% - ROI_between_plot_26_38
% - ROI_between_plot_matrix

ROI_between_plot =  ROI_between_plot_matrix;

figure(), hold on
% iterate over all subjects in database
for subj = 1:length(database)
    
    % create age-vector for plotting
    age_vector(1,subj) = database(subj).age_ses1;
    
    
    % scatterplot
    scatter(repelem(age_vector(1,subj),length(ROI_between_plot)),ROI_between_plot(:,subj))

end

xlim([0 50]),ylim([10 90])
xlabel('age subject')
ylabel('latency in ms')
title('Latency Between Destrieux ROIs: 15, 26 & 38')
hold off

%% Boxplots for between ROI age groups

% list of ROI_between matrices that can be used for plotting:
% - ROI_between_plot_15_26
% - ROI_between_plot_15_38
% - ROI_between_plot_26_38
% - ROI_between_plot_matrix

ROI_between_plot =  ROI_between_plot_matrix;

figure(), hold on

for gg = 1:length(age_vector)
    if age_vector(gg) <= 18
        age_group(gg) = 1;
    elseif age_vector(gg) > 18
        age_group(gg) = 2;
        
    end
end
   
boxplot(ROI_between_plot,age_group)

xlabel('age subject')
ylabel('latency in ms')
title('age-effect on latency two groups, 18- and 18+')
hold off


summed_cceps = (sum(~isnan(ROI_between_plot)));

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

