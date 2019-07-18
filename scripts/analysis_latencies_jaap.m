%% preprocessing before analysis - removing irrelevant channels

subjects = {'RESP0621', 'RESP0706', 'RESP0733', 'RESP0768', 'RESP0294', 'RESP0306', ...
    'RESP0315', 'RESP0348', 'RESP0368', 'RESP0369', 'RESP0400', 'RESP0401', 'RESP0405', ...
    'RESP0435', 'RESP0449', 'RESP0450', 'RESP0458', 'RESP0467', 'RESP0468', 'RESP0477', ...
    'RESP0478', 'RESP0501', 'RESP0502', 'RESP0703', 'RESP0724'};

% RESP0621 
database(1).c_grid = (1:64);
database(1).total_grid = (1:64); 
database(1).age_ses1 = 18;

% RESP0706 
database(2).c_grid = [];
database(2).total_grid = (1:48);
database(2).age_ses1 = 50;

% RESP0733 
database(3).c_grid = (1:48);
database(3).total_grid = (1:64);
database(3).age_ses1 = 8;

% RESP0768 
database(4).c_grid = (1:32);
database(4).total_grid = (1:64);
database(4).age_ses1 = 6;

% RESP0294
database(5).c_grid = [];
database(5).total_grid = (1:64);
database(5).age_ses1 = 35;

% RESP0306 - MKR @ 65
database(6).c_grid = [];
database(6).total_grid = (1:97);
database(6).age_ses1 = 17;

% RESP0315 - MKR @ 65
database(7).c_grid = []; 
database(7).total_grid = (1:89);
database(7).age_ses1 = 10;

% RESP0348 - MKR @ 65
database(8).c_grid = []; 
database(8).total_grid = (1:81);
database(8).age_ses1 = 15;

% RESP0368 - MKR @ 65
database(9).c_grid = []; 
database(9).total_grid = (1:81);
database(9).age_ses1 = 14;

% RESP0369 
database(10).c_grid = (33:48); 
database(10).total_grid = (1:64);
database(10).age_ses1 = 45;

% RESP0400 - MKR @ 65
database(11).c_grid = []; 
database(11).total_grid = (1:73);
database(11).age_ses1 = 10;

% RESP0401
database(12).c_grid = []; 
database(12).total_grid = (1:64);
database(12).age_ses1 = 15;

% RESP0405 - MKR @ 65
database(13).c_grid = []; 
database(13).total_grid = (1:89);
database(13).age_ses1 = 42;

% RESP0435 
database(14).c_grid = []; 
database(14).total_grid = (1:64);
database(14).age_ses1 = 15;

% RESP0449 - MKR @ 65
database(15).c_grid = []; 
database(15).total_grid = (1:89);
database(15).age_ses1 = 11;

% RESP0450 
database(16).c_grid = []; 
database(16).total_grid = (1:56);
database(16).age_ses1 = 12;

% RESP0458 - MKR @ 65
database(17).c_grid = (66:73); 
database(17).total_grid = (1:81);
database(17).age_ses1 = 9;

% RESP0467 - MKR @ 65
database(18).c_grid = []; 
database(18).total_grid = (1:81);
database(18).age_ses1 = 17;

% RESP0468 - MKR @ 65
database(19).c_grid = []; 
database(19).total_grid = (1:81);
database(19).age_ses1 = 16;

% RESP0477 
database(20).c_grid = []; 
database(20).total_grid = (1:64);
database(20).age_ses1 = 11;

% RESP0478
database(21).c_grid = (1:32); 
database(21).total_grid = (1:56);
database(21).age_ses1 = 13;

% RESP0501 - MKR 65
database(22).c_grid = []; 
database(22).total_grid = (1:89);
database(22).age_ses1 = 14;

% RESP0502
database(23).c_grid = []; 
database(23).total_grid = (1:64);
database(23).age_ses1 = 41;

% RESP0703 - MKR 65
database(24).c_grid = []; 
database(24).total_grid = (1:89);
database(24).age_ses1 = 15;

% RESP0724 - MKR 65
database(25).c_grid = []; 
database(25).total_grid = (1:89);
database(25).age_ses1 = 15;

%% averaged latency per run per subj

% iterate over all subjects in database
for subj = 1:length(database)
    % iterate over all their runs
    for runs = 1:length(database(subj).metadata)
        
        database(subj).metadata(runs).avg_latency = nanmean(database(subj).metadata(runs).n1_peak_sample(:));
        
    end
end



%% averaged latency and amplitudes per subj - SPESclin only - with skipping marker

% iterate over all subjects in database
for subj = 1:length(database)
    
    % create new clean matrix per subj, that will be filled with the
    % results of all different SPESclin runs
    full_latency_matrix = [];
    full_amplitude_matrix = [];
    
    % iterate over all their runs
    for runs = 1:length(database(subj).metadata)
        if strcmp(database(subj).metadata(runs).task, 'SPESclin')
            
            % if total grid has even length, there is no marker at 65
            if mod(length(database(subj).total_grid),2) == 0   
                
                full_latency_matrix = [full_latency_matrix, ...
                    database(subj).metadata(runs).n1_peak_sample(1:(length(database(subj).total_grid)),:)];
                full_amplitude_matrix = [full_amplitude_matrix, ...
                    database(subj).metadata(runs).n1_peak_amplitude(1:(length(database(subj).total_grid)),:)];
                
            % if total grid has uneven length, there is a marker at 65
            elseif mod(length(database(subj).total_grid),2) == 1
                
                full_latency_matrix = [full_latency_matrix, ...
                    database(subj).metadata(runs).n1_peak_sample([1:64,66:(length(database(subj).total_grid))],:)];
                full_amplitude_matrix = [full_amplitude_matrix, ...
                    database(subj).metadata(runs).n1_peak_amplitude([1:64,66:(length(database(subj).total_grid))],:)];
                
            end
            
        end
        
        % save full latency and amplitude matrices of all SPESclin
        % first recalculate the samples to ms for latency
        latency_ms_matrix = (full_latency_matrix - 5120) / database(subj).metadata(runs).data_hdr.Fs * 1000;
        database(subj).all_spesclin_latency = latency_ms_matrix;
        database(subj).all_spesclin_amplitude = full_amplitude_matrix;
        
        % calculate averages and save in structure
        database(subj).avg_latency = nanmean(latency_ms_matrix(:));
        database(subj).avg_amplitude = nanmean(full_amplitude_matrix(:));
        
    end
end


% might want to check when the same pair is stimulated in a new run. What
% happens then to data? If effecting: first put data together and then
% preprocess data might be a solution. 

% some check for the latencies because sometimes they all have same number

% warning: currently database only consists of grids and strips but does
% have the stimulations of the depth electrodes

%% age histogram


% iterate over all subjects in database
for subj = 1:length(database)
    
    age_vector_hist(1,subj) = database(subj).age_ses1;
end

figure()
hist(age_vector_hist,50)
hold on
xlim([0 50]), set(gca,'XTick',[0:50])


%% Age grouped latency

for gg = 1:length(age_vector)
    if age_vector(gg) <= 18
        age_group(gg) = 1;
    elseif age_vector(gg) > 18
        age_group(gg) = 2;

    end
end

figure(2) 
boxplot(matrix_reshape_all,age_group)
hold on
xlabel('age subject')
ylabel('latency in ms')
title('age-effect on latency two groups, 18- and 18+')
hold off



%% boxplot for all electrodes and stims

matrix_reshape_all = NaN(12000,25);
% iterate over all subjects in database
for subj = 1:length(database)
    
    matrix_reshape_all(1:(size(database(subj).all_spesclin_latency,1)*size(database(subj).all_spesclin_latency,2)),subj) = reshape(database(subj).all_spesclin_latency,...
        [1,(size(database(subj).all_spesclin_latency,1)*size(database(subj).all_spesclin_latency,2))]);
    
   age_vector(1,subj) = database(subj).age_ses1;
    
end

figure(3)
boxplot(matrix_reshape_all,age_vector)
hold on
xlabel('age subject')
ylabel('latency in ms')
title('age-effect on latency in all electrodes')
hold off


%% percentage of CCEPs in all electrodes

figure(5), hold on 
% iterate over all subjects in database
for subj = 1:length(database)
    scatter(age_vector(subj),(sum(~isnan(matrix_reshape_all(:,subj)))/ ...
        (size(database(subj).all_spesclin_latency,1)*size(database(subj).all_spesclin_latency,2))),50,[0 0 0])
     
end
xlabel('age subject')
ylabel('percentage CCEPs')
title('age-effect on percentage of CCEPs in all electrodes')
hold off




%% layout for all subjects and only SPESclin runs 

% iterate over all subjects in database
for subj = 1:length(database)
    % iterate over all their runs
    for runs = 1:length(database(subj).metadata)
        if strcmp(database(subj).metadata(runs).task, 'SPESclin')
            
            
        end
    end
end

%% Set data for rendering

addpath('/Fridge/users/jaap/github/ecogBasicCode/render/')
top_path = '/Fridge/users/jaap/ccep/dataBIDS/';
subjects = {'RESP0467'};
hemi_cap = {'R'};

ss = 18; 

% pick a viewing angle:
v_dirs = [90 0]; %;90 0;90 -60;270 -60;0 0];

% set stimulated pair you want to render
stim_pair = 40;

% select significant peaks in the other channels
n1_plot_sample = database(ss).metadata(1).n1_peak_sample(:,stim_pair);
n1_plot_amplitude =  database(ss).metadata(1).n1_peak_amplitude(:,stim_pair);


%% Loop for rendering

for s = 1%:3 %1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(top_path,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
    % electrode locations name:
    dataLocName = dir(fullfile(top_path,['sub-' subj],'ses-1','ieeg',...
        ['sub-' subj '_ses-1_electrodes.tsv']));
    dataLocName = fullfile(dataLocName(1).folder,dataLocName(1).name);
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
        
        % ecog_Label(els,30,12) % add electrode positions
        % add electrodes
        ccep_el_add(els,[0.1 0.1 0.1],20)
        % give stimulated electrodes other color
        ccep_el_add(els(database(ss).metadata(1).stimulated_pairs(stim_pair,1:2),:),[0 0 0],40)
        % set size and color of channels with significant peak 
        % based on sample (from stimulation on, so -5120) and the amplitude
        % color = latency, size = amplitude
        ccep_el_add_size_and_color(els,n1_plot_amplitude,(n1_plot_sample-5120),500,200)

        set(gcf,'PaperPositionMode','auto')
        % print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
        % ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end

%% EVERYTHING FOR CONNECTIVITY MATRICES 

%% Add electrode atlases (without the MKR1+! and unneccesary channels) to structure

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
                % to the brain regions

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

title('Connectivity matrix of possible connections in database (25 subj)')
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
title('Connectivity matrix of possible connections in database (25 subj)')
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
rec_stim_regions = [rec_region stim_region];

%% find ROIs with > 1000 epochs and more then 1 electrode to measure within ROI

% find in which patients those areas are represented
% iterate over all subjects in database
% in atlas ROI is 1 less then in the labels, so add 1 in code 
ROI_destrieux = 4; 

for subj = 1:length(database)

   
    if sum(ismember(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label),ROI_destrieux + 1)) > 1
        disp(['subject: ' database(subj).subject ' has coverage: ' ]); database(subj).age_ses1; 
        
        
        % next step: which electrodes are this within the subject
        ROI_elec_loc = find(str2double(database(subj).metadata(runs).atlases_electrodes.Destrieux_label)== ROI_destrieux + 1);
        ROI_elec_name = database(subj).metadata(runs).atlases_electrodes.name(ROI_elec_loc);
        
        % write names of electrodes to the database struct
        database(subj).ROI_within_coverage = ROI_elec_name;
        
        
    end
    
  
end

%% analyse latency ROIs within 

% for every patient that has >= 2 coverage of ROI

% for every elec in the cell 

% find latency in every ROI 
%% 

tic;
save('dbstruct160719','database', '-v7.3')
time3 = toc;
