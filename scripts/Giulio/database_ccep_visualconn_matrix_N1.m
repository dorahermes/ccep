%This script created and plotd a matrix in which possible
%connections are shown. Afterwards another matrix shows the areas that show
%an early response after a stimulation in the source areas. 
%The first part of the script is useful to integrate Benson and
%Wang & Kastner visual atlases.
% D. Hermes and G. Castegnaro, 2019, UMC Utrecht

%for subjects 0315, 0405, 0306 we need to change the class of the data from
%cell array to double 
for subj = [1 4 5 6];  
    for runs = 1:length(database(subj).metadata);
        if iscell(database(subj).metadata(runs).electrodes.Benson_label); 
            database(subj).metadata(runs).electrodes.Benson_label = str2double(database(subj).metadata(runs).electrodes.Benson_label);
        end 
        if iscell(database(subj).metadata(runs).electrodes.Benson_polarangle); 
           database(subj).metadata(runs).electrodes.Benson_polarangle = str2double(database(subj).metadata(runs).electrodes.Benson_polarangle);
        end
        if iscell(database(subj).metadata(runs).electrodes.Wang_label); 
           database(subj).metadata(runs).electrodes.Wang_label= str2double(database(subj).metadata(runs).electrodes.Wang_label);
        end
    end 
end

runs = 1; %[1:length(database(subj).metadata)] %(this is partially right, because it loops around the rows of the metadata, not specifically on the runs columns

% select subjects/iterate over all subjects in database

for subj = 1
    % we need to make a matrix of size Atlas areas X Atlas areas
    % integrate Wang and Benson maps to gte th best of both
    new_label = NaN(size(database(subj).metadata(runs).electrodes.name,1),1);
        for kk = 1:length(new_label)
            % is an electrode Benson V1-2-3?
            if ismember(database(subj).metadata(runs).electrodes.Benson_label(kk),[1 2 3]) && database(subj).metadata(runs).electrodes.Wang_label(kk) ~= 0  ;
                % check ventral/dorsal?
                %if database(subj).metadata(runs).electrodes.Benson_polarangle(kk) < 90 % upper hemifield -->ventral electrode
                    %new_label(kk) = database(subj).metadata(runs).electrodes.Benson_label(kk)*2-1;
                %elseif database(subj).metadata(runs).electrodes.Benson_polarangle(kk)>=90
                    %new_label(kk) = database(subj).metadata(runs).electrodes.Benson_label(kk)*2;
                %polar angle is not added to the electrode file, so we have
                %to find another way, for now we simply take the wang areas 
                new_label(kk) = database(subj).metadata(runs).electrodes.Wang_label(kk);
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 4 %hV4
                new_label(kk) = 7;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 5 %VO1
                new_label(kk) = 8;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 6 %Vo2
                new_label(kk) = 9;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 7 %LO1
                new_label(kk) = 15;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 8 %LO2
                new_label(kk) = 14;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 9 %TO1
                new_label(kk) = 13;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 10 %TO2
                new_label(kk) = 12;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 11 %V3B
                new_label(kk) = 16;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 12 %V3A
                new_label(kk) = 17;
            elseif database(subj).metadata(runs).electrodes.Benson_label(kk)== 0 && ismember(database(subj).metadata(runs).electrodes.Wang_label(kk),[1:25])  % areas: IPS/PHC/SPL/FEF
                new_label(kk) = database(subj).metadata(runs).electrodes.Wang_label(kk);
            end
        end

    % we keep the Wang ROI names, but add larger area estimates from Benson
    Wang_ROI_Names = {...
        'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
        'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
        'IPS5' 'SPL1' 'FEF'};
    
    %% plot how many electrodes are on every area 
    elec_on_area = (1:length(Wang_ROI_Names))';
    elec_on_area(:,2) = histc(new_label(:),elec_on_area);
    elec_on_area_array = num2cell(elec_on_area(:,1:2)); 
    elec_on_area_array(1:25,1) = Wang_ROI_Names;
    database(subj).metadata(runs).areas_elec = elec_on_area_array; % adding the array to the database


    %Highlighting the VOF (hV4-VO1 ventrally, V3A/B dorsally) 
    Wang_ROI_Names{7} = ['\color{red}' Wang_ROI_Names{7}];
    Wang_ROI_Names{8} = ['\color{red}' Wang_ROI_Names{8}];
    Wang_ROI_Names{17} = ['\color{green}' Wang_ROI_Names{17}];
    Wang_ROI_Names{16} = ['\color{green}' Wang_ROI_Names{16}];
    %% make a matrix with the connections electrodes can measure 


    % 1) we want to have the labels of the measured channels (each
    % corresponds to an electrode)
    channel_labels = zeros(size(database(subj).metadata(runs).channels.name,1),1);
    for kk = 1:size(database(subj).metadata(runs).channels.name,1)
        % to which electrode does this channel correspond
        [a,b] = ismember(database(subj).metadata(runs).channels.name{kk},database(subj).metadata(runs).electrodes.name);
        if a==1 % this channel is an electrode
            channel_labels(kk) = new_label(b);
        end
    end

    % 2) now we walk through stimulated pairs (in database(subj).metadata(runs).stimulated_pairs) and add in
    % measured connection matrix
    measured_connection_matrix = zeros(length(Wang_ROI_Names),length(Wang_ROI_Names));

    for kk = 1:size(database(subj).metadata(runs).stimulated_pairs,1)

        % which target channels are measured?
        measured_channel_labels = channel_labels;
        measured_channel_labels(database(subj).metadata(runs).stimulated_pairs(kk,:)) = []; % remove the stimulated channels

        % get the label of the first stimulated electrode
        stim_el1_label = channel_labels(database(subj).metadata(runs).stimulated_pairs(kk,1));

        if stim_el1_label>0 % this electrode has a label
            [nn,~] = hist(measured_channel_labels(measured_channel_labels>0),[1:25]);
            % add the number of measured channels
            measured_connection_matrix(stim_el1_label,:) = ...
                measured_connection_matrix(stim_el1_label,:) + nn;
        end

        % get the label of the second stimulated electrode
        stim_el2_label = channel_labels(database(subj).metadata(runs).stimulated_pairs(kk,2));

        if stim_el2_label>0 % this electrode has a label
            [nn,~] = hist(measured_channel_labels(measured_channel_labels>0),[1:25]);
            % add the number of measured channels
            measured_connection_matrix(stim_el2_label,:) = ...
                 measured_connection_matrix(stim_el2_label,:) + nn;
        end
    end

    %% Plot the possible connections matrix 
    database(subj).metadata(runs).possible_mat = measured_connection_matrix;
    figure('Position',[0 0 800 800])
    imagesc(measured_connection_matrix,[0 max(measured_connection_matrix(:))])
    axis square
    xlabel('Target Visual Area'),ylabel('Source Visual Area')
    set(gca,'XTick',[1:25],'YTick',[1:25],'FontName','Lato','FontSize',12,...
        'XDir','normal', 'XTickLabel', Wang_ROI_Names, 'YTickLabel', Wang_ROI_Names);
    xtickangle(90)
    cm = jet(max(measured_connection_matrix(:)));
    cm = [1 1 1; cm];
    colormap(cm);
    hcb = colorbar;
    title(hcb,'# of stimulations per area')
    set(gcf,'PaperPositionMode','auto')
    % print('-dpng','-r300','/home/giulio/figures/measured_connections_matrix')
    title(['POSSIBLE MATRIX Subject: ' database(subj).subject ' - Run: ' database(subj).metadata(runs).run ])
    %% create a matrix starting from the ccep data (ccep_detect_onlyN1)
    %on the y axis all the stimulation pairs, on the x axis all the target
    %areas (electrodes). If CCEP is significant then 1, if not then 0. 

    %check which areas were stimulated % e.g. V1d avd V2d
    %check which area had an early response % e.g. V3d

    %contains a 1 for the electrode and the stimulated pair that produce an ER
    conn_plot = database(subj).metadata(runs).n1_peak_amplitude < 0; % negative N1 amplitude

    % remove measured N1 in stimulated electrodes (no N1 is possible) 
    %only  removes the first one 
    for kk = 1:size(database(subj).metadata(runs).stimulated_pairs,1)
        for ss = 1:size(database(subj).metadata(runs).channels.name,1) 
            if database(subj).metadata(runs).stimulated_pairs(kk,1) == ss'
                conn_plot(ss,kk) = 0; % set stimulated channels to zero
            end
        end 
    end 

    % necessary to assign a label to every stimulation pair and match it with
    % the channels that do have a label (channel_labels)
    labeled_stimsets = zeros(size(database(subj).metadata(runs).stimulated_pairs));
    for kk = 1:size(database(subj).metadata(runs).stimulated_pairs,1)
        labeled_stimsets(kk,1:2) = channel_labels(database(subj).metadata(runs).stimulated_pairs(kk,1:2));
    end

    % we have: 
    % output_ER:        channels X Nstimpairs X amplitude & latency
    % conn_plot:        1 if significant output_ER, channels X Nstimpairs 
    % database(subj).metadata(runs).stimulated_pairs:      stimulated channels, stimpairs X 2
    % labeled_stimsets: labels of stimulated pairs, stimpairs X 2
    % channel_labels:   labels of the channels, channesl x 1

    % adding a 1 in the matrix when stimulated electrode produce ER in other
    % electrodes 
    % Matrix coming from data, 1 if there is a ccep, 0 if no ccep is detected 
    measured_N1_matrix = zeros(25,25);
    for ss = 1:size(conn_plot,2) % stimulated channels/stimpairs
        for ee = 1:size(conn_plot,1) % measured channels
            if labeled_stimsets(ss,1) ~= 0 && ~isnan(labeled_stimsets(ss,1)) % does the first stimulated electrode have a label
                % then add 1s for measured CCEPs
                    % does the measured channel have a label?
                    if channel_labels(ee) ~= 0 && ~isnan(channel_labels(ee)) 
                        measured_N1_matrix(channel_labels(ee),labeled_stimsets(ss,1)) = ...
                        conn_plot(ee,ss)+ measured_N1_matrix(channel_labels(ee),labeled_stimsets(ss,1));
                    end
            end
        end 

            if labeled_stimsets(ss,2) ~= 0 && ~isnan(labeled_stimsets(ss,2)) % does the second stimulated electrode have a label
                % then add 1s for measured CCEPs
                % does the measured channel have a label?
                    if channel_labels(ee) ~= 0 && ~isnan(channel_labels(ee))
                    measured_N1_matrix(channel_labels(ee),labeled_stimsets(ss,2)) = ...
                        conn_plot(ee,ss)+ measured_N1_matrix(channel_labels(ee),labeled_stimsets(ss,2));
                    end
            end
    end


    % adding a -1 on areas that are not covered
    measured_N1_matrix(database(subj).metadata(runs).possible_mat==0) = -1;

    %measured_N1_matrix(isnan(measured_N1_matrix)) = -1 %making nans become -1 (Non covered areas) 

    %     %Blocking the script if the absolute matrix does not show any ccep
    %     if isempty(measured_N1_matrix(measured_N1_matrix > 0))
    %         writ = ['The absolute matrix for subject ' database(subj).subject '-Run' database(subj).metadata(runs).run ' run ' num2str(runs) ' does not show any activation for none of the selected areas.'];
    %         disp(writ)
    %         break 
    %     end 



    %% plot absolute ccep connections matrix 
    database(subj).metadata(runs).absolute_mat = measured_N1_matrix; % adding the matrix to the database

    figure('Position',[0 0 800 800])
    imagesc(measured_N1_matrix,[-1 max(measured_N1_matrix(:))])
    axis square
    xlabel('Recorded Visual Area'),ylabel('Stimulated Visual Area')
    set(gca,'XTick',[1:25],'YTick',[1:25],'FontName','Lato','FontSize',12,...
        'XDir','normal', 'XTickLabel', Wang_ROI_Names, 'YTickLabel', Wang_ROI_Names);
    xtickangle(90)

    cm = parula(max(measured_N1_matrix(:)));
    cm = [0 0 0; cm]; % connection showing 0 ccep white
    cm = [1 1 1; cm]; % make area not stimulated as white 

    colormap(cm)
    hcb = colorbar;%('TicksMode','manual', 'Ticks',[-1 0 1],'TickLabels', {'not recorded', 'no response', 'Early response(CCEP)'});
    title(hcb,'# CCEP')
    set(gcf,'PaperPositionMode','auto')
    % print('-dpng','-r300','/home/giulio/figures/measured_N1_matrix')
    title(['ABSOLUTE MATRIX Subject: ' database(subj).subject ' - Run: ' database(subj).metadata(runs).run ])
    database(subj).metadata(runs).absolute_mat = measured_N1_matrix; % adding the matrix to the database

    %% plot relative ccep connection matrix 
    % we could have measured all connections in measured_connection_matrix
    % which percentage of these potential connections did we find?
    relative_N1_matrix = measured_N1_matrix./measured_connection_matrix;
    relative_N1_matrix(database(subj).metadata(runs).absolute_mat==-1) = -1; % making nans become -1 (Non covered areas)

    figure('Position',[0 0 800 800])
    imagesc(relative_N1_matrix,[-1 1])
    axis square

    xlabel('Recording Visual Area'),ylabel('Stimulated Visual Area')
    set(gca,'XTick',[1:25],'YTick',[1:25],'FontName','Lato','FontSize',12,...
        'XDir','normal', 'XTickLabel', Wang_ROI_Names, 'YTickLabel', Wang_ROI_Names);
    xtickangle(90)
    cm = summer(21);
    A = ones(10,3)
    cm(1:10,:) = A  % making -1 black 
    cm(11,:) = [0 0 0] %making 0 white (not showing any ccep) 
    colormap(cm)
    hcb = colorbar('TicksMode','manual', 'Ticks',[-1 0 1],'TickLabels', {'not recorded', 'no response', 'Early response(CCEP)'});
    %hcb.TickLabels = arrayfun( @(x) [num2str(x) '%'], hcb.Ticks * 100, 'UniformOutput', false );
    set(gcf,'PaperPositionMode','auto')
    %print('-dpng','-r300','/home/giulio/figures/relative_N1_matrix')
    %   saveas(relative_N1_matrix,sprintf('relative_matrix %d.png',subj)); % will create FIG1, FIG2,...
    %   hold on 
    title(['RELATIVE MATRIX Subject: ' database(subj).subject ' - Run: ' database(subj).metadata(runs).run '(' num2str(runs) ')' ])
    %add the matrix to the database for plotting together 
    database(subj).metadata(runs).relative_mat = relative_N1_matrix; % adding the matrix to the database

end 
    
%% Assessing the reciprocity of the network
% Reciprocity index = reciprocal connections/all connections 


for subj = 1%:length(database)
    database(subj).metadata(runs).all_conn = [];
    database(subj).metadata(runs).reciprocal_conn = [];
    for xx = 1:size(database(subj).metadata(runs).relative_mat,1) 
        for yy = 1:size(database(subj).metadata(runs).relative_mat,2)
            if database(subj).metadata(runs).relative_mat(xx,yy) > 0.01 
                database(subj).metadata(runs).all_conn(end+1,1:2) = ([xx yy]);
            end
        end
    end

%creating an array with reciprocal conections
    for xx = 1:size(database(subj).metadata(runs).relative_mat,1) 
        for yy = 1:size(database(subj).metadata(runs).relative_mat,2) 
            if database(subj).metadata(runs).relative_mat(xx,yy)>0 && database(subj).metadata(runs).relative_mat(yy,xx) && xx ~= yy 
                database(subj).metadata(runs).reciprocal_conn(end+1,1:2) = ([xx yy]);
            end 
        end 
    end 
database(subj).metadata(runs).reciprocity_index = size(database(subj).metadata(runs).reciprocal_conn,1)./size(database(subj).metadata(runs).all_conn,1);
end 

all_reciprocities = ['RESP0315'; 'RESP0751'; 'RESP0401'; 'RESP0405'; 'RESP0306'; {'RESP0703'}];
for subj = 1%:length(database)
    all_reciprocities(subj,2) = num2cell(database(subj).metadata(runs).reciprocity_index); 
end 

