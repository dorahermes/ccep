%% 
%This script has been created to create and plot a matrix in which possible
%connections are shown. Afterwards another matrix shows the areas that show
%an early response after a stimulation in the source areas. 
%The first part of the script is useful to integrate Benson and
%Wang&Kastner visual atlases.
% D. Hermes and G. Castegnaro, 2019, UMC Utrecht

%% we need to make a matrix of size Atlas areas X Atlas areas

% integrate Wang and Benson maps to gte th best of both
new_label = NaN(size(electrodes_table.name,1),1);
for kk = 1:length(new_label)
    % is an electrode Benson V1-2-3?
    if ismember(electrodes_table.Benson_label(kk),[1 2 3])
        % check ventral/dorsal?
        if electrodes_table.Benson_polarangle(kk)<90 % upper hemifield -->ventral electrode
            new_label(kk) = electrodes_table.Benson_label(kk)*2-1;
        elseif electrodes_table.Benson_polarangle(kk)>=90
            new_label(kk) = electrodes_table.Benson_label(kk)*2;
        end
    elseif electrodes_table.Benson_label(kk)==4
        new_label(kk) = 7;
    elseif electrodes_table.Benson_label(kk)==5
        new_label(kk) = 8;
    elseif electrodes_table.Benson_label(kk)==6
        new_label(kk) = 9;
    elseif electrodes_table.Benson_label(kk)==7 %LO1
        new_label(kk) = 15;
    elseif electrodes_table.Benson_label(kk)==8 %LO2
        new_label(kk) = 14;
    elseif electrodes_table.Benson_label(kk)==9 %TO1
        new_label(kk) = 13;
    elseif electrodes_table.Benson_label(kk)==10 %TO2
        new_label(kk) = 12;
    elseif electrodes_table.Benson_label(kk)==11 %V3B
        new_label(kk) = 16;
    elseif electrodes_table.Benson_label(kk)==12 %V3A
        new_label(kk) = 17;
    elseif electrodes_table.Benson_label(kk)==0 && electrodes_table.Wang_label(kk) % areas: IPS/PHC/SPL/FEF
        new_label(kk) = electrodes_table.Wang_label(kk)+100;
    end
end

% we keep the Wang ROI names, but add larger area estimates from Benson
Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};

%% make a matrix with the connections electrodes can measure 
% we put 1 when electrodes cover source and target regions (just based on
% atlases and electrode locations)

% 1) we want to have the labels of the measured channels (each
% corresponds to an electrode)
channel_labels = zeros(size(channel_table.name,1),1);
for kk = 1:size(channel_table.name,1)
	% to which electrode does this channel correspond
    [a,b] = ismember(channel_table.name{kk},electrodes_table.name);
    if a==1 % this channel is an electrode
        channel_labels(kk) = new_label(b);
    end
end

% 2) now we walk through stimulated pairs (in cc_stimsets) and add in
% measured connection matrix
measured_connection_matrix = zeros(length(Wang_ROI_Names),length(Wang_ROI_Names));

for kk = 1:size(cc_stimsets,1)
    
    % which target channels are measured?
    measured_channel_labels = channel_labels;
    measured_channel_labels(cc_stimsets(kk,:)) = []; % remove the stimulated channels

    % get the label of the first stimulated electrode
    stim_el1_label = channel_labels(cc_stimsets(kk,1));
    
    if stim_el1_label>0 % this electrode has a label
        [nn,~] = hist(measured_channel_labels(measured_channel_labels>0),[1:25]);
        % add the number of measured channels
        measured_connection_matrix(stim_el1_label,:) = ...
            measured_connection_matrix(stim_el1_label,:) + nn;
    end
    
    % get the label of the second stimulated electrode
    stim_el2_label = channel_labels(cc_stimsets(kk,2));
    if stim_el2_label>0 % this electrode has a label
        [nn,~] = hist(measured_channel_labels(measured_channel_labels>0),[1:25]);
        % add the number of measured channels
        measured_connection_matrix(stim_el2_label,:) = ...
             measured_connection_matrix(stim_el2_label,:) + nn;
         
    end
end

%%
figure('Position',[0 0 400 400])
imagesc(measured_connection_matrix,[0 max(measured_connection_matrix(:))])
axis square
xlabel('Target Visual Area'),ylabel('Source Visual Area')
set(gca,'XTick',[1 25],'YTick',[1 25],'FontName','Arial','FontSize',18,...
    'XDir','normal');
cm = jet(max(measured_connection_matrix(:)));
cm = [0 0 0; cm];
colormap(cm)
hcb = colorbar;
title(hcb,'#')
set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300','/home/giulio/figures/measured_connections_matrix')

%% Matrix coming from data, 1 if there is a ccep, 0 if no ccep is detected 
measured_N1_matrix = zeros(25,25);

%check which areas were stimulated % e.g. V1d avd V2d
channel_labels; %contains the labeled areas of the electrode that are recorded   

%check which area had an early response % e.g. V3d
conn_plot; %contains a 1 for the electrode and the stimulated pairs that produce an ER

% now we have a table in which channels and electrode number match 
conn_plot; % every column contains a one if that stimulation(column) has produced ccep in that electrode(row)


for ee = 1:size(conn_plot,1)
    for ss = 1:size(conn_plot,2)
        if conn_plot(ee,ss) == 1 && channel_labels(ee) ~= 0
            measured_N1_matrix(channel_labels(ee)) = 1  
        end 
    end
end


%--> put a 1 in the matrix % e.g. at position 2,6 and 4,6 % new_label
    
%gives you the position in the matrix
    
%check which areas did not have an early response % e.g. hV4
    
%--> put a -1 in the matrix % e.g. at position 2,7 and 4,7


%to erase nan from channel_labels 
for kk=1:66 
    if isnan(channel_labels(kk))
        channel_labels(kk) = 0 
    end
end

        
         
%%
%script to create a matrix starting from the ccep data (ccep_detect_onlyN1)
%on the y axis all the stimulation pairs, on the x axis all the target
%areas (electrodes). If CCEP is significant then 1, if not then 0. 

% select significant peaks in the other channels
conn_plot = output_ER(:,:,1) > 0.0001; % amplitude larger than 0.0001

stim_pair = 46;
%to see which stim_pair corresponds to which electrodes 
cc_stimsets (stim_pair, :);

%% we make a matrix of which connections show an N1 

measured_N1_matrix = NaN(25,25);

% now we have a table in which channels and electrode number match 
conn_plot; % every column contains a one if that stimulation(column) has produced ccep in that electrode(row)
for 1:size(conn_plot,1)
    if 
% --> put a 1 in the matrix % e.g. at position 2,6 and 4,6 % new_label
    
% gives you the position in the matrix
    
% check which areas did not have an early response % e.g. hV4
    
% --> put a -1 in the matrix % e.g. at position 2,7 and 4,7