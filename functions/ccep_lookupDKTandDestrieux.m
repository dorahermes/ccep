function [electrodes_tableWithlabels] = ...
    ccep_lookupDKTandDestrieux(g,electrodes_tsv,freesurfer_dir,hemi_small,output_file,electrode_to_vertex_dist)
%
% This function looks up several atlas labels for electrodes in a
% _electrodes.tsv file
% for every electrode, it looks up the gifti vertices within mm_distance
% and extracts the labels from the corresponding atlas.
%
% Inputs:
% g: a giti file with vertices in the same space as the electrode positions
% electrodes_tsv: electrodes.tsv file (BIDS format)
% freesurer_dir: directory with freesurfer output + Benson & Kastner maps
% hemi_small: hemisphere to look up labels (smaller case: l or h)
% output_file: file name to save the new _electrodes.tsv file with labels
% electrode_to_vertex_dist: how far from each electrode to search (default 3 mm)
%
% Output:
% electrode_tsv table with labels added for every electrode
%
% D Hermes, G Castegnaro and J van der Aar, UMC Utrecht, 2019

% load the electrodes.tsv file:
loc_info = readtable(electrodes_tsv,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
elecmatrix = [loc_info.x loc_info.y loc_info.z];

%%%% LOAD SURFACE LABELS FOR ALL THE MAPS

% load Destrieux map
surface_labels_Destrieux = fullfile(freesurfer_dir,'label',...
    [hemi_small 'h.aparc.a2009s.annot']);
[vertices, label, colortable_Destrieux] = read_annotation(surface_labels_Destrieux);
vert_label_Destrieux = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable_Destrieux.table,1) %76 are labels
    vert_label_Destrieux(label==colortable_Destrieux.table(kk,5)) = kk;
end
clear vertices label 

% load DKT map
surface_labels_DKT = fullfile(freesurfer_dir,'label',...
[hemi_small 'h.aparc.DKTatlas.annot']);
[vertices, label, colortable_DKT] = read_annotation(surface_labels_DKT);
vert_label_DKT = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable_DKT.table,1) 
    vert_label_DKT(label==colortable_DKT.table(kk,5)) = kk; % I don't get this line of code
end
clear vertices label 

     
%%% DEFINE THE OUTPUT
DKTatlas_label = NaN(size(elecmatrix,1),1);
DKTatlas_label_text = cell(size(elecmatrix,1),1);
Destrieux_label = NaN(size(elecmatrix,1),1);
Destrieux_label_text = cell(size(elecmatrix,1),1);

%%% LOOP THROUGH ELECTRODES AND ASSIGN LABELS

for elec = 1:size(elecmatrix,1) % loop across electrodes
    
    % dist from electrode to vertices
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
 
    %%%% DESTRIEUX:
    % take the mode of the labels within X mm (X=electrode_to_vertex_dist)
    localized_electrodes = mode(vert_label_Destrieux(find(b<electrode_to_vertex_dist))); 
    % put the labels (vert_label) back in the matrix
    if ~isnan(localized_electrodes)
        Destrieux_label(elec,1) =  localized_electrodes;
        Destrieux_label_text{elec,1} = colortable_Destrieux.struct_names{Destrieux_label(elec,1)};
    else
        Destrieux_label(elec,1) =  NaN;
        Destrieux_label_text{elec,1} = NaN;        
    end
    %%%% DKT:
    % take the mode of the labels within X mm
    localized_electrodes = mode(vert_label_DKT(find(b<electrode_to_vertex_dist))); 
    % put the labels (vert_label) back in the matrix
    if localized_electrodes~=0 && ~isnan(localized_electrodes)
        DKTatlas_label(elec,1) =  localized_electrodes;
        DKTatlas_label_text{elec,1} = colortable_DKT.struct_names{DKTatlas_label(elec,1)};
    else
        DKTatlas_label(elec,1) =  NaN;
        DKTatlas_label_text{elec,1} = NaN;
    end
end


%%%% Integrating electrode positions with TSV-file of coordinates

% Make a new table with added variables
t = table(...
    Destrieux_label, Destrieux_label_text,...
    DKTatlas_label,DKTatlas_label_text); 

% concatenate the table to what is already in loc_info
electrodes_tableWithlabels = horzcat(loc_info,t);

if ~exist(output_file,'file')
    disp(['writing output ' output_file])
    writetable(electrodes_tableWithlabels,output_file,...
        'Filetype','text','Delimiter','\t'); 
else
    disp(['ERROR: can not overwrite, output file already exists ' output_file])
end

