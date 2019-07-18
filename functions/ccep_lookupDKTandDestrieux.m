function [electrodes_tableWithlabels] = ...
    ccep_lookupDKTandDestrieux(dataRootPath, subj, ses, freesurfer_dir, hemi_small, hemi_cap)
%
% This function looks up DKT and Destrieux atlas labels for electrodes and 
% writes them to the _electrodes.tsv file
% for every electrode, it looks up the gifti vertices within mm_distance
% and extracts the labels from the corresponding atlas. It takes the most
% represented brain area within 3 mm. 
%
% Inputs:
% dataRootPath
% subj:             subject(s) number
% ses:              sessios(s) number
% freesurfer_dir:   directory with freesurfer output + Benson & Kastner maps
% hemi_small:       hemisphere to look up labels (smaller case: l or h)


% Output:
% electrode_tsv table with labels added for every electrode

% Same function as ccep_lookupDKTandDestrieux, but without the Benson &
% Wang atlases

% D Hermes, G Castegnaro and J van der Aar, UMC Utrecht, 2019

% load gifti file and electrodes.tsv file
g = gifti(fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap '.surf.gii'])); 
electrodes_tsv = [dataRootPath '/sub-' subj '/ses-' ses '/ieeg/sub-' subj '_ses-' ses '_electrodes.tsv'];

% define name output file
output_file = 'electrode_positions_fouratlases.tsv';

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

% define the range in which the code searches for the most represented
% area. So in this case it will look within 3 mm how which brain regions
% are found most, it takes the most represented one
electrode_to_vertex_dist = 3; % in mm

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


% write table
if ~exist(fullfile(dataRootPath,['sub-' subj],['ses-' ses],'ieeg', ...
        ['sub-' subj '_ses-' ses '_' output_file]),'file')
    disp(['writing output ' output_file])
    writetable(electrodes_tableWithlabels, ...
        fullfile(dataRootPath,['sub-' subj],['ses-' ses],'ieeg', ...
        ['sub-' subj '_ses-' ses '_' output_file]),...
        'Filetype','text','Delimiter','\t'); 
else
    disp(['ERROR: can not overwrite, output file already exists ' output_file])
end

