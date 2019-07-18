function [electrodes_tableWithlabels] = ...
    ccep_lookupAtlases(dataRootPath, subj, ses, freesurfer_dir, hemi_small, hemi_cap)
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
%
% Output:
% electrode_tsv table with labels added for every electrode
%
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

%load Wang map
Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};
surface_labels_name = fullfile(freesurfer_dir,'surf',...
    [hemi_small 'h.wang15_mplbl.mgz']);
surface_labels = MRIread(surface_labels_name);
vert_label_Wang = surface_labels.vol(:);

%load Benson map
Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','L01','L02','T01','T02','V3b','V3a'};
surface_labels_name = fullfile(freesurfer_dir,'surf',...
[hemi_small 'h.benson14_varea.mgz']);
surface_labels_B = MRIread(surface_labels_name);
vert_label_Benson = surface_labels_B.vol(:);

% load Benson Eccen
surface_labels_name = fullfile(freesurfer_dir,'surf',...
    [hemi_small 'h.benson14_eccen.mgz']);
surface_labels = MRIread(surface_labels_name);
vert_eccen_label = surface_labels.vol(:);
clear surface_labels surface_labels_name
%load Benson Angle
surface_labels_name = fullfile(freesurfer_dir,'surf',...
    [hemi_small 'h.benson14_angle.mgz']);
surface_labels = MRIread(surface_labels_name);
vert_angle_label = surface_labels.vol(:);
clear surface_labels surface_labels_name
% load Benson Sigma
surface_labels_name = fullfile(freesurfer_dir,'surf',...
    [hemi_small 'h.benson14_sigma.mgz']);
surface_labels = MRIread(surface_labels_name);
vert_sigma_label = surface_labels.vol(:);
clear surface_labels surface_labels_name
     
%%% DEFINE THE OUTPUT
DKTatlas_label = NaN(size(elecmatrix,1),1);
DKTatlas_label_text = cell(size(elecmatrix,1),1);
Destrieux_label = NaN(size(elecmatrix,1),1);
Destrieux_label_text = cell(size(elecmatrix,1),1);
Wang_label = NaN(size(elecmatrix,1),1); 
Wang_label_text = cell(size(elecmatrix,1),1);
Benson_label = NaN(size(elecmatrix,1),1); 
Benson_label_text = cell(size(elecmatrix,1),1);
Benson_eccen = NaN(size(elecmatrix,1),1);
Benson_polarangle = NaN(size(elecmatrix,1),1);
Benson_sigma = NaN(size(elecmatrix,1),1);

%%% LOOP THROUGH ELECTRODES AND ASSIGN LABELS

% define the range in which the code searches for the most represented
% area. So in this case it will look within 3 mm how which brain regions
% are found most, it takes the most represented one
electrode_to_vertex_dist = 3; % in mm


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
        Destrieux_label_text{elec,1} = 'n/a';        
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
        DKTatlas_label_text{elec,1} = 'n/a';
    end

    %%%% WANG:
    % take the mode of the labels within X mm
    area_of_electrode = mode(vert_label_Wang(find(b<electrode_to_vertex_dist))); 
    % put the labels (vert_label) back in the matrix
    Wang_label(elec,1) = area_of_electrode;
    if area_of_electrode > 0
        Wang_label_text{elec,1} = Wang_ROI_Names{area_of_electrode};
    else
        Wang_label_text{elec,1} = 'n/a';
    end
    
    %%%% BENSON AREA:
    % take the mode of the labels within X mm
    area_of_electrode = mode(vert_label_Benson(find(b<electrode_to_vertex_dist)));
    % put the labels (vert_label) back in the matrix
    Benson_label(elec,1) =  area_of_electrode;
    if area_of_electrode>0
        Benson_label_text{elec,1} = Benson_Area_Names{area_of_electrode};
    else
        Benson_label_text{elec,1} = 'n/a';
    end
         
    % Only add eccentricity, angle and size if there is a Benson area
    if Benson_label(elec,1)>0
        %%%% BENSON ECCENTRICITY:
        % take the mean of the eccentricity within 3 mm
        eccen_of_electrode = mean(vert_eccen_label(find(b<electrode_to_vertex_dist)));     
        %put the labels (eccentricity) back in the matrix 
        Benson_eccen(elec,1) = eccen_of_electrode;

        %%%% BENSON ANGLE:
        % take the mean of the polar angle within 3 mm
        angle_of_electrode = mean(vert_angle_label(find(b<electrode_to_vertex_dist)));     
        %put the labels (polar angle) back in the matrix 
        Benson_polarangle(elec,1) = angle_of_electrode; 

        %%%% BENSON SIGMA:
        % take the mean of the sigma within 3 mm
        sigma_of_electrode = mean(vert_sigma_label(find(b<electrode_to_vertex_dist)));     
        %put the labels (sigma) back in the matrix 
        Benson_sigma(elec,1) = sigma_of_electrode; 
    end
end


%%%% Integrating electrode positions with TSV-file of coordinates

% Make a new table with added variables
t = table(...
    Destrieux_label, Destrieux_label_text,...
    DKTatlas_label,DKTatlas_label_text,...
    Wang_label, Wang_label_text,...
    Benson_label,Benson_label_text, Benson_eccen,...
    Benson_polarangle, Benson_sigma); 

% concatenate the table to what is already in loc_info
electrodes_tableWithlabels = horzcat(loc_info,t);

electrodes_tableWithlabels = bids_tsv_nan2na(electrodes_tableWithlabels);

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
