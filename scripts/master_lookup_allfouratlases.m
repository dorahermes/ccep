% this script localizes electrodes based on 4 atlases
% Destrieux Atlas, DKT Atlas, Benson Atlas , Wang/Karsner Atlas

%% Setting right paths and loading data
% Make sure that this toolbox is in the path:
% can be cloned from: https://github.com/dorahermes/ecogBasicCode.git
addpath('/Fridge/users/jaap/github/ecogBasicCode/render/')

% adding VistaSoft path
addpath('/home/jaap/vistasoft/external/freesurfer');

%%%% Setting data root path
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';

%%%% Subject information
subjects = {'joure','chaam'};
sessions = {'01','01'};
hemi_cap = {'L','R'}; 
hemi_small = {'l','r'};
s = 2;

% get correct subject info
subj = subjects{s};
ses_label = sessions{s};

% load electrode positions
dataLocName = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_acq-clinicalprojected_electrodes.tsv'];
loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
elecmatrix = [loc_info.x loc_info.y loc_info.z];

% load gifti file for surface coordinates
dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
g = gifti(dataGiiName);
 

%% Destrieux Atlas loading and calculating which area is covered

%%% LOAD ALL THE MAPS
% load Destrieux map
surface_labels_Destrieux = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'label',...
    [hemi_small{s} 'h.aparc.a2009s.annot']);
[vertices, label, colortable_Destrieux] = read_annotation(surface_labels_Destrieux);
vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable_Destrieux.table,1) % 76 are labels
    vert_label(label==colortable_Destrieux.table(kk,5)) = kk;
end
vert_label_Destrieux = vert_label;
clear vertices label 

% load DKT map
surface_labels_DKT = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'label',...
[hemi_small{s} 'h.aparc.DKTatlas40.annot']);
[vertices, label, colortable_DKT] = read_annotation(surface_labels_DKT);
vert_label_DKT = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable_DKT.table,1) 
    vert_label_DKT(label==colortable_DKT.table(kk,5)) = kk; % I don't get this line of code
end
clear vertices label 

%%% DEFINE THE OUTPUT
DKTatlas_label = zeros(size(elecmatrix,1),1);
DKTatlas_label_text = cell(size(elecmatrix,1),1);
Destrieux_label = zeros(size(elecmatrix,1),1);
Destrieux_label_text = cell(size(elecmatrix,1),1);


%%% LOOP THROUGH ELECTRODES AND ASSIGN LABELS
% for every electrode, look up the most common label within 3 mm
electrode_to_vertex_dist = 3; % millimeters

for elec = 1:size(elecmatrix,1) % loop across electrodes
    
    % dist from electrode to vertices
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
 
    %%%% DESTRIEUX:
    % take the mode of the labels within 3 mm
    localized_electrodes = mode(vert_label_Destrieux(find(b<electrode_to_vertex_dist))); 
    % put the labels (vert_label) back in the matrix
    Destrieux_label(elec,1) =  localized_electrodes;
    Destrieux_label_text{elec,1} = colortable_Destrieux.struct_names{Destrieux_label(elec,1)};
    
    %%%% DKT:
    % take the mode of the labels within 3 mm
    localized_electrodes = mode(vert_label_DKT(find(b<electrode_to_vertex_dist))); 
    % put the labels (vert_label) back in the matrix
    DKTatlas_label(elec,1) =  localized_electrodes;
    DKTatlas_label_text{elec,1} = colortable_DKT.struct_names{DKTatlas_label(elec,1)};
    
end
%% Integrating electrode positions with TSV-file of coordinates
% nameing variables
name = loc_info.name;
x = loc_info.x;
y = loc_info.y;
z = loc_info.z;
size = loc_info.size;
group = loc_info.group;
hemisphere = loc_info.hemisphere;
DKTatlas_label = DKTatlas_label;
DKTatlas_name = DKTatlas_label_text;
Destrieux_label = Destrieux_label; 
Destrieux_name = Destrieux_label_text;
% Benson_label = Benson_label;
% Benson_name = Benson_label_text;
% Benson_eccentricity = 
% Benson_polarangle = 
% Benson_sigma = 
% Wang_label = 
% Wang_name =

t = table(name,x,y,z,size,group,hemisphere,DKTatlas_label,DKTatlas_name,Destrieux_label,Destrieux_name); % add variables in here if adding above

electrodes_tsv_name ='electrode_positions_fouratlases.tsv';

writetable(t,electrodes_tsv_name,'FileType','text','Delimiter','\t'); 

%% Create _coordsystem.json

% root_dir = '../';
% ieeg_project = 'templates';
% ieeg_sub = '01';
% ieeg_ses = '01';
% 
% electrodes_json_name = fullfile(root_dir,ieeg_project,...
%     ['sub-' ieeg_sub ],['ses-' ieeg_ses],'ieeg',...
%     ['sub-' ieeg_sub ...
%     '_ses-' ieeg_ses ...
%     '_coordsystem.json']);


root_dir = dataRootPath;
% deleted ieeg_project
ieeg_sub = subj;
ieeg_ses = ses_label;

electrodes_json_name = fullfile(root_dir,...
    ['sub-' ieeg_sub ],['ses-' ieeg_ses],'ieeg',...
    ['sub-' ieeg_sub ...
    '_ses-' ieeg_ses ...
    '_coordsystem.json']); % deleted ieeg_project

loc_json.iEEGCoordinateSystem  = '';
loc_json.iEEGCoordinateUnits  = '';
loc_json.iEEGCoordinateProcessingDescripton = '';
loc_json.IndendedFor = ''; 
loc_json.iEEGCoordinateProcessingDescription = ''; 
loc_json.iEEGCoordinateProcessingReference = '';

jsonSaveDir = fileparts(electrodes_json_name);
if ~isdir(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist, create: %s \n',jsonSaveDir)
end

% m file, but unreadable
% pathname = fileparts('/input/file');
%use that when you save
matfile = fullfile(jsonSaveDir, 'hello.m');
save(matfile);

% % real json file option
% jsonStr = jsonencode(loc_json);
% fid = fopen('hello2.json', 'w');
% if fid == -1, error('Cannot create JSON file'); end
% fwrite(fid, jsonStr, 'char');
% fclose(fid);

json_options.indent = '    '; % this just makes the json file look prettier 
jsonwrite(electrodes_json_name,loc_json,json_options)

%% Adding new paths for Benson and Wang/Kastner atlases

%%%% Setting data root path
dataRootPath = '/Fridge/users/giulio/ccep/dataBIDS/';

%%%% Subject information
subjects = {'chaam'};
sessions = {'01'};
hemi_cap = {'R'}; 
hemi_small = {'r'};
s = 1;

% get correct subject info
subj = subjects{s};
ses_label = sessions{s};

% load electrode positions
dataLocName = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_acq-clinicalprojected_electrodes.tsv'];
loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
elecmatrix = [loc_info.x loc_info.y loc_info.z];
          
% load gifti file for surface coordinates
dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
g = gifti(dataGiiName);

%% Benson Atlas loading and calculating which area is covered

% Loadiong Benson labels 
surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'surf',...
    [hemi_small{s} 'h.benson14_varea.mgz']);
surface_labels_B = MRIread(surface_labels_name);
vert_area_label = surface_labels_B.vol(:);

% check whether first label is V1 or Unknowns
Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','L01','L02','T01','T02','V3b','V3a'};
Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};

elec_benson_label = NaN(size(elecmatrix,1),5); % numbers
elec_benson_label_text = cell(size(elecmatrix,1),1); % strings
 
for elec = 1:size(elecmatrix,1) % loop across electrodes

    % calculate the distance between electrode and g.vertices (surface coordinates)
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
  
    % take the mode of the labels within 3 mm
    area_of_electrode = mode(vert_area_label(find(b<3))); % here script necessary that eliminates electrodes that are not on map
    
    % put the labels (vert_label) back in the matrix
    elec_benson_label(elec,1) =  area_of_electrode;
    
    % here we look up the text string for the current label number
    if area_of_electrode>0
        elec_benson_label_text{elec,1} = Benson_Area_Names{area_of_electrode};
    end
end

%% Calculating Benson eccentricity, polar angle and sigma 

%Benson eccentricity 
    surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'surf',...
    [hemi_small{s} 'h.benson14_eccen.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_eccen_label = surface_labels.vol(:);

%Benson_Eccen_Names = [1:ceil(max(vert_eccen_label))];
    
for elec=1:size(elecmatrix,1) %loop acros electrodes from 1 to 56 
    
    % calculate the distance between electrode and g.vertices (surface coordinates)
        c = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
        
    % take the mean of the eccentricity within 3 mm
    eccen_of_electrode = mean(vert_eccen_label(find(c<3)));     
       
    %put the labels (eccentricity) back in the matrix 
    elec_benson_label(elec,2) = eccen_of_electrode;
   
end
    
% assigning polar angle to every electrode  
surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'surf',...
    [hemi_small{s} 'h.benson14_angle.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_angle_label = surface_labels.vol(:);

for elec=1:size(elecmatrix,1) 
    % calculate the distance between electrode and g.vertices (surface coordinates)
       b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
       
    % take the mean of the polar angle within 3 mm
    angle_of_electrode = mean(vert_angle_label(find(b<3)));     
       
    %put the labels (polar angle) back in the matrix 
    elec_benson_label(elec,3) = angle_of_electrode;   
    
    
end

% assigning sigma to every electrode

surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'surf',...
    [hemi_small{s} 'h.benson14_sigma.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_sigma_label = surface_labels.vol(:);
    
for elec=1:size(elecmatrix,1) 
    % calculate the distance between electrode and g.vertices (surface coordinates)
       b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
       
    % take the mean of the sigma within 3 mm
    sigma_of_electrode = mean(vert_sigma_label(find(b<3)));     
       
    %put the labels (sigma) back in the matrix 
    elec_benson_label(elec,4) = sigma_of_electrode;   
    
    
end 

%% Wang/Kastner Atlas loading and calculating which area is covered    
  surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'surf',...
    [hemi_small{s} 'h.wang15_mplbl.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_wang_label = surface_labels.vol(:);
    
 
for elec = 1:size(elecmatrix,1) % loop across electrodes

    % calculate the distance between electrode and g.vertices (surface coordinates)
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
  
    % take the mode of the labels within 3 mm
    area_of_electrode = mode(vert_wang_label(find(b<3))); 
    
    % put the labels (vert_label) back in the matrix
    elec_benson_label(elec,5) = area_of_electrode;
    
    % here we look up the text string for the current label number
    if area_of_electrode>0
        elec_benson_label_text{elec,1} = Wang_ROI_Names{area_of_electrode};
    end
end