% this script looks ffor each electrode which area on teh Destrieux map is
% closest.

%%%% Adding toolboxes:
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
 
% load Destrieux map
surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'label',...
[hemi_small{s} 'h.aparc.a2009s.annot']);
[vertices, label, colortable] = read_annotation(surface_labels_name);
vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable.table,1) % 76 are labels
    vert_label(label==colortable.table(kk,5)) = kk;
end


%%

% create your output matrix
elec_destrieux_label = zeros(size(elecmatrix,1),1);

for elec = 1:size(elecmatrix,1) % loop across electrodes
    % gebruik dist
    tic
    a = dist(elecmatrix(elec,:),g.vertices');
    toc
    
    % program your own dist
    tic 
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
    toc
    
    % take the mode of the labels within 3 mm
    
    
    % put the labels back in the matrix
    elec_destrieux_label(elec,1) =  
end

