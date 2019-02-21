% this script looks for each electrode which area on the Destrieux map is
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
[hemi_small{s} 'h.aparc.DKTatlas40.annot']);
[vertices, label, colortable] = read_annotation(surface_labels_name);
vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable.table,1) 
    vert_label(label==colortable.table(kk,5)) = kk; % I don't get this line of code
end

%% calculate which brain area is covered by the electrode

% create your output matrix
elec_DKTatlas_label = zeros(size(elecmatrix,1),1);
elec_DKTatlas_label_text = cell(size(elecmatrix,1),1);

for elec = 1:size(elecmatrix,1) % loop across electrodes
    
    % program your own dist
    bb = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
 
    
    % take the mode of the labels within 3 mm
    localized_electrodes = mode(vert_label(find(bb<3))); % note: cannot find where the labelnames are
  
    
    % put the labels (vert_label) back in the matrix
    elec_DKTatlas_label(elec,1) =  localized_electrodes;
    elec_DKTatlas_label_text{elec,1} = colortable.struct_names{elec_DKTatlas_label(elec,1)};
end

%% visualisation of electrode positions

% how many electrodes in each area:
figure,hist(elec_DKTatlas_label,[1:36])
colormap(hot)
xlim([1 36])
% how many connections to each electrode would be possible 
my_connect = zeros(36,36);
for kk = 1:size(elec_DKTatlas_label,1)
    % which electrode do we have now
    my_connect(elec_DKTatlas_label(kk),elec_DKTatlas_label(setdiff(1:size(elec_DKTatlas_label,1),kk))) = 1;
end

figure,
imagesc(my_connect,[0 1])
colormap(hot)


%% Add label code and name to electrode positions TSV file

% contents = load([dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_acq-clinicalprojected_electrodes1.tsv']);
% 
% contents(:,end+1)= elec_DKTatlas_label;
% 
% % 
% % fid = fopen([dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_acq-clinicalprojected_electrodes.tsv'], 'a');
% % fprintf(fid, '%i', elec_DKTatlas_label);
% % fclose(fid);


