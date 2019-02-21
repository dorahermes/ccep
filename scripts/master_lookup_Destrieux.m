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
[hemi_small{s} 'h.aparc.a2009s.annot']);
[vertices, label, colortable] = read_annotation(surface_labels_name);
vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable.table,1) % 76 are labels
    vert_label(label==colortable.table(kk,5)) = kk;
end


%% calculate which brain area is covered by the electrode

% create your output matrix
elec_destrieux_label = zeros(size(elecmatrix,1),1);
elec_destrieux_label_text = cell(size(elecmatrix,1),1);

for elec = 1:size(elecmatrix,1) % loop across electrodes
    
    % program your own dist
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
 
    
    % take the mode of the labels within 3 mm
    localized_electrodes = mode(vert_label(find(b<3))); % note: cannot find where the labelnames are
  
    
    % put the labels (vert_label) back in the matrix
    elec_destrieux_label(elec,1) =  localized_electrodes;
    elec_destrieux_label_text{elec,1} = colortable.struct_names{elec_destrieux_label(elec,1)};
end

%% visualisation of electrode positions

% how many electrodes in each area:
figure,hist(elec_destrieux_label,[1:74])
colormap(hot)
xlim([1 75])
% how many connections to each electrode would be possible 
my_connect = zeros(74,74);
for kk = 1:size(elec_destrieux_label,1)
    % which electrode do we have now
    my_connect(elec_destrieux_label(kk),elec_destrieux_label(setdiff(1:size(elec_destrieux_label,1),kk))) = 1;
end

figure,
imagesc(my_connect,[0 1])
% set(gca,'XTick',[1:74],'YTick',[1:74],'XTickLabel',[],'YTickLabel',[])
% axis square
colormap(hot)


