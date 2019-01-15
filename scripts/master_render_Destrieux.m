% This script renders a brian surface with Destrieux maps
% Can potentially add electrodes on top
% dhermes % jvanderaar 2019, UMC Utrecht

% Make sure that this toolbox is in the path:
% can be cloned from: https://github.com/dorahermes/ecogBasicCode.git
addpath('/Fridge/users/dora/github/ecogBasicCode/render/')

%% Write gifti file in derivatives/surfaces with original MRI coordinates from freesurfer

%%% Preperation step 1: create the output directory for your surface
    % mkdir dataRootPath/derivatices/surfaces/subjectLabel
%%% Preperation step 2: create a gifti file from the freesurfer pial in Freesurfer coordinates
    % run next line in the terminal
    % mris_convert lh.pial lh.pial.gii 

dataRootPath = '/Fridge/users/dora/ccep/dataBIDS/'; % BIDS dir
subjects = {'chaam'};
hemi_cap = {'R'};
hemi_small = {'rh'};
s = 1;

subj = subjects{s};

% load the Freesurfer gifti
g = gifti([dataRootPath '/derivatives/Freesurfer/sub-' subj '/surf/' hemi_small{s} '.pial.gii']);

% convert from freesurfer space to original space
mri_orig = ([dataRootPath '/derivatives/Freesurfer/sub-' subj '/mri/orig.mgz']);
orig = MRIread(mri_orig);
Torig = orig.tkrvox2ras;
Norig = orig.vox2ras;
freeSurfer2T1 = Norig*inv(Torig);

% convert vertices to original space
vert_mat = double(([g.vertices ones(size(g.vertices,1),1)])');
vert_mat = freeSurfer2T1*vert_mat;
vert_mat(4,:) = [];
vert_mat = vert_mat';
g.vertices = vert_mat; clear vert_mat

% save as a gifti
gifti_name = [dataRootPath '/derivatives/surfaces/sub-' subj '/sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii'];

save(g,gifti_name,'Base64Binary')


%% Render plain with used electrodes
dataRootPath = '/Fridge/users/dora/ccep/dataBIDS/';
subjects = {'chaam'};
hemi_cap = {'R'};

% pick a viewing angle:
v_dirs = [90 0];%;90 0;90 -60;270 -60;0 0];

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
%     % electrode locations name:
%     dataLocName = fullfile(dataRootPath,'soc_bids',['/sub-' subj],'/ses-01/','ieeg',...
%         ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);

%     % load electrode locations
%     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
%     elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
%         % make sure electrodes pop out
%         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
%         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
%         ecog_Label(els,30,12) % add electrode positions
%         el_add(els(els_NBF{s},:),'k',30) % add electrode positions
        
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
%             ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end

%% Render Destrieux with electrodes

clear all
dataRootPath = '/Fridge/users/dora/ccep/dataBIDS/';

subjects = {'chaam','chaam'};
hemi_cap = {'R','L'};
hemi_small = {'r','l'};

v_dirs = [270 0];%;90 0;90 -60;270 -60;0 0];


for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'label',...
        [hemi_small{s} 'h.aparc.a2009s.annot']);
    %surface_labels = MRIread(surface_labels_name);
    [vertices, label, colortable] = read_annotation(surface_labels_name);
    vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
    % mapping labels to colortable
    for kk = 1:size(colortable.table,1) % 76 are labels
        vert_label(label==colortable.table(kk,5)) = kk;
    end
    
    % make a colormap for the labels
    cmap = colortable.table(:,1:3)./256;
    
    % electrode locations name:
%     dataLocName = [dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_acq-corrected_electrodes.tsv'];
%     % load electrode locations
%     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
%     elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
  
        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,colortable.struct_names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
%         % make sure electrodes pop out
%         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
%         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);      
%         ecog_Label(els,30,12) % add electrode positions

        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/Wang_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end



