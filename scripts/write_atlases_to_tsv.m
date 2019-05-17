% This code is used to fill the electrodes.tsv with the four Atlases
% DKT atlas, Destrieux atlas, Benson atlas, Wang atlas

% Jaap van der Aar, UMC Utrecht, 02-2019

% Setting data root path
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';

% Adding paths functions ccep_lookupAtlases and bids_tsv_nan2na
addpath('/Fridge/users/jaap/github/ccep/functions');

%%

% Subject information
subjects = {'RESP0768'};
sessions = {'1'};
hemi_cap = {'R'}; 
hemi_small = {'r'};
s = 1;

% get correct subject info
subj = subjects{s};
ses_label = sessions{s};


%% Write gifti file in derivatives/surfaces with original MRI coordinates from freesurfer

%%% Preperation step 1: create the output directory for your surface
    % mkdir dataRootPath/derivatives/surfaces/sub-<>
%%% Preperation step 2: create a gifti file from the freesurfer pial in Freesurfer coordinates
    % run next line in the terminal
    % mris_convert lh.pial lh.pial.gii 

% load the Freesurfer gifti
g = gifti([dataRootPath 'derivatives/freesurfer/sub-' subj '/surf/' hemi_small{s} 'h.pial.gii']);

% convert from freesurfer space to original space
mri_orig = ([dataRootPath '/derivatives/freesurfer/sub-' subj '/mri/orig.mgz']);
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

disp('converted to original space')
%% Script to add the atlases to the tsv file

g = gifti(fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']));
electrodes_tsv = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_electrodes.tsv'];
freesurfer_dir = fullfile(dataRootPath, 'derivatives', 'freesurfer', ['sub-' subj]);
output_file = 'electrode_positions_fouratlases.tsv';
electrode_to_vertex_dist = 3; % in mm

ccep_lookupAtlases(g,electrodes_tsv,freesurfer_dir,dataRootPath,hemi_small,output_file,electrode_to_vertex_dist, subj, ses_label, s)





