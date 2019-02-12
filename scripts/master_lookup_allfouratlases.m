%%This script localizes electrodes based on 4 atlases
% Destrieux Atlas, DKT Atlas, Benson Atlas , Wang/Kastner Atlas

%% Setting right paths and loading data
% Make sure that this toolbox is in the path:
% can be cloned from: https://github.com/dorahermes/ecogBasicCode.git
addpath('/Fridge/users/jaap/github/ecogBasicCode/render/')

% add CCEP path
addpath('/Fridge/users/giulio/github/ccep/functions/')

% adding VistaSoft path
addpath('/home/jaap/vistasoft/external/freesurfer');

%%%% Setting data root path
dataRootPath = '/Fridge/users/giulio/ccep/dataBIDS/';

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

% load gifti file for surface coordinates
dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
g = gifti(dataGiiName);

freesurfer_dir = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj]);

output_file = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_acq-clinicalprojectedLabels_electrodes.tsv'];

mm_distance = 3;

[electrodes_tableWithlabels] = ...
    ccep_lookupAtlases(g,dataLocName,freesurfer_dir,hemi_small{s},output_file,mm_distance);


