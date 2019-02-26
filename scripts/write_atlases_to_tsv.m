% This code is used to fill the electrodes.tsv with the four Atlases
% DKT atlas, Destrieux atlas, Benson atlas, Wang atlas

% Jaap van der Aar, UMC Utrecht, 02-2019

% Setting data root path
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';

% Adding paths functions ccep_lookupAtlases and bids_tsv_nan2na
addpath('/Fridge/users/jaap/github/ccep/functions');

% Subject information
subjects = {'RESP0706'};
sessions = {'1'};
hemi_cap = {'L'}; 
hemi_small = 'l';
run_label = '041501';
s = 1;

% get correct subject info
subj = subjects{s};
ses_label = sessions{s};

% Load gifti file for surface coordinates
dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
g = gifti(dataGiiName);

% Setting TSV-file
electrodes_tsv = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_task-SPESclin_run-' run_label '_electrodes.tsv'];

% Setting freesurfer directory
freesurfer_dir = [dataRootPath 'derivatives/freesurfer/sub-' subj];

% Naming output file
output_file = ['sub-' subj '_ses-' ses_label '_task-SPESclin_run-' run_label '_atlases' '_electrodes.tsv'];

% Setting how far from each electrode to search (default 3 mm)
electrode_to_vertex_dist = 3;

% Run the function that writes the atlases to the TSV file and assign to t
% t = ccep_lookupAtlases(g,electrodes_tsv,freesurfer_dir,hemi_small,output_file,electrode_to_vertex_dist)
t = ccep_lookupDKTandDestrieux(g,electrodes_tsv,freesurfer_dir,hemi_small,output_file,electrode_to_vertex_dist);

% Run the function that changes the NaN in TSV file to n/a
% Currently it only overwrites the labels and not label text (due to ...
% Destrieux_label(elec,1) =  NaN; vs. Destrieux_label_text{elec,1} = NaN;) 
t = bids_tsv_nan2na(t)

% Save new TSV in the corresponding sub/ses/ieeg/ folder
writetable(t, fullfile(dataRootPath,['sub-' subj],['ses-' ses_label],'ieeg',...
    [output_file]),'FileType','text','Delimiter','\t');



