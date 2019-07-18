% This script can be used as workflow script for the preprocessing of the
% CCEP data in the RESPect database. Please use the README CCEP workflow
% file for supplementary information and steps that do not include MATLAB.
% For this workflow, the BIDS-specifications are used. 

% Script includes workflow for:
%   - Adding meta-information in JSON-format
%   - Localization of electrodes
%   - Writing electrode coordinates to TSV-file. Cortical rendering based
%     on Destrieux, DKT, Benson & Wang Atlases.

% Add following toolboxes to path: 
%   - ecogBasicCode:    github.com/dorahermes/ecogBasicCode.git
%   - JSONio:           github.com/gllmflndn/JSONio

% Jaap van der Aar, Giulio Castegnaro, Dora Hermes, Dorien van Blooijs, 2019 

%% Metadata: fill in yourself

% set rootpath
dataRootPath = fullfile('/Fridge','users','jaap','ccep','dataBIDS');

% add folder with functions and scripts
addpath(fullfile('/Fridge','users','jaap','github','ccep'))

% add subject(s) information
subjects = {'RESP0754'};
sessions = {'1'};
hemi_smalls = {'r'};
hemi_caps = {'R'};
% s gives the possibility to add more then 1 subjects, or multiple sessions
% in once. Wrap functions into [ss = 1:size(s) ... end] to iterate over subject data cells
s = 1;
subj = subjects{s};
ses = sessions{s};
hemi_small = hemi_smalls{s};
hemi_cap = hemi_caps{s};

%% STEP 0: preperations

% add ecogBasicCode to path
addpath(fullfile('/Fridge','users','jaap','github','ecogBasicCode')) 


%% STEP 2: fill in JSON file

% be sure the JSONio toolbox is added to the path
addpath(fullfile('Fridge','users,','jaap','github','JSONio'))

ccep_write_coordsystemJSON(dataRootPath, subj, ses)


%% STEP 5: localization of electrodes

%%%% use ALICE or script Dora 



%% STEP 6.1: write coordinates to electrodes TSV-file - loading

% GUI to select the file containing the projected electrodes matrix
% prints empty electrodes table and elecmatrix, and their size, for next
% step. GUI currently set to username/ccep/dataBIDS/sub-<>/-ses-<>/ieeg/
ccep_write_coordinates_1(dataRootPath, subj, ses)

%% STEP 6.2: write coordinates to electrodes TSV-file - manually adding

% Put the values of the elecmatrix in the right places in the table: 
% The elecmatrix is shorter then the TSV, so fill in NaN’s in the gaps. 
% Also there might be a marker on the 65th position, and it possible that
% electrodes are in another order then the TSV. Therefore, this has to be 
% done manually. 

% add electrode X positions
t.x(1:64) = elecmatrix(1:64,1);
t.x(65:133) = NaN;

% add electrode Y positions
t.y(1:64) = elecmatrix(1:64,2);
t.y(65:133) = NaN;

% add electrode Z positions
t.z(1:64) = elecmatrix(1:64,3);
t.z(65:133) = NaN;

%% STEP 6.3: write coordinates to electrodes TSV-file - writing & saving

% make sure the bids_tsv_nan2na.m function is in the path. This function
% writes NaN's to n/a's. It can be found in the ccep/functions folder. 
% write table to electrodes.tsv, and saves the elecmatrix for later check
ccep_write_coordinates_2(dataRootPath, subj, ses, t, t_empty)

%% STEP 7: write atlases labels to electrodes TSV-file

g = gifti(fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap '.surf.gii'])); 
electrodes_tsv = [dataRootPath 'sub-' subj '/ses-' ses '/ieeg/sub-' subj '_ses-' ses '_electrodes.tsv'];
freesurfer_dir = fullfile(dataRootPath, 'derivatives', 'freesurfer', ['sub-' subj]);
output_file = 'electrode_positions_fouratlases.tsv';
electrode_to_vertex_dist = 3; % in mm

ccep_lookupAtlases(g,electrodes_tsv, subj, ses, freesurfer_dir,dataRootPath,hemi_small,output_file,electrode_to_vertex_dist)



