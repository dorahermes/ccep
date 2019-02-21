%%%% This script writes coordsystem JSON file
%%%% Make sure data matlab to JSON library is added
%%%% This can be found here: https://github.com/gllmflndn/JSONio

%%%% Dora Hermes, Jaap van der Aar, Giulio Castegnaro, 2019
%%
%%%% Setting data root path
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';

%%%% add toolbox
addpath('/Fridge/users/jaap/github/JSONio')
%%
%%%% Add subject and session
subjects = {'RESP0733'};
sessions = {'1a'};
s = 1;

subj = subjects{s};
ses_label = sessions{s};

%%
root_dir = dataRootPath;
% deleted ieeg_project
ieeg_sub = subj;
ieeg_ses = ses_label;

electrodes_json_name = fullfile(root_dir,...
    ['sub-' ieeg_sub ],['ses-' ieeg_ses],'ieeg',...
    ['sub-' ieeg_sub ...
    '_ses-' ieeg_ses ...
    '_coordsystem.json']); % deleted ieeg_project

% This line is to create a variable for the name of the file intendedFor
filename_T1w = fullfile(['sub-' ieeg_sub ],['ses-' ieeg_ses], 'anat', ['sub-' ieeg_sub '_ses-' ieeg_ses '_T1w.nii']);

loc_json.iEEGCoordinateSystem  = 'Other';
loc_json.iEEGCoordinateSystemDescription = 'The origin of the coordinate system is between the ears and the axis are in the RAS direction. The scaling is with respect to the individuals anatomical scan and no scaling or deformation have been applied to the individuals anatomical scan';
loc_json.iEEGCoordinateUnits  = 'mm';
loc_json.iEEGCoordinateProcessingDescripton = 'Surface projection Hermes or Branco';
loc_json.intendedFor = filename_T1w; 
loc_json.iEEGCoordinateProcessingDescription = ''; 
loc_json.iEEGCoordinateProcessingReference = '{Hermes et al., 2010 JNeuroMeth , Branco et al., 2018 JNeuroMeth}';

jsonSaveDir = fileparts(electrodes_json_name);
if ~isdir(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist, create: %s \n',jsonSaveDir)
end

json_options.indent = '    '; % this just makes the json file look prettier 
jsonwrite(electrodes_json_name,loc_json,json_options)
