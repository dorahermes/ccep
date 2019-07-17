function ccep_write_coordsystemJSON(dataRootPath, subj, ses)

%%%% This script writes coordsystem JSON file
%%%% Make sure data matlab to JSON library is added
%%%% This can be found here: https://github.com/gllmflndn/JSONio

%%%% Dora Hermes, Jaap van der Aar, Giulio Castegnaro, 2019

electrodes_json_name = fullfile(dataRootPath,...
    ['sub-' subj ],['ses-' ses],'ieeg',...
    ['sub-' subj '_ses-' ses '_coordsystem.json']); 

% This line is to create a variable for the name of the file intendedFor
filename_T1w = fullfile(['sub-' subj ],['ses-' ses], 'anat', ['sub-' subj '_ses-' ses '_proc-deface_T1w.nii']);

% assign information and methodology
loc_json.iEEGCoordinateSystem  = 'Other';
loc_json.iEEGCoordinateSystemDescription = 'The origin of the coordinate system is between the ears and the axis are in the RAS direction. The scaling is with respect to the individuals anatomical scan and no scaling or deformation have been applied to the individuals anatomical scan';
loc_json.iEEGCoordinateUnits  = 'mm';
loc_json.iEEGCoordinateProcessingDescripton = 'Surface projection Hermes or Branco';
loc_json.intendedFor = filename_T1w; 
loc_json.iEEGCoordinateProcessingDescription = ''; 
loc_json.iEEGCoordinateProcessingReference = '{Hermes et al., 2010 JNeuroMeth , Branco et al., 2018 JNeuroMeth}';

jsonSaveDir = fileparts(electrodes_json_name);

% ensure there is a /ieeg/ folder 
if ~isdir(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist, create: %s \n',jsonSaveDir)
end

json_options.indent = '    '; % this just makes the json file look prettier 

% write JSON file
jsonwrite(electrodes_json_name,loc_json,json_options)
