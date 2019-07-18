function ccep_write_coordinates_2(dataRootPath, subj, ses, t, t_empty, elecmatrix, proj_elec_filename)

%   This function writes changes the NaN's table to n/a's.
%   It then writes and saves the table to the electrodes.tsv file
%   It also saves the elecmatrix, to it can be checked later whether it is
%   done correctly (because the covertion is party manually this is
%   important)

%   By Dora Hermes, Jaap van der Aar, Giulio Castegnaro 2019 (full script)

%   Modified by Jaap van der Aar, adjusted to current pipeline/structure 05-2019

% print output table to double check
if ~isequal(t.x,t_empty.x) 
    disp('electrodes are placed in table')
    disp('double check "t" if elecmatrix could not be place 1-on-1 in table')
    t
end

% Add path of bids_tsv_nan2na.m function and run
% because NaN is not compatible with BIDS, use this function to change
% NaN's to N/a's 
t = bids_tsv_nan2na(t);

% write table to electrodes.tsv file
writetable(t, fullfile(dataRootPath,['sub-' subj],['ses-' ses],'ieeg',...
    ['sub-' subj '_ses-' ses '_electrodes.tsv']),...
    'FileType','text','Delimiter','\t');

% Save file to have the possibility to check whether the elecmatrix is
% filled in the right way. Save m-file as:  RESP<>_ses<>__convert_electrodes_check.mat
save([fullfile(dataRootPath,['sub-' subj],['ses-' ses],'ieeg',...
    ['sub-' subj '_ses-' ses '_convert_electrodes_check.mat'])],...
    't_empty','t','elecmatrix','proj_elec_filename')

disp('TSV-file saved in folder')