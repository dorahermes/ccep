%%  This is a template in which the temporarily created electrodes matrix can be loaded
%   And can be written into the existing electrodes.tsv file. Save the     
%   template (with right x,y,z coordinates) in the corresponding subject folder 

%   By Dora Hermes, Jaap van der Aar, Giulio Castegnaro 02-2019

% This is where the original identified electrode locations are stored in a
% matrix called 'elecmatrix' that has the size electrodes * x,y,z
working_dir = fullfile('/Fridge','users','jaap','ccep','dataBIDS');
CCEP_dir = fullfile('/Fridge','CCEP');

% load the temp file
load(fullfile(working_dir,'sourcedata','temp','temp.mat'))

sub_label = 'RESP0733';
ses_label = '1a';
task_label = 'SPESclin';
run_label = '031149';

% load CCEP data
t = readtable(fullfile(CCEP_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_electrodes.tsv']),...
    'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% add electrode X positions
t.x(1:64) = elecmatrix(1:64,1);
% t.x(49:68) = NaN;
% add electrode Y positions
t.y(1:64) = elecmatrix(1:64,2);
% t.y(49:68) = NaN;
% add electrode Z positions
t.z(1:64) = elecmatrix(1:64,3);
% t.z(49:68) = NaN;
%% Add path of bids_tsv_nan2na.m function and run 

addpath('/Fridge/users/jaap/github/ccep/functions')
t = bids_tsv_nan2na(t);

%% write electrode file in working dir - check before moving to CCEP

writetable(t, fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_electrodes.tsv']),...
    'FileType','text','Delimiter','\t');



