%%  This script writes the projected electrodes to the existing, but 'empty', electrodes.tsv
% Save the template (with right x,y,z coordinates) in the corresponding subject folder 

%   By Dora Hermes, Jaap van der Aar, Giulio Castegnaro 02-2019

%   Modified by Jaap van der Aar, adjust to current pipeline/structure 05-2019

% This is where the original identified electrode locations are stored in a
% matrix called 'elecmatrix' that has the size electrodes * x,y,z
working_dir = fullfile('/Fridge','users','jaap','ccep','dataBIDS');
% CCEP folder
CCEP_dir = fullfile('/Fridge','CCEP');


% Insert RESP and session
sub_label = 'RESP0701';
ses_label = '1';


% select the .mat file that contains the 'elecmatrix' (because it sometimes has 
% different names, and not in BIDS, this is best done by GUI)
load(fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    (uigetfile('*.mat','Select *.mat file',...
    [fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg')]))));

% CHANGED SCRIPT TILL HERE, T-FUNCTION HERE IS GOOD
t = readtable(fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']),...
    'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'})
%%
% read the electrodes.tsv file that does not contain positions yet
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

% because NaN is not compatible with BIDS, use this function to change
% NaN's to N/a's 
addpath('/Fridge/users/jaap/github/ccep/functions')
t = bids_tsv_nan2na(t);

%% write electrode file in working dir - check before moving to CCEP

writetable(t, fullfile(working_dir,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label '_ses-' ses_label '_task-' task_label '_run-' run_label '_electrodes.tsv']),...
    'FileType','text','Delimiter','\t');
