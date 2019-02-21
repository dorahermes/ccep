% this script looks for each electrode which area on the Benson map is
% closest.

%%%% Adding toolboxes:
% Make sure that this toolbox is in the path:
% can be cloned from: https://github.com/dorahermes/ecogBasicCode.git
addpath('/Fridge/users/giulio/github/ecogBasicCode/render/')

% adding VistaSoft path
addpath('/home/giulio/vistasoft/external/freesurfer');

%%%% Setting data root path
dataRootPath = '/Fridge/users/giulio/ccep/dataBIDS/';

%%%% Subject information
subjects = {'chaam'};
sessions = {'01'};
hemi_cap = {'R'}; 
hemi_small = {'r'};
s = 1;

% get correct subject info
subj = subjects{s};
ses_label = sessions{s};

% load electrode positions
dataLocName = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_acq-clinicalprojected_electrodes.tsv'];
loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
elecmatrix = [loc_info.x loc_info.y loc_info.z];
          
% load gifti file for surface coordinates
dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
g = gifti(dataGiiName);

% Benson labels 
surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'surf',...
    [hemi_small{s} 'h.benson14_varea.mgz']);
surface_labels_B = MRIread(surface_labels_name);
vert_label = surface_labels_B.vol(:);

% check whether first label is V1 or Unknowns
Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','L01','L02','T01','T02','V3b','V3a'};
Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};
   
%% now we have loaded all the necessary data, we are going to look up the
% labels for every electrode
%% 

% creating a matrix for the electrode labels
elec_benson_label = NaN(size(elecmatrix,1),1); % numbers
elec_benson_label_text = cell(size(elecmatrix,1),1); % strings
 
for elec = 1:size(elecmatrix,1) % loop across electrodes

    % calculate the distance between electrode and g.vertices (surface coordinates)
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
  
    % take the mode of the labels within 3 mm
    area_of_electrode = mode(vert_label(find(b<3))); % here script necessary that eliminates electrodes that are not on map
    
    % put the labels (vert_label) back in the matrix
    elec_benson_label(elec,1) =  area_of_electrode;
    if area_of_electrode>0
        elec_benson_label_text{elec,1} = Benson_Area_Names{area_of_electrode};
    end
end

%%

% how many electrodes in each area:
figure,hist(elec_benson_label,[1:15])

% how many connections to each electrode would be possible 
my_connect = zeros(56,56); %maybe use size of benson/or numb elect size(elecmatrix)
for kk = 1:size(elec_benson_label,1)
    % which electrode do we have now
    my_connect(elec_benson_label(kk),elec_benson_label(setdiff(1:size(elec_benson_label,1),kk))) = 1;
end

figure,
imagesc(my_connect,[0 1])
set(gca,'XTick',[1:15],'YTick',[1:15],'XTickLabel',[],'YTickLabel',[])
axis square
colormap(hot)



