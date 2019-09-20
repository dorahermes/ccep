% This script renders a brian surface with Benson and Kastner maps
% Can potentially add electrodes on top


% Make sure that this toolbox is in the path:
% can be cloned from: https://github.com/dorahermes/ecogBasicCode.git
addpath('/Fridge/users/giulio/github/ecogBasicCode/render/')


%% Render plain with used electrodes
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';
subjects = {'RESP0703'};
hemi_cap = {'L'};

% pick a viewing angle:
v_dirs = [90 0];%;90 0;90 -60;270 -60;0 0];

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
%     % electrode locations name:
     dataLocName = fullfile(dataRootPath,['/sub-' subj],'/ses-1/','ieeg',...
         ['sub-' subj '_ses-1_electrode_positions_fouratlases.tsv']);

%     % load electrode locations
     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
     elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
%         % make sure electrodes pop out
         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
         ecog_Label(els,30,12) % add electrode positions
         el_add(els(els_NBF{s},:),'k',30) % add electrode positions
        
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
%             ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end

%% Render Wang/Kastner with electrodes

%clear all
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';

subjects = {'RESP0703'};
hemi_cap = {'L','R'};
hemi_small = {'l','r'};

v_dirs = [45 0];%;90 0;90 -60;270 -60;0 0];

Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'surf',...
        [hemi_small{s} 'h.wang15_mplbl.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);
    
    % make a colormap for the labelss
    cmap = lines(max(vert_label));
    
    % electrode locations name:
     dataLocName = fullfile(dataRootPath,['/sub-' subj],'/ses-1/','ieeg',...
         ['sub-' subj '_ses-1_electrode_positions_fouratlases.tsv']);    
     % load electrode locations
     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
     elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
  
        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
%         % make sure electrodes pop out
         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);      
         ecog_Label(els,10,6) % add electrode positions

        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/Wang_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end


%% Render Benson Areas with electrodes

clear all
dataRootPath = '/Fridge/users/giulio/ccep/dataBIDS/';

subjects = {'RESP0751'};
hemi_cap = {'R','L'};
hemi_small = {'r','l'};

v_dirs = [270 0;0 0];

Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','L01','L02','T01','T02','V3b','V3a'};

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'surf',...
        [hemi_small{s} 'h.benson14_varea.mgz']);
    surface_labels_B = MRIread(surface_labels_name);
    vert_label = surface_labels_B.vol(:);

    % create a colormap for the labels
    cmap = lines(max(vert_label));
    cmap(1,:) = [0.4,0.447,0.7410];
    cmap(3,:) = [0.7,0.694,0.125];
    
%     % electrode locations name:
     dataLocName = fullfile(dataRootPath,['/sub-' subj],'/ses-1/','ieeg',...
         ['sub-' subj '_ses-1_electrode_positions_fouratlases.tsv']);
     %     % load electrode locations
     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
     elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Area_Names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   

%         % make sure electrodes pop out
         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
         ecog_Label(els,30,12) % add electrode positions

        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/BensonAreas_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end
   


%%
%% LEFT OF HERE, ADJUST THIS CODE FOR BIDS LATER
%%
%% Render Benson Eccentricity with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24,9];
hemi = {'L','L','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

for s = 1%1:length(subjects)
    % subject code
    if subjects(s)<10
        subj = ['0' num2str(subjects(s))];
    else
        subj = num2str(subjects(s));
    end
    
    % electrode locations name:
    dataLocName = [dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_acq-corrected_electrodes.tsv'];
    % gifti file name:
    dataGiiName = [dataRootPath '/sub-' subj '/ses-01/anat/sub-' subj '_T1w_pial.' hemi{s} '.surf.gii'];
    % first data run - to get good channels:
    dataName = dir([dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_task-*_run-*_ieeg_preproc.mat']);
    % surface labels 
    surface_labels_name = [dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/surf/' hemi_s{s} 'h.template_eccen.mgz'];
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);

    % cmap = 'lines';
    cmap = hsv(ceil(max(vert_label)));
    Benson_Eccen_Names = [1:ceil(max(vert_label))];
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);

    % load channel info - do not plot bad channels in rendering
    load([dataRootPath '/sub-' subj '/ses-01/ieeg/' dataName(1).name]);

    % remove bad channels (replace with NaN)
    elecmatrix(exclude_channels,:) = NaN;
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Eccen_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',['./figures/render/BensonEccen_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
        close all
    end
end


%% Render Benson Angle with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24,9];
hemi = {'L','L','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

for s = 1%1:length(subjects)
    % subject code
    if subjects(s)<10
        subj = ['0' num2str(subjects(s))];
    else
        subj = num2str(subjects(s));
    end
    
    % electrode locations name:
    dataLocName = [dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_acq-corrected_electrodes.tsv'];
    % gifti file name:
    dataGiiName = [dataRootPath '/sub-' subj '/ses-01/anat/sub-' subj '_T1w_pial.' hemi{s} '.surf.gii'];
    % first data run - to get good channels:
    dataName = dir([dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_task-*_run-*_ieeg_preproc.mat']);
    % surface labels 
    surface_labels_name = [dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/surf/' hemi_s{s} 'h.template_angle.mgz'];
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);

    % cmap = 'lines';
    cmap = hsv(ceil(max(vert_label)));
    Benson_Angle_Names = [1:ceil(max(vert_label))];
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);

    % load channel info - do not plot bad channels in rendering
    load([dataRootPath '/sub-' subj '/ses-01/ieeg/' dataName(1).name]);

    % remove bad channels (replace with NaN)
    elecmatrix(exclude_channels,:) = NaN;
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Angle_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',['./figures/render/BensonAngle_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
        close all
    end
end


