clear all

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('/Users/dora/Documents/git/ecogBasicCode/render/')
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));

%% write gifti

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [19,23,24,9];
hemi = {'L','L','R','R'};
hemi_2 = {'lh','lh','rh','rh'};
s = 2;

if subjects(s)<10
    subj = ['0' num2str(subjects(s))];
else
    subj = num2str(subjects(s));
end

% %%% create a gifti file from a flywheel obj:
[vertex,face] = read_obj([dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/' hemi_2{s} '.pial.obj']);
g.vertices = vertex';
g.faces = face';
g.mat = eye(4,4);
g = gifti(g);
% %%% create a gifti file from the freesurfer pial:
% or run mris_convert lh.pial lh.pial.gii in the terminal and read it with gifti(lh.pial.gii)

% covert to original space:
mri_orig = ([dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/mri/orig.mgz']);
orig = MRIread(mri_orig);
Torig = orig.tkrvox2ras;
Norig = orig.vox2ras;
freeSurfer2T1 = Norig*inv(Torig);

% convert vertices to original space
vert_mat = double(([g.vertices ones(size(g.vertices,1),1)])');
vert_mat = freeSurfer2T1*vert_mat;
vert_mat(4,:) = [];
vert_mat = vert_mat';
g.vertices = vert_mat; clear vert_mat

% save as a gifti
gifti_name = [dataRootPath '/sub-' subj '/ses-01/anat/sub-' subj '_T1w_pial.' hemi{s} '.surf.gii'];

save(g,gifti_name,'Base64Binary')


%% Render plain with used electrodes

subjects = [19,23,24,9];
hemi = {'L','L','R','R'};
els_NBF = {[107 108 109 115 120 121],'',[45 46],''};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

for s = 3%1:length(subjects)
    % subject code
    if subjects(s)<10
        subj = ['0' num2str(subjects(s))];
    else
        subj = num2str(subjects(s));
    end
    
    % electrode locations name:
    dataLocName = fullfile(dataDir,'soc_bids',['/sub-' subj],'/ses-01/','ieeg',...
        ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);
    % gifti file name:
    dataGiiName = fullfile(dataDir,'soc_bids',['/sub-' subj],'/ses-01/','anat',...
        ['/sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);
    % first data run - to get good channels:
    dataName = fullfile(dataDir,'soc_bids',['/sub-' subj],'/ses-01/','ieeg',...
        ['/sub-' subj '_ses-01_task-soc_run-01_ieeg_preproc.mat']);
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);

    % load channel info - do not plot bad channels in rendering
    load(dataName);

    % remove bad channels (replace with NaN)
    elecmatrix(exclude_channels,:) = NaN;
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

        figure
        ecog_RenderGifti(g) % render
%         ecog_Label(els,30,12) % add electrode positions
        el_add(els(els_NBF{s},:),'k',30) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','render',...
            ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end


%% Render plain with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24,9];
hemi = {'L','L','R','R'};

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
    for k = 1%:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

        figure
        ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['../figures/render/subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
        close all
    end
end

%% Render Wang/Kastner with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24,9];
hemi = {'L','L','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};

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
    surface_labels_name = [dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/surf/' hemi_s{s} 'h.wang2015_atlas.mgz'];
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);

    % cmap = 'lines';
    cmap = lines(max(vert_label));
    
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
        ecog_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',['./figures/render/Wang_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
        close all
    end
end


%% Render Benson Areas with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24,9];
hemi = {'L','L','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

Benson_Area_Names = {'V1' 'V2' 'V3'};

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
    surface_labels_name = [dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/surf/' hemi_s{s} 'h.template_areas.mgz'];
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);

    % cmap = 'lines';
    cmap = lines(max(vert_label));
    
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
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Area_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',['./figures/render/BensonAreas_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
        close all
    end
end
   
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


