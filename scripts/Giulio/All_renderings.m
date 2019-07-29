%% Render Wang/Kastner with electrodes
%clear all
% Make sure that this toolbox is in the path:
% can be cloned from: https://github.com/dorahermes/ecogBasicCode.git
addpath('/Fridge/users/giulio/github/ecogBasicCode/render/')

dataRootPath = '/Fridge/users/giulio/ccep/dataBIDS/';

v_dirs = [90 0;0 0;-12 4];
    
subjects = {'RESP0315','RESP0751','RESP0401','RESP0405','RESP0306'}; % 'RESP0703' cannot plot 0703 because of the gifti file
hemi_cap = {'L','R', 'R', 'L', 'R', 'L'};
hemi_small = {'l','r', 'r', 'l', 'r', 'l'};

v_dirs = [270 0];%;90 0;90 -60;270 -60;0 0];

Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};

for s = [2]
    
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
    cmap1 = [0.2,0.447,0.7410;
            0.650,0.3250,0.098;
            0.29,0.6940,0.1250;
            0.5940,0.1840,0.550;
            0.750,0.674,0.1880;
            0.201,0.745,0.9330;
            0.5350,0.0780,0.184;
            0,0.87,0.741;
            0.750,0.3250,0.0980;
            0.92,0.754,0.1250;
            0.8940,0.1840,0.5560;
            0.2360,0.674 ,1 ;
            0.1010,0.245 ,0.933;
            0.63500,0.289  ,0.184;
            0.5,0.447,0.741  ;
            0.350,0.325 ,0.098;
            0.750,0.6940,0.125;
            0.59400,0.1840,0.556;
            0.43300,0.89400,0.1880;
            0.1001,0.7450,0.93;
            0.6350,0.078,0.1840;
            0,0.29,0.941;
            0.81,0.3250,0.198;
            0.1,0.3,0.125;
            0.794,0.7940,0.896];
    
    
    
    % electrode locations name:
     dataLocName = fullfile(dataRootPath,['/sub-' subj],'/ses-1/','ieeg',...
         ['sub-' subj '_ses-1_electrode_positions_fouratlases.tsv']);    
     % load electrode locations
     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}, 'ReadVariableNames', true);
     elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
  
        figure('Position',[0 0 3000 3000])

        ecog_RenderGiftiLabels(g,vert_label,cmap1,Wang_ROI_Names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
%         % make sure electrodes pop out
         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);      
         ecog_Label(els,50) % add electrode positions
         set(gcf,'PaperPositionMode','auto')
         title([ subjects{s}]);
         %print('-dpng','-r300',[dataRootPath, '/figures/render/Wang_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
         %close all

            

    end
end
