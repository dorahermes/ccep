%% Data for ROC-curves

working_dir = fullfile('/Fridge','users','jaap','ccep','dataBIDS');

load(fullfile(working_dir,'validation', 'sub-RESP0621_ses-1_run-021147_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0621_ses-1_run-021147_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0706_ses-1_run-041501_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0706_ses-1_run-041501_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0733_ses-1b_run-050941_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0733_ses-1b_run-050941_validation_matrix.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0768_ses-1_run-021704_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'sub-RESP0768_ses-1_run-021704_validation_matrix.mat'))

% load variables/parameters that are used 

thresh = [1:.2:5];   
minSD = 50;   

n1_samples_start = find(tt>0.01,1);
n1_samples_end = [find(tt>0.04,1), find(tt>0.05,1), find(tt>0.06,1), find(tt>0.07,1),...
    find(tt>0.08,1), find(tt>0.09,1), find(tt>0.1,1), find(tt>0.11,1)]; 



%% figure with multiple ROC-curves for multiple datasets

fig_individual_rocs = figure('Position',[0 0 1500 900]);

% could use sgtitle()in latest versions of matlab (2018b) 

% ROC RESP0621
subplot(2,5,1:2)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(parameters_optimalize_mat_0621,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(parameters_optimalize_mat_0621,2)

    sens_plot = parameters_optimalize_mat_0621(:,time_th,1);
    spes_plot = parameters_optimalize_mat_0621(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1 - specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
axis square
title('ROC-curves for RESP0621 - 21yr')
set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',16)
% ROC RESP0706

subplot(2,5,3:4)
% plot change level line
 plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(parameters_optimalize_mat_0706,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(parameters_optimalize_mat_0706,2)

    sens_plot = parameters_optimalize_mat_0706(:,time_th,1);
    spes_plot = parameters_optimalize_mat_0706(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1 - specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
axis square
title('ROC-curves for RESP0706 - 50yr')
set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',16)

% ROC RESP0733

subplot(2,5,6:7)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(parameters_optimalize_mat_0733,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(parameters_optimalize_mat_0733,2)

    sens_plot = parameters_optimalize_mat_0733(:,time_th,1);
    spes_plot = parameters_optimalize_mat_0733(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1 - specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
axis square
title('ROC-curves for RESP0733 - 9jr')
set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',16)
% ROC RESP0768

subplot(2,5,8:9)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(parameters_optimalize_mat_0768,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(parameters_optimalize_mat_0768,2)

    sens_plot = parameters_optimalize_mat_0768(:,time_th,1);
    spes_plot = parameters_optimalize_mat_0768(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1 - specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
axis square
title('ROC-curves for RESP0768 - 6yr')
set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',16)

% plot legenda for all ROCs
time_end = round(tt(n1_samples_end),2) * 1000;

subplot(2,5,[5 10]),hold on
for time_th = 1:size(parameters_optimalize_mat_0621,2)
    for clrbar=1:75
        
        plot(clrbar,time_end(time_th),'.','MarkerSize',20,'Color',my_colors(time_th,:))
    end
end
xlim([-4 80]), set(gca,'XTick',[], 'FontName', 'arial', 'FontSize', 16)
ylim([35,115]), ylabel('time end (ms)')
ttl = title('endpoint N1 range'); 
% ttl_pos = get(ttl,'position'); ttl_pos(2) = 120; set( ttl , 'position' , ttl_pos);
box on

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/home/jaap/figures/ROCcurves_individual')


%% load data for averaged ROC-curves
load(fullfile(working_dir,'validation', 'averaged_parameters_optimalization.mat'))
load(fullfile(working_dir,'validation', 'averaged_parameters_optimalization_corrected.mat'))
%% figure with averaged ROC-curves, with and without RESP0768

fig_avg_rocs = figure('Position',[0 0 1700 500]);

% plot averaged ROCs
subplot(1,5,1:2)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(averaged_parameter_scores,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(averaged_parameter_scores,2)

    sens_plot = averaged_parameter_scores(:,time_th,1);
    spes_plot = averaged_parameter_scores(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1 - specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
axis square
title('Averaged ROC-curves')
set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',16)

% plot averaged ROCs minus RESP0768
subplot(1,5,3:4)
% plot change level line
plot([0 1],[0 1],'k'),hold on

% use jet as colors
my_colors = jet(size(averaged_parameter_scores_corrected,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(averaged_parameter_scores_corrected,2)

    sens_plot = averaged_parameter_scores_corrected(:,time_th,1);
    spes_plot = averaged_parameter_scores_corrected(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end
xlabel('1 - specificity')
ylabel('sensitivity')
xlim([0 1]),ylim([0 1])
axis square
title('Averaged ROC-curves minus RESP0768')
set(gca,'XTick',[0:.5:1],'YTick',[0:.5:1],'FontName','arial','FontSize',16)

% plot legenda

time_end = round(tt(n1_samples_end),2) * 1000;

subplot(1,5,5),hold on
for time_th = 1:size(averaged_parameter_scores,2)
    for clrbar=1:75
        
        plot(clrbar,time_end(time_th),'.','MarkerSize',20,'Color',my_colors(time_th,:))
    end
end
xlim([-4 80]), set(gca,'XTick',[])
ylim([35,115]), ylabel('time end (ms)')
ttl = title('endpoint N1 range'); 
% ttl_pos = get(ttl,'position'); ttl_pos(2) = 120; set( ttl , 'position' , ttl_pos);
box on
set(gca,'FontName','arial','FontSize',16)
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/home/jaap/figures/ROCcurves_averages')


%% figure for zoomed in on corrected average

figure()
plot([0 1],[0 1],'k'),hold on
% use jet as colors
my_colors = jet(size(averaged_parameter_scores,2));

% for all different ranges plot a ROC- curve
for time_th = 1:size(averaged_parameter_scores,2)

    sens_plot = averaged_parameter_scores(:,time_th,1);
    spes_plot = averaged_parameter_scores(:,time_th,2);

    plot(1-spes_plot,sens_plot,'Color',my_colors(time_th,:))
end

%plot(1-averaged_parameter_scores_corrected(optimum_mat_x,optimum_mat_y,2),...
 %   averaged_parameter_scores_corrected(optimum_mat_x,optimum_mat_y,1),'k*','MarkerSize',10)
% title('Averaged ROC-curves')
% 
% xlabel('1 - specificity')
% ylabel('sensitivity')
xlim([0 .3]),ylim([0.7 1])
axis square

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/home/jaap/figures/ROCcurves_averages_corrected_unzoomed_optimum')

%% ROC and Confusion for all data with default settings

class_1_visual = reshape(validation_matrix(:,:,1),1,size(cc_epoch_sorted_avg,1)*size(cc_epoch_sorted_avg,2));

class_2_code = reshape(validation_matrix(:,:,2),1,size(cc_epoch_sorted_avg,1)*size(cc_epoch_sorted_avg,2));

TP_array = find(class_1_visual == 1 & class_2_code == 1);
TN_array = find(class_1_visual == 0 & class_2_code == 0);
FP_array = find(class_1_visual == 0 & class_2_code == 1);
FN_array = find(class_1_visual == 1 & class_2_code == 0);

TP = length(find(class_1_visual == 1 & class_2_code == 1));
TN = length(find(class_1_visual == 0 & class_2_code == 0));
FP = length(find(class_1_visual == 0 & class_2_code == 1));
FN = length(find(class_1_visual == 1 & class_2_code == 0));

sensitivity = TP/ (TP + FN);
specificity = TN / (TN + FP);

figure()
[c,cm,ind,per] = confusion(class_1_visual,class_2_code);
plotconfusion(class_1_visual,class_2_code,'Performance N1 Peak Detector')
xlabel('Visual Assessment')
ylabel('Code Performance')

figure()
plotroc(class_1_visual,class_2_code)
title('ROC curve with default parameters')

%% plotting confusion matrix + ROC for default for all sets

% create a huge valiation matrix for both the visual as code assessment

% create empty 200 x 200 to ensure enough space to put all data
validation_matrix_all_visual = NaN(200,200);
% fill upper left with resp0768 - if wanted
% validation_matrix_all_visual(1:size(validation_matrix_0768,1),1:size(validation_matrix_0768,2)) = validation_matrix_0768(:,:,1);
% fill bottom left with resp0733
validation_matrix_all_visual(101:100+size(validation_matrix_0733,1),1:size(validation_matrix_0733,2)) = validation_matrix_0733(:,:,1);
% fill upper right with resp0706
validation_matrix_all_visual(1:size(validation_matrix_0706,1),101:100+size(validation_matrix_0706,2)) = validation_matrix_0706(:,:,1);
% fill bottom right with resp0621
validation_matrix_all_visual(101:100+size(validation_matrix_0621,1),101:100+size(validation_matrix_0621,2)) = validation_matrix_0621(:,:,1);

% repeat for code assessment 
% create empty 200 x 200 to ensure enough space to put all data
validation_matrix_all_code = NaN(200,200);
% fill upper left with resp0768 - if wanted
% validation_matrix_all_code(1:size(validation_matrix_0768,1),1:size(validation_matrix_0768,2)) = validation_matrix_0768(:,:,2);
% fill bottom left with resp0733
validation_matrix_all_code(101:100+size(validation_matrix_0733,1),1:size(validation_matrix_0733,2)) = validation_matrix_0733(:,:,2);
% fill upper right with resp0706
validation_matrix_all_code(1:size(validation_matrix_0706,1),101:100+size(validation_matrix_0706,2)) = validation_matrix_0706(:,:,2);
% fill bottom right with resp0621
validation_matrix_all_code(101:100+size(validation_matrix_0621,1),101:100+size(validation_matrix_0621,2)) = validation_matrix_0621(:,:,2);


% reshape to vector for plotting
class_1_visual = reshape(validation_matrix_all_visual,1,size(validation_matrix_all_visual,1)*size(validation_matrix_all_visual,2));

class_2_code = reshape(validation_matrix_all_code,1,size(validation_matrix_all_code,1)*size(validation_matrix_all_code,2));

figure()
[c,cm,ind,per] = confusion(class_1_visual,class_2_code);
plotconfusion(class_1_visual,class_2_code,'Performance N1 Peak Detector')
xlabel('Visual Assessment')
ylabel('Code Performance')

figure()
plotroc(class_1_visual,class_2_code)
title('ROC curve with default parameters')

%% Render responses in other electrodes (ampl + latency)
addpath('/Fridge/users/jaap/github/ecogBasicCode/render/')
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';
subjects = {'RESP0621'};
hemi_cap = {'R'};

% pick a viewing angle:
v_dirs = [90 0]; %;90 0;90 -60;270 -60;0 0];

% number of subjects 
s = 1;


% set stimulated pairs you want to render
for stim_pair = 19%:3 %1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
    % electrode locations name:
    dataLocName = dir(fullfile(dataRootPath,['sub-' subj],'ses-1','ieeg',...
        ['sub-' subj '_ses-1_electrodes.tsv']));
    dataLocName = fullfile(dataLocName(1).folder,dataLocName(1).name);
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % select significant peaks in the other channels
    n1_plot = squeeze(output_ER_all(:,stim_pair,3:4)); % measured electrodes X latency/ampl of N1 peak
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
        
        % ecog_Label(els,30,12) % add electrode positions
        % add electrodes
        ccep_el_add(els,[0.1 0.1 0.1],20)
        % give stimulated electrodes other color
        ccep_el_add(els(cc_stimsets(stim_pair,1:2),:),[0 0 0],40)
        % set size and color of channels with significant peak 
        % based on sample (from stimulation on, so -5120) and the amplitude
        % color = latency, size = amplitude
        ccep_el_add_size_and_color(els,n1_plot(:,2),(n1_plot(:,1)-5121),500,160)
        title(['Responses for stimulation of ' data_hdr.label{cc_stimsets(stim_pair,1)} ' and ' data_hdr.label{cc_stimsets(stim_pair,2)} ])
        set(gcf,'PaperPositionMode','auto')
        % print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
        % ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end

%% Render brain with Destrieux atlas

dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';
% add vistasoft for read_annotation
addpath('/home/jaap/vistasoft/external/freesurfer');

subjects = {'RESP0621','RESP0768'};
sessions = {'1','1'};
hemi_cap = {'R','R'}; 
hemi_small = {'r','r'};

v_dirs = [90 0];%;90 0;90 -60;270 -60;0 0];


for s = 1:2%1:length(subjects)
    % subject code 
    subj = subjects{s};
    ses_label = sessions{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'label',...
    [hemi_small{s} 'h.aparc.a2009s.annot']);
    % surface_labels = MRIread(surface_labels_name);
    [vertices, label, colortable] = read_annotation(surface_labels_name);
    vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
    % mapping labels to colortable
    for kk = 1:size(colortable.table,1) % 76 are labels
        vert_label(label==colortable.table(kk,5)) = kk;
    end
    
    % make a colormap for the labels
    cmap = colortable.table(:,1:3)./256;
    
    % electrode locations name:
    dataLocName = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_electrodes.tsv'];
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
  
        figure('Color',[1 1 1],'Position',[0 0 1300 800])
        ecog_RenderGiftiLabels(g,vert_label,cmap,colortable.struct_names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);      
        ecog_Label(els,30,12) % add electrode positions
        
        % title
        title(['Brain rendering of ' subj ' using Destrieux Atlas '])
        
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/Wang_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end

%% plot averaged epochs
stim_pair_nr = 1;
ccep_elec = 3;
figure('Position',[0 0 300 200]),hold on
plot(tt(tt>-1 & tt<1),zeros(size(tt(tt>-1 & tt<1))),'Color',[.5 .5 .5])
plot([0 0],[-1200 1200],'Color',[.5 .5 .5])
%plot(tt,squeeze(cc_epoch_sorted_avg(ccep_elec,stim_pair_nr,:)))
plot(tt(tt>-1 & tt<1),squeeze(cc_epoch_sorted_avg(ccep_elec,stim_pair_nr,(tt>-1 & tt<1))),'Color',[0 0 0],'LineWidth',2)
xlabel('time(s)')
% set(gca,'Ydir','reverse')
ylabel('amplitude(uV)')
title(['elec ' data_hdr.label{ccep_elec} ' for stimulation of ' data_hdr.label{cc_stimsets(stim_pair_nr,1)} ' and ' data_hdr.label{cc_stimsets(stim_pair_nr,2)} ])
set(gca,'XTick',[0:.2:.4],'YTick',[-1000 0 1000],'FontName','arial','FontSize',14)
box off
xlim([-.1 .6]),ylim([-1200 1200])
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','/home/jaap/figures/example_ccep')

%% Plot to visualise parameters of detect ER algorithm

ii = 3;
jj = 1;

baseline_tt = tt>-1 & tt<-.01;
signal_median = median(cc_epoch_sorted_avg(ii,jj,baseline_tt),3);

% subtract median baseline from signal
new_signal = squeeze(cc_epoch_sorted_avg(ii,jj,:)) - signal_median;
% testplot new signal: plot(tt,squeeze(new_signal))

% take area before the stimulation of the new signal and calculate its SD
pre_stim_sd_orig = std(new_signal(baseline_tt));

% if the pre_stim_sd is smaller that the minimally needed SD, use this
% the minSD as pre_stim_sd
if pre_stim_sd_orig < minSD
    pre_stim_sd = minSD;
else 
    pre_stim_sd =  pre_stim_sd_orig;
end

% display figure in specific size
figure('Position',[0 0 500 500]), hold on


% plot relevant part of the signal (1s before till 500ms after stimulation)

plot(tt(tt>-1 & tt<.5),squeeze(new_signal(tt>-1 & tt<.5)),'Color',[0 0 0],'LineWidth',2)

xlabel('time(s)')
ylabel('amplitude(uV)')
set(gca, 'FontName','arial','FontSize',20)

ylim([-1200 1200])
xlim([-.5 .5])
box off

hold on
% plot stimulation artifact (first 20 samples, +/- 10ms)
plot(tt(5120:5140),squeeze(new_signal(5120:5140)),'r','LineWidth',2)
% plot found peaks and onsets
% plot(tt(output_ER_all(ii,jj,1)),output_ER_all(ii,jj,2),'b*')
% plot(tt(output_ER_all(ii,jj,3)),output_ER_all(ii,jj,4),'r*')
% plot(tt(output_ER_all(ii,jj,5)),output_ER_all(ii,jj,6),'g*')
% plot(tt(output_ER_all(ii,jj,7)),output_ER_all(ii,jj,8),'y*')

% plot calculated baseline standard deviation
plot(tt(baseline_tt), pre_stim_sd_orig+zeros(size(tt(baseline_tt))), 'r-','LineWidth',2)
plot(tt(baseline_tt), -pre_stim_sd_orig+zeros(size(tt(baseline_tt))), 'r-','LineWidth',2)

% plot adjusted baseline (when calculated < minSD)
plot(tt(baseline_tt), (pre_stim_sd*2.5)+zeros(size(tt(baseline_tt))), 'g-','LineWidth',2)
plot(tt(baseline_tt), (-pre_stim_sd*2.5)+zeros(size(tt(baseline_tt))), 'g-','LineWidth',2)

% plot lines with timeframe of n1 peak
plot([tt(n1_samples_start) tt(n1_samples_start)],[-2000 2000], 'c-','LineWidth',2)
plot([tt(n1_samples_end) tt(n1_samples_end)],[-2000 2000], 'c-','LineWidth',2)

hold off


%% EVERYTHING FOR CONNECTIVITY MATRICES 

% load electrodes file labeled with atlases
electrodes_file_atlases = fullfile(top_path,['sub-' database(subj).metadata(runs).subject], ...
    ['ses-' database(subj).metadata(runs).session],'ieeg',...
    ['sub-' database(subj).metadata(runs).subject '_ses-'  ...
    database(subj).metadata(runs).session '_electrode_positions_fouratlases.tsv']);

electrodes_table_atlases = readtable(electrodes_file_atlases,'Filetype','text','Delimiter','\t');


%% Plot figure that displays the amount of electrodes on an area

% extract labels from table and select only the ones that are not NaN
destrieux_labels = str2double(electrodes_table_atlases.Destrieux_label);
destrieux_labels = destrieux_labels(~isnan(destrieux_labels));

% plot histogram of represented electrodes and their amount
figure,hist(destrieux_labels,[1:74])
colormap(hot)
xlim([1 75])

%% Plot connectivity matrix

% how many connections to each electrode would be possible 
possible_connect = zeros(74,74);
for kk = 1:size(destrieux_labels,1)

    % add every connection the the my_connect matrix
    possible_connect(destrieux_labels(kk),destrieux_labels(setdiff(1:size(destrieux_labels,1),kk))) = ...
        possible_connect(destrieux_labels(kk),destrieux_labels(setdiff(1:size(destrieux_labels,1),kk))) + 1;
end

figure,
imagesc(possible_connect,[0 10])


axis square
grid on

set(gca,'XMinorTick','on','GridColor', [0.9, 0.9, 0.9])
% grid on
colormap(hot)
