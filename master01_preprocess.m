
function_path = '/Volumes/DoraBigDrive/data/SPES/m-files/';
addpath(function_path);

clear all


%identify the list of subjects
subjects = {'79'};
s = 1; % subject to analyse
s_nr = subjects{s}; 

%% Load data

s_info      = subj_info_ccep(s_nr);
dataStruct  = load(s_info.data,'data','times','srate');
load(s_info.labels_Wang)
load(s_info.cortex);
load(s_info.els);

if isequal(s_nr,'31')
    elecmatrix = elecmatrix(1:96,:); % last electrodes are NaN anyways, and not in data
elseif isequal(s_nr,'88')
    elecmatrix = elecmatrix(1:80,:); % last electrodes are NaN anyways, and not in data
end

data    = dataStruct.data;
t       = dataStruct.times;
srate   = dataStruct.srate;

% get the stimulation info for subject s
eval(['ccep = subj_info_SPES_' s_nr '()']);

clear dataStruct

disp('data loaded')

%% Detect outliers based on jumps and voltage threshold works well
outlier_matrix = ccep_detectOutliers(data,ccep,t,1);

%% Remove the outliers
% outliers are replaced with NaNs

data = ccep_clearOut(data,outlier_matrix);
disp('outliers removed')

%% Weighted common average reference

dataCAR = zeros(size(data,1),size(data,2),size(data,3));
for k=1:size(data,3)
    if mod(k,5)==0, disp(['epoch ' int2str(k) ' of ' int2str(size(data,3))]),end
    % select the channels to include in the common average
    chans2inc = s_info.chans;   
    % detect the channels that were stimulated (much different mean signal)
    SignalMean = nanmean(abs(data(:,:,k)),2);   
    threshold2exclude = nanmean(SignalMean)+3*nanstd(SignalMean);
    chans2inc(SignalMean>threshold2exclude)=[];
    
    dataCAR(:,:,k) = carreg_Filt(data(:,:,k), chans2inc);
end
clear SignalMean threshold2exclude

%% Correction for an offset
% choose baseline or median correction

% Baseline correction:
% subtracts for every single trial, the average from a baseline period
data = ccep_baseline_corr(data,t);

% Median correction
% corrects for a constant difference with the median response
% data = ccep_median_corr(data,t,ccep,[10 600]);


%% Plot electrode grid
figure
ctmr_gauss_plot(cortex,[0 0 0],0)
label_add(elecmatrix)


%% %%%%%%%
%% %%%%%%% Done with the preprocessing
%% %%%%%%% Look at some raw CCEPs
%% %%%%%%%


%% Plot all CCEPs for chosen electrode

el = 71; % choose electrode for start
data2plot = nanmean(data(:,:,ccep(el).epochs),3); % take the mean across epochs
t_plot = t(t>10 & t<500); % pick the times in ms that I want to plot
data_plot = data2plot(:,t>10 & t<500); % pick the same times for the data

% make a figure
figure('Position', [0 0 800 600]);

for k=1:size(data_plot,1) % loop across electrodes
    subplot(ceil(size(data2plot,1)/10),10,k),hold on % make a subplot, we need 8x8 for 64 electrodes, more if we have more electrodes
    plot(t_plot,zeros(size(t_plot)),'k:') % plot a zero-line
    plot(t_plot,data_plot(k,:)) % plot the data as a function of time
    ylim([-700 300])
    set(gca,'YTick',[],'XTick',[])
    ylabel([int2str(k)])
end


%% Plot CCEP for a specific pair
el = 71;
elm = 69;

figure('Position', [0 0 300 300])
set(gcf,'PaperPositionMode', 'auto')

subplot(2,1,1),hold on

cm = colormap(jet);
for k=1:length(ccep(el).epochs)
    plot(t,data(elm,:,ccep(el).epochs(k)),...
        'LineWidth',1,...
        'Color',cm(round(k*size(cm,1)/length(ccep(el).epochs)),:))
end
plot(t,squeeze(mean(data(elm,:,ccep(el).epochs),3)),'k',...
'LineWidth',2)

title({['CCEP measured el ' int2str(elm)],...
 [' stim el ' int2str(el)]})

% ylim([-900 300])  % max voltage amps differ 
xlim([-100 500])

xlabel('ms')
ylabel('microV')

subplot(2,1,2),hold on

cm = colormap(jet);
for k=1:length(ccep(el).epochs)
    plot(t,dataCAR(elm,:,ccep(el).epochs(k)),...
        'LineWidth',1,...
        'Color',cm(round(k*size(cm,1)/length(ccep(el).epochs)),:))
end
plot(t,squeeze(mean(dataCAR(elm,:,ccep(el).epochs),3)),'k',...
'LineWidth',2)

title({['CAR CCEP measured el ' int2str(elm)],...
 [' stim el ' int2str(el)]})

% ylim([-900 300])  % max voltage amps differ 
xlim([-100 500])

xlabel('ms')
ylabel('microV')

