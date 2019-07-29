function outlier_matrix = ccep_detectOutliers(data,ccep,times,varargin)
% 
% Detect outliers in terms of variance, standard deviations from the
% median, voltage threshold and voltage jumps. Checks between 25 ms and
% 550 ms after stimulation.
% 
% Call as detect_outliers(data,ccep,times,s), where:
%
%   'data'      is the name of the matrix containing the data
%   'ccep'      is the name of the struct that contains stimulation data
%   'times'     is the time vector
% 
% optional input
%   'makeFig'   if you want to make a figure
%   's'         if you want to save the figures with the subjects name, add s
%
% Hermes and Vassileva, 2016, UMC Utrecht

if ~isempty(varargin)
    makeFig = varargin{1};
    if length(varargin)>1
        s = varargin{2};
    end
else
    makeFig = 0;
end


%% Detect outliers
%%%%%%%%% Set thresholds to reject outliers: %%%%%%%%%%
t_start = 25; % time to start check
t_end = 550; % time to end check
volt_threshold = 2000 ;% voltage jump threshold (any direction)
jumps_threshold = 150; % jump threshold in microV

% var_th_num = 3; % variance for the entire trace
% std_th_num = 3; % standard deviation of each time point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outlier_matrix = zeros(size(data,1),size(data,3)); % preallocation
for el=1:size(data,1) % defines the stimulated electrodes as outliers already
    outlier_matrix(ccep(el).els,ccep(el).epochs)=1;
end

% c_med = 0;
% c_var = 0;
c_volt = 0;
c_jump = 0;


for el=1:size(data,1) % Check for outliers for every electrode

%     %%%%%% outliers in variance %%%%%% PROBLEM: also detects cceps
%     var_per_epoch = squeeze(var(data(el,:,:),[],2));
%     var_threshold = var_th_num*mean(var_per_epoch(outlier_matrix(el,:)==0));
%     outlier_matrix(el,:) = var_per_epoch>var_threshold;
%     c_var = c_var + sum(outlier_matrix(el,:),2);
% 
%     %%%%%% outlier threshold median %%%%%% PROBLEM: also detects cceps
%     up_th_med=squeeze(median(data(el,times>t_start & times<t_end,:),3)+std_th_num*std(data(el,times>t_start & times<t_end,outlier_matrix(el,:)==0),[],3));
%     low_th_med=squeeze(median(data(el,times>t_start & times<t_end,:),3)-std_th_num*std(data(el,times>t_start & times<t_end,outlier_matrix(el,:)==0),[],3));
%     % check > median+3std
%     a = bsxfun(@minus,data(el,times>t_start & times<t_end,:),up_th_med);
%     a = a>0; a = squeeze(a); 
%     a = sum(a,1); % sum across time: number of points where data are larger than the median + 3std 
%     outlier_matrix(el,a>0)=2;
%     c_med = c_med + sum(outlier_matrix(el,:)==2,2);
%     % check < median-3std
%     a = bsxfun(@minus,data(el,times>t_start & times<t_end,:),low_th_med);
%     a = a<0; a = squeeze(a); 
%     a = sum(a,1); % sum across time: number of points where data are smaller than the median - 3std 
%     outlier_matrix(el,a<0)=3;
%     c_med = c_med + sum(outlier_matrix(el,:)==3,2);
%     clear a

    %%%%%% outlier threshold jumps %%%%%%
    abs_diff_data = squeeze(abs(diff(data(el,times>t_start & times<t_end,:),1,2)));
    nr_jumps_per_epoch = sum(abs_diff_data > jumps_threshold,1);
    outlier_matrix(el,nr_jumps_per_epoch>0) = 4;
    c_jump =  c_jump + sum(outlier_matrix(el,:)==4,2);
    clear abs_diff_data nr_jumps_per_epoch
    
    %%%%%% outlier threshold voltage %%%%%%
    volt_thr_data = squeeze(abs(data(el,times>t_start & times<t_end,:)));
    nr_volt_out = sum(volt_thr_data >volt_threshold,1);
    outlier_matrix(el,nr_volt_out>0) = 5;
    c_volt =  c_volt + sum(outlier_matrix(el,:)==5,2);
    clear volt_thr_data nr_volt_out 
    
end
%     disp([int2str(c_var) ' epochs exceed variance threshold of ' int2str(var_th_num) ' SDs'])
%     disp([int2str(c_med) ' epochs exceed ' int2str(std_th_num) ' standard deviations from the median'])

    disp([int2str(c_volt) ' epochs exceed threshold of ' int2str(volt_threshold) ' microV'])
    disp([int2str(c_jump) ' epochs with jumps larger than ' int2str(jumps_threshold) ' microV'])

%% Print outlier grid
if makeFig==1
    figure
    set(gcf,'PaperPositionMode', 'auto');
    imagesc(outlier_matrix,[0 5])
    xlabel('epochs'),ylabel('electrodes')
    if exist('s','var')
        print(gcf,'-djpeg', '-r0', strcat(['figures/outliers/subj' s '/subject_' s '_outliers']))
    end
end

%% Plots for each electrode
% close all

% screenSize = get( groot, 'Screensize' );
% test_fig = figure('Position', screenSize);
% set(test_fig,'PaperPositionMode', 'auto');
%     for el = 1:size(outlier_matrix,1)
%         hold on
%         title({['Stimulated electrode: ' int2str(el)];...
%              ['Upper voltage threshold: ' int2str(volt_threshold)]})
% 
%         plot(times,squeeze(data(elm,:,:)), 'g')
%         plot(times,squeeze(data(elm,:,outlier_matrix(el,:)>0)), 'r')
%         xlim([-100 600])
%         ylim([-2000 2000])  
%         
%         hold off
%         print(test_fig,'-djpeg', '-r0', strcat(['figures/outliers/subj' s '/subejct_' s '_electrode_' int2str(el) '_outliers']))
%         clf(test_fig,'reset');  
%     end
% close all
end