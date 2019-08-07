


%%

% why does this give all nan's? - should work for matrices as well 
% [B,BINT,R,RINT,STATS] = regress(age_vector',ROI_within_plot_matrix_all') % ones(length(ROI_within_plot_matrix_all),1)]);



%% Analysis for regressions
% list of ROI_within matrices that can be used for plotting:
% - ROI_plot_matrix_15
% - ROI_plot_matrix_26
% - ROI_plot_matrix_38
% - ROI_within_plot_matrix_all
% or for between ROI analysis:
% - ROI_between_plot_matrix
% or for all found CCEPs
% - matrix_reshape_all

% FILL IN THE ROI HERE - CHOOSE FROM ABOVE
scatter_ROI = matrix_reshape_all;

% calculate means
trendline_ROI = nanmean(scatter_ROI);

% do regress analysis using age vector and the nanmeans latency
% stats: R2 statistic, the F-statistic and its p-value, and an estimate of the error variance
[B,BINT,R,RINT,STATS] = regress(age_vector',[trendline_ROI' ones(length(trendline_ROI),1)]);

% Durbin-Watson test for independent errors, significant due to set up data
% to do DW test, remove nan's from data
trendline_ROI_nonan = trendline_ROI(~isnan(trendline_ROI));  
R_nonan = R(~isnan(R)); 
% test only significant for 'all' 
[p_dw, d_dw] = dwtest(R_nonan,[trendline_ROI_nonan' ones(length(trendline_ROI_nonan),1)])

% degrees of freedom for regression
df = [1  (length(trendline_ROI_nonan)-2)];   


% can be used for all ROIs together 
% [b,stats]  = robustfit(trendline_ROI',age_vector')

% plot all data and draw regression line

figure(), hold on
% iterate over all subjects in database
for subj = 1:length(database)
    
    % create age-vector for plotting
    age_vector(1,subj) = database(subj).age_ses1;
    
    % scatterplot
    scatter(repelem(age_vector(1,subj),length(scatter_ROI)),scatter_ROI(:,subj))
    
end

xlim([0 55]),ylim([0 100])
xlabel('age subject')
ylabel('latency in ms')
title('latency between the three ROIs')

x = [6:1:50];
y = B(1) * x + B(2); % for robust fit: y = b(2) * x + b(1);
plot(x,y)
hold off




%% group level analysis 

% list of ROI_within matrices that can be used for plotting:
% - ROI_plot_matrix_15
% - ROI_plot_matrix_26
% - ROI_plot_matrix_38
% - ROI_within_plot_matrix_all
% or for between ROI analysis:
% - ROI_between_plot_matrix
% or for all found CCEPs
% - matrix_reshape_all

% FILL IN THE ROI HERE - CHOOSE FROM ABOVE
grouped_scatter_ROI = matrix_reshape_all;

grouped_age_ROI = nanmean(grouped_scatter_ROI);


% normality tests - ALL NOT NORMALLY DISTRIBUTED
H = histfit(grouped_age_ROI(:));
[h_ks,p_ks,ksstat,cv] = kstest(grouped_age_ROI);

for gg = 1:length(age_vector)
    if age_vector(gg) <= 18
        age_group(gg) = 1;
    elseif age_vector(gg) > 18
        age_group(gg) = 2;
        
    end
end

% levene test for heterogeneity of variance - no sig results on means
p_levene = vartestn(grouped_age_ROI',age_group','TestType','LeveneAbsolute')

% create vector with data for both groups
for qq = 1:length(age_group)

    if age_group(qq) == 1
    
        group_1_vec(qq) = grouped_age_ROI(qq);
        group_2_vec(qq) = NaN;
     
    elseif age_group(qq) == 2
        
        group_1_vec(qq) = NaN;
        group_2_vec(qq) = grouped_age_ROI(qq);
        
    end
    
end

% delete nan's just to be sure
group_1_vec = group_1_vec(~isnan(group_1_vec))
group_2_vec = group_2_vec(~isnan(group_2_vec))

% https://nl.mathworks.com/help/stats/ranksum.html - still have to try one
% tailed, but all significant except for 26. 
% Wilcoxon rank-sum/Mann-Whitney
[p,h,stats] = ranksum(group_1_vec',group_2_vec','tail', 'right', 'method', 'approximate')

% calculate median + effect size for report
median(group_1_vec)
median(group_2_vec)
r = stats.zval / (length(group_1_vec) + length(group_2_vec))
% robust bootstrapping - also option 
