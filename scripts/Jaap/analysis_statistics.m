
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


% FILL IN THE ROI HERE (as variable for scatter_ROI) - CHOOSE FROM ABOVE
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
plot(x,y, 'r')
hold off


%%%% Anovan

y = trendline_ROI';
g1 = age_group; 
g2 = age_vector; 

p = anovan(y,{g2,g1},'varnames',{'age','subject nr'},'nested',[0 0; 1 0]) 

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

% FILL IN THE ROI HERE (as variable for grouped_ROI) - CHOOSE FROM ABOVE
grouped_ROI = ROI_between_plot_matrix;

grouped_age_ROI = nanmean(grouped_ROI);


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
r = stats.zval / sqrt((length(group_1_vec) + length(group_2_vec)))



% save('matfiles_for_dora','ROI_plot_matrix_15', 'ROI_plot_matrix_26','ROI_plot_matrix_38','ROI_within_plot_matrix_all','ROI_between_plot_matrix' ...
%     ,'matrix_reshape_all','age_vector','-v7.3')


%% ROI analysis 

% create nan matrix to fill in single point calculated AUCs 
AUC_singlepoint = nan(size(averaged_parameter_scores,1),size(averaged_parameter_scores,2));

% for every row (amplitude threshold)
for aa = 1:size(averaged_parameter_scores,1)
    
    % for every column (endpoints)
    for bb = 1:size(averaged_parameter_scores,2)
        
        % calculate the singlepoint based AUC and put into AUC_singlepoint 
        % matrix by calculating squared part + lower left triangle + upper 
        % right triangle, using with 1 - specificity for x-coordinate
        AUC_singlepoint(aa,bb) = (averaged_parameter_scores(aa,bb,1) * averaged_parameter_scores(aa,bb,2))  +  ...
            (((1 - averaged_parameter_scores(aa,bb,2)) * averaged_parameter_scores(aa,bb,1)) / 2)  + ...
            ((averaged_parameter_scores(aa,bb,2) * (1 - averaged_parameter_scores(aa,bb,1))) / 2);

    end
end

% not normally distributed - which makes sense for ROC curve data
for zz = 1:size(averaged_parameter_scores,2)
    
    [h_ks,p_ks,ksstat,cv] = kstest(AUC_singlepoint(:))
    
end    

% first Friendman ANOVA to test whether there are differences,
% there are differences between ROI curves, so with follow-up wilcoxon 
% to test which curves differ significantly
[p,tbl,stats] = friedman(AUC_singlepoint,1)

% choose the column with the highest meanrank, which is colum 7, endpoint
% of 100 ms

% use Wilcoxon Signed Rank Test to see if this parameter is significantly
% better than others. Now compared with the neighbouring parameters (90ms
% and 110ms, but also significant for all other comparisons)
[p,h,stats] = signrank(AUC_singlepoint(:,7)',AUC_singlepoint(:,8)','tail','right')

r = stats.zval / sqrt((length(AUC_singlepoint) * 2))

[p,h,stats] = signrank(AUC_singlepoint(:,7)',AUC_singlepoint(:,6)','tail','right')

r = stats.zval / sqrt(length(AUC_singlepoint) * 2)

%% analyzing the percentages of CCEPs - grouped for within and between

rel_perc_cceps_within % percentages of CCEPs within ROIs
rel_perc_cceps_between % percentages of CCEPs between ROIs

% violation of normality assumptions
H = histfit(rel_perc_cceps_within(:),10,'kernel');
[h_ks,p_ks,ksstat,cv] = kstest(rel_perc_cceps_within)
[h_ks,p_ks,ksstat,cv] = kstest(rel_perc_cceps_between)

% no violation of heterogeinity
p_levene = vartestn([rel_perc_cceps_within rel_perc_cceps_between]',[ones(1,27) ones(1,27)+1]','TestType','LeveneAbsolute')

% rank-sum test
[p,h,stats] = ranksum(rel_perc_cceps_within',rel_perc_cceps_between','tail', 'right' )

% effect size
r = stats.zval / sqrt((length(rel_perc_cceps_within) + length(rel_perc_cceps_between)))

median(rel_perc_cceps_within,'omitnan') % check on other medians
median(rel_perc_cceps_between,'omitnan') % check on other medians

%% analyzing relative CCEPs - regression

% do some recalculation to add within data and between data 
for subj = 1:length(subjects)
    
     if ~isempty(database(subj).ROI_between_all)
         
         total_relative_cceps(subj) = ((database(subj).amount_cceps) + (sum(~isnan(database(subj).ROI_between_all(:,3))))) / ...
             ((database(subj).total_stims) + (length(database(subj).ROI_between_all(:,3))));

     elseif isempty(database(subj).ROI_between_all) 
         
         total_relative_cceps(subj) = (database(subj).amount_cceps) / (database(subj).total_stims);
     end
end

total_relative_cceps

% do regress analysis using age vector and the nanmeans latency
% stats: R2 statistic, the F-statistic and its p-value, and an estimate of the error variance
[B,BINT,R,RINT,STATS] = regress(age_vector',[total_relative_cceps' ones(length(total_relative_cceps),1)]);

% Durbin-Watson test for independent errors, significant due to set up data
% to do DW test, remove nan's from data
total_relative_cceps_nonan = total_relative_cceps(~isnan(total_relative_cceps));  
R_nonan = R(~isnan(R)); 
% not significant test dw
[p_dw, d_dw] = dwtest(R_nonan,[total_relative_cceps_nonan' ones(length(total_relative_cceps_nonan),1)])

% degrees of freedom for regression
df = [1  (length(total_relative_cceps_nonan)-2)];   



figure(),hold on

for subj = 1:length(database)
    
    scatter(age_vector(1,subj),rel_perc_cceps_within(subj),50,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1])
    
    
end


for subj = 1:length(subjects)
    
    scatter(age_vector(1,subj),rel_perc_cceps_between(subj),50,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
    
end
xlim([0 50]), ylim([0 1])
xlabel('age subject')
ylabel('latency in ms')
title('relative percentage cceps, 18- and 18+, within and between all ROIs')

hold off


%% relative CCEPs - grouped 18+ and 18- 

total_relative_cceps


% create vector with data for both groups
for qq = 1:length(age_group)

    if age_group(qq) == 1
    
        group_1_vec(qq) = total_relative_cceps(qq);
        group_2_vec(qq) = NaN;
     
    elseif age_group(qq) == 2
        
        group_1_vec(qq) = NaN;
        group_2_vec(qq) = total_relative_cceps(qq);
        
    end
    
end

% delete nan's just to be sure
group_1_vec = group_1_vec(~isnan(group_1_vec))
group_2_vec = group_2_vec(~isnan(group_2_vec))

[p,h,stats] = ranksum(group_1_vec',group_2_vec','tail', 'right')

% calculate median + effect size for report
median(group_1_vec)
median(group_2_vec)
r = stats.zval / sqrt((length(group_1_vec) + length(group_2_vec)))