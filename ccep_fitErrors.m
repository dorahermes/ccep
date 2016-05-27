%% Analyze rRMSE errors for all fits
errorData = struct('subject', [], 'errors', [], 'total', [],...
    'good', [], 'gperc', [], 'ok', [], 'okperc', [], 'bad', [], 'bperc', []);

for s_nr = 1:5
    errorData(s_nr).subject = s_nr;
    
    switch s_nr
        case 1
            goodness_mat = goodness_mat_99;
%             clear goodness_mat_99
        case 2
            goodness_mat = goodness_mat_31;
%             clear goodness_mat_31
        case 3
            goodness_mat = goodness_mat_78;
%             clear goodness_mat_78
        case 4
            goodness_mat = goodness_mat_81;
%             clear goodness_mat_81
        case 5
            goodness_mat = goodness_mat_88;
%             clear goodness_mat_88
    end
    
    
    % (1)Get error data
    errorData(s_nr).errors = squeeze(goodness_mat(:,:,1));

    % (2)Total number of fits:
    errorData(s_nr).total = numel(find(errorData(s_nr).errors() > 0)); 
    
    % (3)Number of relRMSE<1 for each channel:
    errorData(s_nr).good = numel(find(errorData(s_nr).errors > 0 & errorData(s_nr).errors < 1)); 
    % (4)Percent good fits
    errorData(s_nr).gperc = 100./(errorData(s_nr).total).*(errorData(s_nr).good);
    
    % (5)Number of relRMSE > 1 & < 2.5
    errorData(s_nr).ok = numel(find(errorData(s_nr).errors > 1 & errorData(s_nr).errors < 2.5)); 
    % (6)Percent ok fits
    errorData(s_nr).okperc = 100./(errorData(s_nr).total).*(errorData(s_nr).ok);
    
    % (7)Number of relRMSE > 2.5. This is a condition to classify an N1.
    % Those bad fits have not been analyzed further
    errorData(s_nr).bad = numel(find(errorData(s_nr).errors > 2.5)); 
    % (8)Percent bad fits
    errorData(s_nr).bperc = 100./(errorData(s_nr).total).*(errorData(s_nr).bad);
    
    %%% Plot errors in a matrix %%%
    
    % First get subject el data
    subjects = {'99','31','78','81','88'};
    s = subjects{s_nr};
    [IFG, SMG, AG, STG, MTG, els2check, esmPlus, ~] = ccep_anat_locs(s);
    
    % plot only for els2check:
    firstEl = min(els2check); 
    lastEl = max(els2check); 
    
    figure
    imagesc(errorData(s_nr).errors, [0 3])
    colorbar('vertical')
    colormap(gray)
    title(['rRMSE subject ' int2str(s_nr)])
    xlabel('measured el')
    ylabel('stimulated el')
    axHands = findobj(get(gcf, 'Children'), 'flat', 'type', 'axes');
    axis(axHands, 'square', [firstEl lastEl firstEl lastEl])
    set(axHands, 'XAxisLocation', 'origin')    
    print(gcf, '-dpng', strcat(['figures/pres_figures/errormat' int2str(s_nr)]));
    %%%
    close all
    
    %%% Also plot error pie charts
    piemap = [.5 .9 .5; .5 .7 .9; .9 .5 .5];
     
    data2pie = cat(2,errorData(s_nr).good, errorData(s_nr).ok, errorData(s_nr).bad);

    text2plot = {strcat([int2str(errorData(s_nr).good)]),...
                 strcat([int2str(errorData(s_nr).ok)]),...
                 strcat([int2str(errorData(s_nr).bad)]),...
                };
    figure('Position', [680 521 511 573])
    pieHandle = pie(data2pie, text2plot);
    for iHandle = 2:2:2*numel(text2plot)
        pieHandle(iHandle).Position = 0.7*pieHandle(iHandle).Position;
    end
    
    colormap(piemap)
    colorbar('horizontal', 'Ticks', [1.3, 2, 2.65], 'Ticklabels', {'rRMSE < 1', '1 < rRMSE < 2.5', 'rRMSE > 2.5'})
    title('rRMSE in number of fitted CCEPs')
    print(gcf, '-dpng', strcat(['figures/pres_figures/errorPie' int2str(s_nr)]));
    
    close all
end

%%
errors2plot = zeros(1,1);
for s_nr = 1:5
    errors2plot = cat(2,errors2plot,reshape(errorData(s_nr).errors,1,numel(errorData(s_nr).errors)));
end
%%


data2plot = zeros(4,5);

for s_nr = 1:5
    data2plot(1,s_nr) = errorData(s_nr).total;
    data2plot(2,s_nr) = errorData(s_nr).good;
    data2plot(3,s_nr) = errorData(s_nr).ok;
    data2plot(4,s_nr) = errorData(s_nr).bad;
end
%%
errormap = [.8 .8 .8; .5 .9 .5; .5 .7 .9; .9 .5 .5];

figure
bar(data2plot', 1)
colormap(errormap)
ylabel('number of CCEP fits')
xlabel('subject')
title('Model performance - rRMSE')
axis tight

%%
figure
a = [1:4];
imagesc(a)
colormap(errormap)
set(gca,'XTick',1:5,...
        'XTickLabel',{'CCEPs total', 'rRMSE < 1', '1 < rRMSE < 2.5', 'rRMSE > 2.5'},...
        'YTickLabel', [],...
        'TickLength',[0 0]);


%% Pairwise ttests
[~, good2ok_p, good2ok_ci, good2ok_stats] = ttest(gperc,okperc,0.5);

[~, ok2bad_p, ok2bad_ci, ok2bad_stats] = ttest(okperc,bperc,0.5);

[~, good2bad_p, good2bad_ci, good2bad_stats] = ttest(gperc,bperc,0.5);

 %% Plot Spearman's correlation matrices
subjects = {'99','31','78','81','88'};
total_corrs = 0;
s_low = 0;
s_moderate = 0;
s_high = 0;
 for s_nr = 1:5
s = subjects{s_nr}; % specify subject
s_info = subj_info_ccep(s); % get subject info

eval(['ccep = subj_info_SPES_' s]); % make CCEP structure
 switch s_nr
     case 1
         goodness_mat = goodness_mat_99;
         t_mat = NaN(size(goodness_mat,1),size(goodness_mat,1));
         for k = els(s_nr).els
             if isempty(ccep(k).els)==0
              t_mat(k,:) = n1_99(k).ttest(:,1);
             end
         end
     case 2
         goodness_mat = goodness_mat_31;
         t_mat = NaN(size(goodness_mat,1),size(goodness_mat,1));
         for k = els(s_nr).els
             if isempty(ccep(k).els)==0
               t_mat(k,:) = n1_31(k).ttest(:,1);
             end
         end
     case 3
         goodness_mat = goodness_mat_78;
         t_mat = NaN(size(goodness_mat,1),size(goodness_mat,1));
         for k = els(s_nr).els
             if isempty(ccep(k).els)==0
               t_mat(k,:) = n1_78(k).ttest(:,1);
             end
         end
     case 4
         goodness_mat = goodness_mat_81;
         t_mat = NaN(size(goodness_mat,1),size(goodness_mat,1));
         for k = els(s_nr).els
             if isempty(ccep(k).els)==0
               t_mat(k,:) = n1_81(k).ttest(:,1);
             end
         end
     case 6
         goodness_mat = goodness_mat_88;
         t_mat = NaN(size(goodness_mat,1),size(goodness_mat,1));
         for k = els(s_nr).els
             if isempty(ccep(k).els)==0
               t_mat(k,:) = n1_88(k).ttest(:,1);
             end
         end
 end

total_corrs =  total_corrs + numel(find(~isnan(goodness_mat(:,:,2))));
s_low = s_low + numel(find(goodness_mat(:,:,2)<0.4));
s_moderate = s_moderate + numel(find(goodness_mat(:,:,2) >= 0.4 & goodness_mat(:,:,2) < 0.7));
s_high = s_high + numel(find(goodness_mat(:,:,2) >= 0.7));

firstEl = min(els(s_nr).els); 
lastEl = max(els(s_nr).els); 

figure
imagesc(goodness_mat(:,:,2), [-1 1])
colormap(gray)
colorbar('vertical')
title(['Correlation matrix subject ' int2str(s_nr)])
xlabel('measured el')
ylabel('stimulated el')
axHands = findobj(get(gcf, 'Children'), 'flat', 'type', 'axes');
axis(axHands, 'square', [firstEl lastEl firstEl lastEl])
set(axHands, 'XAxisLocation', 'origin')

print(gcf, '-dpng', strcat(['figures/pres_figures/corrMats' int2str(s_nr)]))
close all
figure
imagesc(t_mat)
colormap(gray)
colorbar('vertical')
title(['T-values for response amplitude ' int2str(s_nr)])
xlabel('measured el')
ylabel('stimulated el')
axHands = findobj(get(gcf, 'Children'), 'flat', 'type', 'axes');
axis(axHands, 'square', [firstEl lastEl firstEl lastEl])
set(axHands, 'XAxisLocation', 'origin')

print(gcf, '-dpng', strcat(['figures/pres_figures/tMats' int2str(s_nr)]))
end
 
 
 
