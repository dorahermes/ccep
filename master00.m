%% Calls complete pipeline for all subjects.
%% First preallocate the net degree coefficient arrays. Those are needed for the results table at the end.
tic

netDegifg = NaN(5,1);
netDegsmg = NaN(5,1);
netDegag = NaN(5,1);
netDegstg = NaN(5,1);
netDegmtg = NaN(5,1);
netDegpla = NaN(5,1);

elsifg = NaN(5,1);
elspla = NaN(5,1);
clusters = NaN(5,5);

%% Pipeline
% aka 6 months of work in 5 lines of code...
for s_nr = 1:5
    master01_preprocess;
    master02_fit;
    master03_connectivity;
    clusters(s_nr,:) = cat(2,clusterifg, clustersmg, clusterag, clusterstg, clustermtg);
end

%% Log10 the coefficients and make a table


logDegrAL = log10(netDegifg);
logDegrPL = log10(netDegpla);
logDegrSMG = log10(netDegsmg);
logDegrAG = log10(netDegag);
logDegrSTG = log10(netDegstg);
logDegrMTG = log10(netDegmtg);


netDegrCoeff = table(logDegrAL, elsifg, logDegrPL, elspla, logDegrSMG, logDegrAG, logDegrSTG, logDegrMTG,... 
    'VariableNames', {'IFG', 'els_IFG', 'PLA', 'els_PLA', 'SMG', 'AG', 'STG', 'MTG'},...
    'RowNames', subjects ) %#ok  
writetable(netDegrCoeff, 'degreeCoefficient.xlsx');
%%
save('m-files/detectN1/data/net_degrees',...
    'netDegifg', 'netDegpla', 'netDegsmg', 'netDegag', 'netDegstg', 'netDegmtg', 'clusters')
toc

%% Effective connectivity mapping

for s_nr = 1:5
    ccep_effectiveConn;
%     plot_adjMats;
end

ampDistr = table(amp_ifg2smg', amp_smg2ifg', amp_ifg2ag', amp_ag2ifg',...
    amp_ifg2stg', amp_stg2ifg', amp_ifg2mtg', amp_mtg2ifg',...
    amp_ifg2ifg', amp_smg2smg', amp_ag2ag', amp_stg2stg', amp_mtg2mtg',...
    'VariableNames', {'ifg2smg', 'smg2ifg', 'ifg2ag', 'ag2ifg', 'ifg2stg', 'stg2ifg', 'ifg2mtg', 'mtg2ifg'...
    'ifg2ifg', 'smg2smg', 'ag2ag', 'stg2stg', 'mtg2mtg'},...
    'RowNames', subjects) %#ok
writetable(ampDistr, 'ampDistrMED.xlsx');

widthDistr = table(width_ifg2smg', width_smg2ifg', width_ifg2ag', width_ag2ifg',...
    width_ifg2stg', width_stg2ifg', width_ifg2mtg', width_mtg2ifg',...
    width_ifg2ifg', width_smg2smg', width_ag2ag', width_stg2stg', width_mtg2mtg',...
    'VariableNames', {'ifg2smg', 'smg2ifg', 'ifg2ag', 'ag2ifg', 'ifg2stg', 'stg2ifg', 'ifg2mtg', 'mtg2ifg'...
    'ifg2ifg', 'smg2smg', 'ag2ag', 'stg2stg', 'mtg2mtg'},...
    'RowNames', subjects) %#ok
writetable(widthDistr, 'widthDistrMED.xlsx');

latDistr = table(lat_ifg2smg', lat_smg2ifg', lat_ifg2ag', lat_ag2ifg',...
    lat_ifg2stg', lat_stg2ifg', lat_ifg2mtg', lat_mtg2ifg',...
    lat_ifg2ifg', lat_smg2smg', lat_ag2ag', lat_stg2stg', lat_mtg2mtg',...
    'VariableNames', {'ifg2smg', 'smg2ifg', 'ifg2ag', 'ag2ifg', 'ifg2stg', 'stg2ifg', 'ifg2mtg', 'mtg2ifg'...
    'ifg2ifg', 'smg2smg', 'ag2ag', 'stg2stg', 'mtg2mtg'},...
    'RowNames', subjects) %#ok
writetable(latDistr, 'latDistrMED.xlsx');

%% How often is the model better at predicting the data, compared to the actual data it was trained on?

errors = NaN(5,3);

load('/Fridge/users/albena/m-files/detectN1/data/goodness31.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness78.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness81.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness88.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness99.mat')

for s_nr = 1:5
    switch s_nr
        case 1
            goodness_mat = goodness_mat_99;
            clear goodness_mat_99
        case 2
            goodness_mat = goodness_mat_31;
            clear goodness_mat_31
        case 3
            goodness_mat = goodness_mat_78;
            clear goodness_mat_78
        case 4
            goodness_mat = goodness_mat_81;
            clear goodness_mat_81
        case 5
            goodness_mat = goodness_mat_88;
            clear goodness_mat_88
    end
    
    [IFG, SMG, AG, STG, MTG, els2check, esmPlus, ~] = ccep_anat_locs(s);
    matcheck = squeeze(goodness_mat(:,:,1));
    matcheck(matcheck>1) = 0;
    matcheck(matcheck<1 & matcheck>0) = 1;
    
    errors(s_nr,1) = nansum(nansum(matcheck,2)); % number of relRMSE<1 for each channel
    
    matcheck(isnan(matcheck)==0) = 1;
    errors(s_nr,2) = nansum(nansum(matcheck,2)); % total number of valid fits 
    errors(s_nr,3) = 100./errors(s_nr,2).*errors(s_nr,1);
    clear matcheck
end
%% Plot some classified CCEPs on the cortex

[IFG, SMG, AG, STG, MTG, els2check, esmPlus, ~] = ccep_anat_locs(s);

% get a nice alpha value string with no zeros before the decimal point
alphaVal = 0.05/numel(els2check);
alphaVal = round(alphaVal,3);
dispAlpha = num2str(alphaVal);
dispAlpha = strsplit(dispAlpha, '.');
dispAlpha = dispAlpha{2};
dispAlpha = strcat(['.' dispAlpha]);

% pick el
el = 1;

a = min(elecmatrix(:,1,1));
elmat = elecmatrix;
elmat(:,1) = a-1;
clear a

figure;
set(gcf,'PaperPositionMode', 'auto')
n1_p = n1_data(el).ttest(:,2);
% amplitude
n1_amp = zeros(size(n1_data,1),1);
n1_lat = zeros(size(n1_data,1),1);
n1_amp(els2check,1) = median(n1_data(el).fit(els2check,:,5),2);
% latency
n1_lat(els2check,1) = mean(n1_data(el).fit(els2check,:,2),2);
% alpha Bonferroni corrected
n1_amp(n1_p>=alphaVal) = 0; % if its not significant, we're setting it to zero
% set as 0 if we're not fitting a real n1
n1_fit = goodnessOfFitClassified(el,:);
n1_amp((n1_fit~=1),1) = NaN;

% render
ctmr_gauss_plot(cortex,[0 0 0],0)
el_add(elecmatrix, 'k', 7)
el_add(elmat(ccep(el).els,:),[.9 .9 .9],20)
el_add_sc(elmat(1:numel(n1_data),:),n1_amp,n1_lat,200,40)

title({['significant CCEPs when ' int2str(el) ' was stimulated']...
    ['p < ' dispAlpha ' (Bonferroni)']})

%% Do a chi-square test across subject to see if the projection targets and sources of influence are random.
% a matrix to store the data
% 1st dim - subject number
% 2nd dim - (1) - projection column
%           (2) - integration column
% 3rd dim - (1) - chi-statistic
%           (2) - p-value
% 4th dim - (1) - ifg
%           (2) - smg
%           (3) - ag
%           (4) - stg
%           (5) - mtg
chiData = NaN(5,2,2,5);
connectionsData = NaN(5,2,5);
% 1 subject number; 2 - projs and integr, 3 - areas
tot_ifg = zeros(2,5);
tot_smg = zeros(2,5);
tot_ag = zeros(2,5);
tot_stg = zeros(2,5);
tot_mtg = zeros(2,5);

for s_nr = 1:5
    ccep_connections_chi;
    chiData(s_nr,:,:,1) = ifg_chi;
    chiData(s_nr,:,:,2) = smg_chi;
    chiData(s_nr,:,:,3) = ag_chi;
    chiData(s_nr,:,:,4) = stg_chi;
    chiData(s_nr,:,:,5) = mtg_chi;
    
    tot_ifg = bsxfun(@plus, tot_ifg, ifg_vec);
    tot_smg = bsxfun(@plus, tot_smg, smg_vec);
    tot_ag = bsxfun(@plus, tot_ag, ag_vec);
    tot_stg = bsxfun(@plus, tot_stg, smg_vec);
    tot_mtg = bsxfun(@plus, tot_mtg, mtg_vec);
end

% chi-square test across subjects
% tot_ifg_chi(1,:) = ccep_chiTest(tot_ifg(1,:));
% tot_ifg_chi(2,:) = ccep_chiTest(tot_ifg(2,:));
% 
% tot_smg_chi(1,:) = ccep_chiTest(tot_smg(1,:));
% tot_smg_chi(2,:) = ccep_chiTest(tot_smg(2,:));
% 
% tot_ag_chi(1,:) = ccep_chiTest(tot_ag(1,:));
% tot_ag_chi(2,:) = ccep_chiTest(tot_ag(2,:));
% 
% tot_stg_chi(1,:) = ccep_chiTest(tot_stg(1,:));
% tot_stg_chi(2,:) = ccep_chiTest(tot_stg(2,:));
% 
% tot_mtg_chi(1,:) = ccep_chiTest(tot_mtg(1,:));
% tot_mtg_chi(2,:) = ccep_chiTest(tot_mtg(2,:));

%%

%% Colormap for areas:

% This part is also used for the cortex plots
fcifg = [0 0 0; .9 .9 .3];
fcsmg = [0 0 0; .9 .3 .3];
fcag = [0 0 0; .3 .9 .3];
fcstg = [0 0 0; .3 .3 .9];
fcmtg = [0 0 0; .7 .3 .7];
cxm = [fcifg(2,:); fcsmg(2,:); fcag(2,:); fcstg(2,:); fcmtg(2,:)];
% Matrix to plot the colormap:
a = [1:5]';


figure
subplot(1,5,1)
bar(tot_ifg)
set(gca,'XTickLabel', {'projections', 'sources of influence'})
ylabel('number of projections')
title(['IFG projections'])

subplot(1,5,3)
bar(tot_ag)
set(gca,'XTickLabel', {'projections', 'sources of influence'})
ylabel('number of projections')
title(['AG projections'])

subplot(1,5,2)
bar(tot_smg)
set(gca,'XTickLabel', {'projections', 'sources of influence'})
ylabel('number of projections')
title(['SMG projections'])

subplot(1,5,4)
bar(tot_stg)
set(gca,'XTickLabel', {'projections', 'sources of influence'})
ylabel('number of projections')
title(['STG projections'])

subplot(1,5,5)
bar(tot_mtg)
set(gca,'XTickLabel', {'projections', 'sources of influence'})
ylabel('number of projections')
title(['MTG projections'])
















