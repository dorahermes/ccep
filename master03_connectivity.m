%% Get directed connectivity measures
%% Pick subject and load data

% pick subject:
%  s_nr = 3;

subjects = {'99','31','78','81','88'};

s = subjects{s_nr};
eval(['ccep = subj_info_SPES_' s]); % make CCEP structure
s_info = subj_info_ccep(s); % get subject info
load(s_info.cortex) % get cortex
load(s_info.els) 

% Load data:
load('/Fridge/users/albena/m-files/detectN1/data/n1_31.mat')
load('/Fridge/users/albena/m-files/detectN1/data/n1_78.mat')
load('/Fridge/users/albena/m-files/detectN1/data/n1_81.mat')
load('/Fridge/users/albena/m-files/detectN1/data/n1_88.mat')
load('/Fridge/users/albena/m-files/detectN1/data/n1_99.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness31.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness78.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness81.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness88.mat')
load('/Fridge/users/albena/m-files/detectN1/data/goodness99.mat')

%% Get data
switch s
    case '99'
        goodnessOfFitClassified = goodnessClassified_99;
        n1_data = n1_99;
    case '31'
        goodnessOfFitClassified = goodnessClassified_31;
        n1_data = n1_31;
    case '78'
        goodnessOfFitClassified = goodnessClassified_78;
        n1_data = n1_78;
    case '81'
        goodnessOfFitClassified = goodnessClassified_81;
        n1_data = n1_81;
    case '88'
        goodnessOfFitClassified = goodnessClassified_88;
        n1_data = n1_88;
    case '50'
        goodnessOfFitClassified = goodnessClassified_50;
        n1_data = n1_50;
end


%% Get anatomical locations and electrode numbers:
[IFG, SMG, STG, els2check] = ccep_anat_locs(s);

% Set alpha 
alphaVal = 0.05;

%% Create matrices for amplitude, width and latency.

% Preallocate matrix for response values
responseAmp = NaN(size(n1_data,2),size(n1_data,2));
responseWidth = NaN(size(n1_data,2),size(n1_data,2));
responseLat = NaN(size(n1_data,2),size(n1_data,2));

% Fill matrices with data for each el-elm pair
for el = 1:size(n1_data,2) % loop through stimulated electrodes
    if isempty(n1_data(el).fit())==1 % fill with NaNs if there's no fitting data when that electrode was stimulated
        responseAmp(el,:) = NaN;
        responseWidth(el,:) = NaN;
        responseLat(el,:) = NaN;
        
    elseif isempty(n1_data(el).fit())==0
        for elm = 1:size(n1_data,2) % loop through measured electrodes 
            if n1_data(el).ttest(elm,2) < alphaVal ... fill with data if t-test p-value is smaller than the preset alpha
                    && goodnessOfFitClassified(el,elm)==1
                responseAmp(el,elm) = nanmean(squeeze(n1_data(el).fit(elm,:,5)),2);
                responseWidth(el,elm) = nanmean(squeeze(n1_data(el).fit(elm,:,6)),2);
                responseLat(el,elm) = nanmean(squeeze(n1_data(el).fit(elm,:,2)),2);
            end
        end
    end
end

clear el elm

% figure
% imagesc(responseAmp)

%% Plot significant responses on cortex, scaled for amplitude
% el = 50;
% elcolor = 'w';
% 
% figure
% ctmr_gauss_plot(cortex,[0 0 0],0)
% el_add_size(elecmatrix,responseAmp, 150, elcolor)
% el_add(elecmatrix((ccep(el).els),:),[.8 0 .2],20);
% close all

%% make sparse adjacency mats

% Matrix for amplitude
adjAmp = ones(size(responseAmp));
adjAmp(isnan(responseAmp)) = 0;
adjAmpProd = adjAmp + adjAmp'; % This shows reciprocal connections as 2s. This one is symmetric

% figure
% imagesc(adjAmp)

%% Calculate number and % projections to each area for each stimulated electrode of interest

projData = struct('el', [], 'ifg',[], 'smg', [], 'stg', [], 'tot', []);
for k = 1:size(n1_data,2)
    projData(k).el = k;
end

for el = els2check
    if isempty(ccep(el).epochs) == 0
        % total number of projections to the els we're interested in:
        projData(el).tot = sum(adjAmp(el,els2check),2);
        
        %%% Projections to: %%%
        % Inferior frontal gyrus
        projData(el).ifg(1) = sum(adjAmp(el,IFG),2);
        projData(el).ifg(2) = 100./(projData(el).tot).*projData(el).ifg(1);
        % Supramarginal gyrus
        projData(el).smg(1) = sum(adjAmp(el,SMG),2);
        projData(el).smg(2) = 100./(projData(el).tot).*projData(el).smg(1);
        % Superior temporal gyrus
        projData(el).stg(1) = sum(adjAmp(el,STG),2);
        projData(el).stg(2) = 100./(projData(el).tot).*projData(el).stg(1);
    else
        continue
    end
end



%% Inputs

recData = struct('el', [], 'ifg',[], 'smg', [], 'stg', [], 'tot', []);
for k = 1:size(n1_data,2)
    recData(k).el = k;
end

for el = els2check
%     if isempty(ccep(el).epochs) == 0
        % total number of projections:
        recData(el).tot = sum(adjAmp(els2check,el),1);
        
        %%% Projections to: %%%
        % Inferior frontal gyrus
        recData(el).ifg(1) = sum(adjAmp(IFG,el),1);
        recData(el).ifg(2) = 100./(recData(el).tot).*recData(el).ifg(1);
        % Supramarginal gyrus
        recData(el).smg(1) = sum(adjAmp(SMG,el),1);
        recData(el).smg(2) = 100./(recData(el).tot).*recData(el).smg(1);
        % Superior temporal gyrus
        recData(el).stg(1) = sum(adjAmp(STG,el),1);
        recData(el).stg(2) = 100./(recData(el).tot).*recData(el).stg(1);
%     else
%         continue
%     end
end

%% Inputs from LOCAL sites

locInputs = NaN(length(n1_data),2);
for el = IFG
    locInputs(el,1) = (sum(adjAmp(IFG,el))./(length(IFG)-1)).*100;
    locInputs(el,2) = 100 - locInputs(el,1);
end
for el = SMG
    locInputs(el,1) = (sum(adjAmp(SMG,el))./(length(SMG)-1)).*100;
    locInputs(el,2) = 100 - locInputs(el,1);
end
for el = STG   
    locInputs(el,1) = (sum(adjAmp(STG,el))./(length(STG)-1)).*100;
    locInputs(el,2) = 100 - locInputs(el,1);
end

%% Save data
switch s
    case '99'
        responseAmp_99 = responseAmp;
        responseWidth_99 = responseWidth;
        responseLat_99 = responseLat;
        locInputs_99 = locInputs;
        recData_99 = recData;
        projData_99 = projData;
        adjAmp_99 = adjAmp;
        save('m-files/detectN1/data/connectivity_99', 'responseAmp_99', 'responseWidth_99',...
            'responseLat_99', 'locInputs_99', 'recData_99', 'projData_99', 'adjAmp_99');
    case '31'
        responseAmp_31 = responseAmp;
        responseWidth_31 = responseWidth;
        responseLat_31 = responseLat;
        locInputs_31 = locInputs;
        recData_31 = recData;
        projData_31 = projData;
        adjAmp_31 = adjAmp;
        save('m-files/detectN1/data/connectivity_31', 'responseAmp_31', 'responseWidth_31',...
            'responseLat_31', 'locInputs_31', 'recData_31', 'projData_31', 'adjAmp_31');
    case '78'
        responseAmp_78 = responseAmp;
        responseWidth_78 = responseWidth;
        responseLat_78 = responseLat;
        locInputs_78 = locInputs;
        recData_78 = recData;
        projData_78 = projData;
        adjAmp_78 = adjAmp;
        save('m-files/detectN1/data/connectivity_78', 'responseAmp_78', 'responseWidth_78',...
            'responseLat_78', 'locInputs_78', 'recData_78', 'projData_78', 'adjAmp_78');
    case '81'
        responseAmp_81 = responseAmp;
        responseWidth_81 = responseWidth;
        responseLat_81 = responseLat;
        locInputs_81 = locInputs;
        recData_81 = recData;
        projData_81 = projData;
        adjAmp_81 = adjAmp;
        save('m-files/detectN1/data/connectivity_81', 'responseAmp_81', 'responseWidth_81',...
            'responseLat_81', 'locInputs_81', 'recData_81', 'projData_81', 'adjAmp_81');
    case '88'
        responseAmp_88 = responseAmp;
        responseWidth_88 = responseWidth;
        responseLat_88 = responseLat;
        locInputs_88 = locInputs;
        recData_88 = recData;
        projData_88 = projData;
        adjAmp_88 = adjAmp;
        save('m-files/detectN1/data/connectivity_88', 'responseAmp_88', 'responseWidth_88',...
            'responseLat_88', 'locInputs_88', 'recData_88', 'projData_88', 'adjAmp_88');
    case '50'
        responseAmp_50 = responseAmp;
        responseWidth_50 = responseWidth;
        responseLat_50 = responseLat;
        locInputs_50 = locInputs;
        recData_50 = recData;
        projData_50 = projData;
        adjAmp_50 = adjAmp;
        save('m-files/detectN1/data/connectivity_50', 'responseAmp_50', 'responseWidth_50',...
            'responseLat_50', 'locInputs_50', 'recData_50', 'projData_50', 'adjAmp_50');
end

disp('done')

%% Plot matrices

% To only plot data between the first and the last electrodes of interest:
firstEl = min(els2check); 
lastEl = max(els2check); 

response2plot = sqrt(responseAmp.^2);

screenSize = get( groot, 'Screensize' );
h = figure('Position', screenSize, 'Color', 'k')                           %#ok
set(h, 'PaperPositionMode', 'auto')
subplot(2,2,1)
imagesc(response2plot, [0 300])
colorbar('horizontal')
colormap(pinco)
title('\color{white}Amp')
xlabel('\color{white}measured')
ylabel('\color{white}stimulated')

subplot(2,2,2)
imagesc(responseWidth)
colorbar('horizontal')
colormap(pinco)
title('\color{white}Width')
xlabel('\color{white}measured')
ylabel('\color{white}stimulated')

subplot(2,2,3)
imagesc(responseLat)
colorbar('horizontal')
colormap(pinco)
title('\color{white}Latency')
xlabel('\color{white}measured')
ylabel('\color{white}stimulated')

subplot(2,2,4)
imagesc(adjAmpProd)
colorbar('horizontal')
colormap(pinco)
title('\color{white}Reciprocal connections')
xlabel('\color{white}measured')
ylabel('\color{white}stimulated')

axHands = findobj(get(h, 'Children'), 'flat', 'type', 'axes');
axis(axHands, 'square', [firstEl lastEl firstEl lastEl])
set(axHands, 'XAxisLocation', 'origin')

% text('\bf bla bla', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

set(gcf, 'InvertHardCopy', 'off');
print(gcf, '-dpng', strcat(['figures/connectivity_matrices_subj' s]))

%% Plot directed connectivity on cortex

plot_cxDirected;

close all

%% Connectivity net degree coefficient for each area
% 
ccep_degree;

%% Clustering analysis
ccep_clustering;













