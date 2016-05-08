%% Identify clusters
% Calculate the clustering coefficient (for directed graph) for each
% electrode.
% The directed graph local clustering coefficient is the ration of the
% total number of existing connections in the neighbourhood of an
% electrode. The neighbourhood are all the electrodes of that electrode to
% which it is connected in any direction. For the total number of
% connections, 1 to 2 and 2 to 1 are counted as separate unique
% connections, because the graph is directed.
%
% A.G.Vassileva & D.Hermes, UMC Utrecht, 2016

%% Load subject data
subjects = {'99','31','78','81','88'};
s = subjects{s_nr};
eval(['ccep = subj_info_SPES_' s]); % make CCEP structure
s_info = subj_info_ccep(s); % get subject info
load(s_info.cortex) % get cortex
load(s_info.els) 
[IFG, SMG, STG, els2check] = ccep_anat_locs(s);
%% Get connectivity data
switch s_nr
    case 1
         adjAmp = adjAmp_99;
    case 2
         adjAmp = adjAmp_31;
    case 3
         adjAmp = adjAmp_78;
    case 4
         adjAmp = adjAmp_81;
    case 5
         adjAmp = adjAmp_88;
end

%% Calculate clustering coefficients for each electrode
clustering = struct('el',[],'neighbours',[], 'neighboursConns', [], 'cluster', []);
clustermat = NaN(1,length(adjAmp));
for k = 1:size(adjAmp,2)
    clustering(k).el = k;
end

for el = els2check
    % Get neighbourhood of el
    [~, a] = find(adjAmp(el,:)==1);
    [~, b] = find(adjAmp(:,el)'==1);
    neighbourhood = [a,setdiff(b,a)];
    % continue only if there are more than 2 els in the neighbourhood
    if length(neighbourhood) > 2
        clustering(el).neighbours = neighbourhood;
        c = adjAmp(neighbourhood,neighbourhood);
        clustering(el).neighboursConns = sum(sum(c,2));
        clustering(el).cluster = clustering(el).neighboursConns/(numel(neighbourhood)*(numel(neighbourhood) - 1));
        clustermat(el) = clustering(el).cluster;
    end
    clear a b c
end

clusterifg = nanmean(clustermat(IFG));
clustersmg = nanmean(clustermat(SMG));
clusterstg = nanmean(clustermat(STG));

%% Save data
switch s_nr
    case 1
         clustermat_99 = clustermat;
         save('clustermat_99', 'clustermat_99')
    case 2
         clustermat_31 = clustermat;
         save('clustermat_99', 'clustermat_99')
    case 3
         clustermat_78 = clustermat;
         save('clustermat_99', 'clustermat_99')
    case 4
         clustermat_81 = clustermat;
         save('clustermat_99', 'clustermat_99')
    case 5
         clustermat_88 = clustermat;
         save('clustermat_99', 'clustermat_99')
end

%% Plot on cortex
elclustermat = elecmatrix(els2check,:);
a = min(elecmatrix(:,1));
elclustermat(:,1) = a;
maxcluster = max(clustermat);
figure('Color', [0 0 0])
ctmr_gauss_plot(cortex,[0 0 0],0)
el_add_size(elclustermat,clustermat(els2check),1, [.9 .9 .9])
title(['\color{white}cluster distribution patient ' int2str(s_nr)])
set(gcf, 'InvertHardCopy', 'off');
print(gcf,'-r300', '-dpng', strcat(['figures/pres_figures/clusterDistribution' s]))

clear a


