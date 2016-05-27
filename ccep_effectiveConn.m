%% Describe effective connectivity measures between areas
% First z-scores the amplitudes in the adjacency matrix. then 
%
%% Get data
get_connectivity_data;
[IFG, SMG, AG, STG, MTG, els2check, esmPlus, ~] = ccep_anat_locs(s);
PL = [SMG, AG, STG, MTG];
%% 
firstEl = min(els2check); 
lastEl = max(els2check); 

% Get z-scored amplitudes of responses
mat2test = sqrt(responseAmp.^2);
mat2test(isnan(mat2test)==1) = 0;
[respStrength, mu, sigma] = zscore(mat2test);
respStrength(adjAmp==0) = NaN;

%% Plot matrix of response z-scores

h = figure;
imagesc(respStrength)
colormap(gray)
colorbar('vertical')
title('Z-score (amp)')
xlabel('measured')
ylabel('stimulated')

axHands = findobj(get(h, 'Children'), 'flat', 'type', 'axes');
axis(axHands, 'square', [firstEl lastEl firstEl lastEl])
set(axHands, 'XAxisLocation', 'origin')
print(gcf, '-painters', '-dpng', '-r0',  strcat(['figures/ampZsubj_' s]))

%% amplitude z-score distribution

amp_ifg2smg(s_nr) = nanmedian(nanmedian(respStrength(IFG,SMG),2));
amp_smg2ifg(s_nr) = nanmedian(nanmedian(respStrength(SMG,IFG),2));

amp_ifg2ag(s_nr) = nanmedian(nanmedian(respStrength(IFG,AG),2));
amp_ag2ifg(s_nr) = nanmedian(nanmedian(respStrength(AG,IFG),2));

amp_ifg2stg(s_nr) = nanmedian(nanmedian(respStrength(IFG,STG),2));
amp_stg2ifg(s_nr) = nanmedian(nanmedian(respStrength(STG,IFG),2));

amp_ifg2mtg(s_nr) = nanmedian(nanmedian(respStrength(IFG,MTG),2));
amp_mtg2ifg(s_nr) = nanmedian(nanmedian(respStrength(MTG,IFG),2));

amp_ifg2ifg(s_nr) = nanmedian(nanmedian(respStrength(IFG,IFG),2));
amp_smg2smg(s_nr) = nanmedian(nanmedian(respStrength(SMG,SMG),2));
amp_ag2ag(s_nr) = nanmedian(nanmedian(respStrength(AG,AG),2));
amp_stg2stg(s_nr) = nanmedian(nanmedian(respStrength(STG,STG),2));
amp_mtg2mtg(s_nr) = nanmedian(nanmedian(respStrength(MTG,MTG),2));

%% width distribution

width_ifg2smg(s_nr) = nanmedian(nanmedian(responseWidth(IFG,SMG),2));
width_smg2ifg(s_nr) = nanmedian(nanmedian(responseWidth(SMG,IFG),2));

width_ifg2ag(s_nr) = nanmedian(nanmedian(responseWidth(IFG,AG),2));
width_ag2ifg(s_nr) = nanmedian(nanmedian(responseWidth(AG,IFG),2));

width_ifg2stg(s_nr) = nanmedian(nanmedian(responseWidth(IFG,STG),2));
width_stg2ifg(s_nr) = nanmedian(nanmedian(responseWidth(STG,IFG),2));

width_ifg2mtg(s_nr) = nanmedian(nanmedian(responseWidth(IFG,MTG),2));
width_mtg2ifg(s_nr) = nanmedian(nanmedian(responseWidth(MTG,IFG),2));

width_ifg2ifg(s_nr) = nanmedian(nanmedian(responseWidth(IFG,IFG),2));
width_smg2smg(s_nr) = nanmedian(nanmedian(responseWidth(SMG,SMG),2));
width_ag2ag(s_nr) = nanmedian(nanmedian(responseWidth(AG,AG),2));
width_stg2stg(s_nr) = nanmedian(nanmedian(responseWidth(STG,STG),2));
width_mtg2mtg(s_nr) = nanmedian(nanmedian(responseWidth(MTG,MTG),2));

%% latency distribution

lat_ifg2smg(s_nr) = nanmedian(nanmedian(responseLat(IFG,SMG),2));
lat_smg2ifg(s_nr) = nanmedian(nanmedian(responseLat(SMG,IFG),2));

lat_ifg2ag(s_nr) = nanmedian(nanmedian(responseLat(IFG,AG),2));
lat_ag2ifg(s_nr) = nanmedian(nanmedian(responseLat(AG,IFG),2));

lat_ifg2stg(s_nr) = nanmedian(nanmedian(responseLat(IFG,STG),2));
lat_stg2ifg(s_nr) = nanmedian(nanmedian(responseLat(STG,IFG),2));

lat_ifg2mtg(s_nr) = nanmedian(nanmedian(responseLat(IFG,MTG),2));
lat_mtg2ifg(s_nr) = nanmedian(nanmedian(responseLat(MTG,IFG),2));

lat_ifg2ifg(s_nr) = nanmedian(nanmedian(responseLat(IFG,IFG),2));
lat_smg2smg(s_nr) = nanmedian(nanmedian(responseLat(SMG,SMG),2));
lat_ag2ag(s_nr) = nanmedian(nanmedian(responseLat(AG,AG),2));
lat_stg2stg(s_nr) = nanmedian(nanmedian(responseLat(STG,STG),2));
lat_mtg2mtg(s_nr) = nanmedian(nanmedian(responseLat(MTG,MTG),2));

