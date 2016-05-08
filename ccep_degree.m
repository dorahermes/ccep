%% The connectivity net degree coefficient is the ratio between the outdegree and the
% indegree for an area. If the coefficient is larger that one, that area
% acts mostly as a projector, rather than an integrator. If it's smaller
% than one, that area is mostly an integrator. 
%
% A.G. Vassileva & D. Hermes, UMC Utrecht, 2016
%% Subject info
subjects = {'99','31','78','81','88'};
s = subjects{s_nr};

%% Define electrodes, areas and get adjacency matrix.
[IFG, SMG, STG, els2check] = ccep_anat_locs(s);
pla = [SMG STG]; % posterior language area electrodes. IFG is anterior language area

switch s
    case '99'
        adjAmp = adjAmp_99;
        elsifg(1) = length(IFG);
        elspla(1) = length(pla);
    case '31'
        adjAmp = adjAmp_31;
        elsifg(2) = length(IFG);
        elspla(2) = length(pla);
    case '78'
        adjAmp = adjAmp_78;
        elsifg(3) = length(IFG);
        elspla(3) = length(pla);
    case '81'
        adjAmp = adjAmp_81;
        elsifg(4) = length(IFG);
        elspla(4) = length(pla);
    case '88'
        adjAmp = adjAmp_88;
        elsifg(5) = length(IFG);
        elspla(5) = length(pla);
end

%% Outdegree, indegree and net degree in number of connections
%% Annterior language area (i.e. inferior frontal gyrus)
% Outdegree is the total number of unique connections from that area.

outConnifg = sum(adjAmp(IFG,pla),1); % outdegree in number of connections
outConnifg(outConnifg==2) = 1;
outConnifg = sum(outConnifg);
inConnifg = sum(adjAmp(pla,IFG),2); % outdegree in number of connections
inConnifg(inConnifg==2) = 1;
inConnifg = sum(inConnifg);

% Net degree coefficient
switch s
    case '99'
        netDegifg(1) = outConnifg/inConnifg;
    case '31'
        netDegifg(2) = outConnifg/inConnifg;
    case '78'
        netDegifg(3) = outConnifg/inConnifg;
    case '81'
        netDegifg(4) = outConnifg/inConnifg;
    case '88'
        netDegifg(5) = outConnifg/inConnifg;
end


%% SMG:
notsmg = [IFG STG]; % els that are not in area of interest
outConnsmg = sum(adjAmp(SMG,notsmg),1); % outdegree in number of connections
outConnsmg(outConnsmg==2) = 1;
outConnsmg = sum(outConnsmg);
inConnsmg = sum(adjAmp(notsmg,SMG),2); % outdegree in number of connections
inConnsmg(inConnsmg==2) = 1;
inConnsmg = sum(inConnsmg);

% Net degree coefficient
switch s
    case '99'
        netDegsmg(1) = outConnsmg/inConnsmg;
    case '31'
        netDegsmg(2) = outConnsmg/inConnsmg;
    case '78'
        netDegsmg(3) = outConnsmg/inConnsmg;
    case '81'
        netDegsmg(4) = outConnsmg/inConnsmg;
    case '88'
        netDegsmg(5) = outConnsmg/inConnsmg;
end

% %% AG:
% notag = [IFG SMG STG MTG]; % els that are not in area of interest
% outConnag = sum(adjAmp(AG,notag),1); % outdegree in number of connections
% outConnag(outConnag==2) = 1;
% outConnag = sum(outConnag);
% inConnag = sum(adjAmp(notag,AG),2); % outdegree in number of connections
% inConnag(inConnag==2) = 1;
% inConnag = sum(inConnag);
% 
% % Net degree coefficient
% switch s
%     case '99'
%         netDegag(1) = outConnag/inConnag;
%     case '31'
%         netDegag(2) = outConnag/inConnag;
%     case '78'
%         netDegag(3) = outConnag/inConnag;
%     case '81'
%         netDegag(4) = outConnag/inConnag;
%     case '88'
%         netDegag(5) = outConnag/inConnag;
% end

%% STG:
notstg = [IFG SMG]; % els that are not in area of interest
outConnstg = sum(adjAmp(STG,notstg),1); % outdegree in number of connections
outConnstg(outConnstg==2) = 1;
outConnstg = sum(outConnstg);
inConnstg = sum(adjAmp(notstg,STG),2); % outdegree in number of connections
inConnstg(inConnstg==2) = 1;
inConnstg = sum(inConnstg);

% Net degree coefficient
switch s
    case '99'
        netDegstg(1) = outConnstg/inConnstg;
    case '31'
        netDegstg(2) = outConnstg/inConnstg;
    case '78'
        netDegstg(3) = outConnstg/inConnstg;
    case '81'
        netDegstg(4) = outConnstg/inConnstg;
    case '88'
        netDegstg(5) = outConnstg/inConnstg;
end

%% MTG:
% notmtg = [IFG SMG AG STG]; % els that are not in area of interest
% outConnmtg = sum(adjAmp(MTG,notmtg),1); % outdegree in number of connections
% outConnmtg(outConnmtg==2) = 1;
% outConnmtg = sum(outConnmtg);
% inConnmtg = sum(adjAmp(notmtg,MTG),2); % outdegree in number of connections
% inConnmtg(inConnmtg==2) = 1;
% inConnmtg = sum(inConnmtg);
% 
% % Net degree coefficient
% switch s
%     case '99'
%         netDegmtg(1) = outConnmtg/inConnmtg;
%     case '31'
%         netDegmtg(2) = outConnmtg/inConnmtg;
%     case '78'
%         netDegmtg(3) = outConnmtg/inConnmtg;
%     case '81'
%         netDegmtg(4) = outConnmtg/inConnmtg;
%     case '88'
%         netDegmtg(5) = outConnmtg/inConnmtg;
% end

%% Posterior language area

outConnpla = sum(adjAmp(pla,IFG),1); % outdegree in number of connections
outConnpla(outConnpla==2) = 1;
outConnpla = sum(outConnpla);
inConnpla = sum(adjAmp(IFG,pla),2); % outdegree in number of connections
inConnpla(inConnpla==2) = 1;
inConnpla = sum(inConnpla);

% Net degree coefficient
switch s
    case '99'
        netDegpla(1) = outConnpla/inConnpla;
    case '31'
        netDegpla(2) = outConnpla/inConnpla;
    case '78'
        netDegpla(3) = outConnpla/inConnpla;
    case '81'
        netDegpla(4) = outConnpla/inConnpla;
    case '88'
        netDegpla(5) = outConnpla/inConnpla;
end








