function ccep_el_add_size_and_color(els,r2_size,r2_color,r2_sizemax,r2_colormax)

%   els = coordinates matrix
%   r2_size = input amplitudes
%   r2_color = input latencies
%   r2_sizemax = scale for amplitudes
%   r2_colormax = scale for latencies

%     Copyright (C) 2006  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%     Modified for own use by J. van der Aar, april 2019, UMC Utrecht
%     whereas color represents latency and size represents amplitude

hold on

% red-green color scale
cm1=[repmat([0 0 0],80,1)];
cm1(1:40,2)=[1];
cm1(41:80,1)=[1];
cm1(41:80,2)=flip([0:1/39:1]);
cm1(1:40,1) = [0:1/39:1];

cm2=[repmat([.5 .5 .5],100,1)];

elsize=[15:(45-15)/(100-1):60];

% set the max to 100
r2_size=100*r2_size/r2_sizemax;
r2_color=80*r2_color/r2_colormax;
r2_size = round(r2_size);
r2_color = round(r2_color);

% set everything larger than 100 to 100
r2_size(r2_size>100) = 100;
r2_color(r2_color>80) = 100;
r2_size(r2_size<-100) = -100;
r2_color(r2_color<-80) = -100;


% electrode with r2:
for k=1:size(els,1)
    if ~isnan(r2_size(k))
        elsize_r2=elsize(abs(r2_size(k)));
        elcol_r2=cm1(r2_color(k),:); 
        plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize_r2)
        plot3(els(k,1),els(k,2),els(k,3),'.','Color',elcol_r2,'MarkerSize',elsize_r2-5)
    else
        plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize(1))
    end

end

%% Legendas

%%%% Run these for the legendas of the size(amplitude) and color(latency)

%% red-green color scale

cm1=[repmat([0 0 0],80,1)];
cm1(1:40,2)=[1];
cm1(41:80,1)=[1];
cm1(41:80,2)=flip([0:1/39:1]);
cm1(1:40,1) = [0:1/39:1];

% elsize=[15:(45-15)/(100-1):45];
elsize = ones(1,100)*50;

figure('Color',[1 1 1],'Position',[30 50 500 150]),hold on
for k=1:80
    plot(k,1,'.','MarkerSize',elsize(k),'Color',cm1(k,:))
end
title('Electrode Color')
xlim([10 70])
set(gca,'YColor','w')
set(gca,'XTick', [10:10:70],'YTick',[], 'FontName','arial','FontSize',18)
xlabel('N1 peak latency (ms)')
set(gcf, 'PaperPositionMode', 'auto');

%% grey scale for amplitude

cm2=[repmat([.5 .5 .5],100,1)];

elsize=[15:(100-15)/(100-1):100];



figure('Color',[1 1 1],'Position',[30 50 500 150]),hold on
for k=1:200
    plot((k*10),1,'.','MarkerSize',elsize(ceil(k/2)),'Color',cm2(ceil(k/2),:))
end
title('Electrode Size')
set(gca,'YColor','w')
xlim([0 2000])
set(gca,'XTick', [0:500:2000],'YTick',[], 'FontName','arial','FontSize',18)
xlabel('amplitude peak (uV)')
set(gcf, 'PaperPositionMode', 'auto');
