function el_add_sizecolor(els,r2_size,r2_color,r2_sizemax,r2_colormax)

% function el_add_withr2size(els,r2)
% 
% els: rows = electrodes, columns = xyz
% r2: for size
% varargin:
% {1} : maximum for scale (if {''} absmax)
% example: el_add_withr2size(els,r2,1)

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

hold on

%gray2green
% cm1=[repmat([0 0 0],100,1)];
% cm1(1:10,2)=[0:0.7/9:0.7]';
% cm1(10:100,2)=[0.7:(1-0.7)/90:1]';
% cm1(1:10,1)=[0]';
% cm1(1:10,3)=[0]';
% cm1(20:100,1)=[0:1/80:1]';
% cm1=[repmat([0 1 0],100,1)];

% cm1(1:100,2)=[0:1/99:1]';
% cm1=[repmat([0 0 0],100,1)];
% cm1(1:10,1)=[0:0.7/9:0.7]';
% cm1(10:100,1)=[0.7:(1-0.7)/90:1]';
% cm1(1:10,2)=[0]';
% cm1(1:10,3)=[0]';
% cm1(20:100,2)=[0:1/80:1]';

cm1 = jet(64);
cm1 = resample(cm1,100,64);
cm1(cm1>1) = 1;
cm1(cm1<0) = 0;

cm2=[repmat([0 0 0],100,1)];
cm2(1:10,3)=[0:0.7/9:0.7]';
cm2(10:100,3)=[0.7:(1-0.7)/90:1]';
cm2(1:10,2)=[0]';
cm2(1:10,1)=[0]';
cm2(20:100,2)=[0:1/80:1]';

elsize=[15:(45-15)/(100-1):60];

% set the max to 100
r2_size=100*r2_size/r2_sizemax;
r2_color=100*r2_color/r2_colormax;
r2_size = round(r2_size);
r2_color = round(r2_color);

% set everything larger than 100 to 100
r2_size(r2_size>100) = 100;
r2_color(r2_color>100) = 100;
r2_size(r2_size<-100) = -100;
r2_color(r2_color<-100) = -100;


% electrode with r2:
for k=1:size(els,1)
    if ~isnan(r2_size(k))
        if abs(r2_size(k))>0.01
            
            elsize_r2=elsize(abs(r2_size(k)));
            if r2_color(k)>0.01
                elcol_r2=cm1(r2_color(k),:); 
                plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize_r2)
                plot3(els(k,1),els(k,2),els(k,3),'.','Color',elcol_r2,'MarkerSize',elsize_r2-5)
            elseif r2_color(k)<0.01
                elcol_r2=cm2(abs(r2_color(k)),:);
                plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize_r2)
                plot3(els(k,1),els(k,2),els(k,3),'.','Color',elcol_r2,'MarkerSize',elsize_r2-5)
            end
        else
            plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize(1))
        end
    end
end


%% example code to plot color scale

 
% cm1(1:100,2)=[0:1/99:1]';
% cm1=[repmat([0 0 0],100,1)];
% cm1(1:10,1)=[0:0.7/9:0.7]';
% cm1(10:100,1)=[0.7:(1-0.7)/90:1]';
% cm1(1:10,2)=[0]';
% cm1(1:10,3)=[0]';
% cm1(20:100,2)=[0:1/80:1]';
% 
% cm2=[repmat([0 0 0],100,1)];
% cm2(1:10,3)=[0:0.7/9:0.7]';
% cm2(10:100,3)=[0.7:(1-0.7)/90:1]';
% cm2(1:10,2)=[0]';
% cm2(1:10,1)=[0]';
% cm2(20:100,2)=[0:1/80:1]';
% 
% % cm2=cm2(end:-1:1,:);
% % cm=[cm2; cm1];
% elsize=[15:(45-15)/(100-1):45];

% figure('Color',[1 1 1],'Position',[30 50 50 300]),hold on
% for k=1:100
%     plot(1,k,'.','MarkerSize',elsize(k),'Color',cm1(k,:))
% end
%     
% for k=1:100
%     plot(1,-k,'.','MarkerSize',elsize(k),'Color',cm2(k,:))
% end

% ylim([-120 120])
% set(gcf, 'PaperPositionMode', 'auto');
% print('-painters','-r300','-dpng',strcat(['./figures/colorscale_el_add_sizable']));
% print('-painters','-r300','-depsc',strcat(['./figures/colorscale_el_add_sizable']));


% %% plot electrode color
% figure('Color',[1 1 1],'Position',[30 50 50 300]),hold on
% for k=1:100
%     plot(1,k,'.','MarkerSize',40,'Color',cm1(k,:))
% end
% set(gcf, 'PaperPositionMode', 'auto');
% print('-painters','-r300','-dpng',strcat(['./figures/el_add_sizecolor_colorscale']));
% print('-painters','-r300','-depsc',strcat(['./figures/el_add_sizecolor_colorscale']));
% 
% %% plot electrode size
% 
% elsize=[15:(45-15)/(100-1):45];
% 
% figure('Color',[1 1 1],'Position',[30 50 50 300]),hold on
% for k=1:100
%     plot(1,k,'.','MarkerSize',elsize(k),'Color',[.5 .5 .5])
% end
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print('-painters','-r300','-dpng',strcat(['./figures/el_add_sizecolor_sizescale']));
% print('-painters','-r300','-depsc',strcat(['./figures/el_add_sizecolor_sizescale']));

