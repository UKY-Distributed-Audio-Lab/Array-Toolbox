function mv = mkmov3d(tpos,fratein,frateout,fname,sfactor)
%
% This function creates a movie out of 3D SSL data organized
% in an ASA file.
%
%       mv = mkmov3d(tpos,fratein,frateout,fname,sfactor)
%
% Inputs:
%       TPOS     => 3D SSL data points
%                       col 1    - Segment number
%                       col 2    - 
%                       col 3    - Intensity
%                       cols 4-6 - 3D Coordinates
%                       col 7    - Frame Number
%                       col 8    - Detection Indicator
%       FRATEIN  => Input frame rate
%       FRATEOUT => Output frame rate
%       FNAME    => String of AVI file name where movie is to be stored
%       SFACTOR  => Scale factor for size of mark based on detection statistic
%                   (100 is recommended)
%
% View code viewmovie.m
%
% Written by Harikrishnan Unnikrishnan      05/27/2008
%                             Modified      06/23/2008
%                             Modified      08/09/2012
%                        - Updated TPOS



close all;
pos = [210   343   690   572];
mv = avifile(fname);         % Opening avi file   
mv.Fps = frateout;
%mv.compression = 'None';
crat = (fratein/frateout);   % Conversion frame rate

fig = figure;
% Setting figure properties
set(fig,'DoubleBuffer','off','Position',pos);
set(gca,'NextPlot','replace','Visible','off');


% Getting easily differentiable colors to assign for each segment
colormap lines;   
    c = colormap;
    c = unique(c,'rows');
    colormap(c);

% Making the size of the point proportional to the intensity
% size of pt = 1 + round[(intensity*scalefactor)/maxintensity]
siz = 1 + round(tpos(:,3)*sfactor/(max(tpos(:,3))));

% Assigning colors in a cycle. (If there are more than 7 segments)
col = mod(tpos(:,1),7)+1;
fnum = 1;
afnum = 0;
crat = (fratein/frateout);  % Conversion frame rate

% Fetching Frames
% while fnum  <  max(tpos(:,6))-ceil(crat)
while fnum  <  max(tpos(:,7))-ceil(crat)
    fnum = round(afnum) + 1;
    afnum = afnum + crat;
    tmp = tpos(tpos(:,7)>=fnum & tpos(:,7)<round(fnum+crat+1),4:6);
    tsiz = siz(tpos(:,7)>=fnum & tpos(:,7)<round(fnum+crat+1));
    tcol = col(tpos(:,7)>=fnum & tpos(:,7)<round(fnum+crat+1));
    
    
    s3 = scatter3((tmp(:,1)),(tmp(:,2)),(tmp(:,3)),tsiz,tcol,'filled','o');
    xlim([0 3.8]); ylim([0 3.8]); zlim([0 2.3])
    hold on; stem3((tmp(:,1)),(tmp(:,2)),(tmp(:,3)),'+');
%    set(gcf,'Position', pos)
    set(gca,'DataAspectRatio', [1.6522    1.6522    1.0000])
%     set(gca,'CameraPosition', [-23.4471   20.6559    6.8520])  %  For camera A 
    set(gca,'CameraPosition', [-17.1261  -22.8952    7.3878])  %  For camera B 
%     d = get(s3,'parent');   % Getting parent object to set camera position
%     set(d,sobj);
    xlabel('X Meters');ylabel('Y Meters');zlabel('Z Meters');
%     title(['time from start = ' num2str((fnum-1)/fratein)])
    mmv = getframe(fig);
    mv = addframe(mv,mmv);
    pause(.1)
    hold off;
end;

mv = close(mv);


