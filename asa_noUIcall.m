function asa_noUIcall(pn,fn,pmp,fmp)
%
%  This function computes an Auditory Scene Analysis (ASA) based on a 
%  recording from a distributed microphone system.  The ASA consists of a 
%  text file identifying an active sound source, its location, and its 
%  connection to sound sources at other times.  The function prompts for a 
%  multichannel wave file from a distributed microphone recording.  It then
%  prompts for a .txt file containing the microphone array information.  
%  It then prompts through a series of dialog boxes for pressure, 
%  temperature, and humidity (for sound speed calculation), and FOV regions.
%  Users enter the range in the Z direction and then selects either square 
%  or rectangular regions using mouse clicks.
%
%  It then performs Sound Source Location (SSL) using a Coherent Steered 
%  Response Power (CSRP) algorithm. 
%  CFAR threshold (1e-5) is applied to identify sound sources. A text file 
%  is then created to identify detected sources, strength of their detection 
%  statistics, and their locations in time and space. For each detection, 
%  a row is created with the following format:
%  strms1 = [timestamp, detection statistic in excess of threshold, x, y, z postions in meters]
%
%  For now this information is written as a variable in a .mat file along
%  with parameters in a data structure FOVP, such as sample frequency (fs) 
%  and speed of sound (c).
%
% Required Functions:
%           signifslist.m
%   
% Written by Kevin D. Donohue (donohue@engr.uky.edu)    May 2012
% Modified by Kirstin Brangers                          July 2012
%


fovp.graphicout_flag = 0;   %  Set to 1 for collapsed acoustic images for each frame while processing
fovp.fc = 400;              %  High-pass cutoff freq. filter for signal conditioning
fovp.fs = 16000;            %  Sampling rate in Hz for resampling audio file
fovp.beta = .85;            %  Beta for whitening
% fovp.shape = 1.0;         %  Shape parameter for probability of false alarm computation 

%  Time window for Frequency Domain Block Processing
fovp.trez = 40e-3;          %  In seconds
fovp.tinc = fovp.trez/2;    %  Time increment in seconds

% CFAR Threshold
fovp.pcfar = 1e-3;        %  False alarm probability

%  Room Resolution: Step through Cartesian grid for mic and sound source plane
rad = [.42, .42, .2];       %  Threshold neighborhood for CFAR detection [x,y,z] in meters
rez = [.04, .04, .08];       %  Grid Resolution in Meters [x,y,z] in meters
 
% Compute Region Extension to Compensate for Threshold Neighborhood
fovp.threshneigh = round(rad ./ rez);       %  Convert neighborhood in meters to samples

% Prompt for wave file:
% [fn, pn] = uigetfile('*.wav','Select Multichannel Recording.');
% if isnumeric(fn)                %  If menu cancelled, quit
%     error('file error')
% end


% Prompt for microphone position file
% [fmp, pmp] = uigetfile('*.txt','Select text file with microphone positions.');
% if isnumeric(fn)                %  If menu cancelled, quit
%     error('file error')
% end
fovp.mp = load([pmp, fmp]);     % Load file with mic positions

% % Prompt for Environment Parameters
% prompt={'Temperature (degrees C):','Pressure (mmHg):','Relative Humidity (percent):'};
% name='Parameters to Determine Sound Speed';
% numlines=1;
% % defaultanswer={'24.0','29.06', '61.5'};
% defaultanswer={'25.6','29.095', '48.4'};
% answer=inputdlg(prompt,name,numlines,defaultanswer,'on');
% tmp = str2double(answer{1});        % Temperature
% pres = str2double(answer{2});       % Pressure
% rh = str2double(answer{3});         % Relative Humidity

tmp = 25.6;
pres = 29.095;
rh = 48.4;

% Compute Speed of Sound with Entered Parameters
fovp.c = SpeedOfSound(tmp,rh,pres);

% Compute Limits for Initial Plots
minx = min(fovp.mp(1,:));
maxx = max(fovp.mp(1,:));
miny = min(fovp.mp(2,:));
maxy = max(fovp.mp(2,:));
minz = min(fovp.mp(3,:));
minz = min([minz,0]);
maxz = max(fovp.mp(3,:));

%  Plot Microphone Positions
figure(1)
plot3(fovp.mp(1,:),fovp.mp(2,:),fovp.mp(3,:),'xb') % Plot mic postions in blue with 'x'
hold on
[rw cl] = size(fovp.mp);
for k=1:cl
    text(fovp.mp(1,k)+.25,fovp.mp(2,k)+.25,fovp.mp(3,k)+.25,int2str(k))
end
hold off
xlabel('X Meters');
ylabel('Y Meters');
zlabel('Z Meters');
title(['Microphone Positions with Sound Speed = ' num2str(fovp.c) ' m/s'])
% axis([minx(1)-1, maxx(1)+1, miny(1)-1, maxy(1)+1, minz(1),maxz(1)+.1])
axis([0, 4, 0, 4, 0, 3])
grid on
pause(1)

% % Prompt User for Height Coordinates
% % Default coords for Z position - range from .6 meters (~2ft) to 2 meters (~6.5ft)
% % NOTE: This range is used for all regions selected in FOV
% prompt={'Enter first Z coord of FOV (meters):', 'Enter second Z coord of FOV (meters):'};
% name='Parameters to Determine Height of FOV Region';
% numlines=1;
% defaultanswer={'0.6','2.0'};        % Default Z range
% answer=inputdlg(prompt,name,numlines,defaultanswer,'on');
% z1 = str2double(answer{1});         % First coord for Z range
% z2 = str2double(answer{2});         % Second coord for Z range

% % 2D Plot So X,Y Coordinates Can Be Selected with Mouse Click
% figure(1)
% plot(fovp.mp(1,:),fovp.mp(2,:),'xb')    % Plot mic postions in blue with 'x'
% xlabel('X Meters');
% ylabel('Y Meters');
% title(['Microphone Positions with Sound Speed = ' num2str(fovp.c) ' m/s'])
% % axis([minx(1)-1, maxx(1)+1, miny(1)-1, maxy(1)+1])
% axis([0, 4, 0, 4])
% grid on
% pause(1)

% Prompt User for Region(s) of Interest
% fovflag = 'Yes';            % Turn on Fov Flag to prompt user for initial region
fovcnt = 0;                 % Initialize counter for number of FOV lines
% % xyposS = 1;                 % Initialize counter to hold X,Y coordinates for selecting source
% while fovflag(1) == 'Y'
%    % Prompt user to either select a source with radius or region of interest
%    roi = questdlg('Select Source or Region.', 'Region of Interest', 'Source' ,'Region', 'Source');
for i=1:4    
    z1 = [1.4, 1.2, 1.3, 1.6];
    z2 = [1.6, 1.4, 1.5, 1.8];
    
   % Prompt User for Height Coordinates
   % Default coords for Z position - range from .6 meters (~2ft) to 2 meters (~6.5ft)
%    prompt={'Enter first Z coord of FOV (meters):', 'Enter second Z coord of FOV (meters):'};
%    name='Parameters to Determine Height of FOV Region';
%    numlines=1;
%    defaultanswer={'1.4','1.8'};        % Default Z range
%    answer=inputdlg(prompt,name,numlines,defaultanswer,'on');
%    z1 = str2double(answer{1});         % First coord for Z range
%    z2 = str2double(answer{2});         % Second coord for Z range
%    fovm(fovcnt+1,3)=z1;
%    fovm(fovcnt+2,3)=z2;
    
    fovm(fovcnt+1,3)=z1(i);
    fovm(fovcnt+2,3)=z2(i);

% %    if roi(1) == 'S'                            % If user selects 'Source'
%         % Select Source of Interest - One Mouse Click
%         h = helpdlg('Select Position of Sound Source with Crosshairs','Select Source Position.');
%         uiwait(h);                             % Delete dialog box when user selects OK
%         [xCoord(xyposS),yCoord(xyposS)] = ginput(1);    % Get X,Y position of mouse click
        
%         % Ask user for radius around source
%         rdef = 0.0;                            % Default radius
%         prompt={'Enter radius around source (meters):'};
%         name='Parameters to Determine Radius of FOV Region';
%         numlines=1;
%         defaultanswer={num2str(rdef,3)}; 
%         answer=inputdlg(prompt,name,numlines,defaultanswer);
%         radxy = str2double(answer{1});         % Radius entered by user
%         % If radius is less than default radius, change to default radius
%         if radxy < rdef                        
%           radxy = rdef;
%         end
%         % Radius cannot be negative
%         if radxy < 0                        
%           radxy = 0;
%         end
     
    radxy = 0.4;
    xCoord = [0.5, 0.5, 3.079, 3.079];
    yCoord = [0.5, 3.079, 3.079, 0.5];
%        % Add radius to selected source position to obtain opposite corners
%        fovm(fovcnt+1,1)= xCoord(xyposS)-radxy;      % Corner 1
%        fovm(fovcnt+1,2)= yCoord(xyposS)-radxy;  
%        fovm(fovcnt+2,1)= xCoord(xyposS)+radxy;      % Corner 2
%        fovm(fovcnt+2,2)= yCoord(xyposS)+radxy;
    fovm(fovcnt+1,1)= xCoord(i)-radxy;      % Corner 1
    fovm(fovcnt+1,2)= yCoord(i)-radxy;  
    fovm(fovcnt+2,1)= xCoord(i)+radxy;      % Corner 2
    fovm(fovcnt+2,2)= yCoord(i)+radxy;
       
%        % Get values for graphing purposes
%        x = min(fovm(fovcnt+1,1),fovm(fovcnt+2,1));      % X position
%        y = min(fovm(fovcnt+1,2),fovm(fovcnt+2,2));      % Y position    
%        w = abs((fovm(fovcnt+1,1)-fovm(fovcnt+2,1)));    % Width of rect
%        h = abs((fovm(fovcnt+1,2)-fovm(fovcnt+2,2)));    % Height of rect
       
       % Plot selected source including radius
%        figure(1)
%        hold on
%        h = rectangle('Position',[x, y, w, h],'EdgeColor','r');
%        hold off
       
%    else                                     % If user selects 'Region'
%        % Select Region of Interest - Two Mouse Clicks
%        h = helpdlg('Select Region of Interest with Crosshairs (opposite corners)','Select Source Corners.');
%        uiwait(h);
%        [xCoordA,yCoordA] = ginput(1);       % First selected corner
%        
%        % Show first selected corner on plot
%        figure(1)
%        hold on
%        h = plot(xCoordA,yCoordA,'r.-');    % Plot selected corner in red with '.'
%        hold off
%        
%        % Get second corner
%        [xCoordB,yCoordB] = ginput(1);
%        
%        % Opposite corners
%        fovm(fovcnt+1,1)= xCoordA;               % Corner 1  X1
%        fovm(fovcnt+1,2)= yCoordA;               %           Y1               
%        fovm(fovcnt+2,1)= xCoordB;               % Corner 2  X2
%        fovm(fovcnt+2,2)= yCoordB;               %           Y2
%        
%        % Get smaller value for graphing purposes
%        x = min(xCoordA,xCoordB);                % X position
%        y = min(yCoordA,yCoordB);                % Y position
%        w = abs(xCoordA-xCoordB);                % Width of rectangle
%        h = abs(yCoordA-yCoordB);                % Height of rectangle
%        
%        % Plot selected region
%        figure(1)
%        hold on
%        h = rectangle('Position',[x, y, w, h],'EdgeColor','r');
%        hold off
%    end
   
   % Increment Counters
%    xyposS = xyposS+1;
   fovcnt = fovcnt+2;
%    fovflag = questdlg('Add another block of interest?', 'Field of View'); 
%    if fovflag(1) == 'C'                     %  If cancel is selected, quit
%     error('User exited program.')
%    end
end

% Delete previously chosen regions on 2d plot
%     delete(h);
    
% Display Selected Corners in Matlab Command Window
    lenfovcnt = fovcnt;
    fovcnt = 0;
    region = 1;
    xyposS = 1;
    fprintf('The selected regions of interest: \n');
    while lenfovcnt ~= 0
        fprintf(['The Source of Interest for Region ',num2str(region), ' is at corners:  ( ',num2str(fovm(fovcnt+1,:)), ' ) and ( ' ,num2str(fovm(fovcnt+2,:)),' ) \n \n']);
        lenfovcnt = lenfovcnt-2;
        fovcnt = fovcnt+2;
        region = region+1;
        xyposS = xyposS+1;
    end

%  Use FOV points to create grid for analysis
%  Plot user selected points only - no CFAR Threshold
for k=1:2:fovcnt
    for kc=1:3
        xmm = min(fovm(k:k+1,kc));
        xmx = max(fovm(k:k+1,kc));
        fov(kc,1) = xmm(1);
        fov(kc,2) = xmx(1);
    end
    vertx{fix((k-1)/2+1)} = [fov(:,1), [fov(1:2,1); fov(3,2)], [fov(1,1); fov(2,2); fov(3,1)], [fov(1,1); fov(2:3,2)], ...
                              [fov(1,2); fov(2:3,1)], [fov(1,2); fov(2,1); fov(3,2)], [fov(1:2,2); fov(3,1)], fov(:,2)];
    fovp.sgrid{fix((k-1)/2+1)} = {[fov(1,1):rez(1):fov(1,2)], [fov(2,1):rez(2):fov(2,2)], [fov(3,1):rez(3):fov(3,2)]};
end
fnum = length(vertx);       %  Number of FOV region blocks

%  Plot mics with FOV regions
figure(1)
plot3(fovp.mp(1,:),fovp.mp(2,:),fovp.mp(3,:),'xb')
hold on
for k=1:cl
    text(fovp.mp(1,k)+.25,fovp.mp(2,k)+.25,fovp.mp(3,k)+.25,int2str(k))
end

xlabel('X Meters');
ylabel('Y Meters');
zlabel('Z Meters');
title(['Microphone Positions with Sound Speed = ' num2str(fovp.c) ' m/s'])
grid on
pause(1)


%  Skech layout of microphone and FOV regions before adding CFAR Threshold
for kv=1:fnum
    %  Y lines
    plot3(vertx{kv}(1,1:2),vertx{kv}(2,1:2),vertx{kv}(3,1:2),'k.-')
    plot3(vertx{kv}(1,3:4),vertx{kv}(2,3:4),vertx{kv}(3,3:4),'k.-')
    plot3(vertx{kv}(1,5:6),vertx{kv}(2,5:6),vertx{kv}(3,5:6),'k.-')
    plot3(vertx{kv}(1,7:8),vertx{kv}(2,7:8),vertx{kv}(3,7:8),'k.-')
    %  Z lines
    plot3(vertx{kv}(1,1:2:3),vertx{kv}(2,1:2:3),vertx{kv}(3,1:2:3),'k.-')
    plot3(vertx{kv}(1,2:2:4),vertx{kv}(2,2:2:4),vertx{kv}(3,2:2:4),'k.-')
    plot3(vertx{kv}(1,5:2:7),vertx{kv}(2,5:2:7),vertx{kv}(3,5:2:7),'k.-')
    plot3(vertx{kv}(1,6:2:8),vertx{kv}(2,6:2:8),vertx{kv}(3,6:2:8),'k.-')
    %  X linesk
    plot3(vertx{kv}(1,1:4:5),vertx{kv}(2,1:4:5),vertx{kv}(3,1:4:5),'k.-')
    plot3(vertx{kv}(1,2:4:6),vertx{kv}(2,2:4:6),vertx{kv}(3,2:4:6),'k.-')
    plot3(vertx{kv}(1,3:4:7),vertx{kv}(2,3:4:7),vertx{kv}(3,3:4:7),'k.-')
    plot3(vertx{kv}(1,4:4:8),vertx{kv}(2,4:4:8),vertx{kv}(3,4:4:8),'k.-')
end
hold off
pause(.1)  %  Pause to let graphics operate

% Add on CFAR Threshold Detection
for k=1:2:fovcnt
    for kc=1:3
        xmm = min(fovm(k:k+1,kc));
        xmx = max(fovm(k:k+1,kc));
        fov(kc,1) = xmm(1)-fovp.threshneigh(kc)*rez(kc);  %  Add extra to account for CFAR averaging
        fov(kc,2) = xmx(1)+fovp.threshneigh(kc)*rez(kc);  %  Add extra to account for CFAR averaging 

    end
    vertx{fix((k-1)/2+1)} = [fov(:,1), [fov(1:2,1); fov(3,2)], [fov(1,1); fov(2,2); fov(3,1)], [fov(1,1); fov(2:3,2)], ...
                              [fov(1,2); fov(2:3,1)], [fov(1,2); fov(2,1); fov(3,2)], [fov(1:2,2); fov(3,1)], fov(:,2)];
    fovp.sgrid{fix((k-1)/2+1)} = {[fov(1,1):rez(1):fov(1,2)], [fov(2,1):rez(2):fov(2,2)], [fov(3,1):rez(3):fov(3,2)]};
end

% Assign wavefile recording filename to input data structure
fovp.fn = fn;   %  file name
fovp.pn = pn;   %  path make empty [] if in current working directory 

%  Detection significant sound sources in each frame
strms1 = signifslist(fovp);              % Find sources in selected regions  
save([fn(1:end-4) '.mat'],'strms1', 'fovp') % Save scipt of all detected sources