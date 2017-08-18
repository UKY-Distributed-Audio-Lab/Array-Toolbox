% This script shows an example of simulating a cocktail party recording 
% from individually recorded signals.  The microphone positions and the 
% source positions for each individual recording can be set by the 
% simulation. The script reads in 6 wave files, generates a planar 
% distribution of microphones in a 3.6 by 3.6 by 2.2 meter room with mild 
% room reverberation, plots the mic and source positions, and plays 4 
% random mic recordings of the simulated cocktail party.  Other options 
% exist for the simulation and playback:
%
% Users can set:
%   1) reverbEn - room reverberation enable
%   2) pflag - play back 4 random (example) signals
%   3) sflag - save all channels
%   4) fsre - (optional) resampling rate
%   5) spos - source positions
%   6) mpos - mic positions
%   7) (optional) other room and recording parameters in data structure recinfo
%      SEE simarraysig.m and simarraysigim.m
%
% Required files from Array Toolbox:
%   1) cocktailp.m
%   2) delayt.m         
%   3) imagesim.m 
%   4) regmicsplane.m 
%   5) roomimpres.m     
%   6) simarraysig.m    
%   7) simarraysigim.m  
%   8) wav2sig.m
%   9) SpeedOfSound.m
%  10) atmAtten.m

%
%   Required Data files:  man1.wav, man2.wav, man3.wav, woman1.wav,
%                         woman2.wav,woman3.wav
%
%
% Type 'runCocktailp' at the MATLAB prompt or Debug->Run
%
% Written by Satoru Tagawa (staga2@uky.edu) and Kevin D. Donohue 6/23/2008
% Updated 6/4/2013


clear

% SET BY USER *************************************************************
% With reverberation (1) or without (0) room reverberation
reverbEn = 1;  %  Note simulation take longer with room reverberations

% Play flag. Set to 1 to play back 4 (or all if N less than 4) channels of sigout at
% end of script
pflag = 1;

% Save flag. Set to 1 to save multichannel recording wave file: sigout.wav
sflag = 1;  

% Sampling frequency to be used throughout (optional)
%   otherwise the defaults is that of the lowest sampled wavefile
fsre = 16000;  % IN Hz COMMENT THIS LINE TO Use SAMPLING RATE OF ORIGINAL FILES

% Single speaker wavefiles for creating cocktail party sound
wavefiles = {'man1.wav';'woman1.wav'; 'woman3.wav'};
nsource = length(wavefiles);
% Read in wavefiles and resample signals if requested
if exist('fsre','var')==1
  sigin = wav2sig(wavefiles,fsre,[0 10]);  %  Read in over requested range (first 10 seconds)
  fs = fsre;  %  Reassign sampling rate
else
  [sigin, fs] = wav2sig(wavefiles,[0 10]);  %  Read in over requested range (first 10 seconds)
end

corners = [0 0 0; 3.6 3.6 2.2]';  %  opposite corners of room coordinates in meters (x,y,z)
% Random location generation of speaker/source positions in the same horizontal plane
spos = [(corners(1,1)+ (corners(1,2)-corners(1,1))*rand(1,nsource)); ...
        (corners(2,1)+ (corners(2,2)-corners(2,1))*rand(1,nsource)); ...
        (corners(3,1)+ (corners(3,2)-corners(3,1))*rand(1,nsource))];
disp('Speaker/Source Positions:')
disp(spos)

% Set position of microphones
% Mic position cannot be at the exact same position as the source position
%  Use function to automatically generate a planar geometry in ceiling 
micplane  =[0 0 2.1; 0 3.6 2.1; 3.6 3.6 2.1]';  %  Define 3 points for mic plane
micspacing = 1.2;
mpos = regmicsplane(micplane, micspacing);  % generate mics rectilinearly in a plane
[dimnum, micnum] = size(mpos);  %  Get number of mics in micnum
if reverbEn == 1
    % Set parameters for simarraysigim here
    reflecc = [.9, .8, .95, .8, .7, .4];  %  Reflection coefficients of walls, ceiling and floor
    recinfo = struct ('fs',fs,'corners',corners,'reflecc',reflecc);
else
    % Set parameters for simarraysig here
    recinfo = struct ('fs',fs);
end

dlmwrite('mpos.txt',mpos); % write mic postions to a .txt file

%  Set enviornmental parameters affecting attenuation and speed of sound
%recinfo.temp = 22; % in centigrade
%recinfo.press = 29.92; % in mm of Hg
%recinfo.hum = 38; % percent humidity

% Prompt for Environment Parameters
prompt={'Temperature (degrees C):','Pressure (mmHg):','Relative Humidity (percent):'};
name='Parameters to Determine Sound Speed';
numlines=1;
defaultanswer={'25.6','29.095', '48.4'};
answer=inputdlg(prompt,name,numlines,defaultanswer,'on');
recinfo.temp = str2double(answer{1});        % Temperature
recinfo.press = str2double(answer{2});       % Pressure
recinfo.hum = str2double(answer{3});         % Relative Humidity

% write temp, press, and hum info to .txt files
dlmwrite('temp.txt',recinfo.temp);
dlmwrite('press.txt',recinfo.press);
dlmwrite('hum.txt',recinfo.hum);
% ************************************************************************

% Display whether room reverb is on or off
if reverbEn == 1
    disp('ROOM REVERB IS TURNED ON');
else
    disp('ROOM REVERB IS TURNED OFF');
end

% Output source positions to screen for obseration
[nRs,nCs] = size(spos);
figure(1)
plot3(mpos(1,:),mpos(2,:),mpos(3,:),'bo')
hold on
plot3(spos(1,:),spos(2,:),spos(3,:),'rx')
hold off
axis([min(corners(1,:)) max(corners(1,:)) min(corners(2,:)) max(corners(2,:)) min(corners(3,:)) max(corners(3,:))])
title('Mic position denote by blue Os, and source positions denoted by red Xs')
xlabel('X-Dimension in meters')
ylabel('Y-Dimension in meters')
zlabel('Z-Dimension in meters')
grid on
pause(.5);  % Pause to let Matlab generate Figure before simulation

%  Simulate cocktail party recording
sigout = cocktailp(sigin, spos, mpos, recinfo);

% Play (up to) 3 random channels of sigout
if pflag == 1
    
    if micnum < 3   % Number of mics may be less than 4
        mn = 1:micnum;
    else    % Create 4 random numbers for mic playback
        [x,mn] = sort(rand(1,micnum));
        mn = mn(1:3);
    end
    
    % Play the simulated input to the mics
    for i=1:length(mn)
        disp(['Playing the input of the mic located at (', ...
            num2str(mpos(1,mn(i))), ', ', ...
            num2str(mpos(2,mn(i))), ', ', ...
            num2str(mpos(3,mn(i))), ')']);
        soundsc(sigout(:,mn(i)),fs);
        pause(length(sigout(:,mn(i)))/fs);
    end
    
end % Play

% Save all channels of sigout
if sflag == 1
    wavwrite(sigout,fs,'sigout.wav');
end