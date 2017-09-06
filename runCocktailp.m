% This script shows an example of simulating a cocktail party recording 
% from individually recorded signals (single channel wav files).  The
% wave files, microphone and source positions, and room parameters can be
% for the simulation by changing lines in this script. 
% For the default example, this script reads in 6 wave files (identified by
% file name strings in a cell array). A regular planar microphone distribution
% is generaged over a 3.6 by 3.6 meter ceiling (2.1 meters high).
% Wall refelction coeficient are set to create mild room reverberation (<.8).
% Pressure, temperate, and humidity are set to result in a given sound speed.
% These setting and other factors can be changed by the user through studying
% the script and reading comments.  Plots of the mic and source positions
% are generated to confirm mic placments, and randomly selected mic
% recordings are played to assess the mix and reverb effects.
% The following simulation options can be changed with relative ease:
%
% Users can set:
%   1) reverbEn - enable/disable room reverberation and set wall reflection
%                 coefficients
%   2) pflag - enable/disable playback of 4 randomly selected mic signals
%   3) saveVarsMat - flag for option on how to save results
%   4) fsre - (optional) resampling rate (down or up sample from original rate)
%   5) spos - source positions (x,y,z coordinates in colums of matrix, where
%             each column corresponds to wave file name in cell array).
%   6) mpos - mic positions
%   7) (optional) other room and recording parameters in data structure recinfo
%      SEE simarraysig.m and simarraysigim.m
%
%  Outputs are a multichannel wave file (sigout.wav) and depending one
%  the setting for saveVarsMat, either a mat file with critical or sperate
%  text files with processing parameters information (mic positions, speaker
%  locations, speed of sound, and labels for each sound source) 
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
% Written by Satoru Tagawa (staga2@uky.edu) and Kevin D. Donohue 6/23/2008
% Updated 6/4/2013
% Updated 7/31/2017 by Grant Cox and Kevin D. Donohue

clear
% SET BY USER *************************************************************
% With reverberation (1) or without (0) room reverberation
reverbEn = 1;  %  Note simulation will take longer with room reverberations
% If reverb enabled, set reflection coeffients for rectangular room in
%  vector below where [Wall1, Wall2, Wall3, Wall4, ceiling, floor]
reflecc = [.75, .75, .75, .75, .5, .5]; %  Reflection coefficients of walls, ceiling and floor

% If running this simulation for use with the cocktail party toolbox beamforming
% analysis script (partyanalysis), it will be useful to save the critical results
% to a .mat file that can be quicky opened by the script. Setting the saveVarsMat
% flag to 1 will automatically save the array recording, source postions, and 
% microphone postions to a single mat file that is open by the partyanalysis
% script, which by default is set to open the mat file.  If not set the
% array recording will be save to a multichannel wavefile and source and
% mic postion separate text files. You must enable the loading of this file 
% inpartyanalysis to use this feature.
saveVarsMat = 1;  % Set to 1 to save all output variable in mat file

% Play flag. Set to 1 to play back 3 (or all if N less than 3) channels of sigout at
% end of script
pflag = 1;

% Sampling frequency set through resampling the originals (optional).
% Otherwise, if commented out the defaults is that of the lowest sampled wavefile
fsre = 16000;  % forced sampling frequency in Hz

% Single speaker wavefiles for creating cocktail party sound.  If they
% were not all recorded at the same sampling rate, you must set parameter
% "fsre" to resample all to the same rate.
wavefiles = {'man1.wav', 'man2.wav', 'man3.wav', 'woman1.wav', ...
             'woman2.wav', 'woman3.wav'};    

Secondsin = 10;  % Number of seconds to read in data (suggest making this 
%                 short (~5 seconds) while debugging or adjusting parameters
%                 to speed up the modify and test cycle)

%  Size of room for limit on scanning for sources 
corners = [0 0 0; 3.6 3.6 2.2]';  %  opposite corners of room coordinates in meters (x,y,z)
% Set position of each speaker/source in meters [x, y, z].  Note matrix is
% transposed, so each line below becomes a column in SPOS.  Number of lines
% should equal the number of wave files loaded.
spos = [1.1  1.4  1.5; ...
        0.9  1.2  1.3; ...
        2.1  0.6  1.1; ...
        2.6  0.8  1.3; ...
        3.1  3.0  0.9; ...
        3.4  3.5  1.4]';
    
%--------------------------SET MICROPHONE POSITIONS------------------------
% Set position of microphones (default uses toolbox funtion to generate
% regular planar patterns, MPOS matrix can also be entered manually, where
% each column is the X,Y, Z postion of each microphone )
% Mic position cannot be at a source position
% Use toolbox function to automatically generate a planar geometry in ceiling 
micplane  =[0 0 2.1; 0 3.6 2.1; 3.6 3.6 2.1]';  %  Define 3 points for mic plane
micspacing = 1.3;  % Spacing for rectangular grid in meters
mpos = regmicsplane(micplane, micspacing);  % generate mics rectilinearly in a plane


%--------------------------READ IN WAVEFILES-------------------------------
% Read in wavefiles to matrix columns and resample signals if requested and
% trim to shortest length signal or requested end point.
if exist('fsre','var')==1
  sigin = wav2sig(wavefiles,fsre,[0 Secondsin]);  %  Read in over requested range (first 10 seconds)
  fs = fsre;  %  Reassign sampling rate
else
  [sigin, fs] = wav2sig(wavefiles,[0 Secondsin]);  %  Read in over requested range (first 10 seconds)
end

%----------------------------------REVERB----------------------------------
if reverbEn == 1
    % Set parameters for simarraysigim here
    recinfo = struct ('fs',fs,'corners',corners,'reflecc',reflecc);
else
    % Set parameters for simarraysig here without reverb (set refectivity of walls to
    % 0)
    reflecc = [0, 0, 0, 0, 0, 0];  %  Reflection coefficients of walls, ceiling and floor
    recinfo = struct ('fs',fs,'corners',corners,'reflecc',reflecc);
end

[dimnum, micnum] = size(mpos);  %  Get number of mics in micnum
if reverbEn == 1
    % Set parameters for simarraysigim here
    reflecc = [.8, .8, .8, .8, .7, .4]; % [.95, .95, .95, .95, .5, .5];  %  Reflection coefficients of walls, ceiling and floor
    recinfo = struct ('fs',fs,'corners',corners,'reflecc',reflecc);
else
    % Set parameters for simarraysig here without reverb (set refectivity of walls to
    % 0)
    reflecc = [0, 0, 0, 0, 0, 0];  %  Reflection coefficients of walls, ceiling and floor
    recinfo = struct ('fs',fs,'corners',corners,'reflecc',reflecc);
end

%  Set environmental parameters affecting attenuation and speed of sound
recinfo.temp = 22; % in centigrade
recinfo.press = 29.92; % in inches of Hg
recinfo.hum = 38; % percent humidity
c = SpeedOfSound(recinfo.temp,recinfo.hum,recinfo.press); % Sound Speed


%---------------------------DISPLAY REVERB SETTING-------------------------
if reverbEn == 1
    disp('ROOM REVERB IS TURNED ON');
else
    disp('ROOM REVERB IS TURNED OFF');
end

%--------------------------DISPLAY SOURCE POSITIONS------------------------
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


%---------------------------RUN SIMULATION---------------------------------
sigout = cocktailp(sigin, spos, mpos, recinfo);


%--------------------------SAVE VARIABLE DATA------------------------------
% Save variables parameters need for use in beamforming in analysis script
% partyanalysis.m if flag is set to 1.
if saveVarsMat == 1
    save('CocktailPartySimulation.mat','wavefiles','spos','mpos','c');
else % If not save multichannel recording, speaker position, and mic postion
     % and list of original wavefile source names in separate files
    fileID = fopen('wavefiles.txt','w');  % Files with wavefile names for labels
    for n=1:length(wavefiles)
        curr = char(wavefiles(1,n));
        fprintf(fileID,'%s',curr);
        if n ~= length(wavefiles)
            fprintf(fileID,',');
        end
    end
    % Create/open file to save simulated speaker/source positions - text file
    fileID = fopen('spos.txt','w');
    sposSize = size(spos);
    for n=1:sposSize(1)
        for j=1:sposSize(2)
            fprintf(fileID,'%f',spos(n,j));
            if j ~= sposSize(2)
                fprintf(fileID,',');
            else
                fprintf(fileID,'\n');
            end
        end
    end
    fclose(fileID);
    % Create/open file to save simulated mic positions - text file
    fileID = fopen('mpos.txt','w');
    mposSize = size(mpos);
    for n=1:mposSize(1)
        for j=1:mposSize(2)
            fprintf(fileID,'%f',mpos(n,j));
            if j ~= mposSize(2)
                fprintf(fileID,',');
            else
                fprintf(fileID,'\n');
            end
        end
    end
    fclose(fileID);    
    % Create/open file to save speed of sound
    fileID = fopen('csim.txt','w');
    fprintf(fileID,'%f',c);
    fclose(fileID);
end

% Save multichannel wave file from simulation
energy = std(sigout);  % Get RMS values from array recording
max_energy = max(energy);  % Use max RMS value for normalizing before saving
sigout = sigout / (5 * (max_energy(1) + eps)); % scale down to limit clipping
audiowrite('sigout.wav', sigout,fs);  % Save multichannel wave file


%-------------------------------PLAYBACK-----------------------------------
% Play (up to) 3 random channels of sigout
if pflag == 1
    
    if micnum < 3   % Number of mics may be less than 4
        mn = 1:micnum;
    else    % Create 3 random numbers for mic playback
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
        pause(Secondsin);
    end
    
end % Play