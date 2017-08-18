%  This script reads in an audio file segment of a phrase about a scare crow
%  and uses the simarraysig function to simulate this word spoken over an
%  array.  It compares the recording with and without multipath scatterers
%  and also for the near and far mic.  This script requires the programs
%  "simarraysig.m," "regmicsline.m," "roomimpres.m," and wave file scare_crow.wav
%
%    Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005 (updated
%    June 2008)

clear
[y,fs]= audioread('scare_crow.wav');  % Read in sample sound source to simulate
                                    % array recording
%  Frequency dependent Attenuation
temp = 28; % Temperature centigrade
press = 29.92; % pressure inHg
hum = 80;  % humidity in percent
dis = 1;  %  Distance in meters (normalized to 1 meter)
prop.freq = fs/2*[0:200]/200;  %  Create 201 point frequency axis
prop.atten =  atmAtten(temp, press, hum, dis, prop.freq);  %  Attenuation vector
% Set speed of sound in propagation data structure for simarraysig.m
prop.c = SpeedOfSound(temp,hum,press);


%  Create mic array, 5 mics (1 meter spacing) on the x-axis (2=y and 0=z)
%  starting at -5 meters with a spacing of 1 meter in the positive
%  x direction
fom = [-5 2 1; 0 2 1]';              %  End coordinates of array
mpos = regmicsline(fom,1);      %  Generate linear microphone array
[dimnum, micnum] = size(mpos);  %  Get number of dimensions and mics
%  Set signal position near far negative end of array and 1 meter in front of
%  the array element
sigpos = [-7; 1; 1.5];
%  In the first case do not use multi-path scatterers (no "mpath" field in
%  data structure "prop")
[sigout, tax] = simarraysig(y, fs, sigpos, mpos, prop);

% Plot the signals received by the first and last mic in array
%compute delay between the first and last mic
d1 = sqrt(sum((mpos(:,1)-sigpos(:)).^2))/prop.c;
dend = sqrt(sum((mpos(:,micnum)-sigpos(:)).^2))/prop.c;
delay = dend-d1;  %  Actual delay based on distances from source
figure(1)
plot(tax,sigout(:,1), 'r', tax, sigout(:,micnum), 'k')
xlabel('seconds')
ylabel('signal amplitudes')
title(['First mic signal (red) and last mic signal (black) delay of ' num2str(delay*1000) 'millisecs'])



%  Now use multi-path scatterers
mpath = mpos+[0; 0.40; 0]*ones(1,micnum);  % for each mic place a satterer behind it (relative to speaker)
                                         % at 40 cm to simulate a wall/scatterer behind the mic array
[rw, cl] = size(mpos);
prop.mpath = [.8*ones(1,cl); mpath];  % give each scatterer a .8 reflection coefficient
[sigoutmp, tax] = simarraysig(y, fs, sigpos, mpos, prop);

%compute delay between the first and last mic
figure(2)
plot(tax,sigoutmp(:,1), 'r', tax, sigoutmp(:,micnum), 'k')
xlabel('seconds')
ylabel('signal amplitudes')
title('Multipath scatterers: First mic signal (red) and last (black)')

% Play sounds at first and last mic without multipath and
% then at first and last mic with multpath 
soundsc([sigout(:,1); sigout(:,micnum); sigoutmp(:,1); sigoutmp(:,micnum)],fs)
