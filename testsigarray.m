%  This script reads in an audio file sement of the word "five"
%  in then uses the simarraysig function to Simulate this word spoken over an
%  array.  This script requires the programs "simarraysig.m," "regmicsline.m"
%  and wave file fivefemale.wav
%
%    Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005

[y,fs]= audioread('fivefemale.wav');
%  Frequency dependent Attenuation
temp = 28; % Temperature centigrade
press = 29.92; % pressure inHg
hum = 80;  % humidity in percent
dis = 1;  %  Distance in meters (normalized to 1 meter)
prop.freq = fs/2*[0:100]/100;  %  Create 100 point frequency axis
prop.airatten =  atmAtten(temp, press, hum, dis, prop.freq);  %  Attenuation vector
prop.c = SpeedOfSound(temp,hum,press);
% Set speed of sound
c = prop.c;

%  Create mic array, 5 mics (1 meter spacing) on the x-axis (0 y and 0 z)
%  starting at -2 meters with a spacing of 1 meter in the positive
%  direction
fom = [-2 2; 0 0];
mpos = regmicsline(fom,1);
[dimnum, micnum] = size(mpos);  %  Get number of dimensions and mics
%  Set signal position at far negative end of array and 1 meter in front of
%  the array
sigpos = [-2; 1];
%  in the first case do not use multi-path scatterers
[sigout, tax] = simarraysig(y, fs, sigpos, mpos, prop);

% Plot the signals recieve by the first and last mic in array

%compute delay between the first and last mic
d1 = sqrt(sum((mpos(:,1)-sigpos(:)).^2))/c;
d5 = sqrt(sum((mpos(:,5)-sigpos(:)).^2))/c;
delay = d5-d1;  %  Actuall delay based on distances from source
figure(1)
plot(tax,sigout(:,1), 'r', tax, sigout(:,5), 'k')
xlabel('seconds')
ylabel('signal amplitudes')
title(['First mic signal (red) and last mic signal (black) delay of ' num2str(delay*1000) 'millisecs'])


%  Now use multi-path scatterers
mpath = mpos-[0; -0.20]*ones(1,micnum);  % for each mic place a satterer behind it 20 cm to simulate a wall/scatterer behind the mic array
[rw, cl] = size(mpos);
mpath = [.7*ones(1,cl); mpath];  % give each scatterer a .7 refection coefficient
[sigout, tax] = simarraysig(y, fs, sigpos, mpos, prop);

%compute delay between the first and 5th mic
figure(2)
plot(tax,sigout(:,1), 'r', tax, sigout(:,5), 'k')
xlabel('seconds')
ylabel('signal amplitudes')
title(['Multipath scatterers: First mic signal (red) and last (black)'])
