%  This script reads in an audio file segment of a phrase about a scare crow
%  and uses the simarraysigim function to simulate this word spoken over an
%  array in a room with reverberation.  It compares the recordings between
%  a near and far mic to compare the impact of the room response (the
%  further the mic from the speaker the more room effects can be heard).
%  This script requires the programs "simarraysigim.m," "regmicsline.m,"
%  "roomimpres.m," and wave file scare_crow.wav
%
%    Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005 (updated
%    June 2008)

[y,fs]= audioread('scare_crow.wav');  % Read in sample sound source to simulate
                                    % array recording
                                    
%  Reflection coefficients of 4 walls, floor, and ceiling.
refcoef = [.90, .95 .90 .95 .2 .4];
%  Coordinates of opposite corners of a rectangular room
foroom = [-7 -3 0; 1 3 1.24]'; 
%  Frequency dependent Attenuation
temp = 28; % Temperature centigrade
press = 29.92; % pressure inHg
hum = 80;  % humidity in percent
dis = 1;  %  Distance in meters (normalized to 1 meter)
prop.freq = fs/2*[0:200]/200;  %  Create 201 point frequency axis
prop.atten =  atmAtten(temp, press, hum, dis, prop.freq);  %  Attenuation vector
prop.c = SpeedOfSound(temp,hum,press);
c = prop.c;


%  Create mic array, 5 mics (1 meter spacing) on the x-axis (2=y and 1=z)
%  starting at -5 meters with a spacing of 1 meter in the positive
%  x direction
fom = [-5 2 1; 0 2 1]';              %  End coordinates of array
mpos = regmicsline(fom,1);      %  Generate linear microphone array
[dimnum, micnum] = size(mpos);  %  Get number of dimensions and mics
%  Set signal position near far negative end of array and 1 meter in front of
%  the array element
sigpos = [-5; 1; 0];
%  Simulate array recording
[sigout, tax] = simarraysigim(y, fs, sigpos, mpos, foroom, refcoef, prop);

% Plot the signals received by the first and last mic in array
%compute delay between the first and last mic
d1 = sqrt(sum((mpos(:,1)-sigpos).^2))/prop.c;
dend = sqrt(sum((mpos(:,micnum)-sigpos).^2))/prop.c;
delay = dend-d1;  %  Actual delay based on distances from source
figure(1)
plot([0:length(y)-1]/fs,y/10,'g-', 'Linewidth', 1)
hold on
plot(tax,sigout(:,1), 'r--',  'Linewidth', 2)
plot(tax, sigout(:,micnum), 'k:', 'Linewidth', 2)
hold off
legend({'Original', 'First Mic', 'Last Mic'})
xlabel('seconds')
ylabel('signal amplitudes')
title(['Original and Microphone signals (first and last) with a delay between them of ' num2str(delay*1000) 'ms'])

% Play original sound, recording at first, and finally from the last mic 
soundsc([y/10; sigout(:,1); sigout(:,micnum)],fs)