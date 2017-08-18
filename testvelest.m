%  This script loads a measured data file recorded from a colinear
%  sound source and 2 microphones separated by .45 meters.  A white noise
%  burst was played through the sound source for aproximately 1 second. 

micdist = .45;  % Microphone distance in meters
[y, fs] = audioread('sound_speed_p45m.wav');  %  2 channel recorded sound with collinear configuration
%  Plot example waveform
figure(1); plot([0:length(y(:,1))-1]/fs, y(:,1))
title('Recorded White Noise Burst on Channel 1')
xlabel('Seconds')
ylabel('Amplitude')
%  Estimate velocities
[vest, cst]  = velest(y(:,1:2),fs,micdist);
% Identify good points that exceed correlation coefficient threshold 
hh = find(~isnan(vest));
%  Average all valid points to estimate the velocity 
v = sum(cst(hh).*vest(hh))/sum(cst(hh))
%  compute 95% confidence limits
con95lim = tinv(.975,length(hh)-1)*std(vest(hh))/length(hh)
%   Display results in text
disp(['Estimated sound velocity is ' num2str(v) ' m/s'])
disp(['with 95% confidence limits of ' num2str(v-con95lim) ' to ' num2str(v-con95lim)])

%  Plot all valid estimates
figure(2); plot(vest,'xr')
title('Etimated Velocities for each window with a .4 or greater correlation')
xlabel('Index of Sliding Window')
ylabel('m/s')

% Plot correlation coefficients
figure(3); plot(cst,'ob')
title('Correlation coefficients for each window')
xlabel('Index of Sliding Window')
ylabel('Correlation Coefficient')