
f.fn = 'shortcp14.wav';
f.pn = [];
tmp = mean([24.3, 24.2]);  %  Average Temperature measurement in Centigrade
rh = mean([63.2, 59.8]);  %  Average Relative Humidity measurement in percent
pres= 29.06;               %  Air pressure measurement in mmHg
f.c = SpeedOfSound(tmp,rh,pres);  % Speed of sound
f.fc = 300; %  High-pass filter cut-off for signal conditioning
f.fs = 8000;  %  Sampling frequency (will be resampled if not the original sampling rate)
load('mp28.mat');
mp = mp(:,1:2:end);
f.mp = mp;

%  Opposite corner points for room, also walls are locations for noise source placement
froom = [0 0 0; 3.57 3.57 2.26]';  % [x1, y1, z1; x2 y2 z2]'
%  Field of view for reconstructing SRP image (opposite corner points)
fov = [0 0 .57; 3.6 3.6 1.7]';
f.bata = .85;  % Beta for whitening
%  Time window for Frequency domain block processing
f.trez = 100e-3;  %  In seconds
%  Room Resolution: Step through cartesion grid for mic and sound source
%  plane
f.tinc = 50e-3;  %  Time increment in seconds
f.rez = .1;  %  In meters
%  Compute grid axis for pixel of the SRP image (reduce resolution in the
%  Z-direction)
f.sgrid = {{[fov(1,1):f.rez:fov(1,2)];[fov(2,1):f.rez:fov(2,2)]; [fov(3,1):3*f.rez:fov(3,2)]}}; 
f.graphicout_flag = 1;
f.pcfar = 1e-1;
strms1 = signifslist(f); 
%  Convert detection file format to one used by Hari's program
x = convtpos(strms1,f);
f.tinc = f.trez/2;
%  Now apply auditory stream identification using the proximity thresholds
timethr = 8;        % Time in seconds of tolerable nondetection (quiet period in speech) to link streams (1.1)
spacethr = .5;      %  Distance in meters where one source is allow(smaller detections are ignored) (.9)
toosmall = 0.15;       % Time in seconds to drop any stream shorter than this

tr = round(timethr/f.tinc) +1;       % Convert time to frame numbers
tsfrm = round(toosmall/f.tinc) +1;   % Convert time to frame numbers
stseg = imsegm(x,spacethr,tr,tsfrm);     % Hari's program to create true ASA matrix

strmnum = max(stseg(:,1));
% Create beamformed streams and save in separate files
[scriptout len] = bfsegment(f.fn,f.pn,[1:strmnum(1)],stseg,f);

if moviemake_flag == 1
    %  Create avi file to for movie of detected streams (must create a unique
    %  avi file name each time.)  Will not overwrite
    fratein = 1/(fovp.trez/2);                  % Frames per second on input
    frateout = 15;                              % Output frame rate
    fname = [basenam '.avi'];                   % File name of movie
    sfactor = 200;                              %  Scale detection stat value to marker size
    mv = mkmov3d(stseg,fratein,frateout,fname,sfactor);
end

