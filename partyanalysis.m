% Cocktail Party Analysis Script
% This test script demonstates how to use the three beamforming
% functions dsb.m, gjbf.m, and time-frequecy masking tfmask.m.  To run, make
% sure that the following files can be found:
%
% * dsb.m
% * gjbf.m
% * bmlms.m
% * ccafbounds.m
% * tfmask.m
% * getmask.m
% * micpos.dat
% * srcpos.dat
% * cocktail11k.wav
%
% This script shows how to implement four beamforming algorithms:
%
% * Delay-Sum Beamform (DSB)
% * Traditional Griffiths-Jim Beamformer (GJBF)
% * Robust Generalised Sidelobe Canceller (Improved GJBF)
% * TF masking (post processing, our variant). 
%
% Note that there are four speakers in the room and Mike will be the
% target, and is talking about hockey. Phil is the primary interferer and he
% even closer to the Microphone closest to Mike and talking louder.
%
% Written by Phil Townsend (jptown0@engr.uky.edu) 6/2/08
% Updated by Kevin D Donohue August 2014 to remove correlation, apply taper
% window for overlap and add, apply the invese distance BF option to delay and summ
% and added time frequency masking


% Setup
% Initialize variables, open sound, define constants
clear
c = 345.6;  %  Speed of sound for recording
tWin = 80e-3;  % Window size for block processing
iw = 0;  %  Set flag to 1 for inverse distance weighting on delay and sum BF
fName = 'cocktail11k.wav';  % source file of cocktail party data
%  Get paramaters from wave file
g = audioinfo(fName);
sigSize = [g.TotalSamples, g.NumChannels];
fs = g.SampleRate;
N = sigSize(1);  % Total samples in each track
nWin = round(tWin*fs);  % Audio window size in samples
if nWin/2 ~= fix(nWin/2)  % Ensure samples are even for overlap and add
    nWin = nWin+1;
end
nInc = round(nWin/2);  % Window increment %50 overlap
M = sigSize(2);  % Number of microphones

% Load mic and speaker positions
m = load('micpos.dat')/100;  % cm -> m
sRaw = load('srcpos.dat')/100;  % speaker positions, cm -> m
people = {'mike' 'kate' 'phil' 'donohue'};  % all speakers
for p=1:length(people)  % Iterate over everyone in the party
    s.(people{p}) = sRaw(:,p);  % Set speaker location
end

% Output source positions to screen for observation
figure(1)
plot3(m(1,:),m(2,:),m(3,:),'bo')
hold on
plot3(s.mike(1),s.mike(2),s.mike(3),'rx', 'MarkerSize', 14)
plot3(s.kate(1),s.kate(2),s.kate(3),'g>', 'MarkerSize', 14)
plot3(s.phil(1),s.phil(2),s.phil(3),'g>', 'MarkerSize', 14)
plot3(s.donohue(1),s.donohue(2),s.donohue(3),'g>', 'MarkerSize', 14)

hold off
%axis([min(corners(1,:)) max(corners(1,:)) min(corners(2,:)) max(corners(2,:)) min(corners(3,:)) max(corners(3,:))])
title({'Mic position denote by blue Os, and source of interest denoted by red Xs', ...
        'and interfering source denoted by green triangles'})
xlabel('X-Dimension in meters')
ylabel('Y-Dimension in meters')
zlabel('Z-Dimension in meters')
grid on



%  Find closest mic to speaker of interest for comparison
mdd = sum((m - s.mike*ones(1,M)).^2);
[dum, cindex] = min(mdd);

[yin, fs] = audioread(fName);
yclose = yin(:,cindex(1));

hwin = hann(nWin+1);  %  Tappering window for overlap and add
hwin = hwin(1:end-1);  % Make adjustment so even windows align
% DSB
x = zeros(nWin, M);  % Current window of audio data
yDsb = zeros(N, 1); % delay-sum beamformer output
hwt = waitbar(0,'Beamformer: DSB');
for n=1:nInc:N-nWin  % iterate over 20ms windows
    % For each 20ms window save the previous and open the new
    xPrev = x;  x = yin(n:n+nWin-1,:);
    dum = dsb(x, xPrev, fs, s.mike, m, c,iw); % run
    yDsb(n:n+nWin-1)  = yDsb(n:n+nWin-1) + dum.*hwin;
    waitbar(n/N,hwt);
end
close(hwt)

% Masking
x = zeros(nWin, M);  % Current window of audio data
yDsb1 = zeros(N, 1); % delay-sum beamformer output
hwt = waitbar(0,'Beamformer: DSB with TF Masking');
for n=1:nInc:N-nWin  % iterate over 20ms windows
    % For each 20ms window save the previous and open the new
    xPrev = x;  x = yin(n:n+nWin-1,:);
    dum = dsb(x, xPrev, fs, s.mike, m, c); % run
    yDsb1(n:n+nWin-1)  = yDsb1(n:n+nWin-1) + dum.*hwin;
    waitbar(n/(4*N),hwt);
end
x = zeros(nWin, M);  % Current window of audio data
yDsb2 = zeros(N, 1); % delay-sum beamformer output
for n=1:nInc:N-nWin  % iterate over 20ms windows
    % For each 20ms window save the previous and open the new
    xPrev = x;  x = yin(n:n+nWin-1,:);
    dum = dsb(x, xPrev, fs, s.kate, m, c); % run
    yDsb2(n:n+nWin-1)  = yDsb2(n:n+nWin-1) + dum.*hwin;
    waitbar(1/4+n/(4*N),hwt);
end
x = zeros(nWin, M);  % Current window of audio data
yDsb3 = zeros(N, 1); % delay-sum beamformer output
for n=1:nInc:N-nWin  % iterate over 20ms windows
    % For each 20ms window save the previous and open the new
    xPrev = x;  x = yin(n:n+nWin-1,:);
    dum = dsb(x, xPrev, fs, s.phil, m, c); % run
    yDsb3(n:n+nWin-1)  = yDsb3(n:n+nWin-1) + dum.*hwin;
    waitbar(2/4+n/(4*N),hwt);
end
x = zeros(nWin, M);  % Current window of audio data
yDsb4 = zeros(N, 1); % delay-sum beamformer output
for n=1:nInc:N-nWin  % iterate over 20ms windows
    % For each 20ms window save the previous and open the new
    xPrev = x;  x = yin(n:n+nWin-1,:);
    dum = dsb(x, xPrev, fs, s.donohue, m, c); % run
    yDsb4(n:n+nWin-1)  = yDsb4(n:n+nWin-1) + dum.*hwin;
    waitbar(3/4+n/(4*N),hwt);
end
yMask = tfmask([yDsb1, yDsb2, yDsb3, yDsb4],1,tWin,4,fs);

close(hwt)




% Traditional GJBF
% Notice that several parameters must be saved and recycled between
% iterations to ensure that the final conditions from one audio
% window become the initial conditions for the next.
mu = .1;  order = 20;  beta = .9;  % LMS filter parameters
p = 0;  q = 0;  % don't need these now, set to zero
phi = [];  psi = [];  % CCAF bounds not needed
K = [];  % MC NLMS norm threshold not needed
bmWForce = [];  mcWForce = [];  % not locking taps right now
snrThresh = [];  snrRate = [];  % no SNR thresholding right now
snrInit = [];
x = zeros(nInc, M);  % Current window of audio data
b = zeros(nInc, 1);  % embedded DSB output
z = zeros(nInc, M-1);  % BM output
bmWall = [];  % BM LMS taps (not needed here still but need [])
mcWall = [];  % MC LMS taps (initialize to [])
yGjbf = zeros(N, 1); % traditional GJBF output

hwt = waitbar(0,'Beamformer: GJBF');
for n=1:nInc:N-nWin  % iterate over 20ms windows     
    xPrev = x;  x = yin(n:n+nWin-1,:); % load audio
    % Call beamforming function for this window
    [dum, bmWall, mcWall, snrAll, b, z] = ...
        gjbf(x, fs, s.mike, m, c, p, q, mu, order, beta, ...
             phi, psi, K, xPrev, b, z, bmWall(:,:,end), ...
             mcWall(:,:,end), snrThresh, snrRate, snrInit, ...
             bmWForce, mcWForce);
         yGjbf(n:n+nWin-1) = yGjbf(n:n+nWin-1) + dum.*hwin;
         waitbar(n/N,hwt);
end
close(hwt)


% Robust GJBF (Hoshuyama)
% Need to recycle initial/final values between windows here, too
p = 5;  q =10; % guess number of signal propagation samples across array
[phi, psi] = ccafbounds(m, fs, c, p, order);  % calculate CCAF bounds
K = .01; % MC LMS norm constaint
snrThresh = -10;  snrRate = 2;  % -10dB threshold, check twice/window
snrAll = -Inf*ones(1,M);  % saved, recycled SNR values
x = zeros(nInc, M);  % Current window of audio data
b = zeros(nInc, 1);  % embedded DSB output
z = zeros(nInc, M);  % BM output
bmWall = [];  % BM LMS taps (initialize to [])
mcWall = [];  % MC LMS taps (initialize to [])
yRobustGsc = zeros(N, 1); % robust GSC output
hwt = waitbar(0,'Beamformer: Robust GSC');
for n=1:nInc:N-nWin  % iterate over 20ms windows
    xPrev = x;  x = yin(n:n+nWin-1,:); % load audio
    % Call beamforming function for this window
    [dum, bmWall, mcWall, snrAll, b, z] = ...
	gjbf(x, fs, s.mike, m, c, p, q, mu, order, beta, phi, ...
	     psi, K, xPrev, b, z, bmWall(:,:,end), ...
	     mcWall(:,:,end), snrThresh, snrRate, snrAll(end, :), ...
	     bmWForce, mcWForce);
     yRobustGsc(n:n+nWin-1) = yRobustGsc(n:n+nWin-1) + dum.*hwin;
     waitbar(n/N,hwt);
end
close(hwt)

% Uncomment next 5 lines to save results for later listening
audiowrite('yClose.wav', yclose/(sqrt(2)*max(abs(yclose))), fs);
audiowrite('yDsb.wav', yDsb/(sqrt(2)*max(abs(yDsb))), fs);
audiowrite('yGjbf.wav', yGjbf/(sqrt(2)*max(abs(yGjbf))), fs);
audiowrite('yRobustGsc.wav', yRobustGsc/(sqrt(2)*max(abs(yRobustGsc))), fs);
audiowrite('yMask.wav', yMask/(sqrt(2)*max(abs(yMask))), fs);


%Play mic closest to speaker of interest
hwt = waitbar(.1,{'Playing closest mic recording', ' Male speaker of interest talking about hockey being', ...
                   'talked over by a speaker discussing snowboarding loudly'});
sound(yclose/(std(yclose)*10),fs)
pause(length(yclose)/fs)
%  Play DSB with inverse weighting result
if iw == 1
   waitbar(.3,hwt,{'Playing DSB (inverse distance weighting) result', ' Male Speaker of Interest Talking About Hockey'});soundsc(yDsb,fs)
else
   waitbar(.3,hwt,{'Playing DSB result', ' Male Speaker of Interest Talking About Hockey'});soundsc(yDsb,fs)
end    
sound(yDsb/(std(yDsb)*10),fs)
pause(length(yDsb)/fs)
%  Play Griffiths Jim result
waitbar(.5,hwt,{'Playing GJBF result', ' Male Speaker of Interest Talking About Hockey'});
sound(yGjbf/(std(yGjbf)*10),fs)
pause(length(yGjbf)/fs)
%  Play Robust GSC result
waitbar(.7,hwt,{'Playing Robust GSC result', ' Male Speaker of Interest Talking About Hockey'});
sound(yRobustGsc/(std(yRobustGsc)*10),fs)
pause(length(yRobustGsc)/fs)
%  Play Masking Result
waitbar(.9,hwt,{'Playing TF Maksing result', ' Male Speaker of Interest Talking About Hockey'});
sound(yMask/(std(yMask)*10),fs)
pause(length(yMask)/fs)
close(hwt)