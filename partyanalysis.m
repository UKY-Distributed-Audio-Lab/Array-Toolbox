% edited by Grant Cox, 08/18/2017
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


% -------------------------------Setup------------------------------------
% Initialize variables, open sound, define constants
clear;

%--set loadVarsMat flag to 1 if you set saveVars in runCocktailp.m to 1. This
%will load the variables from the .mat file into the workspace for use 
%in this analysis script.
%--set loadVarsMat flag to 0 if you would like to import the parameters
% source postions, mic positions, sound speed and source name data from .txt
% files.  This is useful if you have experimental data and can create these
% files manually from the experimental recording . The load statements are
% in the "Load mic and speaker positions" section
loadVarsMat = 1;

% SET SPEAKER OF INTEREST. Choose one speaker from the wavefiles cell array. 
% The integer corresponds to a speaker in the array.
speaker_of_interest_index = 3;


tWin = 80e-3;  % Window size for block processing
iw = 0;  %  Set flag to 1 for inverse distance weighting on delay and sum BF


%----------------------Load mic and speaker positions----------------------
% The load function brings in 3 variables saved from running
% runCocktailp.m: wavefiles(cell array of strings), mpos(3 x Y matrix of
% doubles), and spos (3 x Z matrix of doubles)
fName = 'sigout.wav';  % output from runCocktailp.m simulation
source_info = audioinfo(fName);
sigSize = [source_info.TotalSamples, source_info.NumChannels];
fs = source_info.SampleRate;
Num_samples = sigSize(1);  % Total samples in each track
num_mics = sigSize(2);  % Number of microphones (should be same as channels)

if loadVarsMat == 1
    load('CocktailPartySimulation.mat');
else
    wavefilesraw = fileread('wavefiles.txt');
    wavefiles = strsplit(wavefilesraw,',');
    c = str2num(fileread('csim.txt')); % Sound speed
    spos = load('spos.txt'); % source positions
    mpos = load('mpos.txt'); % mic positions
end

%  Convert window size and increment to samples
nWin = round(tWin*fs);  % Audio window size in samples
if nWin/2 ~= fix(nWin/2)  % Ensure samples are even for overlap and add
    nWin = nWin+1;
end
nInc = round(nWin/2);  % Window increment %50 overlap

for p=1:length(wavefiles)  % Iterate over every wavefile name
      speakers.wavefiles{p} = spos(:,p);  % Set speaker location
end

speaker_of_interest = wavefiles{speaker_of_interest_index};

%--------------Output source positions to screen for observation-----------
figure(1)
plot3(mpos(1,:),mpos(2,:),mpos(3,:),'bo')
hold on
for n=1:length(wavefiles)
    %PLOT SOURCE POSITIONS ON FIGURE
    plot3(speakers.wavefiles{n}(1),speakers.wavefiles{n}(2),speakers.wavefiles{n}(3),...
        'g>','MarkerSize',14);
end

%--------------------------Add labels to figure----------------------------
hold off
title({'Mic position denote by blue Os, and source of interest denoted by red Xs', ...
        'and interfering source denoted by green triangles'})
xlabel('X-Dimension in meters')
ylabel('Y-Dimension in meters')
zlabel('Z-Dimension in meters')
grid on


%---------Find closest mic to speaker of interest for comparison-----------
micDist = sum((mpos - speakers.wavefiles{speaker_of_interest_index} * ones(1,num_mics)).^2);
[dum, cindex] = min(micDist);

%--------------------------Read source .wav file---------------------------
[yin, fs] = audioread(fName);
yclose = yin(:,cindex(1));


%window editing
hwin = hann(nWin+1);  %  Tappering window for overlap and add
hwin = hwin(1:end-1);  % Make adjustment so even windows align

%---------------------------------DSB--------------------------------------
x = zeros(nWin, num_mics);  % Current window of audio data
yDsb = zeros(Num_samples, 1); % delay-sum beamformer output
hwt = waitbar(0,'Beamformer: DSB');
for n=1:nInc:Num_samples-nWin  % iterate over 20ms windows
    % For each 20ms window save the previous and open the new
    xPrev = x;  x = yin(n:n+nWin-1,:);
    dum = dsb(x, xPrev, fs, speakers.wavefiles{speaker_of_interest_index},...
        mpos, c,iw); % run
    yDsb(n:n+nWin-1)  = yDsb(n:n+nWin-1) + dum.*hwin;
    waitbar(n/Num_samples,hwt);
end
close(hwt)



%--------------------------------Masking-----------------------------------
hwt = waitbar(0,'Beamformer: DSB with TF Masking');
yDsbMask = [];
for p=1:length(wavefiles)
    x = zeros(nWin, num_mics);  % Current window of audio data
    yt = zeros(Num_samples, 1); % delay-sum beamformer output
    for n=1:nInc:Num_samples-nWin  % iterate over 20ms windows
        % For each 20ms window save the previous and open the new
        xPrev = x;  x = yin(n:n+nWin-1,:);
        %DSB to one source per loop (the for loop with variable p)
        dum = dsb(x, xPrev, fs, speakers.wavefiles{p}, mpos, c); % run
        yt(n:n+nWin-1) = yt(n:n+nWin-1) + dum.*hwin;
        waitbar((n-1)/4 + n/(4*Num_samples),hwt);
    end
    %concactenate
    yDsbMask = [yDsbMask,yt];
end
yMask = tfmask(yDsbMask,speaker_of_interest_index,tWin,4,fs);
close(hwt);


%----------------------------Traditional GJBF------------------------------

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
x = zeros(nInc, num_mics);  % Current window of audio data
b = zeros(nInc, 1);  % embedded DSB output
z = zeros(nInc, num_mics-1);  % BM output
bmWall = [];  % BM LMS taps (not needed here still but need [])
mcWall = [];  % MC LMS taps (initialize to [])
yGjbf = zeros(Num_samples, 1); % traditional GJBF output

hwt = waitbar(0,'Beamformer: GJBF');
for n=1:nInc:Num_samples-nWin  % iterate over 20ms windows     
    xPrev = x;  x = yin(n:n+nWin-1,:); % load audio
    % Call beamforming function for this window
    [dum, bmWall, mcWall, snrAll, b, z] = ...
        gjbf(x, fs, speakers.wavefiles{speaker_of_interest_index}, ...
             mpos, c, p, q, mu, order, beta, ...
             phi, psi, K, xPrev, b, z, bmWall(:,:,end), ...
             mcWall(:,:,end), snrThresh, snrRate, snrInit, ...
             bmWForce, mcWForce);
         yGjbf(n:n+nWin-1) = yGjbf(n:n+nWin-1) + dum.*hwin;
         waitbar(n/Num_samples,hwt);
end
close(hwt)




%-------------------------Robust GJBF (Hoshuyama)--------------------------
% Need to recycle initial/final values between windows here, too
p = 5;  q =10; % guess number of signal propagation samples across array
[phi, psi] = ccafbounds(mpos, fs, c, p, order);  % calculate CCAF bounds
K = .01; % MC LMS norm constaint
snrThresh = -10;  snrRate = 2;  % -10dB threshold, check twice/window
snrAll = -Inf*ones(1,num_mics);  % saved, recycled SNR values
x = zeros(nInc, num_mics);  % Current window of audio data
b = zeros(nInc, 1);  % embedded DSB output
z = zeros(nInc, num_mics);  % BM output
bmWall = [];  % BM LMS taps (initialize to [])
mcWall = [];  % MC LMS taps (initialize to [])
yRobustGsc = zeros(Num_samples, 1); % robust GSC output
hwt = waitbar(0,'Beamformer: Robust GSC');
for n=1:nInc:Num_samples-nWin  % iterate over 20ms windows
    xPrev = x;  x = yin(n:n+nWin-1,:); % load audio
    % Call beamforming function for this window
    [dum, bmWall, mcWall, snrAll, b, z] = ...
	gjbf(x, fs, speakers.wavefiles{speaker_of_interest_index}, ...
         mpos, c, p, q, mu, order, beta, phi, ...
	     psi, K, xPrev, b, z, bmWall(:,:,end), ...
	     mcWall(:,:,end), snrThresh, snrRate, snrAll(end, :), ...
	     bmWForce, mcWForce);
     yRobustGsc(n:n+nWin-1) = yRobustGsc(n:n+nWin-1) + dum.*hwin;
     waitbar(n/Num_samples,hwt);
end
close(hwt)




%---------------------------------PLAYBACK---------------------------------
%--------------------------------------------------------------------------
% Uncomment next 5 lines to save results for later listening
audiowrite('yClose.wav', yclose/(sqrt(2)*max(abs(yclose))), fs);
audiowrite('yDsb.wav', yDsb(1)/(sqrt(2)*max(abs(yDsb))), fs);
audiowrite('yGjbf.wav', yGjbf/(sqrt(2)*max(abs(yGjbf))), fs);
audiowrite('yRobustGsc.wav', yRobustGsc/(sqrt(2)*max(abs(yRobustGsc))), fs);
audiowrite('yMask.wav', yMask/(sqrt(2)*max(abs(yMask))), fs);

%------------------Play mic closest to speaker of interest-----------------
hwt = waitbar(.1,{'Playing closest mic recording to ', speaker_of_interest});
sound(yclose/(std(yclose)*10),fs)
pause(length(yclose)/fs)

%-------------------Play DSB with inverse weighting result-----------------
if iw == 1
   waitbar(.3,hwt,{'Playing DSB (inverse distance weighting) result for ', speaker_of_interest});
   soundsc(yDsb,fs)
else
   waitbar(.3,hwt,{'Playing DSB result for ', speaker_of_interest});soundsc(yDsb,fs)
end    
sound(yDsb/(std(yDsb)*10),fs)
pause(length(yDsb)/fs+.5)

%-------------------------Play Griffiths Jim result------------------------
waitbar(.5,hwt,{'Playing GJBF result for ',speaker_of_interest});
sound(yGjbf/(std(yGjbf)*10),fs)
pause(length(yGjbf)/fs+.5)

%---------------------------Play Robust GSC result-------------------------
waitbar(.7,hwt,{'Playing Robust GSC result for ', speaker_of_interest});
sound(yRobustGsc/(std(yRobustGsc)*10),fs)
pause(length(yRobustGsc)/fs+.5)

%----------------------------Play Masking Result---------------------------
waitbar(.9,hwt,{'Playing TF Maksing result for ', speaker_of_interest});
sound(yMask/(std(yMask)*10),fs)
pause(length(yMask)/fs+.5)
close(hwt)

%clear workspace
clearvars