function [rt60, sigpow, nospow] = rt60est(snd,fs, fhc)
%  This function estimates the RT 60 time of a signal containing a recorded
%  sound vector, SND, with sampling rate FS.  In order to get accurate RT60 value
%  a broadband continous sound (i.e. white noise) must be played loud enough and
%  long enough so that the diffuse sound in the room has reached steady state.
%  the source should be about 1.5 to 3 meters away from the measurment mic so the 
%  direct path does not dominate the recording.  Then the sound must 
%  aburptly stop and recording contine until the sound falls below the noise floor.
%  Input vector SND must be a segment of the recorded signal where at least the first 20%
%  of the signal is the excitation sound and at least the last 20% of the signal must be the
%  noise floor.  The program uses the begining and ending parts of the signal to estimate
%  the signal power and noise floor power.  The roll-off of sound from the room reverberation
%  is found based on these 2 estimates.  The slope of the roll-off is estimate in dB per second
%  and extend to find the amount of time for a 60dB drop in sound.
%
%          [rt60, sigpow, nospow] = rt60est(snd, fs, fhc)
%
%  The inputs include the sound vector (SND), the sampling frequency (FS), and an
%  optional parameter FHC, which is the cut off frequency for a high-pass filter applied
%  to the data before the estimate.  The outputs are the RT60 time (RT60) in seconds, the
%   peak signal power (SIGPOW) in dB and noise floor power (NOSPOW) in dB.
%   THe estimate will always plot the lines fitted to the signal envelope to
%   ensure that reasonable fit to the data was achieved.
%
%
%   Kevin D. Donohue (donohue@engr.uky.edu) May 29, 2006
%


%  If optional parameter present, apply high-pass filter
if nargin ==3
    [b,a] = butter(4, 2*fhc/fs, 'High');
    snd = filtfilt(b,a,snd);
end
len = length(snd);
endest = .2*len;  % Number of seconds from end of segment from which to estimate
%  the signal and noise powers.

% Compute squared envelope of sound file
envsnd = abs(hilbert(snd)).^2;
%  Remove spike and valleys for a better envelop estimate with a 5ms
%  3rd quartile filter (75th percentile value) 
nwind = round(5e-3*fs);  % integer window length
rnk = max([round(.75*nwind),1]);  % Integer Rank of 75th percentile value
fenv = osfilt(envsnd,nwind,rnk);  %  Get order statistic for 75th percentile
                                  %  estimate in sliding window mode.
                                  %  This OSFILT function is included at
                                  %  the end of this program
fenv = fenv(nwind:end-nwind);  %  Trim off beging and ending window to limit edge effects         
len = length(fenv);
tp = [0:len-1]/fs;   % Compute time axis for trimed segments
lfenv = 10*log10(fenv); % Convert envelope to dB
% Get first 10% of segment to estimate signal power
sigon = lfenv(1:round(endest));
tsigon = tp(1:length(sigon)); 
%  Fit line to signal part of envelop
p = lfm(tsigon,sigon);  %  Line fit program included at the end of this program
%  Signal power estimate is the mean of the line over the signal portion
sigpow = mean(p(1)*tp(1:length(sigon))+p(2));
%  Trim the power envelop from the peak signal point on
%  get last 10% seconds of segment to estimate noise power
tsigoff = tp(len-round(endest):len);
sigoff = lfenv(len-round(endest):len);
%  Fit line to noise part of envelop
poff = lfm(tsigoff,sigoff);
%  Noise power estimate is the mean of the line over the noise portion
nospow = mean(poff(1)*tp(len-round(endest):len)+poff(2));
%  Find envelop values 5 db Below signal power 
%   and 5dB above noise power
ulim = sigpow - 5;
llim = nospow + 5;
%  Find all point 5 dB below the maximum power level
gg = find(lfenv < ulim);
%  find first point to reach 5 dB above the noise floor
trngg = find(lfenv(gg) < llim);
gg = gg(1:trngg(1));
trollon = tp(gg);
%  Fit line to roll-off
prt = lfm(trollon,lfenv(gg));
%  Based on slope estimate determine time
%  for sound to drop 60 dB after sound is shut off
%rt60 =(sigpow-60-prt(2))/prt(1) - trollon(1);
rt60 = -60/prt(1);

%  plot envelop and line fit for signal power
figure; plot(tp,lfenv,'g',tsigon, p(1)*tsigon+p(2)); hold on
plot(tsigoff, poff(1)*tsigoff+poff(2))
plot(tp(gg), prt(1)*trollon+prt(2),'k')
hold off
xlabel( 'seconds')
ylabel('dB')
title(['RT60 Estimate = ' num2str(rt60) 's based on line fit to envelope in dB'])


%  Assuming the point at which the RT60 decay falls below the noise floor 
%  is less than half the distance to the end, the median
%  would be a good estmate of the noise floor power


%  Line fit program that include trimming of peaks and other outlier
%  features
function p = lfm(x1,y1)
%  p = itlnfte(x,y)
%  This program will fit a straight line to the data
%  defined by abscisa point x and ordinant points y.
%  This function operate iteratively by fitting a line using
%  all the data points and then by eliminating 50% of the points
%  with the highest signed error values and redoing it.
%
%
%   Written by Kevin D. Donohue 6-25-1997
%   donohue@engr.uky.edu


pkf = polyfit(x1, y1, 1);           % First fit to get initial error
difl = (pkf(1)*x1 + pkf(2) - y1);  % Find error at each point
uind = find(abs(difl) < median(abs(difl))); % Find all points with signed error
x1 = x1(uind);
y1 = y1(uind);

if (isempty(x1) || isempty(y1))
    p(1)=0;p(2)=0;
else
    for k = 1:3
        pkf = polyfit(x1, y1, 1);  %  Fit line to censored set of points
        difl = (pkf(1)*x1 + pkf(2) - y1);  % Find error at each point
        %  less than median error (outliers with positive error)
        uind = find((difl) > median((difl))); % Find all points with absolute error
        %  greater than median error
        x1 = x1(uind);
        y1 = y1(uind);

        %uind = find(-difl < median(-difl)); % Find all points with signed error
        %  less than median error (outliers with negative error)
        %p = polyfit(x1(uind), y1(uind), 1);  %  Fit line to censored set of points
    end
end
p = pkf;


function gout = osfilt(x,n,r)
%  This function processes a signal with an order statistic (OS) filter.  The input
%  signal vector (X) is processed with a sliding window of sample length N and the 
%  sample with rank R (where 1 is the minimum and N is the maximum).
%
%        gout = osfilt(x,n,r)
%
%  Output:  GOUT is a vector of the same size as input X containing the order statistic
%  values.  If N is odd, symmetrical windowing is applied with respect to the ouput point.
%  If N even, symetry is not possible, so 2 consectutive outputs are averaged together 
%  and sent to the output.
% 
%  written by Kevin D. Donohue (donohue@engr.uky.edu) May 2006



len = length(x);  %  Get length of input signal
hn = fix(n/2);   %  Determine half the processing window size for symetric placement
%  Check for even or odd sliding window size
if hn ~= n/2     %   If window length is odd,
                 % a natural symetric between input and output exists
    gt = zeros(size(x));
    %  Initial window placement over first sample
    %   up to the point where window is within signal
    for k=1:hn
        wind = x(1:hn+k);  %  Get window from the begin
        gt = sort(wind);   %  rank
        rank = max([round((r/n)*(hn+k)),1]);  %  find proportinal rank
        gout(k) = gt(rank(1));  % Assign to output
    end
    %  Sample where window fits completely in signal
    rank = r;
    for k = hn+1:len-hn
        wind = x(k-hn:k+hn);
        gt = sort(wind);
        gout(k) = gt(r);
    end
    %  Final window placement where it extends beyond the end
    for k = len-hn+1:len
        wind = x(k:len);
        gt = sort(wind);
        rank = max([round((r/n)*(len-k+1)),1]);
        gout(k) = gt(rank);
    end
else   %  If N is even, average 2 consectutive ouputs to restore symetry
    gt = zeros(size(x));
    %  Get initial segemnt
    for k=1:hn
        wind = x(1:hn+k-1);
        gt = sort(wind);
        rank = max([round((r/n)*(hn+k-1)),1]);
        gout(k) = gt(rank(1));
    end
    rank = r;
    for k = hn+1:len-hn
        wind = x(k-hn:k+hn-1);
        gt = sort(wind);
        gout(k) = gt(r);
    end
    for k = len-hn+1:len+1
        wind = x(k-hn:len);
        gt = sort(wind);
        rank = max([round((r/n)*(len-k+hn+1)),1]);
        gout(k) = gt(rank);
    end
    [r,c] = size(gout);
    if r == 1
        kern = ones(1,2)/2;
    else
        kern = ones(2,1)/2;
    end
    gout = conv2(gout,kern,'valid');
end