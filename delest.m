function [d, cv] = delest(sig1,sig2,fs)
%  This function inputs 2 signal segments (SIG1 and SIG2) sampled
%  over the same time duration with sampling rate (FS):
%
%       [d, cv] = delest(sig1,sig2,fs)
%
%  The output of the function is:
%  D => the relative time delay of SIG1 with respect to SIG2.
%  CV => normalizeed value of peak (correlation coefficient)
%  This function uses a center clip at 20% of signal power to remove
%  low level noise, performs a cross correlation, finds the maximum
%  peak, and then uses a sinc interpolation to find effective time
%  delay between the signals  (at a resolution of 10 times the
%  sampling rate).
%
%    written by Kevin D. Donohue, April 5, 2006 (donohue@engr.uky.edu)

clpper = 20;  %  Set clipping threshold
upsamp = 10;  %  Increase in resolution for time lag peak
interpint = 10;   % Number of samples to either side of the max point
                  %  in low resolution cross correlation
                  

%  Detrend signals to remove linear and DC offsets
sig1 = detrend(sig1); 
sig2 = detrend(sig2);
%  Compute signal RMS values 
rms1 = std(sig1);
rms2 = std(sig2);
%  Compute thresholds for clipping
threshd1 = rms1*sqrt(clpper/100);
threshd2 = rms2*sqrt(clpper/100);
%  Zero out signal energy below threshold
xcind = find(abs(sig1) < rms1*threshd1);
sig1(xcind) = 0;
xcind = find(abs(sig2) < rms2*threshd2);
sig2(xcind) = 0;

%  Compute cross correlation between signals
[xc, lags] = xcorr(sig1,sig2);
[mv, iv] = max(xc);  %  find max value in low resolution cross-correlation
mxpt = iv(1);        %  find position of max value
%  Extract 10 points around maximum for upsampling and getting a higher
%  resolution delay estimate
st = max([1,(mxpt-10)]);
ed = min([length(xc),(mxpt+10)]);
segest = xc(st:ed);
% Apply a tapering window to limit edge effects
[r,c] = size(segest);
taperwin = hamming(length(segest));
if r>1   %  if column vector
    segest = segest.*taperwin;
else    % if row vector
    segest = segest.*taperwin';
end
%  Resample at higher rate
hrseg = resample(segest,upsamp,1);
[cv,iv] = max(hrseg);  %  find max on finer grid
cv = cv/(eps+(sqrt(sum(sig1.^2)*sum(sig2.^2))));  % Compute correlation coefficient
% add to beginging of extracted segment for actual delay
d = (lags(st) + (iv(1)-1)/upsamp)/fs;
%figure(1);
%plot([0:length(sig1)-1], sig1, 'r', [0:length(sig1)-1], sig2, 'b');
%ds=d*fs
%title(num2str(ds))
%pause