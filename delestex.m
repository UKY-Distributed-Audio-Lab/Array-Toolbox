%  This script demonstrates the delay estimator (DELAYESTTM) toolbox program and
%  computes its performance as a function of SNR.  It performs a Monte Carlo simulation
%  of a pulse in noise.  The signal and simulation parameters can be set by editing the
%  initial lines of this script (see comments for details).  The signal
%  processing options for the delay estimator can also be set (see help on
%  the DELAYESTTM function.
%  
%    Script written by Kevin D. Donohue  (May 30, 2006)
%

runs = 10;  %  Number of simulations estimates.  Not too high if ploting in each loop pass
%  Signal parameters
fs = 44.1e3;   %  Sampling frequency
dl = .00371;  %  Simulated Delay in seconds
%  Generate target pulse bandwidths
f12p = 3800;  %  Corresponding upper frequency limit
f11p = 400;  %  Lower frequency limit
xsp = 20; %  Percent for stopband interval
f11s = f11p*(1-xsp/100);  %  Lower stopband
f12s = f12p*(1+xsp/100);  %  Upper stopband
%  White noise snr
wgnsnr = 20;  %  SNR in db
sclnos = 10^(-wgnsnr/20);
%  Processing paraemters (these can be changed to observe thier effects
%  on performance, see help on the delayesttm function
winlen = 16*dl;  %  Set window length in seconds (must be bigger than delay)
%  Create signal length long enough to include signal
len = 2*(4/(f12p-f11p) + 4/f11p);
td = winlen;
simsiglen = td+len+winlen;
%  Signal processing parameters
%  Center clipping is useful for removing contributions of low level noise
%   especially important of a long segment is used and the signal of
%   interest has relatively short support
pro.cliptype = 'cnsc';  %  Select clipping type
pro.cliplevel = .3;     %  Select clipping level (percentage of maximum between 0 and 1)
%  Detrending is important for removing DC and low frequency artifact that
%  may dominate the correlation peak
pro.detrendord = 1; %  Detrending by removing best line fit
%  The peak rarely fall on a sample point.  Additional accuracy or grid resolution can be
%  obtained by interpolating.  Sinc function interpolation is used
pro.interpord = 5; %  Number of actual data points to each side of the interpolated point
pro.interpfac = 10; %  Factor to increase samping rate (i.e. 10 interpolate samples between
                    %  every original sample point and peak taken from finer grid

%  Create pulse samples
[target1] = simimp(f11p,f12p,f11s,f12s,td,fs,simsiglen);
target1 = target1 / max(abs(target1));  %  Normalize by maximum amplitude
target2 = delayf(target1,fs,dl,simsiglen);  %  Create delayed version
for k=1:runs
    s1 = target1 + sclnos*randn(size(target1));  %  Add Gaussian noise to 1st pulse
    s2 = target2 + sclnos*randn(size(target2));  %  Add Gaussian noise to 2nd pulse
    %  Randomly place pulse in processin window
    stpoint = max([td - (rand(1))*winlen/2,0]);  %  Starting point of processing window  
    edpoint = min([winlen+stpoint,simsiglen]);     % Ending point of processing window
    stind = round(stpoint*fs) + 1;  %  Index of starting point
    edind = round(edpoint*fs) +1;   %  Index of ending point
    %  Apply delay estimator
    [d(k), cv(k), pro] = delayesttm(s2(stind:edind),s1(stind:edind),fs,pro);
    %  Comment out next 6 lines to avoid the plotting and stoping in the
    %  loop
    figure(1); plot(1000*[stind:edind]/fs, s2(stind:edind),'r');
    hold on; plot(1000*[stind:edind]/fs, s1(stind:edind),'g'); hold off
    title(['delay estimate = ' num2str(1000*d(k)) 'ms  and % error = ' num2str(100*(d(k) - dl)/dl) ])
    xlabel('milliseconds')
    disp('Hit any key to continue')
    pause
end
%  Compute percent error
pe = 100*(d - dl)/dl;
pe = sort(pe);
%  Drop outliers from error computation (i.e. 10% of the largest positive and
%  negative errors)
mpe = mean(pe(round(runs*.1)+1:round(runs*.9)+1));  %  Trimmed mean error in percent
spe = std(pe(round(runs*.1)+1:round(runs*.9)+1));   %  Trimmed standard deviation in percent
disp(['All done! The average error (droping 10% outliers) is: ' num2str(mpe) '%'])
disp(['The corresponding standard deviation is: ' num2str(spe)])

