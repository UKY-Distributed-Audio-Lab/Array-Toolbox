%  This script demonstrates the delay estimators using a time domain approach 
%  (function DELAYESTTM), and a frequency domain approach (function
%  DELAYSETFR).  The signal and simulation parameters can be set by editing the
%  initial lines of this script (see comments for details).  The signal
%  processing options for the delay estimator can also be set (see help on
%  the DELAYESTTM and DELAYESTFR functions) .  Signal simulation values include
%  the delay, SNR, sampling frequency, and bandlimits for the test signals. It performs
%  a Monte Carlo simulation of the signal in noise, where the noise pattern and
%  signal relationship to the window is varied in each run.  The figure for the
%  signal in noise are presented in each pass of the loop along with the estimation
%  errors for both approaches presented in the titles of the figure.
%  The program pauses with each plot.  These line can be commented out for
%  larger Monte Carlo runs
%  
%    Script written by Kevin D. Donohue  (May 30, 2006)
%

runs = 20;  %  Number of simulations estimates.  Not too high if plotting in each loop pass
%  Signal parameters
fs = 44.1e3;   %  Sampling frequency
dl = .0025;  %  Simulated Delay in seconds
%  Generate target signal pulse bandwidths
f12p = 10000;  %  Corresponding upper frequency limit
f11p = 500;  %  Lower frequency limit
xsp = 20; %  Percent for transition band interval
%  White noise snr
wgnsnr = 30;  %  SNR in dB
%  Ratio of signal window length to delay (should be at least 2 times the delay length
% to include the duration of the signal and the delay between them)
winlenrat = 6;

%  Processing parameters (these can be changed to observe thier effects
%  on performance, see help on the delayesttm  and delayestfr functions.

%  Signal processing parameters for time domain approach
%  Center clipping is good for removing low-level noise distribute over the
%  windowed segment
prot.cliptype = 'none';  %  Select clipping type
prot.cliplevel = 0;    %  Select clipping level (percentage of maximum between 0 and 1)
%  Detrending is important for removing DC and low frequency artifacts that
%  may dominate a correlation peak
prot.detrendord = 1;  %  Select type of detrending, this detrending removes the best line fit
%  Time domain method interpolates to locate a peak on a higher resolution 
%  grid spacing allowed for by the original sampling rate
prot.interpord = 4;  %  Half the Sinc function interpolation order
prot.interpfac = 10; %  Number of interpolated points between samples (factor to decrease grid spacing)

% Signal Processing parameters for frequency domain approach
%  It is best to limit the frquency range used by the estimator to those
%  where the signal energy is known to be high
prof.flow = 400;  %  Low frequency limit to exclude lower frequency from estimate
prof.fhigh = 3800;    %  High frequency limit to exclude higher frequencies from estimate 
%  Detrending is important for removing DC and low frequency artifacts that
%  may dominate a correlation peak
prof.detrendord = 1;  %  Select type of detrending, this detrending removes the best line fit  

% Compute from above parameters 
sclnos = 10^(-wgnsnr/20);   %  Scaling ratio for noise and normalized signal
f11s = f11p*(1-xsp/100);  %  Lower stopband for simulated signal
f12s = f12p*(1+xsp/100);  %  Upper stopband for simulated signal
winlen = winlenrat*dl;  %  Window length in seconds (must be bigger than delay)

simsiglen = 2*winlen + 2*(4/(f12p-f11p) + 4/f11p);  %  length of simulated signal to ensure full signal fit in window 

%  Create pulse samples
[target1] = simimp(f11p,f12p,f11s,f12s,winlen,fs,simsiglen);
target1 = target1 / max(abs(target1));  %  Normalize by maximum amplitude
target2 = delayf(target1,fs,dl,simsiglen);  %  Create delayed version
for k=1:runs
    s1 = target1 + sclnos*randn(size(target1));  %  Add Gaussian noise to 1st pulse
    s2 = target2 + sclnos*randn(size(target2));  %  Add Gaussian noise to 2nd pulse
    %  Randomly place pulse in processing window
    stpoint = max([winlen - (rand(1))*winlen/2,0]);  %  Starting point of processing window  
    edpoint = min([winlen+stpoint,simsiglen]);     % Ending point of processing window
    stind = round(stpoint*fs) + 1;  %  Index of starting point
    edind = round(edpoint*fs) +1;   %  Index of ending point
    
    %  Apply delay estimators
    [df(k), stds(k), prof] = delayestfr(s2(stind:edind),s1(stind:edind),fs,prof);
    [dt(k), cv(k), prot] = delayesttm(s2(stind:edind),s1(stind:edind),fs, prot);
    %  Comment out next 7 lines to avoid the plotting and stopping in the
    %  loop (comment up to and including the pause statement)
    figure(1); plot(1000*[stind:edind]/fs, s2(stind:edind),'r');
    hold on; plot(1000*[stind:edind]/fs, s1(stind:edind),'g'); hold off
    title1 = ['delay estimate (frequency domain) = ' num2str(1000*df(k)) 'ms  and % error = ' num2str(100*(df(k) - dl)/dl) ];
    title2 = ['delay estimate (time domain) = ' num2str(1000*dt(k)) 'ms  and % error = ' num2str(100*(dt(k) - dl)/dl) ];
    title({title1;title2}) % cell array to make title 2 lines
    disp('Hit any key to continue')
    pause
end
%  Compute percent error for frequency domain approach
pef = 100*(df - dl)/dl;
pef = sort(pef);
%  Drop outliers from error computation (i.e. 10% of the largest positive and
%  negative errors)
mpef = mean(pef(floor(runs*.05)+1:floor(runs*.95)+1));  %  Trimmed mean error in percent
spef = std(pef(floor(runs*.05)+1:floor(runs*.95)+1));   %  Trimmed standard deviation in percent
disp(['The average error for frequency domain (droping 10% outliers) is: ' num2str(mpef) '%'])
disp(['The corresponding standard deviation is: ' num2str(spef)])
%  Compute percent error for time domain approach
pet = 100*(dt - dl)/dl;
pet = sort(pet);
%  Drop outliers from error computation (i.e. 10% of the largest positive and
%  negative errors)
mpet = mean(pet(floor(runs*.05)+1:floor(runs*.95)+1));  %  Trimmed mean error in percent
spet = std(pet(floor(runs*.05)+1:floor(runs*.95)+1));   %  Trimmed standard deviation in percent
disp(['All done! The average error for time domain (droping 10% outliers) is: ' num2str(mpet) '%'])
disp(['The corresponding standard deviation is: ' num2str(spet)])
