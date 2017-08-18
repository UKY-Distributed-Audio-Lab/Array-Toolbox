% This script demonstrates the usage of the functions to estimate the
% intelligibility of speech-in-noise conditions. The script reads two
% wavefiles; the first one (man's voice) is the target signal and the
% second one (woman's voice) is the interfering speaker. To contrast the 
% effects of an interference speaker (non-stationary noise) with a white
% noise inference (stationary noise/constant power) a white noise signal
% is also used. These signals are then scaled (varying SNR) to illustrate 
% Good intelligibility (~.6), barely intelligible (~.2), and unintelligible
% (.1).
% 
% The output of the script plots the Speech intelligibility Index for 100 ms
% windows over the speech signal, and displays the computed mean and standard
% deviation of SII (for active speech, silence intervals were excluded).
% Users can change the envelope threshold which is used to remove the
% silence intervals so the mean SII is not as dependent on pauses between
% words (if speech is not active intelligibility is low as expected). They
% can also change the length of the overlaping time window which is used
% to window the signals while estimating the SII.  The scaled target
% signal together with the interfering signal(noise) is played for reference
% purpose.
%
% Functions required from Array Toolbox
%  1. intel.m
%  2. spectrumlevel.m
%  3. sii.m
%  4. rmsilence.m
% Data Files required -> woman1.wav, man1.wav
%
% Written by Arulkumaran Muthukumarasamy (arulkumaran@uky.edu) July 2008

hpc = 100;  %  Highpass cuttoff
wts=[2.0, 0.14, .06]; % weights used to scale the signal for different levels of intelligibility
twin=100e-3; % overlapping time window used to move along the total signal while estimating SII
env_thresh=0.15; % threshold used to remove silence intervals in the signals
% Single speaker wavfiles which are assumed to be the speech signal and masking
% noise respectively
wavfiles={'man1.wav', 'woman1.wav'};
% Readin wavefiles and truncates it to be of equal length
[sigin, fs] = wav2sig(wavfiles,[0 10]);
[b,a] = butter(5,hpc/(fs/2),'high'); % designing 5th order high pass butterworth filter
% Filter the input signals through a high pass filter to remove the room
% noise
sigfil(:,1) = filter(b,a,sigin(:,1)); % filtering the speech signals using the designed high pass filter
sigfil(:,2) = filter(b,a,sigin(:,2)); % filtering the noise signal using the designed high pass filter
sigfil(:,3) = randn(size(sigfil(:,2)))*std(sigfil(:,2))/4;  % white noise with same RMS as interfering speech 
% Remove the silence intervals using the function 'rmsilence.m'
% The silence intervals in interfering signal does not mask the target signal and hence intelligibility
% would be higher at those points. Removing these silence intervals reduces
% the deviation of intelligibility over time
[sig_rm,less]=rmsilence(sigfil,fs,env_thresh);

for rr=1:length(wts)
    ysig = sig_rm(:,1)*wts(rr); % scale the signal using the weights specified
    ynos1 = sig_rm(:,2);  % Interferring speech
    ynos2 = sig_rm(:,3);  % Interferring white noise
    %The Speech Intelligibility Index is computed using the function "intel.m"
    [sii_val,tax]=intel(ysig,ynos2,fs,twin);
    sii_mean=mean(sii_val);% Estimating the mean of SII over time
    sii_std=std(sii_val); % Estimating the standard deviation of SII
    mn=sprintf('The mean Speech Intelligibility Index is %f',sii_mean);
    sd=sprintf('The Standard Deviation of the Speech Intelligibility Index is %f',sii_std);
    % plots the estimated SII over the time axis 
    figure;
    
    plot(tax,sii_val);
    set(gca,'Ylim', [0 1])
    xlabel('time(seconds)');
    ylabel('Speech Intelligibility Index(SII)');
    % displays the mean and the standard deviation of SII 
    title(['White noise background with mean SII = ', num2str(sii_mean)])
    disp(mn);
    disp(sd);
    
    soundsc(sigfil(:,1)*wts(rr)+sigfil(:,3),fs);% playback the signal along with the masking noise
    pause;
    [sii_val,tax]=intel(ysig,ynos1,fs,twin);
    sii_mean=mean(sii_val);% Estimating the mean of SII over time
    sii_std=std(sii_val); % Estimating the standard deviation of SII
    mn=sprintf('The mean Speech Intelligibility Index is %f',sii_mean);
    sd=sprintf('The Standard Deviation of the Speech Intelligibility Index is %f',sii_std);
    % plots the estimated SII over the time axis 
    figure;
    plot(tax,sii_val);
    set(gca,'Ylim', [0 1])
    xlabel('time(seconds)');
    ylabel('Speech Intelligibility Index(SII)');
    % displays the mean and the standard deviation of SII 
    disp(mn);
    disp(sd);
    title(['Interfering speech background with mean SII = ', num2str(sii_mean) '(for man''s voice)'])
    soundsc(sigfil(:,1)*wts(rr)+sigfil(:,2),fs);% playback the signal along with the masking noise
    pause;
end
