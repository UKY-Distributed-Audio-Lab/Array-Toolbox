function [det, detreject] = rollingdetect(filename)
% This script reads through a script to detect when a source emits a sound.
% Audio from FILENAME is read in segments, the size of which can be changed in the script
% through NWLEN. The 'forgetting factor' FF is a way to tweak the influence
% of previous detections/non-detections on the processing of the current
% sample. NOTE: this script will occasionaly throw out true sound source
% segments and still allow some non-detection segments to pass through. DET
% is a matrix which contains the timestamps of all the segments that the
% script marked as detections. DETREJECT contains all the other timestamps,
% and is useful for determining the effectiveness of the script for getting
% rid of non-detection audio segments.

[y,fs]= audioread(filename, [1 5]); %just read to get the sampling frequency
nm= size(y,2); %get the number of mics
nwlen=fs*10e-3;
pfa = 1e-3;  %  Probability of false positive
tscl = -log(pfa);  %  Scale to convert energy estimate to threshold
bseg = 1;  %  Start segment index
eseg = bseg+nwlen-1;  % End segment index
ninc = (nwlen)/2; %by default each segment is half a window length apart
ff = .75;  %  forgetting factor for noise power estimate (>0 and <1)
b = 1;  %  Set high so it will start off not detecting
det=[];
detreject=[];
siz= wavread(filename, 'size');

% setup the parameters for filtering and whitening
[bn,an] = butter(4,[40e3,130e3]/(fs/2)); %set butterworth filter range in brackets
tapwin = flattap(nwlen,20);
winrec = tapwin*ones(1,nm);
beta= .45; % beta factor for use in whitening

% initialize the current number of detections/rejections
cnt = 0;
cntreject= 0;
while eseg<= siz(1,1)
 % Read in chunk of data and filter and whiten
 [a,fs]= audioread(filename, [bseg,eseg]);
 a = filtfilt(bn,an,a);
 a = whiten(a.*winrec, beta, [45e3 60e3]/(fs/2));
 %  Get energy in segment
 gg = (std(a)).^2;
 % Compare to treshold based on noise power history
 if mean(gg > tscl*b) >=0.5 %if at least 50% channels have detection
     cnt=cnt+1;  % increment detection counter
     det(cnt) = bseg/fs;  %  If detected, save position don't update noise power
     %provide visual (optional)
     plot(a);
     pause(.1);
 else
    b = (ff*b + gg)*(1-ff);  %  Updated noise energy
    cntreject= cntreject+1;
    detreject(cntreject)= bseg/fs;
 end
 %  Go to next segment
 bseg = ninc + bseg;
 eseg = bseg+nwlen-1;
end
% flip the matrices if format required by other scripts
det=det';
detreject=detreject';
 
    
