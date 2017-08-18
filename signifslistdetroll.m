function strms1 = signifslistdetroll(filename, nwlen, fovp, fsb,fsresample)
% This function combines the rollingdetect script with signifslist to roll
% through an audio file and only process the segments where rollingdetect
% detects a source. This function is almost the same as rollingdetect,
% except that for each detection, signifslistdetect will be called to
% process the segment. This script also returns a matrix of detection data.
% NOTE: nwlen should be passed through in seconds. It will be converted to
% samples later in the script.

strms1=[];
g = audioinfo(filename); %just read to get the sampling frequency
fs = g.SamplingRate;
nm = g.NumChannels; %get the number of mics
pfa = 1e-3;  %  Probability of false positive
tscl = -log(pfa);  %  Scale to convert energy estimate to threshold
bseg = 1;  %  Start segment index
nwlen= nwlen*fs;
eseg = bseg+nwlen-1;  % End segment index
ninc = (nwlen)/2; %by default each segment is half a window length apart
ff = .75;  %  forgetting factor for noise power estimate (>0 and <1)
b = 1;  %  Set high so it will start off not detecting
det=[];
detreject=[];
g = audioinfo(filename, 'size');
siz = [g.TotalSamples, g.NumChannels];
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
% a = whiten(a.*winrec, beta, [45e3 60e3]/(fs/2));
 %  Get energy in segment
 gg = (std(a)).^2;
 % Compare to treshold based on noise power history
 % also keep b*tscl in a 'good' range (around 0.5e-4)
 % 'bad' range is around 0.5e-5
%  if any( tscl*b< 0.1e-4)
%     b= 0.3e-4/tscl;
%  end
 if mean(gg > tscl*b) >=0.5 %if at least 50% channels have detection
     cnt=cnt+1;  % increment detection counter
     strms1= [strms1; signifslistdet(bseg/fs, filename, nwlen, fovp, fsb, fsresample)];
     %det(cnt) = bseg/fs;  %  If detected, save position don't update noise power
     %provide visual (optional)
     %plot(a);
     %pause(.1);
 else
    b = (ff*b + gg)*(1-ff);  %  Updated noise energy
    cntreject= cntreject+1;
    detreject(cntreject)= bseg/fs;
 end
 %  Go to next segment
 bseg = ninc + bseg;
 eseg = bseg+nwlen-1;
end