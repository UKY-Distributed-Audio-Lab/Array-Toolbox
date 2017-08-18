%  This script tests the RT60 program using a segment of real data
%  recorded in a the Audio Cage in May 2006.  The 12 by 12 by 8.5 foot
%  Cage was plexiglass wall with 1 in acoustic foam.
%  
%  Data was collect by Adam Leedy and program written by Kevin D. Donohue
%  (donohue@engr.uky.edu) May 29 2006


hpc = 100;  %  Set highpass filter limit
%  Load data file
fnam = 'rt60exdata.wav'
[sigseg,fs] = audioread(fnam);
tp = [0:length(sigseg)-1]/fs;  % compute time axis
figure; plot(tp, sigseg)  %  Plot recorded signal
title(['Recording of Noise Burst and Return to Quiet'])
xlabel('Seconds')
[rt60, sigpow, nospow] = rt60est(sigseg, fs, hpc);  %  compute RT 60 time
%  and plot envelope used to make estimates

