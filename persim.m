%  This script will use the SIGSIM function and the TESTPD function to 
%  Examine the performance of variations on 2 DS beamformers used in the
%  SRP algorithms.

%  Set the simulated positions of the sound sources (columns are x,y,z
%  coordinates in meters)
sigpos = [-.6 .9 0 .2; .7, .4 .1 .9];
%  Set the coordinates of the mic postions (colums are x,y,z in meters for each mic);
mpos = [-.8 -.4 0 .4 .8; 0 0 0 0 0];
%  Set limits for the field of view in meters (x,y,z coordinate for upper
%  left and lower right)
flims = [-1 1; 1 0];
%  Speed of sound in medium
c = 345;
%  Positions of scatterers generating multipath echos
mpaths = [.5, .5; 1.2 -.5; -2.2 -.1];
%  Frequency limits of simulated sound
f1 = 200;
f2=4000;
%  Window length for simulated signal
winlen = 100/(f2-f1);   %  Make it 100 times the inverse bandwidth
fs = 44.1e3;            %  Sampling rate
pcttap = 20;            %  taper edges of spectrum of simulated sound by this % of the bandwidth
snr = 20;               %  SNR (peak signal the RMS noise) in dB
%[sigout, tax] = simarraysig(sigit, fs, sigpos, mpos, flims, c);
[sigout, tax, td] = sigsim(f1,f2,winlen,fs,pcttap, sigpos, mpos, flims, c);
nosgain = 10^(-30/20);
sigout = sigout + nosgain*max(max(abs(sigout)))*randn(size(sigout));
