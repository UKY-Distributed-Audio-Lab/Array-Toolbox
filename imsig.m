
f1 = 100;
f2 = 8000;
winlen = 100e-3;
fs = 44100;
pcttap = 40;
sigpos = [.5; .5];
mpos = [0; 0];
flims = [-1 1; 0 1];
c = 345;

[sigit, tax, td] = sigsim(f1,f2,winlen,fs,pcttap, sigpos, mpos, flims, c);