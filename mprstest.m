%  This script is to test mlefreampairs

%Readin data to fool with
fils = {'twospeaker1_sent2_mic1.wav' ...
        'twospeaker1_sent2_mic2.wav' ...
        'twospeaker1_sent2_mic3.wav' ...
        'twospeaker1_sent2_mic4.wav' ...
        'twospeaker1_sent2_mic5.wav' ...
        'twospeaker1_sent2_mic6.wav' ...
        'twospeaker1_sent2_mic7.wav' ...
        'twospeaker1_sent2_mic8.wav'}
          
[y,fs, nbits] = audioread(fils{1}, [1 10]);

start_second = 6;
end_second = start_second + 0.5; % .0232*10;
lfco = 500;  %  Low cut-off frequency 
hfco = 7000;  % High cut-off
[b,a] = butter(4,2*[lfco, hfco]/fs);

tint = [fix(start_second*fs)+1:fix(end_second*fs)];
clear y

for k = 1:length(fils)
    [y,fs, nbits] = audioread(fils{k}, [fix(start_second*fs)+1, fix(end_second*fs)]);
    ba(:,k) = filtfilt(b,a,y);
end
clear y
plot(ba(:,1))
pause
%  set parameters
c = 345.8824;  %  Speed of Sound

% Mic postions [x; y] column for each mic 
mpos = [-140 -100 -60 -20 20 60 100 140; 0 0 0 0 0 0 0 0]/100;  
%  Limits of rectangular field of view (FOV) - diagonal corners
flims = [-110 110; 20 200]/100;
rez = .06;
%flims = [-35 150; 50 200]/100;    
win = 20e-3; %  Correlation winddow
winlen = round(win*fs)+1;
corwin = hamming(winlen);
grd = gridspace(flims(:,1), flims(:,2), rez);

prs = mposanaly(mpos,2);
[g, msrt] = sort(prs(:,3));
out = prs(msrt,1:3);
hh = find(out(:,3) < 10);
prs = out(hh,:);
nrm = max(prs(:,3));
prs(:,3) = prs(:,3)/(nrm+eps);

[rws,clms] = size(ba);
blen = fix((win+.01)*fs);
skip = fix(winlen/2);
imc = 0;
for k=1:skip:rws-blen
    imc=imc+1
    sigar = whiten(ba(k:k+blen,:),.9);
    mleim = srpframepairs(sigar, grd, mpos, prs, fs, c, corwin);
    da(:,:,imc) = localdnoise(mleim,6,3);
    %figure(1); imagesc(grd{1},grd{2}, localdnoise(mleim,6,1), [0 .3]); colorbar, axis('xy');
    %figure(2); plot(sigar(:,1))
    %pause(.01)
end