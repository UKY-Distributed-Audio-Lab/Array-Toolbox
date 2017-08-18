fils = {'s2150cm.wav'}
[y,fs, nbits] = wavread(fils{1}, [1 10]);

start_second = 40;
end_second = start_second + 50; % .0232*10;
lfco = 100;  %  Low-pass cut-off frequency 
hfco = 15000;  % High-pass cut-off
[b,a] = butter(4,2*[lfco, hfco]/fs);

tint = [fix(start_second*fs)+1:fix(end_second*fs)];
clear y

%    [y,fs, nbits] = wavread(fils{1}, [fix(start_second*fs)+1, fix(end_second*fs)]);

    [y,fs, nbits] = wavread(fils{1});

[smp, chns] = size(y); 
for k = 1:chns
    ba(:,k) = filtfilt(b,a,y(:,k));
end
clear y
