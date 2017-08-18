fils = {'twospeaker1_sent1_mic1.wav', ...  
 'twospeaker1_sent1_mic2.wav', ...   
'twospeaker1_sent1_mic3.wav', ...   
'twospeaker1_sent1_mic4.wav', ...   
'twospeaker1_sent1_mic5.wav', ...   
'twospeaker1_sent1_mic6.wav', ...   
'twospeaker1_sent1_mic7.wav', ...   
'twospeaker1_sent1_mic8.wav'};
[y,fs, nbits] = wavread(fils{1}, [1 10]);

start_second = 7.25;
end_second = start_second + 3; % .0232*10;
lfco = 100;  %  Low-pass cut-off frequency 
hfco = 7000;  % High-pass cut-off
[b,a] = butter(4,2*[lfco, hfco]/fs);
ba = [];
tint = [fix(start_second*fs)+1:fix(end_second*fs)];
clear y
for k=1:length(fils)
    [y,fs, nbits] = wavread(fils{k}, [fix(start_second*fs)+1, fix(end_second*fs)]);
    ba(:,k) = filtfilt(b,a,y);
%ba(:,k) = y;
end
%    [y,fs, nbits] = wavread(fils{1});
clear y
