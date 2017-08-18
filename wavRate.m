function [ awavOut ] = wavRate( fwav1, fwav2, twinlen, tend, fs )
% Input two .wav files. Play segments of length TWIN from each, where the
% segments are concatenated in the order they occur in the recording. i.e. a
% segment from the first file will play, then the audio stream will
% continue with a segment from the second file.
% INPUTS: 
%       fwav1-> A character array representing the first .wav file
%       fwav2-> A character array representing the second .wav file
%       twinlen-> The desired window length in seconds
%       tend-> The end timestamp on the file
%       fs-> The sampling frequency of the files (should be the same for
%       both)
% OUTPUT:
%       awavOut-> An output signal array of the concatenated audio stream. 
%
%   Written by Kyle Ihli  October 2014

swinlen= fs*twinlen; %window length in samples
send= fs*tend;
bseg=1;
eseg=swinlen;
numwin= floor(tend/twinlen/2); %calculate number of full windows
ratetab= zeros(numwin, 2);
awavOut= [];

for i=1:numwin
    [wav1, fs1] = audioread(fwav1, [bseg eseg]);
    [wav2, fs2] = audioread(fwav2, [eseg+1 eseg+swinlen]);
    if fs1~=fs2
       disp('ERROR: both input .wav files should have same sampling frequency...');
       break
    end
    fs=fs1;
    disp('File 1 playing. . .');
    soundsc(wav1(:,1), fs);
    %pause(1+twinlen); % 1-second break between clips
    pause(twinlen);
    disp('File 2 playing. . .');
    soundsc(wav2(:,1), fs);
    pause(twinlen);
%     bIn= false; % 'good' input flag
%     while bIn==false
%         fselect= input('Which file has the lesser quality (enter 1 or 2)? ','s');
%         if (fselect(1)==int2str(1)) || (fselect(1)==int2str(2))
%             bIn=true;
%             disp(['File ' fselect ' selected as having lower quality.']);
%         else
%             disp('ERROR: invalid input. Input must be 1 or 2. . .');
%         end
%     end
%     chquality= input('Rate the quality of the file on a scale from 1 to 10, with 10\n being equal in quality to the other file: ','s');
    
%     ratetab(i,:)= [str2double(fselect) str2double(chquality)];
    
    bseg= eseg+swinlen+1;
    eseg= eseg+1+2*swinlen;
    awavOut= [awavOut; wav1(:,1); wav2(:,1)];
end
audiowrite(['wavRate results-' fwav1(1:end-4) '_' fwav2], awavOut,fs);

end

