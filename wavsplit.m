function wavsplit(wavName, N, seglen, wininc)
%
%   This function opens a wavefile and splits it into N partitioned files  
%   with a partial overlap between contiguous partitions for use as an input
%   to another function that would typically apply a sliding window over the 
%   original file(such as a convolution operation). The last partition 
%   may trim the orignal file so that computations end with a full window.
%   
%      wavsplit(wavName, N, seglen, wininc)
%
%   Inputs:
%       1) wavName - string identifying wave file for partitioning
%       2) N - number of files to split the original into (positive integer)
%       3) seglen - Processing Window length for sliding window in terms of
%                   samples (positive integer)
%       4) wininc - Number of samples to increment sliding window(positive integer)
%
%   Outputs:
%       No Return values - creates N wavfiles into a dynamically created
%           directory "out_[wavName]", with each file named
%           part(integer).wav
%       The beginning and the ending indices of where the wav file is split
%           can be outputted to the screen by setting outputEn to 1 (line 86)
%
%   written by Satoru Tagawa (staga2@uky.edu) 5/21/08

% Read file into y
[y,fs] = audioread(wavName);

% Make sure the value of wininc makes sense
if (wininc < 1 || wininc >= seglen - 1)
    error('wininc is too large for an overlap to occur');
end

% Find the approximate output wav length
appWavLen = length(y(:,1))/N;
% Make sure the value of N makes sense
if (appWavLen < seglen)
    error('N is too large for the segment length');
end

% Pre-allocate storage for begin_index and end_index to make MATLAB happy
begin_index = zeros(1,N);
end_index = zeros(1,N);

% Compute the beginning and the ending index for splitting the waveform
begin_index(1) = 1;
%   The file is cut at the the segment that ends right before appWavLen
end_index(1) = seglen;
while (end_index(1) <= appWavLen-wininc)
    end_index(1) = end_index(1) + wininc;
end

% Compute how many points need to be added in the beginning of files 2 to N
    numAddPts = seglen - wininc - 1;

% Find the beginning and the ending indices of where each wav file is cut      
for i=2:N
    begin_index(i) = end_index(i-1) - numAddPts;
    if (i<N)
        % Find end_index(i)
        end_index(i) = end_index(i-1);
        while (end_index(i) <= appWavLen*i - wininc)
            end_index(i) = end_index(i) + wininc;
        end
        
        % Error Check
        if (end_index(i) > length(y(:,1)) || end_index(i) <= end_index(i-1))
            error('end_index(i) should never have this value');
        end
    else 
        % Find end_index(N)
        end_index(i) = end_index(i-1);
        while (end_index(i) <= length(y(:,1)) - wininc)
            end_index(i) = end_index(i) + wininc;
        end
        
        % Error check
        if (end_index(i) > length(y(:,1)) || end_index(i) <= end_index(i-1))
            error('end_index(i) should never have this value');
        end
    end
end

% Output beginning and ending indices
outputEn = 0;   % Set this equal to 1 to enable output
if outputEn == 1
    begin_index
    end_index
end

% Write to files
mkdir(['out_',wavName(1:end-9)]);
for j=1:N
    audiowrite(['out_',wavName(1:end-9),'/part', ...
        num2str(j),'.wav'],y(begin_index(j):end_index(j),:), fs);
end