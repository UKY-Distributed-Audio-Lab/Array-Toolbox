function [scriptout len] = bfsegment(fn,pn,segn,imseg,fovp)
%
%  This function beamforms on detected streams in file FN located at path PN
%  identified by their segment number in vector SEGN.  A delay and sum
%  beamformer is used with inverse distance weighting.  Each beamformed
%  stream will be saved as a wave file and a script create to identify
%  the spatial position and starting time associated with the file.
%
%     [scriptout, len] = bfsegment(fn,pn,segn,imseg,fovp)
%
% Inputs:
%       FN    => String for filename of multichannel wave file recording
%       PN    => String for path name of file (use empty vector if current
%               directory [])
%       SEGN  => Vector of indices denoting segments to extract
%       IMSEG => Text file with script information identifying segment
%               where each row is of the form:
%               [segment #, detection statistic, x, y, z location, frame#,
%                  detection/rejection description]
%       FOVP  => Data structure with recording and processing information:
%     Critical fields of FOVP:
%        mp     =>	Microphone positions, each column corresponds to a microphone [x,y,z]' in meters
%        c      => 	Speed of sound
%        fs     => 	Sampling frequency 
%        trez   =>	Time window length (resolution) in seconds for detecting each block (increment will be 1/2
%                   of time window length) 
%
% Outputs:
%       scriptout   =>  Text file identifying position, start time, and file
%                        name of valid beamformed streams for the form:
%                        [segment#, xs,ys,zs, xe, ye, ze, time start, filename]
%                        Filename will have the same base name of the 
%                        original file with the 2 digit segment number
%                        appended, i.e. base01.wav, base02.wav, ...
%                        Any streams with valid segments not found will not 
%                        be saved and passed on to the masking component.
%       len         =>  Number of valid streams detected
% 
% Required Functions: 
%               dsb.m, arweights.m, intel.m, and sii.m  
%
% Written by Kevin D. Donohue July 23, 2012  (donohue@engr.uky.edu)
% Modified by Kirstin Brangers July 30, 2012
%           Deleted any nonvalid streams before passing on to masking component
% 



% Name of multichannel wavefile with source recordings
sfile = [pn, fn]; 
% Processing parameters
% wp => mic weight distibtion where 0 results in equal weights,
% and 1 results in an inverse distance weighting where they are
% scaled such that the closest mic gets a weight of 1. A positive
% number gives more weight to closer mics, and a negative number
% gives more weight to distant mics.
wp = .8;

% Processing window size in seconds for reading in data and applying 
% the current source position number associated with the corresponding
% time the window includes.
ws = fovp.trez;
if isempty(fovp.tinc)
    tinc= ws/2;
else
    tinc = fovp.tinc;
end

% High-pass cutoff in Hertz for filtering signal before beamforming
fc = 100;

% Load files associated with recordings
% load('ASA1nfOBJ100fps.mat');      %  Source Positions
% stps = find(imsegsw(:,7) > 0);    %  Find points not rejected by secondary threshold
% imt = imsegsw(stps,:);


% Read in dimensions of file
[sizz, fs]= wavread(sfile,'size');  %  Speaker recording file size and sample rate
siglen = sizz(1);
chans = sizz(2);
% Compute High-pass filter coefficients
[b,a] = butter(4,fc/(fs/2),'high');


% Loop to extract each requested segment
numofstrm = length(segn);
scriptout = cell(numofstrm,8) ;     % Initialize output cell array
hwait = waitbar(0,['Processing potential streams ' int2str(0) ' of ' int2str(numofstrm)]);
scnt = 0;  % Initalize valid script count
for k=1:numofstrm
  waitbar(k/numofstrm,hwait,['Processing potential streams ' int2str(k) ' of ' int2str(numofstrm)]);
   
    % Find points of requested stream and not rejected by secondary threshold
    stps = find(imseg(:,1) == segn(k) & imseg(:,8) > 0);
    if isempty(stps)
        disp(['No valid segments found for stream number ', num2str(segn(k))])
        
        %scriptout{k,1} = 0;
        %scriptout{k,2} = 0;
        %scriptout{k,3} = 0;
        %scriptout{k,4} = 0;
        %scriptout{k,5} = 0;
        %scriptout{k,6} = 0;
        %scriptout{k,7} = 0;
        %scriptout{k,8} = fout;
        %yo = zeros(siglen,1); 

    else    % If points found in requested stream, beamform on it to extract
        scnt = scnt +1;  %  Increment script count
        
         % Output file name
        if scnt<10                         % If single digit, append a pre zero
            fout =  [fn(1:end-4), '0' int2str(scnt) '.wav'];
        else
            fout = [fn(1:end-4), int2str(scnt) '.wav'];
        end
        
        imt = imseg(stps,:);  % Group all segments associated with same source together
        [rp,cp] = size(imt);
        scriptout{scnt,1} = imt(1,1);
        scriptout{scnt,2} = imt(1,4);
        scriptout{scnt,3} = imt(1,5);
        scriptout{scnt,4} = imt(end,6);
        scriptout{scnt,5} = imt(end,4);
        scriptout{scnt,6} = imt(end,5);
        scriptout{scnt,7} = imt(end,6);
        scriptout{scnt,8} = imt(7,6)*tinc;
        scriptout{scnt,9} = fout;
        medfilord = 11;  %  Order of median filter to limit anomolous jumps in locatation estimates
        if rp > medfilord
            imt(:,4) = medfilt1(imt(:,4),medfilord);
            imt(:,5) = medfilt1(imt(:,5),medfilord);
            imt(:,6) = medfilt1(imt(:,6),medfilord);
        end
        
        %  Initialize sample indicies to step through recordings one
        %  window at a time
        yo = zeros(siglen,1);       % Initialize ouput arrays for beamforming
        seglen = round(ws*fs);      % Signal time index for end of first window
        wsc = ws;                   % Initialize current time window end point in seconds 
        %tc = tinc*(imt(1,7)-1);     % Initialize window beginning in seconds
        tc = 0;
        segstart = 1; %round(tc*fs)+1;
    %    segstart = max(segstart,1);% Signal time index for beginning of first window
        segend = segstart+seglen-1;
        
        % Sine square window for output
        swin = sin(pi*[0:seglen-1]'/(seglen-1)).^2;

        %  Beamform the data in loop until current segment end sample
        %  is less than total signal length
        segend = min(segend,siglen);
        
        x = wavread(sfile,[segstart segend]);   % Speaker of Interest (SOI)
        [x,xf1] = filter(b,a,x);
        xPrev = zeros(size(x));
        
        %  Find spatial location for current time
        [mv, mp] = min(abs(tc - (imt(:,6)-1)*tinc));
        sloc = imt(mp(1),4:6);
        d = (sloc'*ones(1,chans)-fovp.mp);
        ar = arweights(sqrt(sum(d.^2)));
        
        %  Normalize weight value of channel closest to SOI
        [nscl, lmax] = max(ar);
        
        % Raise weights to power to emphasize/deemphasize distant mics
        arw = (ar/nscl(1)).^wp; 
        
        % Expand out weight in a matrix to multiply with signal
        ww = ones(length(x),1)*arw;
        xp = x.*ww;         % Apply weights to multichannel segement
        xo = dsb(xp, xPrev, fs, sloc', fovp.mp, fovp.c); % Combined
        yo(segstart:segend) = xo.*swin; 
        xPrev = xp;
        tc = tc + tinc;     % Increment to next segment 
        segstart = round(tc*fs)+1;  % Signal time index for beginning of next window
        segend = segstart+seglen-1;

        while segend <= siglen
            % Read in short segment of data form
            x = wavread(sfile,[segstart segend]);   % Speaker of Interest (SOI)
            [x,xf1] = filter(b,a,x,xf1);
            
            % Find spatial location for current time
            [mv, mp] = min(abs(tc - imt(:,6)*tinc));
            sloc = imt(mp(1),4:6);
            d = (sloc'*ones(1,chans)-fovp.mp);
            ar = arweights(sqrt(sum(d.^2)));
            
            % Normalize weight value of channel closest to SOI
            [nscl, lmax] = max(ar);
            
            % Raise weights to power to emphasize/deemphasize distant mics
            arw = (ar/nscl(1)).^wp; 
            
            % Expand out weight in a matrix to multiply with signal
            ww = ones(length(x),1)*arw;
            xp = x.*ww;         % Apply weights to multichannel segement
            xo = dsb(xp, xPrev, fs, sloc', fovp.mp, fovp.c); % Combined
            yo(segstart:segend) = xo.*swin + yo(segstart:segend);  %  Accumulate with overlap 
            xPrev = xp;
            tc = tc + tinc;     % Increment to next segment 
            segstart = round(tc*fs)+1;  % Signal time index for beginning of next window
            segend = segstart+seglen-1;

        end
    end
  
    wavwrite(yo/4,fs,[pn,fout])
    
end
close(hwait)

% Check streams: If columns 1-7 are zeros, then delete streams
valid = cell(1,9) ;     % Initialize output cell array
temp= zeros(1,8);       % Temp vector to hold values of scriptout
count = 1;              % Initialize stream number 
for i=1:scnt
    % Convert cell array to numeric array
    for n=1:8
        temp(1,n) = cell2mat(scriptout(i,n));
    end
    % If columns 1-7 are not all zeros, copy over to output cell array VALID
    if sum(temp) ~= 0
        valid{count,1} = temp(1,1);
        valid{count,2} = temp(1,2);
        valid{count,3} = temp(1,3);
        valid{count,4} = temp(1,4);
        valid{count,5} = temp(1,5);
        valid{count,6} = temp(1,6);
        valid{count,7} = temp(1,7);
        valid{count,8} = temp(1,8);
        valid{count,9} = scriptout{i,9};
        
        count = count+1;        % Increment stream number
    end
end

scriptout = valid;
[len cols] = size(scriptout);
