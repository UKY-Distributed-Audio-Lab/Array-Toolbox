%TODO Test
%"I would like a function (can adapt wave2sig.m) that will also ask for a
%vector  of weights to scale the relative RMS values of each of the files.
%So you would normalize each file as it is read in and then multiply it by
%the weight given in the input scaling vector. The output matrix should be
%relatively scaled wave files in each column.  If you then do std() along
%each column (std is effectively the same as RMS) you should get numbers
%proportional to the vector of weights.  This way in the simulation we can
%also adjust the relative strengths of the signals for testing (i.e. make
%the the signal weak/strong enough to find the limits of
%detectability)." -Dr. Donohue
%^this is for reference and not inspiration

%added sig = wav2sig(fnames,weights) functionality but other errors have
%come up. need logic to
%continue.---------------------------------------------------------------------------------------=--

function [sig,fs] = wav2sig(fnames, varargin)
% usage_instr = ['sig = wav2sig(fnames)\nsig = wav2sig(fnames,fs)\n'...
%                  'sig = wav2sig(fnames,weights)\nsig = wav2sig(fnames,tInt)\n'...
%                  'sig = wav2sig(fnames,fs,weights)\nsig = wav2sig(fnames,fs,tInt)\n'...
%                  'sig = wav2sig(fnames,fs,tInt,weights)\n' ];
% This function reads in wave files and stores all the information into a 
% single matrix with equal number of rows, with each column representing
% the different wave files.  
% 
%   sig = wav2sig(fnames)
%   sig = wav2sig(fnames,fs)
%   sig = wav2sig(fnames,tInt)
%   sig = wav2sig(fnames,weights)
%   sig = wav2sig(fnames,fs,tInt)
%   sig = wav2sig(fnames,fs,weights)
%   sig = wav2sig(fnames,fs,tInt,weights)
%
%   Inputs:
%   1) fnames - name of the wave files in a cell array. 
%      Format:
%       a) for 1 file: {'filename1.wav'} 
%       b) multiple files: {'filename1.wav';'filename2.wav'}
%       NOTE: 
%        Each wave file should only have one channel.  If more than one 
%        channels are present, only the first channel will be used.
%   2) varargin (optional):
%       a) fs - resample the wave file to this frequency
%       b) tInt - 1x2 vector to specify time interval (in seconds) to 
%          trim down to
%       c) weights - scaling vector of weights multiplied to each input
%          file after it is normailzed. [ w1; w2;...;wn ]
%
%   Output:
%   sig - matrix with the following properties:
%       a) Number of rows is determined by the number of rows in the
%          longest wave file
%          i) shorter wave files will be padded with 0's
%       b) Number of columns is determined by the number of input wave 
%          files
%
%   Written by Satoru Tagawa (staga2@uky.edu) 6/12/08
%   Edited by Grant Cox 8/22/17


% Parameter check ********************************************************
% The function must have at least 1 parameter, and at most 3 parameters
if nargin == 0
    error('wavToSigMat must have at least 1 parameter');
end
if nargin > 4
    error('There can only be a maximum of 4 parameters');
end

% If more than one parameter, check to see what the parameters are
if nargin > 1
    [numR1,numC1] = size(varargin{1});
%     if numR1 ~= 1  % Argument cannot have more than 1 row
%             error('2nd parameter must have the dimension 1x1 or 1x2');
%     end

    if length(varargin) == 1 %---------------------------------2 parameters
%         if isvector(varargin)
        if [numR1 , numC1] == size(fnames)
            weights = varargin{1};
        elseif numC1 == 1
            fs = varargin{1};
        elseif numC1 == 2
            tInt = varargin{1};
        else
            error('tInt parameter must have the dimension 1x1 or 1x2');
        end
    elseif length(varargin) == 2 %-----------------------------3 parameters
        if numC1 == 1
            fs = varargin{1};
        else
            error('fs must have the dimension 1x1');
        end
        
        [numR2,numC2] = size(varargin{2});
        if [numR2 , numC2] == size(fnames)
            weights = varargin;
        end
        
        %TODO FIX THE REST OF THE INPUTS FROM HERE. USE SIMILAR METHOD TO
        %PREVIOUS LINES
        
        if numR2 ~= 1
            error('tInt must have the dimension 1x1 or 1x2');
        end
        
        if numC2 == 2
            tInt = varargin{2};                 
        else
            error('tInt must have the dimension 1x2');
        end
    else         %---------------------------------------------4 parameters
        if numC1 == 1
            fs = varargin{1};
        else
            error('fs must have the dimension 1x1');
        end
        
        [numR2,numC2] = size(varargin{2});
        
        if numR2 ~= 1
            error('tInt must have the dimension 1x1 or 1x2');
        end
        
        if numC2 == 2
            tInt = varargin{2};               
        else
            error('tInt must have the dimension 1x2');
        end
        
        %vector of weights to scale rms vals
        weights = varargin{3};
    end
end
%*************************************************************************


% Read in each wave files and place in cell array "y"
% for fno=1:length(fnames)
%     Read the wave file, determine its sample frequency
%     [y{fno},nfs(fno)]=audioread(fnames{fno});
%     disp(['loop: ',num2str(fno),' pre-processed std: ',num2str(std(y{fno}))])
%     disp([num2str(max(abs(y{fno})))])
%     
%     If more than one channel present, eliminate all but the first channel
%     [nR,nChanOrig] = size(y{fno});
%     if nChanOrig ~= 1
%         y{fno} = y{fno}(:,1);
%     end
%     
%     if a vector of weights to scale the relative RMS values is present,
%     normalize each signal, then multiply it by the scalar weight.
%     if length(varargin) == 3
%         y{fno} = normc(y{fno});     %normalize the column
%         y{fno} = y{fno} * weights{fno};     %multiply by the weight
%         disp(['loop: ',num2str(fno),' post-processed std: ',num2str(std(y{fno}))])
%     end
% end

%TODO********************************************************************************************************************

for fno=1:length(fnames)
    [y{fno},nfs(fno)]=audioread(fnames{fno});
    %If more than one channel present, eliminate all but the first channel
    [nR,nChanOrig] = size(y{fno});
    if nChanOrig ~= 1
        y{fno} = y{fno}(:,1);
    end
end

% if a vector of weights to scale the relative RMS values is present,
% normalize each signal, then multiply it by the scalar weight.
if length(varargin) == 3
    y_std = zeros(1,length(fnames));
    for k=1:length(fnames)
        y_std(1,k) = std(y{k});
        y{k} = y{k} / y_std(k);
        y{k} = y{k} * weights{k};
    end
end

% If fs is not given, down sample to the signal with the lowest fs
if nargin == 1 || nargin == 2 && numC1 == 2
    fs = min(nfs);
end    

% Resample the wave files
for fno=1:length(fnames)   
    yf{fno} = resample(y{fno},fs,nfs(fno));
    
    % Find length of the resampled signal
    siglen(fno) = length(yf{fno});
end

maxSigLen = max(siglen);
% Conform the number of rows to that of the longest signal
for fno=1:length(fnames)
    if siglen(fno) == maxSigLen
        sig(:,fno) = yf{fno};
    elseif siglen(fno) < maxSigLen
        sig(:,fno) = [yf{fno}; zeros(maxSigLen-siglen(fno),1)];
    else
        error(['You can''t have a longer signal than the longest', ...
            'signal...']);
    end
end

% if tInt is passed in, trim down or zero pad 'sig'
if nargin==3 || nargin==2 && numC1 == 2
    if tInt(1) > tInt(2)
        error('2nd column in tInt must be greater than the 1st column');
    elseif tInt(1) < 0 || tInt(2) < 0
        error('tInt cannot contain negative values');
    end
    
    begin_index = fs*tInt(1)+1;
    end_index = fs*tInt(2);
    
    % Warning if begin_index is greater than signal length
    if begin_index > length(sig)
        warning('Did you mean to create a silent signal?');
    end
    
    % zero pad if end_index is greater than the signal length
    if end_index > length(sig)
        sig = [sig; zeros(end_index-length(sig),length(fnames))];
    end
    
    sig = sig(begin_index:end_index,:);
end

% Write out to a wav file
% For debugging purpose
audiowrite('wav2sigout.wav',sig,fs);