%edited by Grant Cox, September 6, 2017

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
%       b) tInt - 1x2 matrix to specify time interval (in seconds) to 
%          trim down to. ex: [0 5] for 0->5 seconds
%       c) weights - scaling vector of weights multiplied to each input
%          file after it is normailzed. { w1; w2;...;wn } n = size(fnames)
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


tInt_flag = false;
weight_flag = false;
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
            weight_flag = true;
        elseif numC1 == 1
            fs = varargin{1};
        elseif numC1 == 2
            tInt = varargin{1};
            tInt_flag = true;
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
            weights = varargin{2};
            weight_flag = true;
        elseif numC2 == 2
            tInt = varargin{2};
            tInt_flag = true;
        else
            error('Check argument specifications.');
        end
        
        %TODO FIX THE REST OF THE INPUTS FROM HERE. USE SIMILAR METHOD TO
        %PREVIOUS LINES
        
%         if numR2 ~= 1
%             error('tInt must have the dimension 1x1 or 1x2');
%         end
        
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
            tInt_flag = true;
        else
            error('tInt must have the dimension 1x2');
        end
        
        %vector of weights to scale rms vals
        weights = varargin{3};
        weight_flag = true;
    end
end
%*************************************************************************

for fno=1:length(fnames)
    [y{fno},nfs(fno)]=audioread(fnames{fno});
    %If more than one channel present, eliminate all but the first channel
    [nR,nChanOrig] = size(y{fno});
    if nChanOrig ~= 1
        y{fno} = y{fno}(:,1);
    end
end

% % if a vector of weights to scale the relative RMS values is present,
% % normalize each signal, then multiply it by the scalar weight.
% if weight_flag
%     y_std = zeros(1,length(fnames));
%     for k=1:length(fnames)
%         y_std(1,k) = std(y{k});
%         y{k} = y{k} / y_std(k);
%         y{k} = y{k} * weights{k};
%     end
% end

% If fs is not given, down sample to the signal with the lowest fs
if exist('fs') == 0
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
if tInt_flag
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

% if a vector of weights to scale the relative RMS values is present,
% normalize each signal, then multiply it by the scalar weight.
if weight_flag
    sig_std = zeros(1,length(fnames));
    for k=1:length(fnames)
        sig_std(k) = std(sig(:,k));
        sig(:,k) = sig(:,k) / sig_std(k);
        sig(:,k) = sig(:,k) * weights{k};
    end
end

% Scale to fit in wavefile and limit clipping
normme = mean(std(sig));
sig = sig / (10*normme);
%9-6-17, Grant Cox: this still clips very often. Perhaps we should find the
%max and divide by it? Or put in some logic to decide which method to use.

% Write out to a wav file
% For debugging purpose
audiowrite('wav2sigout.wav',sig,fs);