function [sig,fs] = wav2sig(fnames,varargin)
% This function reads in wave files and stores all the information into a 
% single matrix with equal number of rows, with each column representing
% the different wave files. 
%   Inputs:
%   1) fnames - name of the wave files in a cell array. 
%      Format:
%       a) for 1 file: {'filename1.wav'} 
%       b) multiple files: {'filename1.wav';'filename2.wav'}
%       NOTE: 
%        Each wave file should only have one channel.  If more than one 
%        channels are present, only the first channel will be used.
%   2) (struct) vars - optional:
%       a) fs - resample the wave file to this frequency
%       b) tInt - 1x2 matrix to specify time interval (in seconds) to 
%          trim down to. ex: [0 5] for 0->5 seconds
%       c) weights - 1xn array of weights multiplied to each input
%          file after it is normailzed. [w1 w2 w3 wn] -> n = size(fnames)
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


if nargin == 0
    error('wavToSigMat must have at least 1 parameter');
end
if nargin > 2
    error('There can only be a maximum of 2 parameters');
end

%***************************Parameter Checking*****************************
tInt_flag = false;
weight_flag = false;

%if the strucutre is provided
if nargin == 2
    vars = varargin{1}; %assign strucutre
    if isfield(vars,'tInt')
       tInt_flag = true;
       tInt = vars.tInt;
    end
    if isfield(vars,'weights')
       weight_flag = true;
       weights = vars.weights;
    end
    if isfield(vars,'fs')
        fs = vars.fs;
    end
end

%*******************************Signal Checks******************************
for fno=1:length(fnames)
    [y{fno},nfs(fno)]=audioread(fnames{fno});
    %If more than one channel present, eliminate all but the first channel
    [nR,nChanOrig] = size(y{fno});
    if nChanOrig ~= 1
        y{fno} = y{fno}(:,1);
    end
end

if exist('fs') == 0
    fs = min(nfs);
end

%********************************Resample**********************************
for fno=1:length(fnames)   
    yf{fno} = resample(y{fno},fs,nfs(fno));
    
    % Find length of the resampled signal
    siglen(fno) = length(yf{fno});
end


%*******************************Check Sig Length***************************
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


%*****************************Trim time interval***************************
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


%*******************************Normalize**********************************
%9-11-2017 GC: I put the normalization here before the weights were
%applied. The effect of the weights was being minimalized by the next two
%lines whenever they were after the weight scaling in the code.
normme = mean(std(sig));
sig = sig / (10*normme);


%******************************Scalar Weights******************************
% if an array of weights to scale the relative RMS values is present,
% normalize each signal, then multiply it by the scalar weight.
if weight_flag
    sig_std = zeros(1,length(fnames));
    for k=1:length(fnames)
        sig_std(k) = std(sig(:,k));
        sig(:,k) = sig(:,k) / sig_std(k);
        sig(:,k) = sig(:,k) * weights(1,k);
    end
end

% Scale to fit in wavefile and limit clipping
% normme = mean(std(sig));
% sig = sig / (10*normme);
%9-6-17, Grant Cox: this still clips very often. Perhaps we should find the
%max and divide by it? Or put in some logic to decide which method to use.

%9-11-17 GC: clipping varies very much with the different combinations of
%files. Dividing sig by 10 * normme seems to cover most cases, but it also
%looks like it reduces the effect of the weight vector on the signal. i.e.
%changing the weights doesn't do a lot to the amplitude of the outputted
%file

% Write out to a wav file
% For debugging purpose
%audiowrite('wav2sigout.wav',sig,fs);