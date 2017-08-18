function [outsig,less,act_indexo] =  rmsilence(sig,fs,threshin)
% This function removes the silence intervals from the input speech based
% on an envelope threshold. The input signal is up-sampled, segmented to
% remove samples that fall below a threshold, and then re-sampled back to
% the original sampling rate, and filtered to smooth out the discontinuities
% where pauses in active speech occured.
% The threshold used is a scaled function of the median of the envlope.
% The default threshold is one-fourth of the median of the envelope. 
%
% [outsig,less,act_index] =  rmsilence(sig,fs,threshin)
%  
% INPUTS
% sig-   Vectorized Speech signal input.  If sig has multiple columns
%        the first signal is considered the reference signal in which
%        the pause segments are detected and these are applied to
%        all other columns in the the matrix.
% fs -   Signal sampling frequency in Hertz
% threshin -  (OPTIONAL) Scale on threshold applied to the envelope for
%              detecting scilence periods.  The actual threshold is 
%              is computed by multiplying 'threshin" by the median
%              value of the envelope samples. The default is (.25), which
%              is 25% of the median envelope value over input 'sig'.
% OUTPUTS
% outsig -   Vectorized Signal Output with silence intervals corresponding to
%            the reference signals removed
% less -     Length of the silenced signal segement removed from the original
%            signal in seconds
% act_index - indices of original signal corresponding to active speech
%
% Written by Arulkumaran Muthukumarasamy (arulkumaran@uky.edu ) April 2008
% and updated by Kevin D. Donohue ( donohue@engr.uky.edu ) August 2008

%Check for Correct Parameters
 if nargin > 3
    error('There can only be a maximum of 3 parameters');
 end

  ufs=48000;  %  Upsample limit for better interpolation and filtering
              %   on scilence removal points
 
 nfs=fs; % sampling frequency of the output signal
 [rr,cc]=size(sig); % size of the input signal
 %  Check to ensure vector is a column vector
 if rr < cc
     sig = sig';
     [rr,cc]=size(sig); % size of the input signal
 end
 % Use only first column as signal and remove signal scilence intervals 
 % from corresponding noise sections.
for sno=1:1
    sig_r=[];
    % upsampling the signal if the sampling frequency of input signal is not 44100 hz 
    if (fs ~=48000)
        sig_r=resample(sig(:,sno),ufs,fs);
    end
    env=abs(hilbert(sig_r));% finding the envelope of the signal
    %assigning the threshold for the envelope as a scaled function of the
    %envelope of the signal
    if (nargin ==3)
        thresh=threshin*(median(env));
    else
    % If threshold is not provided, use a third of the envelope median
        thresh=median(env)/3;   
    end
    indx=find(env>thresh);%finding the indexes of those points in the envelope which are above the threshold
    rmsig=sig_r(indx);% extracting the corresponding points from the signal using the indexes
    % rmsig=medfilt1(rmsig,5);% uses a median filter to smoothen the extracted signal
    extsig=resample(rmsig, fs,ufs);% resampling it to Original rate to reduce the distortions in the extracted signal
    less(sno)=(length(sig_r)/ufs)-(length(extsig)/nfs);
    siglen(sno) = length(extsig);
end
act_indexo=resample((indx-1)*fs/ufs, fs,ufs,0);% resampling it to Original rate to reduce the distortions in the extracted signal
act_indexo = floor(act_indexo+1);

if length(act_indexo) ~= siglen
    errordlg('Odd Resampling!',  'Indecies for extracting trimed signal from original do not match extracted signal length');
end
%  Initalized output matrix of signals
outsig = zeros(length(extsig),cc);
outsig(:,1) = extsig;  %  Assign reference signal back to first column  
%  Apply silence intervals to the other signals in the matrix so
%  signals are still synchronized in time.
for k = 2:cc
    if (fs ~=48000)
        sig_r=resample(sig(:,k),ufs,fs);
    end
    rmsig =sig_r(indx);
    % resampling it to Original rate to reduce the distortions in the
    % extracted signal
    extsig = resample(rmsig, fs,ufs);
    outsig(:,k) = extsig;
end

