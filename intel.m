function [ssi_val,tax] = intel(sig,nos,fs,twin)
% This function is used to estimate the intelligibility of the given
% speech-in-noise condition, where SIG is the signal of interest
% and NOS is the correpsonding noise or interfering signal.  Note 
% that signal and noise samples must be separable for this function
% to provide accurate estimates. Speech Intelligibility can be defined
% as the degree to which the speech can be understood correctly by
% the listener.
% The function calls the function 'sii.m' to calculate the speech
% intelligibility index based on the standard ANSI S3.5 1997 "Methods
% for calculation of the Speech Intelligibility Index" (Refer
% http://www.sii.to/html/programs.html).  Since intelligibility can change
% over hte course of a speech sample, this function breaks this extended
% signal into smaller segement and estimates the intelligibility index
% sequentially over the lenth of the signal:
%
%            [ssi_val,tax] = intel (sig,nos,fs,twin)
%
% The output is a vector of speech intelligibility indices which range
% between zero(completely unintelligible)and one(perfect intelligibility).
% The second output argument is the corresponding time axis in seconds
% where 0 denotes the starting point of the input signal. This function
% also calls 'spectrumlevel.m,' which estimates the spectrum levels at
% different bands of the signal segments.
%
% The input arguments are:
%   SIG => Vectorized Signal segment whose intelligibility has to be estimated in the given speech-in-noise condition.
%   NOS => Noise Signal segment(vector), SIG and NOS must have the same length.
%   fs => sampling frequency in Hertz
%   TWIN => Time interval in seconds which is used to segment the sig and nos
%           into 'twin' overlapping blocks(OPTIONAL)
%Dependencies
% -> spectrumlevel.m
% -> sii.m
% 
% Written by Kevin D. Donohue (donohue@engr.uky.edu)& Arulkumaran Muthukumarasamy(arulkumaran@uky.edu) January 2008

%------- Parameter Check ------
 if nargin > 4
    error('There can only be a maximum of 4 parameters');
 end
 % check for the size of the signal and noise
 if (size(sig)~=size(nos))
    error('The signal and noise must be of same length');
 end
 % Check for the time window parameter which is optional
if nargin==3
    tint=100e-3; %default time window
else
    tint=twin;
end

[r,c]=size(sig);%size of the signal
nwin=round(fs*tint);% window size in samples which is used to move along the signal length and calculate the intelligibility
k=0;
bp=1;% beginnning point of the window
ep=nwin;%end point of the window
tinc=round(nwin/2);%used to increment and move the window along the signal overlapping(50%) the previous window

while(ep<=r)
    k=k+1;
    tax(k)=((ep+bp)/2-1)/fs;%computing the time axis corresponding to sii
    Es=spectrumlevel(sig(bp:ep),fs);% estimating the spectrum level of the input speech signal segment
    En=spectrumlevel(nos(bp:ep),fs); %estimating the spectrum level of the input noise signal  segment
    ssi_val(k)=sii('E',10*log10(Es),'N',10*log10(En),'I',7);% calculating the intelligibility for a particular segment of the signal
    bp=tinc+bp-1; %reset the beginning point of next window to center of the current window
    ep=bp+nwin-1; %reset the end point of next window such that it is 'nwin' samples(window size)apart from the beginning of the window
end
