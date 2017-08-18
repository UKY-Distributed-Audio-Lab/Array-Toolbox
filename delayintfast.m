function sd = delayintfast(s,fs,d,ns)
%
%  This function delays the signal in vector S sampled at FS by an amount
%  D denoted in seconds.  If NS is present in the input argument,
%  the output signal length will be equal to this value in seconds (either by
%  truncating or zero padding).  If not present, the output signal length
%  will vary depending on the delay (it will be the original
%  signal length plus the max delay rounded up to the next sample).
%
%    sd = delayintfast(s,fs,d,ns)
%
%  Delays must be positive and cannot exceed signal length. S must be
%  a column-wise vector or matrix.  The program will not check for these
%  so system will crash if these are not fed in properly.
%  If S is a matrix, the program will assume each column is a signal and
%  the number of elements in D must be equal to the number of columns.
%
%  Written by Kevin D. Donohue (donohue@engr.uky.edu)  July 2012
%

%  Find length and dimension of input signal
[rn,c] = size(s);
%  Determine final output length
slen = ceil(ns*fs);

sdd = zeros(slen,1);    % Initialize dummy vector for storing integer shift
sd = zeros(slen,c);     % Initialize output matrix

%  Loop through each row of signal matrix and apply delay
for k=1:c
    %  Shift integer sample component of delay
    id = round(d(k)*fs)+1;
    sdd(id:(rn+id-1)) = s(1:rn,k);
    sd(:,k) = sdd(1:slen);
end

