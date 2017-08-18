function [sd, tn] = delayf(s,fs,d,ns)
%  This function delays the signal in vector S sampled at FS by an amount
%  D denoted in seconds.  If NS is present in the input arguments, the output
%  signal length will be equal to this value in seconds (either by
%  truncating or zero padding).  If not present, the output signal length
%  will vary depending on the delay (it will be the original
%  signal length plus the max delay rounded up to the next sample).
%
%    [sd, tn] = delay(s,fs,d,ns)
%
%  The delayed version of the signal SD starts at the same time point (D should
%  not be negative otherwise an error will occur).  This version performs the delay
%  in the frequency domain.  IF TN is present in the output argument, a time axis is
%  created for the output array where 0 denotes the original starting point for S.
%  If S is a matrix, the program will assume each column is a signal and
%  the number of elements in D must be equal to the number of columns.
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu)  July 2005
%

%  Get signal processing parameters:
%  Find length and dimension of input signal
[r,c] = size(s);


%  Covert to column vectors for processing
if r == 1;
    s = s.';  %  Ensure it is a column vector
    [rn,c] = size(s);
elseif c == 1
    [rn,c] = size(s); % Otherwise rn = r indicating it came in as a column vector
else   %  If multi column and row then process each column as its own signal
    [rn,c] = size(s); % Otherwise input was columnwise matrix
    if length(d) ~= c
        error('Number of columns in signal matrix must equal number of elements in delay vector')
    end
end

if min(d) < 0;
    error('Delay cannot be negative')
end
% Compute requested delay in sample points rounded up
nd = ceil(d*fs);

%  Determine final output length
if nargin == 4
    slen = ceil(ns*fs);
else
    slen = max(nd)+rn;   %  Take max of all delays in vector to determine final signal length
end
sd = zeros(slen,c);
if slen <= rn;
    sd = s(1:slen,:);
else
    sd(1:rn,:) = s;
end

%  Create frequency shift vector in the frequency domain
nfft = 2^nextpow2(2*slen);
fax = fs*(-nfft/2:nfft/2-1)'/nfft;

%  Loop through each column of signal matrix and apply delay
for k=1:c
    shft = exp(-j*d(k)*2*pi*fax);   % Frequency function for delay
    shft = ifftshift(shft);         % Make axis compatable with numeric FFT
    fsd = fft(sd(:,k),nfft);        % Take FFT
    fsd = fsd.*shft;                %  Apply delay
    dum = ifft(fsd);                %  Return to time domain
    sd(:,k) = dum(1:slen);          %  Trim time domain signal to required length
end

%  Restore dimension of signal vector to original orientation
if rn == r   %  Was already a column vector or multiple signal vector
   if nargout == 2      %   Create time axis if requested
     tn = (0:slen-1)' /fs;
   end
else     %  If it was originally a column vector, transpose
   sd = sd.';
   if nargout == 2      %   Create time axis if requested
     tn = (0:slen-1)/fs;
   end
end
%  If original signal was real, make output real
if isreal(s)
    sd = real(sd);
end
