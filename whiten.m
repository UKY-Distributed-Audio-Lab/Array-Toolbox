function w = whiten(s,pw,flim)
% This function performs a partial whitening of signal contained in vector
% S.  The parameter PW ranges from 0 to 1 where 0 represents no spectral
% whitening and 1 represent complete spectral whitening:
%
%   w = whiten(s,pw,flim)
% 
% The whitening is performed in the frequency domain with parameter PW
% according to:
%
%                        S(f)
%        W(f) = ---------------------------------
%                      (abs(S(f))).^PW
%
%  where the division is only performed for non zero values of S(f), otherwise 
%  W(f) is set to 0.  If optional vector FLIM = [fl, fh] is given, then
%  only the frequencies between FL and FH are whitened.  The rest of the
%  frequencies are set to zero.
%
%  The output W will be the same size as S and is the time domain whitened
%  signal.  If S is a matrix, this program will perform a columnwise
%  whitening of the signal.
%
%      Edited/Updated by Sayed Saghaian(05/12/09)

%  Find length and dimension of input signal

[r,c] = size(s);


if c == 1;
    s = s.';  %  Ensure it is a row vector
    [r,cn] = size(s);
elseif r == 1
    [r,cn] = size(s); % Otherwise cn = c indicating it came in as a row vector
else   %  If multi column and row then process each column as its own signal
    s = s.';
    [r,cn] = size(s); % Otherwise cn = c indicating it came in as a row vector
end
if pw < 0 || pw > 1
    error('Whitening parameter must be between 0 and 1')
end

nfft = 2^nextpow2(2*cn);    %  pad with zeros to prevent circular convolution and make power of 2
w = zeros(r,cn) + 1j*zeros(r,cn);  % Initialize output array

fw = zeros(r,nfft) + 1j*zeros(r,nfft);   %  Initialize frequency array

if nargin == 3
    lowindp = round(flim(1)*(nfft/2))+1;
    highindp = round(flim(2)*(nfft/2))+1;
    lowindn = (nfft)-highindp+2;
    highindn = (nfft)-lowindp+2;
    frq = [lowindp:highindp, lowindn:highindn];
else
    frq = (1:nfft);
end
%  Loop to take FFT of each row and apply whitening operation
%  for all non-zero frequency components
for k=1:r
    fs = fft(s(k,:),nfft); %  Take FFT
    fw(k,frq) = fs(frq);
    nzt = find(abs(fw(k,:)) ~= 0);  %  Find non-zero spectral values
    fw(k,nzt) = fw(k,nzt) ./ abs(fw(k,nzt)).^pw;  %  Normalize magnitudes for non-zero values
    dum = ifft(fw(k,:));   %  return to time domain
    w(k,:) = dum(1:cn);  %  Remove zero padded samples
end

%  Restore dimension of signal vector to original orientation
if cn == c && r == 1
   w = w;
else
   w = (w.');
end
%  If original signal was real, make output real
if isreal(s)
    w = real(w);
end
