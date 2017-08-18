function [E,bcf]=spectrumlevel(sig,fs)
% This function estimates the spectrum power levels of signals in the
% columns of SIG.  The individual bands follow the band based on the 
% one-third-octave-band method. 
%
%            [E,bcf] = spectrumlevel(sig,fs)
%    
% The inputs are the signal in vector 'sig' sampled at a sampling
% frequency 'fs' (multiple signals are stored column-wise in 'sig'.
% The output is the matrix 'E' where each row is the power level of
% the individual bands. Each column of 'sig' corresponds to columns in
% output matrix 'E'. The second output argument 'bcf' is a vector of
% the center frequencies of the spectral bands.
% The one-third-octave-band procedure is a standard procedure, which
% divides the total spectrum into 18 bands, wherein each band is centered
% at the center frequencies 'fc' given below with a bandwidth of one
% third of the center frequency.
%   fc =[160 200 250 315 400 500 630 800 1000 1250 1600 2000, ...
%        2500 3150 4000 5000 6300 8000];
%
% Written by Kevin D. Donohue (donohue@engr.uky.edu)
%  and Arulkumaran Muthukumarasamy(arulkumaran@uky.edu) January 2008

[rr,cc] =size(sig); % size of the input signal
nfft=2*rr; % number of fft points
bcf=[160 200 250 315 400 500 630 800 1000 1250 1600 2000, ...
     2500 3150 4000 5000 6300 8000]; %center frequencies according to one third octave band
 % estimating the lower and upper frequencies from the center freq of the
 % bands individually
 for i=1:length(bcf)
    fc= bcf(i); 
   %Quadradic formula coefficents to calculate the upper and lower band
   %frequencies
    a=1;
    b=(fc/3);
    c=(fc^2);
    fl(i)=round((-b+sqrt(b^2+4*a*c))/2);% lower freq of the band
    fu(i)=round((fc/3)+fl(i));   % upper freq of the band
 end
 % finding the indices corresponding to the upper and lower band frequencies
   fli=1+round((nfft-1)*fl/fs);
   fui=1+round((nfft-1)*fu/fs);
 for k=1:cc
     pf=fft(sig(:,k),nfft);% estimating the fft of the signal
     for n=1:length(fli) 
     E(n,k)=sum(abs(pf(fli(n):fui(n)).^2));%estimating the spectrum level in each of the induvidual bands
     end
 end
 