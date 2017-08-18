function [sigout, tax] = simimp(f1p,f2p,f1s,f2s,td,fs,siglen)
% This function will create an impulse-like signal associated with 
% a delay and bandwidth. The output includes the signal SIGOUT  and its
% corresponding time axis "tax."  The output is the impulse response of 
% a butterworth filter, the stopband and passbands are constrained by
% the first 4 input parameters.  If the constraints are unreasonable
% (resulting in an order greater than 20), the program will end in an
% error.
% 
%  [sigout, tax] = simimp(f1p,f2p,f1s,f2s,td,fs,siglen)
%
% The input arguments are sound source descriptors:
% f1p =>  Starting frequency of 3dB passband of sound source
% f2p =>  Ending frequency of 3dB passband of sound source
% f1s =>  Ending frequency of first 12dB stopband of sound source (before f1p)
% f2s =>  Starting frequency of second 12dB stopband of sound source (after f2p)
% td =>   time delay in seconds with reference to the first point in 
%         output signal array.
% fs =>   Signal sampling frequency
% siglen => optional input, it is the length of output signal in seconds
%           it(should be at least 2/(f2s-f1s) plus the delay, maybe more if
%            lower frequencies are involved).  The default
%            is td+20/(f2s-f1s) + 5/f1p (I love a good ending).
%   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005


%  Ensure a workable frequency range
if ( (f2p-f1p) <= 0) | (abs(f2p-f1p) > fs/2)
    error('Frequency range must have nonzero positive durration and be less than fs/2')
end
if f1s >= f1p | f2s < f2p
    error('Stopband must be non-overlapping with passband')
end
%  If signal length not provided, set to default
if nargin <= 6
    siglen = (td+20/abs(f2s-f1s)+5/f1p);
end
%  Compute required filter order
n = buttord(2*[f1p, f2p]/fs, 2*[f1s, f2s]/fs, 3, 12);
if n > 10
    error(['Be Reasonable! Require fiter order is ' num2str(n) '. Try making stopband point further away from passband.'])
else
    [b,a] = butter(n,2*[f1p, f2p]/fs);
end
% Initialize impulse array
nsiglen = ceil(fs*siglen);
sigin= zeros(nsiglen,1);
sigin(round(fs*td)+1)=1;  %  Create impulse with specified delay
tax = [0:nsiglen-1]/fs;  %  Create output time axis
sigout = filter(b,a,sigin);  % Pass impulse through filter
sigout = sigout(1:nsiglen);  % trim the impulse respose to specified signal length
%  Test if output was stable
test = max(sigout);   %  If max value out of range kill program     
if test > 100 | test == 0
    error(['Something bad happened to the signal, filter probably unstable. Try relaxing constraints'])
end