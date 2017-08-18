function [sd, tn] = delayt(s,fs,d,ns)

%  This function delays the signal in vector S sampled at FS by an amount
%  D denoted in seconds.  If NS is present in the input argument,
%  the output signal length will be equal to this value in seconds (either by
%  truncating or zero padding).  If not present, the output signal length
%  will vary depending on the delay (it will be the original
%  signal length plus the max delay rounded up to the next sample).
%
%    [sd, tn] = delayt(s,fs,d,ns)
%
%  The delayed version of the signal SD starts at the same time point (D should
%  not be negative otherwise an error will occur).  This version performs the delay
%  in the time domain with an FIR shifting (interpolation) filter.  The FIR
%  interpolation filter is a 4-th order weighted SINC function, the default
%  weighting window is the cosine square (one of my favorites).  Comments
%  can be changed inside code to activate different windowing options.
%  IF TN is present in the output, a time axis is created for the output array
%  where 0 denotes the original starting point for S.
%  If S is a matrix, the program will assume each column is a signal and
%  the number of elements in D must be equal to the number of columns.
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu)  July 2005
%   Modified by Kevin D. Donohue September 2014 (corrected scaling issue on
%   coefficients)
%

%  Signal Processing parameters:
ord = 4;    %  Order of interpolation filter
ordh = ceil(ord/2);
sp = (-ordh+1:ordh);  % grid for evaluation of sinc interpolators

%  Find length and dimension of input signal
[r,c] = size(s);

%  Covert to column vecotors for processing
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

sdd = zeros(slen+ord,1);  % Initalize dummy vector for storing integer shift
sd = zeros(slen,c);       % Initialize output matrix

%  Loop through each row of signal matrix and apply delay
for k=1:c
    %  Shift integer sample component of delay
    id = fix(d(k)*fs)+1;
    sdd(id:min([slen,(rn+id-1)])) = s(1:min([slen-id+1, rn]),k);
    %  Interpolate to fractional part of sample part of delay
    fd = d(k)*fs-(id-1);
    t = sp-fd;  % compute noninteger delay axis
%  WINDOW WEIGHT OPTIONS .....
    %  If a different weighted filter is desired, arrange comments
    %  so your desired window FIR filter is uncommented
    h = (cos(0.5*pi*(t)/(max(abs(t))+1)).^2).*sinc(t);  %  Cosine squared windowed sinc
%    h = (1-abs(t)/(max(abs(t))+1)).*sinc(t);   % Triangle windowed sinc
%    h = (.54 -.46*cos(2*pi*(t+min(t))/(max(t)-min(t)))).*sinc(t); % Hamming windowed sinc
    %    For a linear interpolation filter (good for baseband signals with
    %    sharp transitions, comment in the next 2 lines and comment out the
    %    3rd.
%       h = (0:1)-k/(n+1);
%       h = abs(h);
%  END WINDOW WEIGHT OPTIONS .....
    
    dum = filter(h,1,sdd);  % Apply FIR filter
    %  Compensate starting point for integer filter delay
    sd(:,k) = dum(ordh:slen+ordh-1);
end

%  Restore dimension of signal vector to original orientation
if rn == r     % If input was originally a column or multi-signal matrix we are done
   if nargout == 2      %   Create time axis if requested
     tn = [0:slen-1]'/fs;
   end
else      % If input was originally a row vector, take transpose
   sd = (sd.');
   if nargout == 2      %   Create time axis if requested
     tn = [0:slen-1]/fs;
   end
end
