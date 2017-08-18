function w = flattap(len,pct)
%  This function creates a tapering window of length LEN that is 
%  flat except for a percentage toward the begining and the end, which
%  has a hann tapper:
%
%    w = flattap(len,pct)     
%
%  PCT is the percentage (between 0 (boxcar) and 50 (hann)) of tapper
%  at the window ends, and W is a column vector.
%
%     written by Kevin D. Donohue (donohue@engr.uky.edu) August 2005
%

%  Percentage can't be over 50% or less than 0
if pct > 50   %  if it is clip to 50
    pct = 50;
elseif pct < 0   %  if it is clip to 0
    pct = 0;
end
wh1 = round(pct*len/100);  %  Find number of samples corresponding to percentage
tax = [0:len-1]';  %  Compute time axis
w = ones(size(tax));  %  Initialize a flat window
ind1 = find(tax < wh1);   % find point of percentage at beginning of window
ind2 = find(tax > len-1-wh1);   %  Find points of percentage at end of window
tapwin = hann(length(ind1)+length(ind2));  %  Compute both sides of hann window
w(ind1) = tapwin(1:length(ind1));       %  set first half to beginning of window
w(ind2) = tapwin(length(ind2)+1:end);   %  set last half to end of window
