function ar = arweights(rst)
%  This function computes a set of weights for each sensor in an
%  array based on the distances from sources to the sensors.
%
%    ar = arweights(rst)
%
%  input:
%  RST    vector of distances of each sensor to the point in space
%  output:
%  AR     corresponding vector of weights for the array elements
%
%    Written by Kevin D. Donohue August 2005
%
[mt, nmic] = min(rst);  % Find closest mic to point being considered
%  Create shading values to weight mic inputs as function of
%  distance giving closest mic the most weight.
ar = 1 ./ ((rst+eps));
ar = ar/(sum(ar)+eps);