function mpos = regmicsline(fom,sp)
%  This function generates the microphone coordinates for regular
%  line array placement over a line, plane, or cube.
%
%      mpos = regmicsline(fom, sp)
% 
%  The inputs are the field over which the mics are placed, FOM, and the
%  spacing between mics, SP.  FOM is a 2 column matrix, where each column
%  indicates the begining and ending point of the microphones line array.
%  The first mic will be placed at the first point in FOM. The subsequent
%  mics are placed with a constant SP spacing along the line
%  connecting the 2 points in FOM.  The last mic will be placed at the
%  furthest point from the beginning that is less than or equal to the
%  distance to the last point in FOM.
%  The output will be columnwise coordinates of the microphones in matrix MPOS
%  The dimensionality of the output coordinates in MPOS will be the same as
%  those in FOM.
%
%   written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005

[r,c] = size(fom);  %  Dimensions of corner points

%  Use 2-point form for equation of line in 3-D space
%  t = (x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)

%  Find distances from first to last point along each dimension
a = fom(:,2)-fom(:,1);

% Determine parametric increment on t that corresponds to spacing
tsp = sqrt(sp^2 / sum(a.^2));
t = 0;  %  Starting point for t
%  Determine number of mics on array
numinc = fix(1/tsp + 2*eps) + 1;
mpos = zeros(r,numinc);  %  Initialize output array
%  Loop to compute coordinates for mic position along the line
for k=1:numinc 
   mpos(:,k) = a*t+fom(:,1);
   t = t+tsp;  %  Increment to next position
end
