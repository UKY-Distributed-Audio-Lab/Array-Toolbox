function [posx, posy, mag] = peakfind2d(sig, tx, ty)
% This function finds local maxima in the 2-D matrix SIG.
% TX and TY are the spatial axes for SIG. The columns are indexed by TX and
% rows are indexed by TY.  If TX and TY are not present they will be assigned
% integer spacings starting with 0.
%
%    [pos, mag] = peakfind(sig, tx, ty)
%
% The output values are in vecotors POSX, and POSY correspond to the peak positions
% relative to the input vectors TX and TY, and MAG is the value at the corresponding
% local maximum.  Points on the edge of the matrix are not included in the
% peak search.
%
%  written August 2005, by Kevin D. Donohue (donohue@engr.uky.edu)

%  Get size and dimension of input Matrix
[r,c] = size(sig);
%  If axes are not given assign them
if nargin == 1;
    tx = [0:c-1];
    ty = [0:r-1];
end

%  Initalize counter and arrays for outputs
count = 0;
posx=[];
posy=[];
mag=[];
%  First and last points are never selected
%  Step through rest of array
for kr=2:r-1
    for kc=2:c-1
        % Find the difference of point with all adjacent neighbors 
        tprod = sig(kr,kc) - sig(kr-1:kr+1,kc-1:kc+1);

        dummin =  min(min(tprod));  %  Find minimum difference
        dummax = max(max(tprod));   % Find maximum difference
        %  Local maximum will have all positive or 0 differences and 
        %   cannot have all zeros, which would be a flat field
        if dummin(1) >= 0 & dummax(1) > 0  
              %  Check for diagonal change too
              count = count +1;
              posy(count) = ty(kr);
              posx(count) = tx(kc);
              mag(count)  = sig(kr,kc);
        end
    end
end
