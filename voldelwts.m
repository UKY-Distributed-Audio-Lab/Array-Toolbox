function [rstref, at, mpos] = voldelwts(gridax, mpos, c)
% This function computes delays and weights for a SSL problem given
% the FOV points, mic positions, and speed of sound.
%
%   [rstref, at, mpos] = voldelwts(gridax, mpos, c)
%
% Inputs:
%
%  "gridax"       Cell array of grid points for x-axis, y-axis and z-axis.
%                 The number of elements in this cell array can be less than
%                 3.  If so, the dimension of output space will be reduced.
%
%  "mpos"         Matrix indicating the x, y, and z Cartesian positions of each
%                 mic, where the first row contains the x positions, second y
%                 position, and third z positions.  If an axis is missing,
%                 only the x-y plane will be used.  The dimensions should be
%                 consisent with the axes in cell array "gridax" and number of
%                 mics must be the same as columns in "sigar".
%  "c"            The speed of sound in the room. (scalar value)
%
% Required Functions: 
%           arweights.m
%
% Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005
% Modified by Kirstin Brangers                       July 2012


dimfov = max(size(gridax));    %  Get dimensions for FOV 

% X Dimension
dimx = length(gridax{1});
if dimx > 1
   xax = gridax{1}(:);          % Access all elements in the first cell of gridax
else
   xax = gridax{1}(1);          % Access the first element in the first cell of gridax
end
% Y Dimension
if dimfov >= 2                  %  If second dimension present
    dimy = length(gridax{2});
    if dimy > 1
       yax = gridax{2}(:);      % Access all elements in the second cell of gridax
    else
       yax = gridax{2}(1);      % Access the first element in the second cell of gridax
    end
else                            % If not, add a singleton to hold place
    dimy = 1;
    yax = 0;
end
% Z Dimension
if dimfov == 3                  %  If third dimension present
    dimz = length(gridax{3});
    if dimz > 1;
        zax = gridax{3}(:);     % Access all elements in the third cell of gridax
    else
        zax = gridax{3}(1);     % Access the first element in the third cell of gridax
    end
else                            % If not, add a singleton to hold place
    dimz = 1;
    zax = 0;
end

% Obtain mic array information
[mr,mc] = size(mpos);           % Determine number of mics = mc 

% Get signal dimensionality

% Extend mic coordinates to 3 dimensions with zeros if dimension is less than 3
if mr == 1
   mpos(2:3,:) = zeros(2,mc);
elseif mr == 2
   mpos(3,:) = zeros(1,mc);
end

% Initialize matrices
rstref = zeros(dimx,dimy,dimz,mc);  % Correlation signal matrix
at = zeros(dimx,dimy,dimz,mc);      % Shading weights matrix

%  Loop through every point in FOV
for kz=1:dimz                       % Z-Dimension Loop
    for ky=1:dimy                   % Y-Dimension Loop
        for kx = 1:dimx             % X-Dimension Loop
            %  Distance of FOV position from all microphones
            ds = mpos - [xax(kx); yax(ky); zax(kz)]*ones(1,mc);            
            %  Convert distances to time
            rst = sqrt(ds(1,:).^2 + ds(2,:).^2 + ds(3,:).^2)/c;            
            %  Create shading values to weight mic inputs as function of
            %  distance giving closest mic the most weight.
            at(kx,ky,kz,:) = arweights(rst);            
            % Find mic with maximum delay to FOV
            md = max(rst);            
            % Implement all other delays with respect to the furthest mic
            rstref(kx,ky,kz,:) = md - rst;
            
        end                     %  End X-Dimension Loop
    end                         %  End Y-Dimension Loop
end                             %  End Z-Dimension Loop 
