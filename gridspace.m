function grd = gridspace(p1, p2, inc)
% This function computes the coordinates of a rectangular grid space
% where P1 and P2 denote the oposite corners of the rectangle.  INC is
% the increment between grid points.
%
%     grd = gridspace(p1, p2, inc)
% 
%  The output GRD is a cell array with N vectors containing coordinates in the space
%  where N is the same dimension as vectors P1 and P2.  If INC is a scaler
%  then the increment between gridpoints is the same in all dimensions.
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005


% Determine number of dimesnsion for grid.
dim = length(p1);
% Determine nature of increment in each dimension
if length(inc) < dim
    inc = inc(1)*ones(1,dim);  %  make increments the same in all directions
end

% Initalize output cell array
grd = cell(1,dim);
% Create grid vectors from the smaller value corner coordinate to the
% larger one
for k=1:dim
    grd{k} = [min([p1(k),p2(k)]):inc(k):max([p1(k),p2(k)])];
end
