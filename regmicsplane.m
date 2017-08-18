function mpos = regmicsplane(fom,sp,geom)
%  This function generates the microphone coordinates for a regular
%  plane array placement over a plane or cube.
%
%      mpos = regmicsplane(fom, sp, geom)
% 
%  The inputs are:
%  FOM => the field over which the mics are placed. FOM is a 3 by 3
%         matrix where the first and last column define the opposite corners
%         of a rectangular space within which the mics are distributed. The
%         2nd column is an additional point that defines the plane.
%  SP =>  the spacing between mics in meters.  If a scalar, then spacing
%         is the same in all dimensions, otherwise if a 3 element vector
%         it defines the spacing in the x,y, and z directions (x,y,z).
%  GEOM => (optional input) flag indicating that the geometry of the spacing
%           is either rectangular or hexagonal.The options are:
%                GEOM = 'Rect' (default) is a rectilinear arrangment
%                GEOM = 'Hex' is a hexagonal spacing (in this case SP must be a scalar).


%  The first mic is placed at the first point (first column) in FOM. The
%  subsequent mics are placed with a constant spacing along the plane defined
%  by the FOM.  The last mic will be placed at the furthest point from the
%  beginning but still within or on the boundaries defined by FOM. 
%  The output will be columnwise coordinates of the microphones in matrix
%  MPOS.  The dimensionality of the output coordinates in MPOS will be the
%  same as those in FOM.
%    Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005, updated
%    June 2008.
%


[r,c] = size(fom);  %  Get dimension of Mic Field
if c ~= 3
    error('The mic plane needs 3 points to define the distriubtion (3 by 3 matrix for 3-D)')
end
lenr = length(sp);  %  Get number of spacings provided
%  If only one spacing is given, apply it to all dimensions
if lenr == 1
    sp(2) = sp(1);
    sp(3) = sp(1);
end

%  Use 3 vector form for equation of plane in 3-D space
%  OP = OA + s*AB +t*AC

%  Generate 3 vectors for pointing to every point in plane
oa = fom(:,1);       %  From origin to first point in plane
ob = fom(:,2)-fom(:,1);  %  From first point to OB or S coordinate
oc = fom(:,3)-fom(:,2);  %  From OB or S coordinate to OC or t coordinate

%  If 3 optional argument not given set to default - Rectilinear
if nargin == 2
    geom(1) = 'r';
end

%  If a rectilinear array
if geom(1) == 'r' | geom(1) == 'R'
    % Determine parametric increment on t and s that corresponds to spacing
    obsinc = sqrt(sp(1)^2 / sum(ob.^2));  % s spacing
    octinc = sqrt(sp(2)^2 / sum(oc.^2));  % t spacing
    %  Determine number of increments 
    nob = fix(1/obsinc+2*eps)+1;  % Number along vector ob (s)
    noc = fix(1/octinc+2*eps)+1;  % Number along vector oc (t)
    mpos = zeros(r,nob*noc);  %  initialized mic coordinate matrix

    t = 0;  %  Starting point for t;
    mc = 0; %  Initialize mic element counter
    %  Loop for t increments
    for k=1:noc
        s = 0;  %  reset starting point for s;
        %  Loop for s increments
        for m=1:nob
            mc = mc+1;  %  update mic element counter
            mpos(:,mc) = oa+ ob*s + oc*t; % compute mic coordinates
            s = s+obsinc; % increment along s direction
        end
    t = t+octinc;  %  increment along t direction 
    end

%  Else if hexagonal array is chosen
elseif geom(1) == 'h' | geom(1) == 'H'
    % Determine parametric increment on t and s that corresponds to spacing
    obshex = sqrt(sp(1)^2 / sum(ob.^2));  % isotropic spacing along s
    octhex = sqrt((sp(1)^2 -(sp(1)/2)^2) / sum(oc.^2));  % isotropic spacing along t
       %  Determine number of increments 
    nob = fix(1/obshex+2*eps)+1;  % Number along vector ob (s)
    noc = fix(1/octhex+2*eps)+1;  % Number along vector oc (t)
    mpos = zeros(r,nob*noc);  %  initalized mic coordinate matrix

    t = 0;  %  Starting point for t;
    mc = 0; %  Initialize mic element counter
    %  Loop for t increments
    for k=1:noc
        %  Compute odd and even row offset 
        if k/2 == fix(k/2);
            s = sqrt((sp(1)/2)^2 / sum(ob.^2));  % isotropic spacing along s
        else
            s = 0;
        end
        %  Loop for s increments
        for m=1:nob
            if s <= 1
                mc = mc+1;  %  update mic element counter
                mpos(:,mc) = oa+ ob*s + oc*t; % compute mic coordinates
                s = s + obshex; % increment along s direction
            end
        end
    t = t+octhex;  %  increment along t direction 
    end
end
mpos = mpos(:,1:mc);  %  Trim array down to actual size