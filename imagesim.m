function [dlays, scals] = imagesim(forv, s, m, bs, c, dbd)
%  This function generates a sequence of time delays and scale factors associate with
%  the reflections from a sound source and microphone in a rectangular room according
%  to the image method described by J.B. Allen and D.A. Berkley in "Image method for
%  efficiently simulating small-room acoustics," J. Acoust. Soc. Am., Vol. 65, No. 4,
%  Apr. 1979.  The scales and delays are useful for computing.  
%  The function requires the opposite corners of the room in matrix FORV,
%  where each column is the 3D coordinates [x,y,z]' of the room, vector BS
%  contains the wall reflection coefficients, where first 4 elements correspond to
%  the 4 walls ordered clockwise from the first point of FORV (looking down on the
%  room), and then last 2 elements are the reflection coefficient of the floor and
%  ceiling, respectively.  The coordinates of the source are contained in
%  vector S, coordinates of the microphone in M, the speed of sound is given by C,
%  and optional parameter DBD is the number of decibels below which the scale factors
%  will not be included.  If DBD is not provided, the default value is set to 60 dB.
%
%                 [dlays, scals] = imagesim(forv, s, m, bs, c, dbd)
%
%   The output vectors are:
%   DLAYS =>  the time delays in seconds for each source image, including the original
%             source which is the first element in the array.
%   SCALS =>  the scale factors the received sound as a result of spherical wave propagation
%             and reflection losses for multipath signals.
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu)  September 2005
%   (Updated June 2008)


%  If minimum dB value not given, set to default 
if nargin < 6
    dbd = 60;
end

% Test validity of reflection coefficients
if (max(bs) > 1) || (min(bs) < 0)
    error('  !!All reflection coefficient must be in the interval [0 1]!! ')
end

%  Locate lower, bottom, left corner of room (looking down on cartesion
%  room and standard cartesion axes).
refvertx = min(forv(1,:));  %  Left most coordinate
refverty = min(forv(2,:));  %  Bottom most coordinate
refvertz = min(forv(3,:));  %  lower most coordinate
%  Use lower, bottom, left point as the reference  [0 0 0] for new coordinate
%  system
refv = [refvertx(1);refverty(1);refvertz(1)];  % Reference to scale to 0
s = s-refv;   %  Shift source position
m = m-refv;   %  Shift mic position
forv = forv-refv*ones(1,2);  %  Shift room verticies

%  compute room dimensions
a = forv(:,2)-forv(:,1);  %  Room dimensions
rd = abs(a);  %  Absolute distances of room dimensions

%  Get all room verticies
if a(3) > 0 && a(1)*a(2) > 0  %  First point is on the floor, x,y same sign
    v(:,1) = [forv(1,1); forv(2,1); forv(3,1)];
    v(:,2) = [forv(1,1); forv(2,2); forv(3,1)];
    v(:,3) = [forv(1,2); forv(2,2); forv(3,1)];
    v(:,4) = [forv(1,2); forv(2,1); forv(3,1)];
    v(:,5) = [forv(1,1); forv(2,1); forv(3,2)];
    v(:,6) = [forv(1,1); forv(2,2); forv(3,2)];
    v(:,7) = [forv(1,2); forv(2,2); forv(3,2)];
    v(:,8) = [forv(1,2); forv(2,1); forv(3,2)];
elseif  a(3) < 0 && a(1)*a(2) > 0 % First point on the ceiling, x,y same sign
    v(:,1) = [forv(1,1); forv(2,1); forv(3,2)];
    v(:,2) = [forv(1,1); forv(2,2); forv(3,2)];
    v(:,3) = [forv(1,2); forv(2,2); forv(3,2)];
    v(:,4) = [forv(1,2); forv(2,1); forv(3,2)];
    v(:,5) = [forv(1,1); forv(2,1); forv(3,1)];
    v(:,6) = [forv(1,1); forv(2,2); forv(3,1)];
    v(:,7) = [forv(1,2); forv(2,2); forv(3,1)];
    v(:,8) = [forv(1,2); forv(2,1); forv(3,1)];
elseif  a(3) < 0 && a(1)*a(2) < 0 % First point on the celing, x,y different sign
    v(:,1) = [forv(1,1); forv(2,1); forv(3,2)];
    v(:,4) = [forv(1,1); forv(2,2); forv(3,2)];
    v(:,3) = [forv(1,2); forv(2,2); forv(3,2)];
    v(:,2) = [forv(1,2); forv(2,1); forv(3,2)];
    v(:,5) = [forv(1,1); forv(2,1); forv(3,1)];
    v(:,8) = [forv(1,1); forv(2,2); forv(3,1)];
    v(:,7) = [forv(1,2); forv(2,2); forv(3,1)];
    v(:,6) = [forv(1,2); forv(2,1); forv(3,1)];
else   %  First point is on the floor, x,y different sign
    v(:,1) = [forv(1,1); forv(2,1); forv(3,1)];
    v(:,4) = [forv(1,1); forv(2,2); forv(3,1)];
    v(:,3) = [forv(1,2); forv(2,2); forv(3,1)];
    v(:,2) = [forv(1,2); forv(2,1); forv(3,1)];
    v(:,5) = [forv(1,1); forv(2,1); forv(3,2)];
    v(:,8) = [forv(1,1); forv(2,2); forv(3,2)];
    v(:,7) = [forv(1,2); forv(2,2); forv(3,2)];
    v(:,6) = [forv(1,2); forv(2,1); forv(3,2)];
end

% Find location of lower, bottom, left corner vertex so reflection
% coefficients of the walls can be properly assigned
gb = find(sum(v) == 0);
vst = mod(gb(1)-1,4)+1;  %  It does not matter if point is on floor or ceiling

% Assign wall reflection coefficients
if vst == 1    %  if first vertex at the reference point
    by1 = bs(1);
    bx2 = bs(2);
    by2 = bs(3);
    bx1 = bs(4);
elseif vst == 2  %  If first vertex one away from reference point in clockwise direction
    by1 = bs(4);
    bx2 = bs(1);
    by2 = bs(2);
    bx1 = bs(3);
elseif vst == 3 %  If first vertex 2 away from reference point in clockwise direction
    by1 = bs(3);
    bx2 = bs(4);
    by2 = bs(1);
    bx1 = bs(2);
else   %  If first vertex 1 away in counter clockwise direction.
    by1 = bs(2);
    bx2 = bs(3);
    by2 = bs(4);
    bx1 = bs(1);
end
bz1 = bs(5);  %  Floor reflection coefficient
bz2 = bs(6);   %  Ceiling reflection coefficient


%  Get all image permutations for first order image sources
pv = zeros(3,8);
Rp = zeros(3,8);
for p=1:8
    pv(:,p) = [1+(-1)^ceil(p/4); 1+(-1)^ceil(p/2); 1+(-1)^ceil(p)]/2;  %  Binary permutation for indexing
    Rp(:,p) = [s(1)-m(1)+2*m(1)*pv(1,p); s(2)-m(2)+2*m(2)*pv(2,p); s(3)-m(3)+2*m(3)*pv(3,p)]; % image point permutations
end

%  Find threshold for dB down amplitudes
fac = 4*pi;
thr = 10^(-dbd/20);  %  Threshold

%  Set x-dimensions to upper bound
nxmax = 1;
while nxmax < (bx1*bx2)^nxmax/(fac*thr*rd(1));
    nxmax=nxmax+1;
end

%  Set y-dimensions to upper bound
nymax = 1;
while nymax < (by1*by2)^nymax/(fac*thr*rd(2));
    nymax=nymax+1;
end

%  Set z-dimensions to upper bound
nzmax = 1;
while nzmax < (bz1*bz2)^nzmax/(fac*thr*rd(3)+eps);
    nzmax = nzmax+1;
end

%  Set up initial array for outputs based on max number of points computed
dlays = zeros(1,(2*nxmax+1)*(2*nymax+1)*(2*nzmax+1)*8);
scals = zeros(1,(2*nxmax+1)*(2*nymax+1)*(2*nzmax+1)*8);
%  Loop to compute delays and magnitudes for each image
cnt = 0;  %  Initialize counter for output array
Rr=zeros(3,1); % Initialize image position vector
for p=1:8  % First order permutations
    for n=-nxmax:nxmax  % Multiple images in x-direction
        bx=bx1^abs(n-pv(1,p))*bx2^abs(n);
        Rr(1) = 2*n*rd(1);
        for l=-nymax:nymax %  Multiple images in y-direction
            bxy=bx*by1^abs(l-pv(2,p))*by2^abs(l);
            Rr(2) = 2*l*rd(2);
            for m=-nzmax:nzmax  %  Multiple images in z-direction
                Rr(3) = 2*m*rd(3);
                bxyz = bxy*bz1^abs(m-pv(3,p))*bz2^abs(m); %  Compute losses from all reflections
                dst = norm(Rp(:,p)+Rr);  %  Compute total distance traveled for image source
                stst = bxyz/(fac*(dst+1));   %  Combined reflection and diffraction losses
                %  Test to see if dB level of scale factor is above limit
                if stst > thr   %  If so, save delay and scale factor to output array
                    cnt = cnt+1;
                    scals(cnt) = stst;
                    dlays(cnt) = dst/c;
                end
            end
        end
    end
end
%  Trim ouput array of excess zeros
dlays = dlays(1:cnt);
scals = scals(1:cnt);
%  Sort output array according to smallest to largest delay
[dlays, in] = sort(dlays);
scals = scals(in);
