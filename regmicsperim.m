function mpos = regmicsperim(fom,sp,stp)
%  This function generates the microphone coordinates for an 
%  array placement over a perimeter of a rectangle
%
%      mpos = regmicsline(fom, sp, stp)
% 
%  The inputs are the field over which the mics are placed, FOM, and the
%  spacing between mics, SP.  FOM is a 2 column matrix, where each column
%  indicates the opposite corners of the rectangle for the microphone
%  positions.
%  The first mic will be placed at the first point in FOM. The subsequent mics
%  are placed with a constant SP spacing along the perimeter.  If the optional
%  parameter STP is given the first point will be displaced by a distance of SP
%  from the first vertex point toward the next vertex in the clockwise direction.
%  If a mic does not fall evenly on a corner, the distance between the mic on
%  the next edge is taken directly (not along the perimeter) to maintain mutual
%  mic distances. The last mic will be placed to maintain the distance with the prevous
%  mic; however, if it is less than SP to the original starting mic it will
%  be trimmed off.
%  The output will be columnwise coordinates of the microphones in matrix
%  MPOS.  If dimensions of vectors are only 2-D, they will be padded with
%  zeros to make them 3-D.  MPOS will always have 3 dimensions.
%  
%   written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005

[r,c] = size(fom);  %  Dimensions of corner points
%  Create vertex points of rectangle
if r == 2  %  If only 2 dimensions given, expand to a third
   fom = [fom; zeros(1,c)];
end
a = fom(:,2)-fom(:,1);
%  Create all 4 vertices
   v(:,1) = fom(:,1);
   v(:,2) = [fom(1,1); fom(2,2); fom(3,1)+a(3)/2];
   v(:,4) = [fom(1,2); fom(2,1); fom(3,2)-a(3)/2];
   v(:,3) = fom(:,2);

%  loop to step through all edges and place mics
if nargin < 3
   mpos(:,1) = v(:,1);   %  Put first mic at first vertex
else
    dv = v(:,2)-v(:,1);
    nv = dv/norm(dv,2);
    mpos(:,1) = v(:,1)+ stp*nv;
end
kt = 1;   %  Initialize mic counter
for k=1:4
    a = v(:,mod(k,4)+1)-v(:,k);  %  Compute vector between verticies.
    % Determine parametric increment that corresponds to spacing
    %tsp = sqrt(sp^2 / sum(a.^2));
    newa = a / norm(a,2);          %  Normalize vector to get unit direction
    % Compute offset length along new side to maintain mic distance   
    % Apply offset
    if k ~= 1  %  If first time through loop, offset has already been applied           
        cq = sqrt(sp^2-sum(stub.^2));  %  Find length of other triangle leg
        offset = cq*newa;  %  Point in direction of current vertex  
        kt = kt+1;    %  Increment count to place next mic
        mpos(:,kt) = v(:,k)+offset;  %  Position new mic with offset
    end
    %  Once first point established place mics at equal increments along
    %  current perimeter edge
     while norm(mpos(:,kt)-v(:,k),2) <= norm(a,2)
        %  Increment to next mic position
        kt = kt+1;
        mpos(:,kt) = newa*sp+mpos(:,kt-1);
        %plot(mpos(1,:),mpos(2,:),'sr','MarkerSize', 12);
        %pause
     end

% compute value to account for spacing from a vertex
% this is not a multiple of the specifed
% spacing. Step through all 4 sides of rectangle
        stub = v(:,mod(k,4)+1)-mpos(:,kt-1);
        kt = kt-1;  % Call back counter to last valid mic position
end
mpos = mpos(:,1:end-1);
%  If too close to first mic, trim one more point
if norm(mpos(:,1)-mpos(:,end),2) < sp/2
    mpos = mpos(:,1:end-1);
end
%end
