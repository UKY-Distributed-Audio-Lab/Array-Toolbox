function mpos = regmicsperim2(fom,sp,num, stp)
%  This function generates the microphone coordinates for an 
%  array placement over a perimeter of a rectangle
%
%      mpos = regmicsperim2(fom, sp, num, stp)
% 
%  The inputs are::
%  FOM =>  the field over which the mics are place defined by the lower left and upper right
%          verticies of the rectangular plane or room. FOM is a 2 column matrix, where each column
%          indicates the opposite corners or the rectangle for the
%          microphone positions.  If vector in FOM are 2-D they will be
%          padded with zeros to make them 3-D.
%  SP =>   uniform spacing between mics
%  NUM =>  the number of mics.
%  STP = > optional parameter indicating the first mic position along the perimeter.
%
%  Mics will be place in a clockwise direction (looking down on the field  of view) from the initial
%  point, which is either STP if given, the first point in FOM. The subsequent mics are placed with
%  a constant SP spacing along the perimeter in the clockwise direction.  If a mic does not fall evenly
%  on a corner, the distance between the mic on the next edge is taken directly (not along the rectangular
%  perimeter). If the number and mic spacing result in a length greater than the
%  perimeter, it will wrap around.
%
%  The output will be columnwise coordinates of the microphones in matrix MPOS
%  The dimensionality of the output coordinates in MPOS will be the same as
%  those in FOM.
%
%   written by Kevin D. Donohue (donohue@engr.uky.edu) October 2005

[r,c] = size(fom);  %  Dimensions of corner points
%  Create vertex points of rectangle
if r == 2  %  If only 2 dimension given, expand to a third
   fom = [fom; zeros(1,c)];
end
%  Create all 4 vertices
a = fom(:,2)-fom(:,1);
v(:,1) = fom(:,1);
v(:,2) = [fom(1,1); fom(2,2); fom(3,1)+a(3)/2];
v(:,4) = [fom(1,2); fom(2,1); fom(3,2)-a(3)/2];
v(:,3) = fom(:,2);


%  Find where first mic should be placed
if nargin < 4  %  If initial point not given, set first point to first vertex
   mpos(:,1) = v(:,1);   %  Put first mic at first vertex
   dv = v(:,2)-v(:,1);   %  Vector pointing to next vertex
   nv = dv/norm(dv,2);   %  Clockwise unit direction vector
   idst = 1;
else  %  If first point specified
    %  Expand dimension of STP if necessary
    sr = length(stp);
    if sr < 3
        stp = [stp; 0]
    end
    
    %  Find 2 verticies that stp is between
    for kd=1:4  %  Find distance to all verticies
        distv(kd) = norm(v(:,kd)-stp,2);
    end
    [ds,ids] = sort(distv);
    %  Find 2 closest verticies
    if max(ids(1:2)) == 4 & min(ids(1:2)) == 1
        dv = v(:,4) - v(:,1);  %  Point in clockwise direction
        idst = 1;   %  Identify index of vertex behind point
    else
        dv = v(:,min(ids(1:2))) - v(:,max(ids(1:2)));  %  Point in clockwise direction
        idst = max(ids(1:2));   %  Identify index of vertex behind point
    end
    nv = dv/norm(dv,2);  %  Clockwise unit direction vector
    mpos(:,1) =stp;   %  First mic placement at stp
end
   



% step through all 4 sides of rectangle
kt = 1;   %  Initalize mic counter
k = 0;
while k<5 & kt < num
    a = v(:,mod(k+idst-1,4)+1)-v(:,mod(k+idst-2,4)+1);  %  Compute vector between verticies.
     newa = a / norm(a,2);          %  Normalize vector to get unit direction
    %  Check to see if next point exceed the end vertex
    stub = (mpos(:,kt)-v(:,mod(k-1+idst,4)+1));
    dirv = (newa')*stub;
    if dirv <= 0
        kt=kt+1;
        mpos(:,kt) = mpos(:,kt-1) + (newa)*sp;    
    else  %  Turn corner
     k=k+1;
     a = v(:,mod(k-1+idst,4)+1)-v(:,mod(k+idst-2,4)+1);  %  Compute vector between verticies.
     newa = a / norm(a,2);          %  Normalize vector to get unit direction
     stub = newa*sqrt(sp^2-(sp-norm(stub,2))^2);
%     if norm(stub,2) == 0
  %       mpos(:,kt) = v(:,mod(k+idst-2,4)+1) + newa*sp;
 %    else
         mpos(:,kt) = v(:,mod(k+idst-2,4)+1) + stub;
   %  end
    end
    
end

