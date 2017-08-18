function pa = mpos2source(mc,sp,n,rad)
%  This function performs an analysis on microphone subset geometry
%  for mics taken N at a time out of the entire array
%
%      pa = mpos2source(mc,sp,n,rad)
%
%  Coordinates for all microphones are in matrix MC where each
%  column denotes a mic position and each row represents the
%  dimensions of the coordinates. SP is the vector denoting the source
%  position.  The first N columns of the output matrix PA are the
%  indices identifying the mics in the subset, and the indices correspond
%  to columns of MC. The next column is the maximum interdistance between
%  microphones in the subset, then the minimum interdistance, and finaly the standard
%  deviation between all mic, where the spacing to the closest mic used.
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005
%

%  Determine number of mics and dimension of space
[dim,num] = size(mc);

if n > num
    error('number of Microphones in subsets cannot exceed total number of Mics')
end

%  Create all possible combinations of mic subsets
msubsets = nchoosek(1:num,n);  % Get indecies of all mic n-tuples

%  Find number of combinations found
[nsubs, c] = size(msubsets);

%  Initialize output array
pa = zeros(nsubs,n+3);
pa(:,1:n) = msubsets;

%  Compute total number of pairs in subset 
npairs = factorial(n) / (factorial(2)*factorial(n-2));

%  Loop through all pair in the subsets to determin minimum and maximum
%  distances and the standard devaition of all interval between closest
%  pairs
radx = [rad; 0; 0];
rady = [0; rad; 0];
radz = [0; 0; rad];

for k=1:nsubs
    exmn = nchoosek(msubsets(k,:),2);
    %  Look at all mic pair and find the closest and furthest in the set
    d = zeros(1,npairs);
    for pk = 1:npairs
        dp(pk) = sqrt(sum((mc(:,exmn(pk,1))-sp).^2))-sqrt(sum((mc(:,exmn(pk,2))-sp).^2));
        dx(pk) = sqrt(sum((mc(:,exmn(pk,1))-sp+radx).^2))-sqrt(sum((mc(:,exmn(pk,2))-sp+radx).^2));
        dy(pk) = sqrt(sum((mc(:,exmn(pk,1))-sp+rady).^2))-sqrt(sum((mc(:,exmn(pk,2))-sp+rady).^2));
        dz(pk) = sqrt(sum((mc(:,exmn(pk,1))-sp+radz).^2))-sqrt(sum((mc(:,exmn(pk,2))-sp+radz).^2));
    end
    d = [abs(dx-dp), abs(dy-dp), abs(dz-dp)];
    dum = min(d);    %  Find minimum distance
    pa(k,n+1) = dum(1);
    dum = max(d);    %  Find maximum distance
    pa(k,n+2) = dum(1);
    pa(k,n+3) = std(d);
end