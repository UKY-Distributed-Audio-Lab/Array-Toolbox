function pa = mposanaly(mc,n)
%  This function performs an analysis on microphone subset geometry
%  for mics taken N at a time out of the entire array
%
%      pa = mposanaly(mc,n)
%
%  Coordinates for all microphones are in matrix MC where each
%  column denotes a mic position and each row represents the dimensions of
%  the coordinates.  The first N columns of the output matrix PA are the
%  indices identifying the mics in the subset, and the indices correspond
%  to columns of MC. The next column is the maximum distance between
%  microphones in the subset, then the minimum distance, and finally the standard
%  deviation between all mics, where the spacing to the closest mic is used.
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

%  Loop through all pair in the subsets to determine minimum and maximum
%  distances and the standard deviation of all interval between closest
%  pairs
for k=1:nsubs
    exmn = nchoosek(msubsets(k,:),2);
    %  Look at all mic pairs and find the closest and furthest in the set
    d = zeros(1,npairs);
    for pk = 1:npairs
        d(pk) = sqrt(sum((mc(:,exmn(pk,1))-mc(:,exmn(pk,2))).^2));
    end
    dum = min(d);    %  Find minimum distance
    pa(k,n+1) = dum(1);
    dum = max(d);    %  Find maximum distance
    pa(k,n+2) = dum(1);
    
    %  Compute standard deviation between all mic pairs in subset
    if n ==2   %  In only 2 mics in subset there is no variation so
        pa(k,n+3) = 0;
    else
        %  Step through all mics and find sequence of minimum distances
        %  between mic pairs
        mindist = [];   % Initialize array to store the minimum distances (to closest neighbor) for each mic
        %  Loop through all first mics in pair
        for pk=1:n-1
            d =[];  %  Initialize the array to store the distances to all pairs
            %  Loop to compute distance to all pairs
            for nk=pk+1:n
                d(nk-pk) = sqrt(sum((mc(:,msubsets(k,pk))-mc(:,msubsets(k,nk))).^2));
            end
            dum = min(d);
            mndist(pk)=dum(1);  % Find closest mic distance
        end
        pa(k,n+3) = std(mndist);  %  Compute standard deviation of all distances
    end
end


