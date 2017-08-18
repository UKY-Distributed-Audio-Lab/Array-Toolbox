function imseg = imsegmb(tpos,srez,frez,seglen,lim)
% 
% Function to segment audio scene using spatial and time threshold limits
%
%        IMSEG = imsegm(TPOS,SREZ,FREZ,SEGLEN,LIM)
% 
% Inputs:
%   TPOS --->    N X 5 matrix
%                   col1      --> Intensity
%                   col2,3,4  --> x,y,Z co-ordinates  
%                   col5      --> Frame number
%                   Note: If no target is present in the frame, intensity=0
%   SREZ --->   The radius of the spherical area to search for a candidate for
%               any given segment
%   FREZ --->   Number of frames to look for the candidate for the current segment
%   SEGLEN -->  Minimum segment length in terms of frame number. Any segments
%               of shorter duration will be neglected.
%   LIM -->     [min max] specifies the limits for the intensities of the point.
%               Any point with intensities outside this limit is omitted
%
% Output:
%   IMSEG --->  N X 7 matrix
%                  col 1    ---> Segment number.
%                  col 2-6  ---> TPOS
%                  col 7    ---> Detection indicator
%                      0 for no detection
%                      -1 if omitted due to SEGLEN criteria
%                      -2 non maximum detction for the same segment
%                      1 if detected and belongs to a valid segment
%
%   Note:   
%   Though the limits LIM is optional, it's highly recommended to use it.
%
%   Conflict between candidates in the same frame for same segment
%   is resolved by taking the point with max magnitude. Other
%   candidates Detection Indicator (COL 7) is set to -2
%
%   Segments dropped because of segment length threshold are assigned
%   Detection Indicator (COL 7) as -1
%
%   Frames with no detection are marked with one pseudo point with intensity 0
%   Detection Indicator (COL 7) as 0
%	
% Written by Harikrishnan Unnikrishnan  (harikrishnan@uky.edu)   June 2008
% Modified by Kevin D. Donohue          (donohue@engr.uky.edu)   June 2012
% Modified by Harikrishnan Unnikrishnan (harikrishnan@uky.edu) August 2012
%           - Detected false stream using Detection Threshold



% Apply secondary threshold on TPOS  (initial done in dent3dwei.m function) 
 if nargin >4
    if length(lim) ~= 2
        error('LIM must be of the form [min_value max_value]');
    end
    lfram = max(tpos(:,5));
	
	% Ignore frames whose intensity is not between the range in 'LIM'
	tpos =  tpos(tpos(:,1)<=lim(2) & tpos(:,1)>=lim(1),:);    
end;


frame_length = max(tpos(:,5));  % 5th column contains the frame numbers

% Initialize output array
imseg  = zeros(1,8);
 
hwait = waitbar(0,['Processing frame ' int2str(0) ' of ' int2str(frame_length)]);
% Pass One : Assign initial segment number using <SREZ> and <FREZ> conditions
% Store the points and segment number to IMSEG, adding 2 elements to TPOS
%	Column 1:  Segment number
%	Column 2:  Index to TPOS for curPt (current point)
for fnum =1:frame_length;	% Iterate for all points in the frame
   
    % Find timestamps
    frm = find(tpos(:,5)==fnum & tpos(:,1) >0); 
   
    % Find indices to TPOS that correspond to the frame window
    cfram_indx = find(tpos(:,5)>fnum & tpos(:,5)<fnum+frez & tpos(:,1) >0);
    
  
    for kk=1:length(frm)                % Iterate for all points in the frame
		% Index for the current point in the frame
		currPt = frm(kk);    
		
        for jj =1: length(cfram_indx)	% Iterate for all points in the window
			% Index to TPOS for the point in the window
			winPt = cfram_indx(jj);
           
		    % Calculate the distance between currPt and winPt
            d = sqrt(sum((tpos(currPt,2:4)-tpos(winPt,2:4)).^2));
 
			% If currPt and winPt are within <SREZ> of each other
            if d<srez
				% See if currPt already has an imseg_num associated
                lkup =  imseg(imseg (:,2)==currPt,1);

					if (isempty(lkup))  %  If  not create new stream
						imseg_num = 1+max(imseg (:,1));
						imseg  = [imseg ;imseg_num,currPt,tpos(currPt,1:5),1];
                    else  %  If so assign to same stream
						imseg_num = lkup;
                    end                    
                   
					
				% See if winPt already has an imseg_num associated
				lkup =  imseg(imseg (:,2)==winPt,1);
                if (isempty(lkup))
					imseg  = [imseg ;imseg_num,winPt,tpos(winPt,1:5),1];                          
                end    
            end     % if d<srez
                               
        end         % for jj
    end             % for kk   
    waitbar(fnum/frame_length,hwait,['Processing frame ' int2str(fnum) ' of ' int2str(frame_length)]);
end                 % for fnum
close(hwait)



% Pass Two : Check for short invalid segment using <SEGLEN> criteria
% resolving conflicts of multiple points in the same segment and frame
% and reassigning seg numbers


% Fetch a list of segment numbers to reassign smallest possible integer.
% This is done for better quality of 3D movie.

hwait = waitbar(0,['Validating first pass streams for minimum length']);
segnum =1:length(unique(imseg(:,1)));
count =0;


for ii= 1:max(imseg(:,1))
    L = find(imseg(:,1)==ii);  
    if(length(L)<seglen)
      %  imseg(L,1:5)=kron([1 0 -1 -1 -1],ones(length(L),1));
        imseg(L,8) =-1;
    else
        count = count +1;
        imseg(L,1)=segnum(count);  % Assign smallest possible seg num
    end
end

 waitbar(1/2,hwait,['Checking for multiple detections of same stream']);
% Apply max criteria to resolve multiple detection in 
% any given region of interest <SREZ> and <FREZ> criteria
 
    for jj = min(imseg(:,7)):max(imseg(:,7))
        for ii= 1:max(imseg(:,1))
            fn = find(imseg(:,7)==jj & imseg(:,1)==ii);
            if(length(fn)>1)
                temp = imseg(fn,3:6);
                coord= temp(temp(:,1)==max(temp(:,1)),2:4);
                coord = coord(1,:);
                imseg(fn(1),4:6) = coord;
                % imseg(fn(2:end),3:5) = -2;    
                imseg(fn(2:end),8)= -2;  
            end
            clear fn;
        end
    end
     
    %  Check for weak streams that are close to stronger streams
    %  Assess average strength of each stream
    numstrs = max(imseg(:,1));  % number of detected streams
    strnth = zeros(1,numstrs);
    for ii= 1:numstrs
        astrms = find(imseg(:,1)==ii);  %  Find all stream frames
        strnth(ii) = sum(imseg(:,3));
    end
    frameclose = 5;
    %  Check relative strenghts for close streams
    [sst, sstind] = sort(strnth);
    for ii=1:numstrs
        strmindx = find(imseg(:,1)==sstind(ii));  % get rows associated with each stream
        lo1 = mean(
    
    
% Assign a pseudo point for frames with detection below secondary threshold
if nargin>4
    for jj = 1:lfram
        fn = imseg(imseg(:,7)==jj,7);
        if(isempty(fn))
            indx=   find(imseg(:,7)==jj-1,1,'last');
            imseg = [imseg(1:indx,:);[1 0 NaN NaN NaN jj 0];imseg(indx+1:end,:)];          
        end
    end
end

waitbar(1,hwait,['Done']);


imseg = imseg(2:end,:);

% Ensure the data is sorted with respect to frame number (default case)
imseg = sortrows(imseg,7);

close(hwait)
end
  