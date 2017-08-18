function tposout = convtpos(tpos, fovp)
%  This function converts ASA files to a frame based format  
%  Function is of the form:
%
%      tposout = convtpos(tpos, fovp)
%
% Input:
%   tpos = [time stamp, detection statistic, x, y, z]
%       with time stamp in seconds, x,y,z coordinates in meters
%       and detection statistic is excess threshold (0 implies equal to thrshold)
%   fovp is a data structure with file trez.  This is the time window
%       used in the detection program.  Half of this is the time increment
%       between window (and time stamps of consecutive detections).
% 
% Output:
%  tposout = [detection statistic, x, y, z, frame number]
%  
%  Where segment number detection statistic and x,y,z is a repeat of  TPOS
%  and frame number is an integer indicating the sequential time stamps.  
%  For missing time stamps (i.e. no detected targets) a line of zeros is
%  inserted for the first 4 entrees.  So on the output there is at least 
%  one line per frame, sometimes multiple detections in one frame (i.e.
%  multiple spatial locations at the same time).
%
% Written by Kevin D. Donohue (donohue@engr.uky.edu) June 2012
%



[r,c] = size(tpos);
if ~isfield(fovp,'tinc')
    fovp.tinc = fovp.trez/2;
end

frstpnt = min(tpos(:,1));
lstpnt = max(tpos(:,1));

frmlimit = floor((lstpnt(1)-frstpnt(1))/fovp.tinc);
frmnum = 1;
outcnt = 1;
tposindx = round((tpos(:,1)-frstpnt(1))/fovp.tinc)+1;
while frmnum <= frmlimit
    indxid = find(frmnum == tposindx);
    if isempty(indxid)
        tposout(outcnt,:) = [zeros(1,4), frmnum];
        outcnt=outcnt+1;
    else
        for k=1:length(indxid)
        tposout(outcnt,:) = [tpos(indxid(k),2:5), frmnum];
        outcnt = outcnt+1;
        end
    end
    frmnum = frmnum+1;  %  Increment frame number
end
