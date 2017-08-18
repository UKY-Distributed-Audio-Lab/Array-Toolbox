function clipdata = cclip(xin,cliptype,clpercent)
%  This function applies various types of center clipping to
%  the input signal (XIN) based on a percentage of its peak values.
%
%     clipdata = cclip(xin,cliptype,clpercent)
%
%  inputs:
%  XIN       => Vector of input signal values
%  CLIPTYPE  => String indicating type of clipping to apply
%          'cnsc' = Soft Center Clip:
%          C[x(n)]=x(n)-Cl for (x(n)>=Cl), 0 for (|x(n)|<Cl),  and x(n)+Cl for (x(n)<=-Cl)
%          'cnhc' = Hard Center Clip:
%          C[x(n)]=x(n) for (|x(n)|>=Cl), and 0 for (|x(n)|<Cl)
%          'cnhl' = Center Hard Limiter Clip:
%          C[x(n)]= 1 for (|x(n)|>=Cl), 0 for (|x(n)|<Cl), and -1 for (x(n)<=-Cl)
%  CLPERCENT => percentage of maximum absolute signal value to obtain the
%               clipping level.  The percentage value must be given between 0 < cl < 1
%
%   Written by Tim Black, updated by Kevin D. Donohue (donohue@engr.uky.edu) Jan 2006



%  Check for validity of parameters
if clpercent > 1 | clpercent < 0
    error('Precent of signal for clip threshold should be between 0 and 1')
end
if length(cliptype) ~= 4
    error('Clip Type not properly specified!  Use strings cnsc, cnhc, cnhl, for soft, hard, or hard limit center clip.')
end

% The clipping level, CL, is set to the percentage of the maximum magnitude of the
% positive and negative peaks (effective Magnitude)
effmag= max(abs(xin));  %  Effective Magnitude
cl=clpercent*effmag(1);                      %  Clipping level

% Center clip XIN (the current window of xin) according to argument cliptype: 
clipdata=zeros(1,length(xin));  %  Initialize output to zeros
ind1=find(xin>=cl);     %  find all positive points greater than clip threshold 
ind2=find(xin<=-cl);    %  find all negative points with magnitudes greater than clip threshold

if sum(cliptype=='cnsc') == 3      %  Center soft clip
   clipdata(ind1) = xin(ind1)-cl;   %  Reduce amplitude by clip threshold (for positive values)
   clipdata(ind2) = xin(ind2)+cl;   %  Reduce amplitude by clip threshold (for negative values)
elseif sum(cliptype=='cnhc') == 3  %  Center hard clip
   clipdata(ind1) = xin(ind1);  % Pass positive values greater than threshold 
   clipdata(ind2) = xin(ind2); % Pass negative values with magnitudes greater than threshold
elseif sum(cliptype == 'cnhl') == 3   %  Hard limit clip
   clipdata(ind1) = 1;  %  Set postive values to 1;
   clipdata(ind2) = -1; %  Set negative values to -1;
else
    error('Clip Type not specified!  Use strings cnsc, cnhc, cnhl, for soft, hard, or hard limit center clip.')
end
