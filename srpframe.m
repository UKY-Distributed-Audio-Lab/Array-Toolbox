function [mleim, mleimn] = srpframe(sigar, gridax, mpos, fs, c, srchwin, powrez)
% This function processes signals from a microphone array relative to
% measurments in space containing the microphones and the sound sources.
% The mic positions field of view grid point with any uint as long as
% the units for parameters listed below are consistent (I like meters
% and seconds).  This version considers mic pair correlations and combines
% them separates (as opposed to simply computing the coherent power at the
% point in the spatial grid with all pairs).
% The output is the liklihood of the sound source as function of the
% spatial grid.  The number of axes in cell array GRIDAX determines the
% dimesion of the ouput liklihood grid.
%
%   [mleim] = mlescanprs(sigar, gridax, mpos, prs, fs, c, srchwin, powrez)
%
% Inputs
%  "sigar"        Matrix where each column is the rf signal segment from each mic.
%
%  "gridax"       Cell array of grid points for x-axis, y-axis and z-axis.
%                 The number of elements in this cell array can be less than
%                 3.  If so it will do a reduce dimension output space.
%
%  "mpos"         Matrix indicating the x,y, and z cartesion positions of each
%                 mic, where the first row contains the x positions, second y
%                 position, and third z positions.  If an axis is missing
%                 it will be assume to be the z and an array in the x-y
%                 plane will be assumed.  The Dimensions should be
%                 consisent with the axes in cell array "gridax" and number
%                 of mics must be the same as columns in "sigar."
%
%  "prs"          Matrix identifying all mic pairs for processing in each
%                 row.  Columns 1 and 2 contain the indices (corresponding
%                 to column indices of "mpos") for the particular mics in the
%                 pair.  An additional column can be present for a weighting
%                 of the pairs.  If not present, they all have the same
%                 weights
%           
%  "fs"           Sampling frequency of the rf mic signals in "sigar"
%
%  "c"            scalar is the speed of sound in the room.
%  "srchwin"       Search window length for a coherent sound.  This is 
%                 time in milliseconds from the most recent sample on the
%                 furthest mic from the point of interest. It is aligned
%                 such that the last point in the window corresponds to
%                 the last point (most recent time sample) of the non-delayed
%                 version of "sigar." Typically "srchwin" would be smaller than
%                 the length of "sigar" so the delays will result in
%                 measured data samples shifting into the listening window.
%                 If "sigar" is not long enough to have have samples to
%                 shift into the correlation window, zero padding will occur.
%
%  "powrez"       window over which power is computed.  Can be smaller than
%                 "srchwin"; however it cannot be larger
%
%                   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005

%  UNCOMMENT FOR TABLE FIR SHIFTING ....
nint = 8;  %  Number of subinterval points
ord = 4;
tab = subsamplefir(nint+2,ord,'c');
% END UNCOMMENT ....

dimfov = max(size(gridax));    %  Get dimensions for FOV 
dimx = length(gridax{1});
xax = gridax{1}(:);
prezx = (xax(2)-xax(1));
if dimfov >= 2  %  If second dimension present
    dimy = length(gridax{2});
    yax = gridax{2}(:);
    prezy = (yax(2)-yax(1));
else  % If not, add a singleton to hold place
    dimy = 1;
    yax = 0;
    prezy = 0;
end

if dimfov == 3  %  If third dimension present
    dimz = length(gridax{3});
    zax = gridax{3}(:);
    prezz = (zax(2)-zax(1));
else  % If not, add a singleton to hold place
    dimz = 1;
    zax = 0;
    prezz = 0;
end

if nargin < 7
%  Compute symetric distance for Response power estimate
powrez = 4*ceil(fs*sqrt(prezx^2+prezy^2+prezz^2)/(c));
end
% Obtain mic array information
[mr, mc] = size(mpos);  % Determine # number of mics = mc 
% Get signal dimensionality
[sr,sc] = size(sigar);  % Determine the mic signal array size [# of samples, # of mics]
%  Test for consistency of mic and signal array
if mc ~= sc
    error('Number of mic signal not consistent with number of mics')
end

% extend mic coordinates to 3 dimension with zeros, if less than 3
if mr == 1
   mpos(2:3,:) = zeros(2,mc);
elseif mr ==2
   mpos(3,:) = zeros(1,mc);
end
% Find number of samples in listening window
winsamp = round(srchwin*fs) + 1;

temp = zeros(winsamp,sc);  %  Initialze correlation signal matrix
% If listening window longer than input segment, reducing listening window
% to signal length.
if winsamp > sr
    winsamp = sr;
end

%  Loop through every point in FOV
for kz=1:dimz    
    for ky=1:dimy
        for kx = 1:dimx
            %  Distance of position from microphones
            ds = mpos - [xax(kx); yax(ky); zax(kz)]*ones(1,mc);  % Spatial
            rst = sqrt(ds(1,:).^2 + ds(2,:).^2 + ds(3,:).^2) /c;  %  Corresponding Time differences
   %         [mt, nmic] = max(rst);  % Find furthest mic to point being considered
            %  Make all delays relative to furthest mic, setting furthest mic to a 0 delay
            %rstref = mt-rst;
            %  Create shading values to weight mic inputs as function of
            %  distance giving closest mic the most weight.
            at = arweights(rst);
            md = max(rst);
            rstref = md - rst;

%  ARRANGE COMMENTS FOR DESIRED DELAY PROGAM
%            sd = delayt(sigar, fs, rstref, sr/fs);    %  Time domain FIR shifting
%            sd = delayf(sigar,fs,rstref, sr/fs);      %  Frequency Domain shifting
%            sd = delayint(sigar,fs,rstref,sr/fs);      %  Time domain sample interval shifting 
            sd = delaytab(sigar,tab,fs,rstref, sr/fs);  %  Time domain with table, uncomment first 3 lines
%  END ARRAGE COMMENTS
            %  Load up matrix with aligned mic signal corresponding to mic
            %  position and apply shading weights
            for ksig = 1:sc
                 temp(:,ksig) = at(ksig)*sd(sr-winsamp+1:sr,ksig); 
            end

%               figure(1)
%               imagesc(temp); colorbar
%               pause
%          Sum coherent power in array
           dumval = sum(temp,2).^2 - (sum(temp.^2,2));
           %avg = ones(1:(2*powrez)/(2*powrez),1);
           %dumval = conv2(dumval,avg,'same');
           %  Find max point of coherence as most like sound source position 
%            plot(dumval)
%            pause
           [ma, mp] = max((dumval));
           [man, mpn] = min((dumval));
           %nv = find(dumval < 0);
           %sdt = sqrt(mean(dumval(nv).^2))
           %plot([0:length(dumval)-1], sum(temp,2), 'k',[0:length(dumval)-1], (sum(temp.^2,2)),'r')
           %pause
           %  Get range of delay and sum values in the neighborhood 
           %  of max point
           rngp = [max([mp(1)-powrez,1]):min([mp(1)+powrez,winsamp])];
           rngn = [max([mpn(1)-powrez,1]):min([mpn(1)+powrez,winsamp])];
           %  Estimate power in this neighborhood
           mleim(ky,kx,kz) = (mean(dumval(rngp))); % - sum(sum(temp(rngp,:).^2));
           mleimn(ky,kx,kz) = (mean(dumval(rngn))); % %- sum(sum(temp(rngn,:).^2));
           %mleim(ky,kx,kz) = ma(1);
        end
    end
end
mleim = squeeze(mleim);  %  Remove singleton dimensions