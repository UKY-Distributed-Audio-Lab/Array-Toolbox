function mleim = srpframenn(sigar, gridax, mpos, fs, c, srchwin)
% This function processes signals from a microphone array and
% creates an acoustic image according to a Steered Response
% Coherent Power (SRCP) algorithm for sound source location.
% The mic positions and field of view (FOV) grid points can
% be in any units as long as they are consistent with parameters
% listed below (I like meters and seconds).  A delay and sum
% beamformer is applied for each pixel in the FOV and the coherent
% power is computed for the sound eminating from the pixel location.
% Coherent power differs from regular power (squared sum normalized
% by time) in that only the power from signal cross correlation
% between 2 different mics is used in the summation.  So coherent
% power can be negative.  Large positive values indicate a strong
% correlation or coincidence between mic signals.  A large negative
% value represents a strong 180 degrees out of phase coincidence (very
% rare that this number can get large for more than 2 mics). 
% The output is the coherent power at each pixel and represents a
% liklihood of the sound source as a function of the spatial grid.  The
% number of axes in cell array GRIDAX determines the dimension of the
% ouput liklihood grid.  This particular version uses only integer delays
% for the delay and sum beamformer.
%
%   mleim = srpframenn(sigar, gridax, mpos, fs, c, srchwin)
%
% Inputs
%  "sigar"        Matrix where each column is the audio signal segment from
%                 each microphone.
%
%  "gridax"       Cell array of grid points for x-axis, y-axis and z-axis.
%                 The number of elements in this cell array can be less than
%                 3.  If so, it will do a reduce dimension output space.
%
%  "mpos"         Matrix indicating the x,y, and z cartesion positions of each
%                 mic, where the first row contains the x positions, second y
%                 positions, and third z positions.  If an axis is missing,
%                 only the x-y plane will be used.  The Dimensions should be
%                 consisent with the axes in cell array "gridax" and number of
%                 mics must be the same as columns in "sigar."
%           
%  "fs"           Sampling frequency of the audio mic signals in "sigar"
%
%  "c"            scalar is the speed of sound in the room.
%
%  "srchwin"      The window length over which coherent power is computed
%                 and represents the time resolution of acoustic image.
%                 Typically, the signals in "sigar" are longer than "searchwin"
%                 so that as delays are made to each signal, actual data will be
%                 shifted into the window over which the delay and sum operation
%                 is performed.
%                 In this implementation the signal from the furthest mic is not
%                 delayed. All other signals are delayed to sync up with the
%                 max delayed signal.
%                 If "sigar" is not long enough to have samples to
%                 shift into the correlation window, zero padding will
%                 occur.
%
%    This function requires the functions: ARWEIGHTS and DELAYINT from the array
%    toolbox to run.
%
%                   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005


dimfov = max(size(gridax));    %  Get dimensions for FOV 
%  x-dimension
dimx = length(gridax{1});
if dimx > 1
   xax = gridax{1}(:);
else
   xax = gridax{1}(1);
end
% y-dimension
if dimfov >= 2  %  If second dimension present
    dimy = length(gridax{2});
    if dimy > 1
       yax = gridax{2}(:);
    else
        yax = gridax{2}(:);
    end
else  % If not, add a singleton to hold place
    dimy = 1;
    yax = 0;
end
%  Z-axis
if dimfov == 3  %  If third dimension present
    dimz = length(gridax{3});
    if dimz > 1;
        zax = gridax{3}(:);
    else
        zax = gridax{3}(:);
    end
else  % If not, add a singleton to hold place
    dimz = 1;
    zax = 0;
end

% Obtain mic array information
[mr, mc] = size(mpos);  % Determine # number of mics = mc 
% Get signal dimensionality
[sr,sc] = size(sigar);  % Determine the mic signal array size [# of samples, # of mics]
%  Test for consistency of mic and signal array
if mc ~= sc
    error('Number of mic signal not consistent with number of mics')
    mleim = [];
    return
end

% extend mic coordinates to 3 dimension with zeros, if dimension less than 3
if mr == 1
   mpos(2:3,:) = zeros(2,mc);
elseif mr == 2
   mpos(3,:) = zeros(1,mc);
end
% Find number of samples in listening window
winsamp = round(srchwin*fs) + 1;

temp = zeros(winsamp,sc);  %  Initialze correlation signal matrix
% If listening window longer than input segment, reduce listening window
% to signal length.
if winsamp > sr
    winsamp = sr;
end
%  Initialize acoustic image array
mleim = zeros(dimy, dimx, dimz);
%  Loop through every point in FOV
for kz=1:dimz    %  Z-Dimension Loop
    for ky=1:dimy  %  Y-Dimension Loop
        for kx = 1:dimx  % X-Dimension Loop
            %  Distance of FOV position from all microphones
            ds = mpos - [xax(kx); yax(ky); zax(kz)]*ones(1,mc);
            %  Convert distances to time
            rst = sqrt(ds(1,:).^2 + ds(2,:).^2 + ds(3,:).^2) /c;
            %  Create shading values to weight mic inputs as function of
            %  distance, giving closest mic the most weight.
            at = arweights(rst);
            % Find mic with maximum delay to FOV
            md = max(rst);
            % Implement all other delays with respect to the furthest mic.
            rstref = md - rst;
            sd = delayint(sigar,fs,rstref,sr/fs); %  Time domain sample interval shifting 
            %  Load up matrix with aligned mic signal corresponding to mic position and apply
            %  shading weights
            for ksig = 1:sc
                 temp(:,ksig) = at(ksig)*sd(sr-winsamp+1:sr,ksig); 
            end
%          Beamform: Sum coherent power in array
           dumval = sum(temp,2).^2 - (sum(temp.^2,2));
           mleim(ky, kx, kz) = mean(dumval)*fs; 
           
        end  %  End X-Dimension Loop
    end  %  End Y-Dimension Loop
end    %  End Z-Dimension Loop 
mleim = squeeze(mleim);  %  Remove singleton dimensions