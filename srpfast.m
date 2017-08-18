function mleim = srpfast(sigar, gridax, mpos, fs, c, srchwin, rstref, at)
%
% This function processes signals from a microphone array and
% creates an acoustic image according to a Steered Response
% Coherent Power (SRCP) algorithm for Sound Source Location (SSL).
% The mic positions and field of view (FOV) grid points can
% be in any units as long they are consistent for parameters
% listed below (I like meters and seconds).  A Delay and Sum
% Beamformer is applied for each pixel in the FOV and the coherent
% power is computed for the sound emanating from the pixel location.
% Coherent power differs from regular power (squared sum normalized
% by time) in that only the power from signal cross correlation
% between 2 different mics is used in the summation.  Therefore, coherent
% power can be negative.  Large positive values indicate a strong
% correlation or coincidence between mic signals.  A large negative
% value represents a strong 180 degrees out of phase coincidence (very
% rare that this number can get large for more than 2 mics). 
% The output is the coherent power at each pixel and represents a
% likelihood of the sound source as a function of the spatial grid.  The
% number of axes in cell array GRIDAX determines the dimension of the
% ouput likelihood grid.  This particular version uses only integer delays
% for the Delay and Sum Beamformer. 
%
%   mleim = srpfast(sigar, gridax, mpos, fs, c, srchwin, rstref, at)
%
% Inputs
%  "sigar"        Matrix where each column is the audio signal segment from
%                 each microphone.
%
%  "gridax"       Cell array of grid points for x-axis, y-axis and z-axis.
%                 The number of elements in this cell array can be less than
%                 3.  If so, dimension of output space will be reduced.
%
%  "mpos"         Matrix indicating the x, y, and z cartesian positions of each
%                 mic, where the first row contains the x positions, second y
%                 position, and third z positions.  If an axis is missing,
%                 only the x-y plane will be used.  The dimensions should be
%                 consisent with the axes in cell array "gridax" and number of
%                 mics must be the same as number of columns in "sigar".
%           
%  "fs"           Sampling frequency of the audio mic signals in "sigar".
%
%  "c"            The speed of sound in the room. (scalar value)
%
%  "srchwin"      The window length over which coherent power is computed.
%                 Also represents the time resolution of acoustic image.
%                 Typically, the signals in "sigar" are longer than "srchwin"
%                 so that as delays are made to each signal, actual data will be
%                 shifted into the window over which the delay and sum operation
%                 is performed.
%                 In this implementation, the signal from the furthest mic is not
%                 delayed. All other signals are delayed to sync up with this
%                 max delayed signal.
%                 If "sigar" is not long enough to have have samples to
%                 shift into the correlation window, zero padding will
%                 occur.
%   "rstref"      Delays with respect to the furthest mic. Computed in 
%                 VOLDELWTS.M
%   "at"          Weights for each mic as a function of distance. Weights 
%                 are computed in VOLDELWTS.M 
%
% Required Functions:
%           delayintfast.m
%
%    Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005
%    Modified by Kirstin Brangers                       July 2012


dimfov = max(size(gridax));     % Get dimensions for FOV 

% X Dimension
dimx = length(gridax{1});
if dimx > 1
   xax = gridax{1}(:);          % Access all elements in the first cell of gridax
else
   xax = gridax{1}(1);          % Access the first element in the first cell of gridax
end
% Y Dimension
if dimfov >= 2                  % If second dimension present
    dimy = length(gridax{2});
    if dimy > 1
       yax = gridax{2}(:);      % Access all elements in the second cell of gridax
    else
       yax = gridax{2}(1);      % Access the first element in the second cell of gridax
    end
else                            % If not, add a singleton to hold place
    dimy = 1;
    yax = 0;
end
% Z Dimension
if dimfov == 3                  % If third dimension present
    dimz = length(gridax{3});
    if dimz > 1;
        zax = gridax{3}(:);     % Access all elements in the third cell of gridax
    else
        zax = gridax{3}(1);     % Access the first element in the third cell of gridax
    end
else                            % If not, add a singleton to hold place
    dimz = 1;
    zax = 0;
end

% Obtain mic array information
[mr,mc] = size(mpos);           % Determine number of mics = mc 

% Get signal dimensionality
[sr,sc] = size(sigar);          % Determine the mic signal array size [# of samples, # of mics]

%  Test for consistency of mic and signal array
if mc ~= sc
    error('Number of mic signals is not consistent with number of mics.')
    mleim = [];
    return
end

% Extend mic coordinates to 3 dimensions with zeros if dimension is less than 3
if mr == 1
   mpos(2:3,:) = zeros(2,mc);
elseif mr == 2
   mpos(3,:) = zeros(1,mc);
end

% Find number of samples in listening window
winsamp = round(srchwin*fs) + 1;
temp = zeros(winsamp,sc);       %  Initialze correlation signal matrix

% If listening window is longer than input segment, reduce listening window
% to signal length.
if winsamp > sr
    winsamp = sr;
end

%  Initalize acoustic image array
mleim = zeros(dimx, dimy, dimz);
lenrat = sr/fs;                 % Length of extended vector in seconds

%  Loop through every point in FOV
for kz=1:dimz                   % Z-Dimension Loop
    for ky=1:dimy               % Y-Dimension Loop
        for kx = 1:dimx         % X-Dimension Loop
            %  Distance of FOV position from all microphones
            sd = delayintfast(sigar,fs,squeeze(rstref(kx,ky,kz,:))',lenrat); %  Time domain sample interval shifting             
            %  Load up matrix with aligned mic signal corresponding to 
            %  mic position and apply shading weights
            for ksig = 1:sc
                 temp(:,ksig) = at(kx,ky,kz,ksig)*sd(sr-winsamp+1:sr,ksig); 
            end            
           % Beamform: Sum coherent power in array
           dumval = sum(temp,2).^2 - (sum(temp.^2,2));
           mleim(kx, ky, kz) = mean(dumval); 
           
        end                     % End X-Dimension Loop
    end                         % End Y-Dimension Loop
end                             % End Z-Dimension Loop 
mleim = squeeze(mleim)*fs;      % Remove singleton dimensions


