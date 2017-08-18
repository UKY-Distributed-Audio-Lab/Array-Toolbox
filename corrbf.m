%% Correlative GSC Beamformer
% Implementation of a correlative-shift beamformer.  This beamformer
% has a similar structure to the traditional GSC and uses the same
% delay-and-sum beamformer for the Fixed Beamformer and adaptive
% Multiple-Input Canceller.  However, pairs of tracks in the Blocking
% Matrix are shifted to yield the greatest cross correlation over a
% window and are subtracted using an alpha parameter calculated
% from expected values.
%
% *Syntax*
%
%  |[y, b, z, tauFinal, newTrackOrder, mcWall, alpha] 
%     = corrbf(x, fs, s, m, c, mu, order, beta, maxShift, 
% 	     trackOrder, tauInit, xPrev, bPrev, zPrev, mcWinit, 
% 	     mcWForce, forceBmOrder, forceBmLags, forceAlpha)|
%
% *Required Inputs*
%
% * |x| - Matrix of raw audio data, each column a mic array track
% * |fs| - Audio sampling rate (Hertz)
% * |s| - 3x1 source location point (meters)
% * |m| - column matrix of 3D mic positions (meters)
% * |c| - speed of sound (meters per second)
% * |mu| - step size of the LMS filters in [0, 2]
% * |order| - order of the LMS filters
% * |beta| - forgetting factor of the LMS filters in [0, 1]
% * |maxShift| - Size of search window during cross correlation in
%      samples 
% * |trackOrder| - Vector of indices for arrangement of the columns
%      of x from closest to furthest from the speaker source.
%
% *Optional Inputs*
%
% * |tauInit| - Vector of lags to apply to the microphone tracks
%      for target signal alignment.  Leave empty to force
%      calculation of taus in the Fixed Beamformer using the source
%      and mic positions.
% * |xPrev| - Matrix of audio data from a previous window.  This
%      data allows for real data to be padded at the front of each
%      shifted track in the DSB algorithm.  Simply supply zeros if
%      no previous data exists.
% * |bPrev| - Vector of DSB output from a previous audio window.
%      This end of this data is used as initial input to the BM LMS
%      filters to allow seamless transitions when beamforming over
%      consecutive windows of data.  Supply zeros if no previous
%      data is available.
% * |zPrev| - Matrix of BM output for a previous audio window.  Used
%      as initial input to the MC LMS filters.  Supply zeros if no
%      previous data is available.
% * |mcWinit| - Matrix of initial tap weights for the MC LMS
%      filters.  Set to [] if not needed.
% * |mcWForce| - Cubic matrix of tap weights to force the MC LMS
%     filters to use.  Useful for recreating a previous run where taps
%     have already been determined.
% * |forceBmOrder| - Force the blocking matrix track order to these
%      values
% * |forceBmLags| - Force the blocking matrix to use these lags.
% * |forceAlpha| - Force the alphas (difference parameters in the
%      blocking matrix) to these values.  Good for recreating a run.
%
% *Outputs*
%
% * |y| - Beamformed output track
% * |b| - Vector of DSB output
% * |z| - Matrix of output tracks from the blocking matrix.  These
%      tracks should hopefully approximate the noise in the simulation
%      as closely as possible.  This is supplied mostly as a debugging
%      output and can be checked to see how much target signal leakage
%      occurs in the beamformer.
% * |tauFinal| - Vector of lags for target signal alignment
%      calculated from the best cross correlation between adjacent
%      microphone tracks.
% * |newTrackOrder| - Once new taus are calculated using cross
%      correlation, this vector contains a possibly new ordering of
%      indices for the microphone tracks from closest to furthest
%      from the speaker source.
% * |mcWall| - Cubic matrix of taps saved from the MC.  Each matrix
%      (1st and 2nd dims) of the cube will be |order| tall and |M|
%      wide such that each column is the taps for an adaptive filter
%      for an individual channel.  The 3rd dimension specifies the
%      iteration of the taps.  This matrix can be fed into another run
%      as |mcWForce|.
% * |alpha| - Vector of alpha's calculated in the BM
%
% Written by Phil Townsend (jptown0@engr.uky.edu) 7-23-08

%% Function Declaration
function [y, b, z, tauFinal, newTrackOrder, mcWall, alpha] ...
    = corrbf(x, fs, s, m, c, mu, order, beta, maxShift, ...
	     trackOrder, tauInit, xPrev, bPrev, zPrev, mcWinit, ...
	     mcWForce, forceBmOrder, forceBmLags, forceAlpha,iw) 

%% Argument Error Checking
if nargin == 19
    iw = 0;
end
if isempty(x) || ~isreal(x) || ~all(all(isfinite(x)))
    error('x must be a real matrix')
elseif isempty(fs) || ~isreal(fs) || ~isscalar(fs) || ~isfinite(fs) ...
	|| fs <= 0 
    error('fs must be positive real scalar');
elseif ~isreal(s) || ~all(isfinite(s)) || size(s,1) ~= 3 || ...
        size(s,2) ~= 1  || length(size(s)) ~= 2
    error('s must be a real 3x1 vector');
elseif ~isreal(m) || ~all(all(isfinite(m))) || size(m,1) ~= 3
    error('m must be real with 3 rows')
elseif ~isreal(c) || ~isscalar(c) || ~isfinite(c) || c <= 0
    error('c must be positive real scalar');
elseif isempty(mu) || ~isscalar(mu) || ~isreal(mu) || ...
        ~isfinite(mu)
    error('mu must be a real scalar');
elseif isempty(order) || ~isscalar(order) || ~isreal(order) || ...
        ~isfinite(order) || order < 0 || ...
        abs(mod(order,floor(order))) > 0
    error('order must be a positive integer');
elseif isempty(beta) || ~isscalar(beta) || ~isreal(beta) || ...
        ~isfinite(beta) || beta < 0 || beta > 1
    error('beta must be a real scalar in [0,1]');
elseif isempty(xPrev) || ~isreal(xPrev) || ~all(all(isfinite(xPrev)))
    error('xPrev must be a real matrix')
elseif ~isvector(bPrev) || ~isreal(bPrev) || ~all(isfinite(bPrev))
    error('bPrev must be a real vector');
elseif isempty(zPrev) || ~isreal(zPrev) || ~all(all(isfinite(zPrev)))
    error('zPrev must be a real matrix')
elseif ~isempty(mcWinit) && ...
	(~isreal(mcWinit) || ~all(all(isfinite(mcWinit))))
    error('mcWinit must be a real matrix')
end

%% Setup
[N, M] = size(x);  % number of samples in track, mics in array
alpha = zeros(M-1,1);  % optimal constant track diffs


%% Beamform
% The blocking matrix of this beamformer uses a special scaling to
% attempt to remove the correlation between tracks, rather than a
% simple difference in the traditional GSC.  That is, rather than
%
% $$ z_k(n) = x_k(n) \ ^\_ \ x_{k+1}(n) $$
%
% we use
%
% $$ z_k(n) = x_k(n) \ ^\_ \  \alpha x_{k+1}(n) $$
%
% $$ \alpha = \frac{ E[x_k(n)]E[x_{k+1}(n)] }{ E[x_{k+1}^2(n)] } $$
%
% In addition, our delays for target alignment in the blocking matrix
% are calculated by looking at the cross correlation between
% microphone tracks in a small window around the calculated/supplied
% lags for alignment.  Our initial delays are taken from a traditional
% delay-sum beamform, but thereafter are updated if the max cross
% correlation exceeds a threshold P (currently chosen at .8)
%
% $$  \hat{\tau}_k = \textrm{argmax}_n
%  \left( \sum_{\xi = \tau_k^{(i)}-\Delta\tau}
%  ^{\xi = \tau_k^{(i)}+\Delta\tau}
%  \beta^{(i)}_k(\xi)v_k^{(i)}(n+\xi) \right) $$
%
% $$ \tau_k^{(i+1)} = \{
%   \begin{array}{ c c }
%      \hat{\tau}_k & \rho_{max} > P \\
%      \tau_k^{(i)} & \textrm{otherwise}
%   \end{array} 
% $$
%
% For instance, if supplied/computed lag for DSB alignment between
% tracks 1 and 2 were 15 at a sample rate of 8 kHz we would have a
% range of 6 and thus look at the cross correlation for lags 9 through
% 21 and use the peak.
%

% Enforce track order.  Previous and current windows must match
xPrev = xPrev(:, trackOrder);  x = x(:, trackOrder);

% Fixed Beamformer (FBF)
if isempty(tauInit)  % calculate delays from distances
    d = sqrt(sum((m - s*ones(1, length(m(1,:)))).^2, 1))'; % distances
    t = (max(d)-d)/c;  % Relative TDOA
    tau = round(t*fs);  % convert to samples
else
    tau = tauInit;
end

if iw == 1
    wghts = arweights(tau*c);  % inverse distance weighting
else
    wghts = ones(size(tau));  % no weighting
end


xShift = zeros(size(x));  % Initialize shifted matrix
for k=1:M  % align signals, pulling real data in the rear
    xShift(:,k) = wghts(k)*[xPrev(end-tau(k)+1:end, k); ...
		   x(1:end-tau(k), k)];
end
b = sum(xShift, 2);  % add data across columns


% Correlative Blocking Matrix (CBM)
% 1. Calculate correlation shifts and rearrange track order to
%    ensure tracks are ordered from closest to furthest
if isempty(forceBmLags)  % calculate only if not given lags
    trlags = diff(tau);  % convert absolute lags to track relative
    for k = 1:M-1  % iterate over track pairs
	trlags(k) = findoptlag(x(:,k), x(:,k+1), trlags(k), ...
			       maxShift);
    end
    % Convert back to individual lags for each track
    tau = [0; cumsum(trlags)];
    tau = tau + abs(min(tau));  % Assure positive since relative
    
    % Now ensure lags are sorted and enforce the new track order
    [tau, newTrackOrder] = sort(tau); % get sorted vector, indices
    xPrev = xPrev(:, newTrackOrder);  x = x(:, newTrackOrder);
else % use supplied, forced lags and order
    tau = forceBmLags;
    xPrev = xPrev(:, forceBmOrder);  x = x(:, forceBmOrder);
end

% Apply shift, supplying real data at front.
for k = 1:M  % iterate over all tracks
    xShift(:,k) = wghts(k)*[xPrev(end-tau(k)+1:end, k); ...
		   x(1:end-tau(k), k)];
end
tauFinal = tau; % save for output

% 2. Calculate alpha, the scaling parameter that minimizes the
%    difference between tracks for pairwise subtraction
if isempty(forceAlpha)  % not forcing alpha values, so calculate
    for k = 1:M-1  % iterate over all track pairs
	alpha(k) = mean(xShift(:,k).*xShift(:,k+1)) / ...
	    mean(xShift(:,k+1).^2);  % E[.] -> mean(.)
    end
else  % use supplied values for alpha
    alpha = forceAlpha;
end

% 3. Pairwise normalization and subtraction for the blocking
%    matrix.  Apply alpha to each column of xShift
z = xShift(:,1:end-1) - xShift(:,2:end)*diag(alpha);


% Multiple Input Canceller (MC)
% The MC here works the same as in the traditional GSC, so
% simply call the function handling that operation for the GJBF.
mcAdapt = [];  % MC LMS filters always on
K = [];  % No MC norm constraint here

% Call external MC function
[y, mcWall] = mclms(z, b, mu, order, beta, zPrev, mcWinit, ...
                    mcWForce, mcAdapt, K);

end % function corrbf




%% Optimal Cross Correlation Lag Solver
% This private function for the Correlative GSC computes the cross
% correlation between vectors |a| and |b| centered at lag |startInd|
% and extending |maxShift| lags around |startInd|.  We search for
% the largest value of this cross correlation that's a local
% maximum and select its corresponding lag if it exceeds a
% hardcoded threshold (currently .8); otherwise the original
% startInd is returned.  This thresholding is intended to ensure
% that we change the lags between two tracks only if we're
% reasonably sure a correlation exists--in the beamforming sense
% this should mean changing lags only for the target signal and
% interference sources.
%
% *Syntax*
%
%  |[optlag, alpha] = findoptlag(a,b,startInd,maxShift)|
%
% *Inputs*
%
% * |a| - Vector to be delayed
% * |b| - Reference vector for cross correlation
% * |startInd| - Offset from lag zero to center the search window
% * |maxShift| - Maximum number of samples to the left and right of
%      startInd to look for the greatest lag
%
% *Outputs*
%
% * |optlag| - Number of delays applied to a; a positive value implies
%      a forward shift while a negative value implies a backward
%      shift.
% * |alpha| - Value of the maximum cross correlation, which can be
%      used as a normalization factor to minimize the difference
%      between the two signals (also called the Wiener coefficient).
%      (Note that the value derived here isn't working so we don't use
%      it above)
%
% Written by Phil Townsend (jptown0@engr.uky.edu)  3/19/08

%% Function Declaration
function [optlag, alpha] = findoptlag(a, b, startInd, maxShift) 

%% Argument Error Checking
if ~isreal(a) || ~isvector(a) || size(a,2) ~= 1
    error('a must be a real column vector');
elseif ~isreal(b) || ~isvector(b) || size(a,2) ~= 1
    error('b must be a real column vector');
elseif ~isreal(maxShift) || ~isscalar(maxShift) || maxShift <= 0 || ...
        mod(maxShift,1) ~= 0
    error('maxShift must be a positve integer');
elseif ~isreal(startInd) || ~isscalar(startInd) || mod(startInd,1) ~= 0
    error('startInd must be an integer');
end

%% Cross Correlation and Analysis

% Perform cross correlation, taking only as many points as necessary
[crosscor, lags] = xcorr(a, b, abs(startInd)+maxShift, 'coeff');

% Now take a window of just allowable points in the crosscor vector.
% Handle the bounds as robustly as possible
midind = find(lags==0);
ubound = midind + startInd + maxShift;
if ubound > length(crosscor);  ubound = length(crosscor);  end
lbound = midind + startInd - maxShift;
if lbound < midind, lbound = midind; end
ccwin = crosscor(lbound:ubound);

% Get the index of the max inside the window, then convert that
% to an index for the original crosscor vector
alpha = max(ccwin);

% Change the lags ONLY if this max correlation coefficient is above a
% threshold
if alpha > .8
    ccwinmaxind = find(ccwin==alpha, 1);
    optlag = ccwinmaxind+lbound-midind-1;  % points to add
else
    optlag = startInd;
end

end  % function findoptlag
