%% Griffiths-Jim Beamformer
% Beamforms on a target sound source using a delay and sum beamformer
% to enhance the target and a blocking matrix and multiple input
% canceller to suppress noise.  Also known as the Generalized Sidelobe
% Canceller (GSC).
%
% *Syntax*
%
%  |[y, bmWall, mcWall, snrAll, b, z] = 
%     gjbf(x, fs, s, m, c, p, q, mu, order, beta, phi, psi, K, 
% 	 xPrev, bPrev, zPrev, bmWinit, mcWinit, snrThresh, 
% 	 snrRate, snrInit, bmWForce, mcWForce,iw)|
%
% *Required Inputs*
%
% * |x| - Matrix of raw audio data, each column a mic array track
% * |fs| - Audio sampling rate (Hertz)
% * |s| - 3x1 source location point (meters)
% * |m| - column matrix of 3D mic positions (meters)
% * |c| - speed of sound (meters per second)
% * |p| - Constant delay placed on all tracks entering the BM
% * |q| - Constant delay placed on DSB output just before final summer
% * |mu| - step size of the LMS filters in [0, 2]
% * |order| - order of the LMS filters
% * |beta| - forgetting factor of the LMS filters in [0, 1]
%
% *Optional Inputs*
%
% Note that these modifications transform the traditional
% Griffiths-Jim Beamformer into an embellished structure put forth by
% Hoshuyama and Sugiyama (see references).
%
% * |phi| and |psi| - Upper and lower bounds for a CCAF blocking
%     matrix.  Leave blank for the traditional Griffiths-Jim blocking
%     matrix where pairs of adjacent microphone tracks are subtracted
%     from each other.  Must specify both or neither.
% * |K| - Apply a norm constraint to the multiple-input canceller (MC)
%     LMS filters.  Leave blank for no constaint.
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
% * |bmWinit| - Matrix of initial tap weights for the BM LMS
%      filters.  Set to [] if not needed.
% * |mcWinit| - Matrix of initial tap weights for the MC LMS
%      filters.  Set to [] if not needed.
% * |snrThresh| and |snrRate| - Parameters for selective adaptation.
%     |snrRate| times the SNR for each track is estimated as the
%     10*log10 ratio of the power of the fixed beamformer to each of
%     the blocking matrix output tracks.  If these estimated SNR's
%     exceed |snrThresh| then only the BM LMS coefficients are allowed
%     to adapt while the taps of the MC are kept constant; otherwise
%     the BM taps are locked and the MC coefficients adapt.  Must
%     specify both |snrThresh| and |snrRate| or neither.
% * |snrInit| - Assume this vector of SNR values at the start
% * |bmWForce| and |mcWForce| - Cubic matrices of tap weights to force
%     the BM and MC LMS filters to use.  Useful for recreating a
%     previous GJBF run where taps have already been determined.
% * |iw| - (optional) if set to 1, the dsb will apply an inverse distance 
%          weighting on the microphones, useful for near-field/immersive
%          applications.  Any other value or its absense will result in
%          the uniform weight dbs.
%
% *Outputs*
%
% * |y| - Beamformed output audio track.
% * |bmWall| - Cubic matrix of taps saved from the blocking matrix.
%      Each matrix (1st and 2nd dims) of the cube will be |order| tall
%      and |M| wide such that each column is the taps for an adaptive
%      filter for an individual channel.  The 3rd dimension specifies
%      the iteration of the taps.  This matrix can be fed into another
%      run as |bmWForce|.
% * |mcWall| - Cubic matrix of taps saved from the multiple-input
%      canceller, referenced in the same mammer as |bmWall| and
%      applicable in subsequent runs as |mcWForce|.
% * |snrAll| - Matrix of estimated SNR's at every iteration in dB
%      (matrix of size NxM)
% * |b| - Vector of DSB output
% * |z| - Matrix of output tracks from the blocking matrix.  These
%      tracks should hopefully approximate the noise in the simulation
%      as closely as possible.  This is supplied mostly as a debugging
%      output and can be checked to see how much target signal leakage
%      occurs in the beamformer.
%
% *References*
%
% * Hoshuyama, Osamu, and Akihiko Sugiyama. "Robust Adaptive
% Beamforming."  Microphone Arrays : Signal Processing Techniques and
% Applications. Ed. Michael Brandstein and Darren Ward. New York:
% Springer, 2001.
%
% Written by Phil Townsend (jptown0@engr.uky.edu) 6-10-08
%  Updated August 2014 by Kevin D Donohue to include inverse distance
%  weighting option


%% Function Declaration
function [y, bmWall, mcWall, snrAll, b, z] = ...
    gjbf(x, fs, s, m, c, p, q, mu, order, beta, phi, psi, K, ...
	 xPrev, bPrev, zPrev, bmWinit, mcWinit, snrThresh, ...
	 snrRate, snrInit, bmWForce, mcWForce,iw)

%% Argument Error Checking
if nargin == 23
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
elseif isempty(p) || ~isscalar(p) || ~isreal(p) || ...
        ~isfinite(p) || p < 0 || abs(mod(p,floor(p))) > 0
    error('p must be a positive integer');
elseif isempty(q) || ~isscalar(q) || ~isreal(q) || ...
        ~isfinite(q) || q < 0 || abs(mod(q,floor(q))) > 0
    error('q must be a positive integer');
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
elseif ~isempty(phi) && (length(size(phi))~=2 || ...
			 ~isreal(phi) || ~all(all(isfinite(phi))))
    error('phi must be a real matrix')
elseif ~isempty(psi) && (length(size(psi))~=2 || ...
			 ~isreal(psi) || ~all(all(isfinite(psi))))
    error('psi must be a real matrix')
elseif xor(isempty(phi), isempty(psi))
    error('phi and psi must both be specified or both empty')
elseif ~isempty(K) && ...
	(~isscalar(K) || ~isreal(K) ||  ~isfinite(K) || K < 0)
    error('K must be a real positive scalar');
elseif isempty(xPrev) || ~isreal(xPrev) || ~all(all(isfinite(xPrev)))
    error('xPrev must be a real matrix')
elseif ~isvector(bPrev) || ~isreal(bPrev) || ~all(isfinite(bPrev))
    error('bPrev must be a real vector');
elseif isempty(zPrev) || ~isreal(zPrev) || ~all(all(isfinite(zPrev)))
    error('zPrev must be a real matrix')
elseif ~isempty(bmWinit) && ...
	(~isreal(bmWinit) || ~all(all(isfinite(bmWinit))))
    error('bmWinit must be a real matrix')
elseif ~isempty(mcWinit) && ...
	(~isreal(mcWinit) || ~all(all(isfinite(mcWinit))))
    error('mcWinit must be a real matrix')
elseif ~isempty(snrThresh) && ...
    (~isreal(snrThresh) || ...
     ~isscalar(snrThresh) || ~isfinite(snrThresh))
    error('snrThresh must be a positive real scalar')
elseif ~isempty(snrRate) && ...
    (~isreal(snrRate) || ...
     ~isscalar(snrRate) || ~isfinite(snrRate) || snrRate <= 0)
    error('snrRate must be a positive integer')
elseif xor(isempty(snrThresh), isempty(snrRate))
    error(['snrThresh and snrRate must both be specified or' ...
	   ' both empty'])
end

%% Setup
[N, M] = size(x);  % number of samples in track, mics in array

% SNR parameters
if isempty(snrRate)  % no SNR thresholding, so assign dummy vals
    snrAll = [];
else  % apply SNR thresholds
    snrAll = zeros(N, M);  % all estimated SNRs mat
end

% Fixed Beamformer (FBF)
% Find the distances of each mic to the source, calculate the
% time-difference of arrival, and then add shifted waveforms.  Note
% that we don't simply use the DSB function since we want to keep
% access to the aligned tracks and it's not too much code.
d = sqrt(sum((m - s*ones(1, length(m(1,:)))).^2, 1)); % distances 
t = (max(d)-d)/c;  % Relative TDOA
delays = round(t*fs);  % convert time to integer samples
xShift = zeros(size(x));  % Initialize shifted matrix
if iw == 1
    wghts = arweights(d);  % inverse distance weighting
else
    wghts = ones(size(d));  % no weighting
end

for k=1:size(x,2)  % align signals, pulling real data from rear
    xShift(:,k) = wghts(k)*[xPrev(end-delays(k)+1:end, k); ...
		   x(1:end-delays(k), k)];
end
b = sum(xShift, 2);  % add data across columns


% Blocking Matrix (BM)
if isempty(phi)  % traditional GJBF BM
    z = diff(xShift, [], 2); % column differences
    mcAdapt = []; % mc always adapts
    bmWall = [];  % Matlab complains if an output arg isn't set..
else  % call function for adaptive BM
    for k=1:size(x,2)  % Apply p delays to each xShift column
	xShift(:,k) = wghts(k)*[xPrev(end-p+1:end, k); ...
		       x(1:end-p, k)];
    end
    % Call external function to handle BM calculations
    [z, bmWall, mcAdapt, snrAll] =  ...
	bmlms(xShift, b, mu, order, beta, phi, psi, bPrev, ...
	      bmWinit, bmWForce, snrThresh, snrRate, snrInit);
end


% Multiple Input Canceller (MC)
% Apply q delay to the DSB vector before running the MC
b = [bPrev(end-q+1:end); b(1:end-q)];
% Call external function for calculations
[y, mcWall] = mclms(z, b, mu, order, beta, zPrev, mcWinit, ...
                    mcWForce, mcAdapt, K);

end  % function gjbf
