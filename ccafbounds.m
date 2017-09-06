%% CCAF Bound Calculator
% This function calculates the constant upper and lower matrices phi
% and psi for the Coefficient-Constained Adaptive Filters (CCAF's) for
% the Blocking Matrix (BM) of a Robust Generalized Sidelobe Canceller
% (GSC) Beamformer.
%
% *Syntax*
%
% |[phi, psi] = ccafbounds(m, fs, c, p, order)|
%
% *Inputs*
%
% * |m| - 3xM matrix of microphone positions, each column a coordinate
%     in R^3 (meters)
% * |fs| - Audio sample rate (Hertz)
% * |c| - Speed of sound (meters/sec)
% * |p| - Estimated propagation time across the array in samples
% * |order| - Order of the adaptive filters
%
% *Outputs*
%
% * |phi| - Matrix of upper bounds, where each column is a vector of
%      bounds for a the adaptive filter of a single track in the BM.
% * |psi| - Matrix of lower bounds with the same structure as psi.
%
% *Notes*
%
% This code is based on equations derived for a linear microphone
% array (see references) where the effect of the coefficient bounds
% would be to expand the main lobe of the beamformer depending on the
% sine of a parameter delta-theta, the supposed error in steering
% angle of the beamformer.  However, here we consider a beamformer of
% arbitary geometry and thus the parameter "delta-theta" no longer
% makes sense.  Our current adaptation of the algorithm is
% 
% # Consider the center of the array as the centroid of array,
%   calculated as the arithmetic mean of the mic coordinates.
% # Hard-code a value for sin(delta-theta), which we know must be
%   bounded [-1 1].  At present we've selected .05
% 
% Note that these adaptations will still recreate the original
% performance of the algorithm for a linear array, where our selection
% of .05 for sin(delta-theta) corresponds to a steering angle error of
% about +/- 3 degrees.
%
% *References*
%
% * Hoshuyama, Osamu, Akihiko Sugiyama, and Akihiro Hirano. "A Rboust
% Adaptive Beamformer with a Blocking Matrix Using
% Coefficient-Constrained Adaptive Filters." IEICE Trans. Fundamentals
% E82-A (1999): 640-47.
%
% Written by Phil Townsend (jptown0@engr.uky.edu) 8-12-08

%% Function Declaration
function [phi, psi] = ccafbounds(m, fs, c, p, order)

%% Argument Error Checking
narginchk(5,5);
if ~isreal(m) || ~all(all(isfinite(m))) || size(m,1) ~= 3
    error('m must be a real matrix with 3 rows')
elseif ~isscalar(fs) || ~isreal(fs) || ~isfinite(fs) || ...
        fs < 0 || abs(mod(fs, floor(fs))) > 0
    error('fs must be a positive integer')
elseif ~isreal(c) || ~isscalar(c) || ~isfinite(c) || c <= 0
    error('c must be positive real scalar');
elseif ~isscalar(p) || ~isreal(p) || ~isfinite(p) || ...
        p < 0 || abs(mod(p, floor(p))) > 0
    error('p must be a positive integer')
elseif ~isscalar(order) || ~isreal(order) || ~isfinite(order) || ...
        order < 0 || abs(mod(order, floor(order))) > 0
    error('order must be a positive integer')
end

%% Calculate Bounds
% In the original paper, the bound vectors for each adaptive filter are
% calculated as
%
% $$\phi_{m,n} = \frac{1}{\pi\ max(.1, \ 
% (n\ ^\_ \ P)\ ^\_ \ T_m, \ ^\_ (n\ ^\_ \ P)\ ^\_ \ T_m)} \quad
% \psi = \ ^\_ \phi\ \forall\ m, n$$
%
% $$T_m = \frac{b_mf_s}{c}\sin\Delta\theta $$
%
% where bm is the distance of the mth microphone to the center of the
% array.  Remember that we must fudge for sin(delta-theta) in R^3.

M = size(m,2);  % number of microphones in the array
phi = zeros(order, M);  % initialze upper bound matrix for iteration

for mIter = 1:M   % iterate over all microphones
    sinDt = .05;  % kludge for 3-D (see Notes section above)
    arrayCentroid = mean(m,2);  % use centroid as "center" of array
    bm = sqrt(sum((m(:,mIter)-arrayCentroid).^2));  % Get mic distance
                                                    % from centroid
    Tm = bm*fs*sinDt/c;  % Hoshuyama equation
    for nIter = 1:order  % Set bound for each tap of this adaptive filter
        phi(nIter, mIter) = 1/(pi*max([.1, (nIter-p)-Tm, ...
            -(nIter-p)-Tm]));  % directly from Hoshuyama paper
    end
end
psi = -phi;  % psi is simply the opposite of phi

end
