function s = wsincc(t)
%  This function creates a tappered sinc function.  This is useful
%  for generating noncausal FIR interpolation coefficients.
%
%     s = wsincc(t)
%
%  T is the vector the time axis points for the sinc function for a normalized
%  sampling interval (zero nulls occur on all non-zeros integer values).  Since
%  this is designed for interpolation the T axis should always cover the 0 point
%  (i.e. point on either side of 0).
% 
%  A weighting window is applied on the truncated sinc coefficients, which is 
%  simple is a cosine squared window

%   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005

    %  Apply half Cosine squared window
       s = (cos(0.5*pi*(t)/max(abs(t))).^2).*sinc(t);
       s = s/sqrt(sum(s.^2));
