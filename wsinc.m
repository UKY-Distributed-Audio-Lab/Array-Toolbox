function s = wsinc(t,wi)
%  This function creates a tappered sinc function.  This is useful
%  for generating noncausal FIR interpolation coefficients.
%
%     s = sinc(t,wi)
%
%  T is the vector the time axis points for the sinc function for a normalized
%  sampling interval (zero nulls occur on all non-zeros integer values).  Since
%  this is designed for interpolation the T axis should always cover the 0 point
%  (i.e. point on either side of 0).
% 
%  WI is an optional string argument identifing the type of weighting window to
%  apply on the truncated sinc coefficients:
%  WI = 'Linear' is a linear interoplator coefficient (no sinc involved)
%  WI = 'Triangle' is a triangluar window
%  WI = 'Cos' is a cosine squared window
%  WI = 'Hanning' is a hanning window (if you type an string not recognized
%                 this will be applied by default
%  
%
%  If a second argument is not present no weighting occurs.
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005

if nargin == 2
    %  Apply Linear interpolation coefficients
    if wi(1) == 'l' | wi(1) == 'L'
        s = zeros(size(t));
        g1 = find(abs(t) == min(abs(t)));
        rng =[(g1(1)-1):(g1(1)+1)];
        s(rng) = 1-abs(t(rng))/max(abs(t(rng)));
        s = s/(sum(s));
    %  Apply triangle window
    elseif wi(1) == 't' | wi(1) == 'T'
        s = (1-abs(t)/max(abs(t))).*sinc(t);
        s = s/sqrt(sum(s.^2));
    %  Apply half Cosine squared window
    elseif wi(1) == 'c' | wi(1) == 'C'
        s = (cos(0.5*pi*(t)/max(abs(t))).^2).*sinc(t);
        s = s/sqrt(sum(s.^2));
    %  Apply hamming window
    else
        s = (.54 -.46*cos(2*pi*(t+min(t))/(max(t)-min(t)))).*sinc(t);
        %s = s/sqrt(sum(s.^2));
    end
%  Otherwise apply no window
else
    s = sinc(t);
    s = s/sqrt(sum(s.^2));
end
