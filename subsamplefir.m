function tab = subsamplefir(n,ord, wi)
%  This function generates a table of interpolation filter coefficients for N+1
%  equi-spaced points between and including 2 consecutive points in the original
%  sample grid with filter order ORD (order should be even).
%
%    tab = subsamplefir(n,ord, wi)
%
%  The FIR is a noncausal weighted sinc fiter with even order (if ORD is odd it
%  will be rounded up to the next even value).  The ouput can be used
%  directly with function "delaytab" for implementing delays.
%
%  WI is an optional string argument identifying the type of weighting window to
%  apply on the truncated sinc coefficients:
%  WI = 'Triangle' is a triangular window
%  WI = 'Cos' is a cosine squared window
%  WI = 'Hanning' is a hanning window (if you type a string not recognized
%                 this will be applied by default
%  WI = 'Linear' is a linear interpolator coefficient (no sinc involved,
%  always order 2)
%  
%  If a second input argument is not present, no weighting occurs on the sinc function.
%
%   written by Kevin D. Donohue (donohue@engr.uky.edu)
%


%  Signal Processing parameters
ordh = ceil(ord/2);
sp = (-ordh+1:ordh);  % grid for evaluation of sinc interpolators, ensure equal number
                      % of samples on either side of interpolation point.

%  If linear interpolation selected
if nargin == 3 && (wi(1) == 'l' | wi(1) == 'L')
  tab = zeros(n+1,2);   %  Initalize output table, only need second order
    %  Compute Linear interpolation coefficients for each sub resolution point
    for k=1:n+1
        t = (0:1)-(k-1)/n;  %  Shift for fractional weights
        tab(k,:) = abs(t);  %  store in output table
    end
    
%  If another weighting function is selected for the SINC interpolation
elseif nargin == 3
    tab = zeros(n+1,length(sp));  % Initalize output table
    %  Compute FIR coefficients for each sub resolution point
    for k=1:n+1
            t = sp-(k-1)/(n);   %  Fractional shift
            %  Triangle window
            if wi(1) == 't' | wi(1) == 'T'
                s = (1-abs(t)/(max(abs(t))+1)).*sinc(t);
                tab(k,:) = s;
            %  Apply half Cosine squared window
            elseif wi(1) == 'c' | wi(1) == 'C'
                s = (cos(0.5*pi*(t)/(max(abs(t))+1)).^2).*sinc(t);
                tab(k,:) = s; 
            %  Apply hamming window
            else
                s = (.54 -.46*cos(2*pi*(t+min(t))/(max(t)-min(t)))).*sinc(t);
                tab(k,:) = s; 
            end
    end
% If no weighting window was indicated, just use the SINC function weights
else   
     tab = zeros(n+1,length(sp));  % Initalize output table
      for k=1:n+1
        t = sp-(k-1)/(n);  %  Fractional shifts
        s = sinc(t);
        tab(k,:) = s/sqrt(sum(s.^2));  %  Normalize energy
     end
end
