function [mags,pos] = dent3dwei(imi,rad,fprob,b)
%
%  This function finds the magnitudes and positions of SRCP values in
%  3D or 2D array IMI, that cross a CFAR threshold.
%
%       [mags,pos]  = dent3dwei(imi,rad,fprob,b)
%
%
%  The CFAR threshold is based on the Weibull distribution of SRCP noise 
%  values with the optional shape parameter B (Default: 0.8) and the scale
%  parameter calcuated using the negative SRCP values in the given
%  neighborhood of RAD pixels around the point in question, for a specific
%  false alarm probability FPROB.  
%  
%  It steps through every point in the image, but only calculates for those
%  that are positive. RAD can be either scalar or vector.  If scalar, the 
%  neighborhood is extended from the test point in all dimensions by this 
%  amount in pixels.  If vector, the neighborhood is extended by
%  different amounts for each dimension corresponding to the entries in
%  the RAD vector:    [x-neighborhood, y-neighborhood, z-neighborhood]
%  
%  Note: Since we are indexing matrices, first element in array is the row
%  or x, which on plots is height or y.  So x and y dimensions should be
%  switched to fit a basic plotting notion.
% 
%  Inputs: 
%  imi    =>    Acoustic image computed by function SRPFAST.M
%  rad    =>    Neighborhood of pixels around point in question.
%               Can be either scalar or a vector.
%  fprob  =>    False positive detection probability for threshold
%               computation.
%  b      =>    Partial whitening parameter (between 0 and 1).
% 
% 
%  Outputs:
%  mags   =>    Magnitudes of significant sources
%  pos    =>    Positions of significant sources    
%
%
%  Written by Kevin D. Donohue (donohue@engr.uky.edu) April 2010
%  Modified by Kirstin Brangers                       July 2012


if nargin == 3
b = 1.2;         % Set default shape parameter 
end

% Check RAD neighborhood
if rad == 0
	error('rad == 0, Cannot process thresholding algorithm without a neighborhood.');
end

%  Compute scale factor based on shape parameter with unit scaling
scl = log(1/fprob)^(1/b);

% Initalize output arrays
mags=[];                % Excess threshold value
pos=[];                 % Position of detection

dims = size(imi);       % Find dimensions of input image
dimnum = length(dims);  % Find number of dimensions

if dimnum == 3
    dimn = length(rad); % Find dimension of neighborhood values
    if dimn == 1
        radlenx = rad;
        radleny = rad;
        radlenz = rad;
    elseif dimn == 2
        radlenx = rad(1);
        radleny = rad(2);
        radlenz = 1;
    else
       radlenx = rad(1);
       radleny = rad(2);
       radlenz = rad(3);
    end
elseif dimnum == 2
    dimn = length(rad); % Find dimension of neighborhood values
    if dimn == 1
        radlenx = rad;
        radleny = rad;
        radlenz = 1;
    elseif dimn == 2
        radlenx = rad(1);
        radleny = rad(2);
        radlenz = 1;
    end
end

if dimnum == 3
    dimx = dims(1);
    dimy = dims(2);
    dimz = dims(3);
    count=0;

    for kz = radlenz+1:dimz-radlenz             % kz = radlenz+1:dimz-radlenz
        zr = (kz-radlenz):(kz+radlenz);
        for ky = radleny+1:dimy-radleny         % ky = radleny+1:dimy-radleny
            yr = (ky-radleny):(ky+radleny);
            for kx = radlenx+1:dimx-radlenx     % kx= radlenx+1:dimx-radlenx
                xr = (kx-radlenx):(kx+radlenx); 
                
                %  Test if positive and local maximum
                if (imi(kx,ky,kz) > 0 && imi(kx,ky,kz) == max(max(max(imi(xr,yr,zr)))))                    
                    % Vectorize all points in neighborhood
                    neigh = reshape(imi(xr,yr,zr),length(xr)*length(yr)*length(zr),1);
                    % If there are no negative values count as crossing threshold
                    if isempty(neigh(neigh<0))
                        count=count+1;                   % Detection count increment
                        mags(count,1) = imi(kx,ky,kz);   % Save magnitude
                        pos(count,:) = [kx,ky,kz];       % Save position
                    else                                 % Negative values exist
                        % Compute Weibull scale parameter
                        a = (mean(abs(neigh(neigh<0)).^b))^(1/b); 
                        thresh = a*scl ;                 % Compute threshold          
                        detstat = imi(kx,ky,kz)- thresh; % Excess threshold value 
                        if(detstat>0)                    % If threshold exceeded
                            count=count+1;               % Count it!
                            mags(count,1) = detstat;
                            pos(count,:) = [kx,ky,kz];
                        end   
                    end         % If positive peak and there are neg points in neighbhorhood
                end             % If positive peak

            end                 % End for kx
        end                     % End for ky
    end                         % End for kz
    
elseif dimnum ==2
    dimx = dims(1);
    dimy = dims(2);
    count=0;

        for ky = (radleny+1):(dimy-radleny)  
            yr = (ky-radleny):(ky+radleny);
            for kx = (radlenx+1):(dimx-radlenx) 
                xr = (kx-radlenx):(kx+radlenx);
                
                %  Test if positive and local maximum
                if (imi(kx,ky) > 0 && imi(kx,ky) == max(max(max(imi(xr,yr)))))
                    % Vectorize all points in neighborhood
                    neigh = reshape(imi(xr,yr),length(xr)*length(yr),1);
                    % If there are no negative values, count as crossing threshold
                    if isempty(neigh(neigh<0))
                        count=count+1;                % Detection count increment
                        mags(count,1) = imi(kx,ky);   % Save magnitude
                        pos(count,:) = [kx,ky];       % Save position
                    else                              % Negative values exist
                        % Compute Weibull scale parameter
                        a = (mean(abs(neigh(neigh<0)).^b))^(1/b); 
                        thresh = a*scl ;              % Compute threshold          
                        detstat = imi(kx,ky)- thresh; % Excess threshold value 
                        if(detstat>0)                 % If threshold exceeded
                            count=count+1;            % Count it!
                            mags(count,1) = detstat;
                            pos(count,:) = [kx,ky];
                        end   
                    end         % if positive peak and there are neg points in neighbhorhood
                end             % if positive peak
                
            end                 % end for kx
        end                     % end for ky
end                             % end else statement
