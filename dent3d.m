function [mags,pos] = dent3d(imi,rad,farate)
%
%       [mags,pos]  = dent3d(imi,rad,farate)
%
%  This function finds the magnitude and its corresponding position of when
%  a SRCP image 'imi', crosses a certain threshold.
%
%  The threshold is based on the Weibull distribution of the noise SRCP 
%  values with the shape parameter, 'b=.63', given and the scale parameter 
%  calcuated using the negative values in the given neighborhood of 'rad' 
%  around the point in question, for a specific false alarm rate.  
%  
%  It steps through every point in the image, but only calculates for those
%  that are positive.  'farate' can be a scalar or a row vector.  
%  
%  **Note: All values in 'imi' are squared.**

if rad == 0
	error('rad == 0, can not process thresholding algorithm');
end

mags=[];
pos=[];

dims = size(imi);       % Find dimensions of input image
dimnum = length(dims);  % Find number of dimensions
if dimnum ~= 3
    disp('Input dimension not equal to 3')
    return
end

%  Assign parameters for local sliding window
kval = .5;
radlenx = rad;
radleny = rad;
radlenz = rad;

dimx = dims(1);
dimy = dims(2);
dimz = dims(3);

scl = gamthresh(farate,kval);

count=0;
for kz = radlenz+1:dimz-radlenz
    for ky = radleny+1:dimy-radleny
        for kx = radlenx+1:dimx-radlenx
            % [kz,ky,kx]
            if (imi(kx,ky,kz) > 0 && imi(kx,ky,kz) == max(reshape(imi(kx-1:1+kx,ky-1:1+ky,kz-1:1+kz),1,27)))
        	xr = [max([1,kx-radlenx]):min([kx+radlenx,dimx])];
                yr = [max([1,ky-radleny]):min([ky+radleny,dimy])];
                zr = [max([1,kz-radlenz]):min([kz+radlenz,dimz])];

                [ti, tj] = find(imi(xr,yr,zr) < 0);
                % If there are no negative values count as crossing
                % threshold
                if isempty(ti)
                	count=count+1;
                    mags(count) = imi(kx,ky,kz);
                    pos(count,:) = [kx,ky,kz];
                else
                    imio = imi(xr,yr,zr).^(2);
                    tjj = mod(tj-1,length(yr)) + 1;
                    zj = ceil(tj/length(yr));
                    as = zeros(length(ti),1);
                    for kj = 1:length(ti)
                         as(kj) = imio(ti(kj),tjj(kj),zj(kj));
                    end
                    
                    a = mean(as);
                    thresh = a*scl; 
                    detstat = imi(kx,ky,kz)^2 - thresh;
                    
                    if(detstat>0)
                        count=count+1;	
                        mags(count) = detstat;
                        pos(count,:) = [kx,ky,kz];
                    end                        
                end
            end
        end % end for kx
    end % end for ky
end % end for kz

end