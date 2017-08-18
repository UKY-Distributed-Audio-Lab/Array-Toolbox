function [mleim] = mlescanprs(sigar, envpar, procpar)

% This function processes signals from a microphone array relative to
% measurments in space containing the microphones and the sound sources.
% The mic positions field of view grid point with any uint as long as
% the units for parameters listed below are consistent (I like meters
% and seconds).
% The output is the liklihood of the sound source as function of the
% spatial grid.  If only a 2-D grid is input for the Mic positions
% and FOV, then only a 2-D grid at the output is computed.
%
%   [mleim, fovxa, fovya] = mlescan(sigar, envpar, procpar)
%
% Inputs
% "sigar"         matrix where each column is the rf signal from each mic.
% "envpar.mpos"   3-row matrix indicating the x,y,z cartesion positions of each
%                 mic, where the first row contains the x positions, second y
%                 position, and third z positions.  If an axis is missing
%                 it will be assume to be the z and an array in the x-y
%                 plane will be assumed.
%           
%  "envpar.xg"    array of grid points for x-axis (from meshgird)
%  "envpar.yg"    array of grid points for y-axis (from meshgird)
%  "envpar.zg"    array of grid points for z-axis (from meshgird), if this
%                 variable is missing, the function will only do a 2-D scan.
%  "procpar.fs"    scalar is the sampling frequency of the rf mic signals
%  "envpar.c"     scalar is the speed of sound in the room.
%  "procpar.win"   scalar the correlation window in time units.  This should
%                 be less than the segment lenght of "sigar" or else zeros
%                 will be padded to allow for relative shifts between signals.
%  "procpar.mxdelay" Maximimum delay to consider in forming beamfields    
%  "procpar.mprs" Indeces of all mic pairs to be combined in beamfield creation 
%  "procpar.combine" Sting indicating how pairwise mic correlation should be combined
%                    mean => average, min = Minimum, med => Median, max =>
%                    Maximum, absmed = Median of Absolute values



dimfov = size(envpar.xg);  %  Get dimensions for FOV 
ndims = length(dimfov);    %  Check to see how many dimensions were provided
mleim = zeros(dimfov);     %  Initialize output image with Zeros

% Obtain mic array information
[mr, mc] = size(envpar.mpos);  % Determine # number of mics = mc 

[combs, nmics] = size(procpar.mprs);  % Determine # number of mics combinations
%                                       and number of mics per combination.  


[sr,sc] = size(sigar);  % Determine the mic signal array size [# of samples, # of mics]
winsamp = procpar.win*procpar.fs;  %  Convert processing window length to samples
temp = zeros(winsamp,sc);  %  Initialze correlation matrix
dlysamp = procpar.mxdelay*procpar.fs; %  Convert maximum delay length to samples
totwin = ceil(winsamp+dlysamp);  %  Compute maximum data window need to 
                                 %  accomodate maximum delay for
                                 %  beamforming

% If maxdelay requires data not in array, pad with zeros 
if totwin > sr;
    zpadd = totwin-sr+1;
    sigar = [sigar; zeros(zpadd, sc)];
end
nth = 0;
%  If only 2 dimension space is specified use double loop
if ndims == 2
    for kx=1:dimfov(1);
        for ky = 1:dimfov(2);
            %  Distance of position from microphones
            ds = envpar.mpos - [envpar.xg(kx,ky); envpar.yg(kx,ky)]*ones(1,mc);  % Spatial
            rst = sqrt(ds(1,:).^2 + ds(2,:).^2)/envpar.c;  %  Corresponding Time differences
            [mt, nmic] = min(rst);  % find closest mic to point being considered
            %  Make all delays relative to closest mic, setting closest mic to
            %  0 delay
            rstref = rst-mt;
            at = rst;
            at = (min(at) ./ (abs(at)));
            %  Convert initial and points to sample indices for each mic signal
            indbeg = round(rstref*procpar.fs) + 1;
            indend = winsamp + round(rstref*procpar.fs);
            %  Load up matrix with aligned mic signal corresponding to mic
            %  position
            pn = 0;
            for ksig = 1:sc
                 temp(:,ksig) = at(ksig)*(sigar(indbeg(ksig):indend(ksig),ksig));
                 pn = pn + sum(abs(temp(:,ksig)).^2)/2;
                 %temp(:,ksig) = temp(:,ksig)./(abs(temp(:,ksig))+eps);
                 %.*hamming(winsamp);
       %         temp(:,ksig) = temp(:,ksig); %/(sqrt(ksig)*std(temp(:,ksig)));
            end
%             figure(1)
%             imagesc(temp); colorbar
%             pause
           dumval = zeros(1,sc*(sc-1));
           icn = 1; 
           for kmps=1:sc-1
                parray = ((temp(:,kmps)))';
                 for da =kmps+1:sc
                     dumval(icn) = parray*temp(:,da);
                     icn=icn+1;
                 end
           end
             %[h, s] = ttest2(imag(dumval), real(dumval),.01);
             mleim(kx,ky) = (sum((dumval))/(sc-1))/pn;
             %mleim(kx,ky) = -log(s);
             
         %    nth = nth + imag(datmp)^2;

        end
    end
%  mleim  = mleim - 3*sqrt((nth)/prod(dimfov));
           
    
    %  
elseif ndims == 3
    %  Check to see if mics have 3 indeces for position
    mcs = size(envpar.mpos);
    %  If not assume mics are in a reference zero plane
    if mcs < 3;
        envpar.mpos = [envpar.mpos; zeros(1,msc(2))];
    end
    for kx=1:dimfov(1)
        for ky = 1:dimfov(2)
            for kz = 1:dimfov(3)
                %  Distance of position from microphones
                ds = envpar.mpos - [envpar.xg(ky,kx,kz); envpar.yg(ky,kx,kz); envpar.zg(ky,kx,kz);]*ones(1,mc);  % Spatial
                rst = sqrt(ds(1,:).^2 + ds(2,:).^2 + ds(3,:).^2)/envpar.c;  %  Corresponding Time differences
                [mt, nmic] = min(rst);  % find closest mic to point being considered
                %  Make all delays relative to closest mic, setting closest mic to
                %  0 delay
                rstref = rst-mt;
                %  Convert initial and points to sample indices for each mic signal
                indbeg = round(rstref*fs) + 1;
                indend = winsamp + round(rstref*fs) + 1;
                %  Load up matrix with aligned mic signal corresponding to mic
                %  position
                for ksig = 1:sc
                    temp(:,ksig) = sigar(indbeg(ksig):indend(ksig),ksig);
                end
                for kmps=1:length(procpar.mprs)
                     mleim(ky,kx,kz) = mleim(ky,kx,kz)+(temp(:,procpar.mprs(kmps,1))').*temp(:,procpar.mprs(kmps,2));
                end
                
            end
        end
    end

end
