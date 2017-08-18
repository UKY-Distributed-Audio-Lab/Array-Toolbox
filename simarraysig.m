function [sigout, tax] = simarraysig(sigarray, fs, sigpos, mpos, prop)
% This function inputs an array of sampled sound sources, assigns spatial
% position to each, and simulates corresponding signals received
% over a distribution of microphones. All distances and times MUST BE in METERS
% and SECONDS.  This includes the spatial positions of mics, signals, and the
% sound speed. Each microphone corresponds to a column in the output matrix
% SIGOUT.  The second output argument is the corresponding time axis
% TAX, where 0 corresponds to the first time sample in sigarray starting
% from its source position.
% 
%  [sigout, tax] = simarraysig(sigarray, fs, sigpos, mpos, prop)
%
% The input arguments are:
% sigarray => Matrix where each column represents a source signal
%             and each row represents a time sample of that signal
% fs =>       Signal sampling frequency
% sigpos =>   2 or 3 row matrix where each column represents the coordinates of a
%             signal source and each row is the dimension x,y, and z (3D simulation
%             if z row present)
% mpos =>     2 or 3 row matrix where each column represents the coordinates
%             of the mic and each row is the dimension x,y,z (z if 3D included)
% prop =>     A data structure with various fields related to the
%             propagation of sound.
%             The required fields are the speed of sound:
%             "c"        speed of sound in meters/second (FYI, c =
%                        331.4+0.6*Temperature in Centigrade)  This value
%                        can be computed from temperature, pressure, and
%                        humidity using function SpeedOfSound.m
%             and the decible limit at which point to stop simulating
%             wall reflections:
%             "mpath"    3 or 4 row vecotor containing the x,y,z (z if 3D included)
%                        coordinates of multipath scatterers in the last 2 or 3 rows.
%                        The first row is the isotropic scattering coefficient (unitless
%                        and less than 1 for a passive scatterer).  Only first order
%                        multipath is considered.  If not present, multipath is not generated.
%             If field "airatten" is present, it is taken as the frequency dependent
%             scale factor in units of dB per meter-Hz. A typical value is -3.2808399e-5.
%             If set to 0 no attenuation due to the air path is applied (the value must
%             be either 0 or negative). They will be applied directly to the propagating wave.
%                  If the "airatten" field is NOT present and instead "freq"
%             and "atten" are included where "prop.freq" is a vector of 
%             frequency points with "prop.atten" the corresponding
%             attenuation points in dB, the these will be used to implement
%             the frequency dependent attenuation.  See atmAtten.m for
%             generating these points based on temperature pressure and humidity.
%          
%
%  Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005 (updated
%  October 2009, updated August 2014
%


% Obtain mic array information
[mr, mc] = size(mpos);      % Determine number of mics = mc
[pgr, pgc] = size(sigpos);  % Determine number of sound sources = pgc
[sgr, sgc] = size(sigarray); % Determine number of signals for each source = sgc
%  Check for consistency of input information
if sgc ~= pgc
    error('Each signal source must have a position - columns of SIGPOS = columns of SIGARRAY')
end


if mr ~= pgr
    error('ErrorTests:convertTest',['The mic and signal coordinates must be in same dimension:\n rows of SIGPOS = rows MPOS'])
end



%  If optional parameter data structure not given, set default parameters
%
   if ~isfield(prop,'c')  % Set speed of sound
       c = 345;   % if not given
   else
       c = prop.c; %  if given
   end
   if ~isfield(prop,'dbd')  % Set speed of sound
       dbd = 60;   % if not given
   else
       dbd = prop.dbd; %  if given
   end
   if ~isfield(prop,'freq')
      prop.freq = (fs/2)*[0:200]/200;
   end
    %  Set attenuation in air and check for other optional fields
   if ~isfield(prop,'airatten') && ~isfield(prop,'atten')
    %   If no fields for air attenuation are given use default
    %  Create attenuation values with nominal temperature, humidity, and
    %  pressure
    temp = 22; % Temperature centigrade
    press = 29.92; % pressure inHg
    hum = 38;  % humidity in percent
    dis = 1;  %  Distance in meters
    prop.atten =  atmAtten(temp, press, hum, dis, prop.freq);
   elseif ~isfield(prop,'atten') && isfield(prop,'airatten') 
       prop.atten = prop.freq*prop.airatten;  % if given
   end
   
   if ~isfield(prop,'mpath')  %  Set multipath scatterers
       mpath = [];  % if not given
       mgc = 0;
   else
       mpath = prop.mpath;  % if given
       [mgr, mgc]= size(mpath);   % Get dimension and number of scatterers
       if (mgr-1 ~= pgr) 
        error('ErrorTests:convertTest',['mic position and multipath scatterers position', ...
         ' \n must be in same spatial dimension (note the multipath matrix has inital row', ...
         ' \n providing the reflection coefficient of the scatterer, with subsequent rows', ...
         ' \n being the position coordinates']);
       end
   end



                
%  Find largest delay to initialize ouput array with appropiate size
dmax = 0;
for k=1:pgc  %  Every signal position
    for r=1:mc  %  with every mic position
        if mgc > 0    %  and if present ..
            for m = 1:mgc   %  with every multipath scatterer
                d = norm((sigpos(:,k)-mpath(2:end,m)),2) + norm((mpath(2:end,m)-mpos(:,r)),2);
                if d > dmax
                    dmax = d;  % new delay bigger than previous, replace max value
                end
            end
        else   %  If no multipath scatterers specified, just do direct paths
            d = norm((sigpos(:,k)-mpos(:,r)),2);
            if d > dmax
                 dmax = d;
            end
        end
    end
end
%  Compute output array length based on the greatest delay
outlen = ceil(dmax*fs/c) + sgr;
tax = (0:outlen-1)/fs;  %  Create corresponding time axis
sigout = zeros(outlen,mc);  % Initialize output array
%  Compute a number of frequency domain coefficients to describe the 
%  attenuation spectrum
len = round(fs/40);
nfax = (0:len)/len;
fax = (fs/2)*nfax;    %  Create frequency axis for the attenuation

%  Loop for each source
for k = 1:pgc
    % First take care of direct path delay and attenuation from each sound source to mic    
    for r = 1:mc
        md = zeros(1,1+mgc);  %  Initialize delay vector
        scl = zeros(1,1+mgc); %  Initialize scale vector
        md(1) = norm((sigpos(:,k) - mpos(:,r)),2);  %  Compute distance between mic and source
        scl(1) = 1/(4*pi*(1+md(1)));   %  No scattering loss for direct path 
        %  If multipath present
        if mgc > 0    %  and if present ..
            for m = 1:mgc   %  with every multipath scatterer
                %  Compute distance from the source to scatterer to mic
                md(1+m) = norm((sigpos(:,k)-mpath(2:end,m)),2) + norm((mpath(2:end,m)-mpos(:,r)),2);
                scl(1+m) = mpath(1,m)/(4*pi*(1+md(1+m)));
            end
        end
        impres = roomimpres(md/c, scl, c, fs, prop);  %  Room impulse response
        dum = conv(sigarray(:,k),impres');  %  Filter input signal with response
         %  Check to ensure output matrix is large enough to accumulate the new
         %   filtered and delayed mic signal
          [orw,ocl] = size(sigout);     %  Current size of output array
          slen = length(dum);           %   length of current source to mic signal
          %  If source to mic signal is bigger than output array, extend it
          %   with zero padding
          if slen > orw
              sigout = [sigout; zeros(slen-orw,ocl)];
          end
          sigout(1:slen,r) = sigout(1:slen,r)+dum;  %  Accumulate source to mic signal in output array   
    end 
end
[orw,ocl] = size(sigout);     %  Current size of output array
tax = [0:orw-1]/fs;   %  Create time axes
   