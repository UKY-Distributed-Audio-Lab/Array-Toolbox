function [sigout, tax] = simarraysigim(sigarray, fs, sigpos, mpos, forv, bs, prop)
% This function simulates stationary sounds received over an array of microphones in a
% rectangular room. Segments of multiple sound souces are contained in columns of
% SIGARRAY. Each source sound is assigned a position in space (corresponding to
% columns of SIGPOS), projected into the room (omnidirectional), and received over
% a microphone array. All distances and times provided in vectors MPOS, SIGPOS, FORV,
% and speed and attenuation parameters prop.c and prop.airatten must be in METERS, Hz,
% and SECONDS. Each microphone corresponds to a column in the output matrix SIGOUT.
% Multipath is computed using the image method for a rectangular room. The parameters
% associated with the multipath room properties are in matrix FORV and vector BS, in
% addition to optional parameters DBD, and SFAC.  The second output argument is the
% corresponding time axis TAX, where 0 corresponds to the first time sample in sigarray
% starting from its source postion.
% 
%  [sigout, tax] = simarraysigim(sigarray, fs, sigpos, mpos, forv, bs, prop)
%
% The input arguments are:
% sigarray => Matrix where each column represents a source signal and each row represents a time
%             sample of that signal.  All segments must be the same length and should be trimmed
%             to just the sound of interest.  Signals can be padded with zeros to make them all
%             the same length
% fs =>       Sampling frequency of input and output audio signals
% sigpos =>   3 row matrix where each column represents the x,y, and z coordinates of a
%             signal source.
% mpos =>     3 row matrix where each column represents the coordinates
%             of the mic and each row is the dimension x,y, and z.
% forv =>     2 column vector containing vertex indices of opposite corners
%             of rectangular room.
% bs   =>     6 element vector containing the reflection coefficients of the
%             walls, ceiling, and floor.  The first 4 elements correspond to
%             the walls starting at the vertex given in column 1 of FORV
%             and going in a clockwise looking down on the room. BS(5) is the
%             floor and BS(6) is the ceiling reflection coefficient.
% prop =>     A data structure with various fields related to the
%             propagation of sound.
%             The required fields are the speed of sound:
%             "c"        speed of sound in meters/second (FYI, c =
%                        331.4+0.6*Temperature in Centigrade)  This value
%                        can be computed from temperature, pressure, and
%                        humidity using function SpeedOfSound.m
%             and the decible limit at which point to stop simulating
%             wall reflections:
%             "dbd"     Decibles down from first room reflection at which
%             to stop simulating subsequet reflections.  60 is a typical
%             number.  If this is too high the simulation will take a long
%             time as it accounts for more attenuating reflections.
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
%  Output parameters:
%   sigout => matrix with each column corresponding to the signal recieved at each
%             microphone defined in MPOS.
%   tax =>    (optional output)  time axis associated with the rows of
%             SIGOUT
%
%
%  Written by Kevin D. Donohue (donohue@engr.uky.edu) September 2005
%  updated August 2014 to include other ways to implement frequency
%  dependent attenuation.
%

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
   
   if ~isfield(prop,'dbd')  %  Set reverberation threshold
       dbd = 60;  % if not given
   else
       dbd = prop.dbd;  % if given
   end

%  Check for limitations on simulation parameters
if dbd <= 0
    error('dB DOWN image scatterer threshold must be greater than 0')
end

% Obtain mic array information
[mr, mc] = size(mpos);      % Determine number of mics = mc
[pgr, pgc] = size(sigpos);  % Determine number of sound sources = pgc
[sgr, sgc] = size(sigarray); % Determine number of signals for each source = sgc
%  Check for consistency of input information
if sgc ~= pgc
    error('Each signal source must have a position - columns of SIGPOS = columns SIGARRAY')
end

if mr ~= pgr
    error('ErrorTests:convertTest',['The mic and signal coordinates must be in same dimension:\n rows of SIGPOS = rows MPOS'])
end

%  Initialize output array
sigout = zeros(1,mc);
% Loop to process every signal
for k=1:pgc  %  Every signal position
    for r=1:mc  %  with every mic position
        %  Create image scatterer for mic source combination
        [dlays, scals] = imagesim(forv, sigpos(:,k), mpos(:,r), bs, c, dbd);
        %  If scatterers exist, convert them into an impulse response
        if ~isempty(dlays)
          impres = roomimpres(dlays, scals, c, fs, prop);  %  Room impulse response
          dum = conv(sigarray(:,k),impres');  %  Filter input signal with response
          %dum = dum(fix(sgr/2-1):end);
          %  Check to ensure output matrix is large enough to accumulate the new
          %  filtered and delayed mic signal
          [orw,ocl] = size(sigout);     %  Current size of output array
          slen = length(dum);           %  length of current source to mic signal
          %  If source to mic signal is bigger than output array, extend it
          %  with zero padding
          if slen > orw
              sigout = [sigout; zeros(slen-orw,ocl)];
          end
          sigout(1:slen,r) = sigout(1:slen,r)+dum;  %  Accumulate source to mic signal in output array    
        end
    end
end

%  If second output requested compute time axis.
if nargout == 2
    [ro, co] = size(sigout);
    tax = ([0:ro-1]') /fs;
end