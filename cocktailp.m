function sigout = cocktailp(sigin, spos, mpos, recinfo)
% This function simulates placing sound sources at specified positions in
% a room and recording the resulting sound field with spatially distributed
% microphones.  The function syntax is given by:
%
%      sigout = cocktailp(sigin, spos, mpos, recinfo)
%
% The inputs are:
%    sigin  =>  Matrix of sampled sound source files.  Each column
%               represents a different source
%    spos    => Matrix with same number of columns as "sigin" corresponding to 
%               the position of each source. Rows are the X,Y, and Z coordinates.
%    mpos    => Matrix describing the position of the microphones.  Each column
%               corresponds to the X,Y, and Z coordinate of a microphone.
%    recinfo => Data structure with recording environment parameters.  The fields are:
%               "fs"       sampling frequency used to resample the wave file to create "sigin"
%               "corners"  (Optional) a 2 column matrix denoting the coordinates of a
%                          pair of opposite corners of the rectangular room in which
%                          the recording was made.  If not present room reverberation
%                          will NOT be simulated.
%               "reflecc"  (optional) a 6 element vector denoting the reflection
%                          coefficients (between 0 and 1) of each of the room surfaces.
%                          The first 4 elements correspond to the walls going clockwise
%                          around the vertex given in the first column of "corners" (looking
%                          down from the ceiling), the next element is the ceiling and the
%                          last element is the floor.  If "corners" is not given these
%                          parameters will not be used.  The default values are .8 for 
%                          all walls and .6 for the floor and ceiling.
%                "dbd"    (optional) If reberation is simulatef, for multiple reflection this is the
%                          dB down value of wall reflection power from the first reflection
%                          for which is simulation no longer generates the reflection.  
%                          Default value is 60. 
%                "temp"    (optional) Temperature in Centigrade, default 22
%               "press"    (optional) Ambient air pressure in mmHg, default 29.92
%                 "hum"    (optional) Realative Humidity percent, default 38
%                "freq"    (optional) Frequencies at which to compute air
%                          attenuation function, default (fs/2)*[0:200]/200;
% 
% Written by Arul Kumaran Muthukumarasamy, Kevin D. Donohue (donohue@engr.uky.edu), Satoru Tagawa 06/25/2008
% Updated by Kevin D. Donohue June 4, 2013, Updated by Kevin D. Donohue
% August 5, 2014

%  Check for sampling frequency
if ~isfield(recinfo,'fs')
    errordlg('Fourth argument must be a data structure with field, fs, denoting sampling frequency if input signal array!', 'Missing Required Field')
    return
else
    fs = recinfo.fs;
end

%  Set attenuation parameters to default if not defined in input structure
if ~isfield(recinfo,'temp')
    recinfo.temp = 22; % Default
end
if ~isfield(recinfo,'press')
    recinfo.press = 29.92; % Default
end
if ~isfield(recinfo,'hum')
    recinfo.hum = 38; % Default
end
%  If frequencies for attenuation function computation not set, 
%  use default
if ~isfield(recinfo,'freq')
   recinfo.freq = (fs/2)*[0:200]/200;  % Default
end
dis = 1;  %  Distance in meters (normalized = 1)
[recinfo.atten] =  atmAtten(recinfo.temp, recinfo.press, recinfo.hum, dis, recinfo.freq);
recinfo.c = SpeedOfSound(recinfo.temp,recinfo.hum,recinfo.press);


% Check dimensionality of inputs
[samps, sourcesig] = size(sigin);
[dims, sourcepos] = size(spos);
if sourcesig ~= sourcepos
    errordlg('Number of given source signals must match source locations')
    return
end




% Check for room reverberation parameter
if isfield(recinfo,'corners')
    iflag = 1;  %  Set flag to perform image method in simulation
    
    forv = recinfo.corners;  %  Set opposite corners of rectangular room
    if isfield(recinfo,'reflecc')
        bs = recinfo.reflecc;  % Set room reflection coefficients
    else
        bs = [.8 .8 .8 .8 .6 .6]; %  Set default reflection coefficients
    end
else  % If room dimensions not given
    iflag = 0;  % set flag to not do the image method simulation
end

               % dB down level at which program stops simulating reverberations
if ~isfield(recinfo,'dbd')
    prop.dbd = 60; % Default
else
    prop.dbd = recinfo.dbd; 
end

% simulate the cocktail party signals over the microphones ****************

sigout = simarraysigim(sigin, fs, spos, mpos, forv, bs, recinfo);

