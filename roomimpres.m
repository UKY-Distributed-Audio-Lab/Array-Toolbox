function [impres, t] = roomimpres(dlays, scls, c, fs, sfac)
%  This function computes the impulse response of a source-microphone pair in a room
%  modeled by a direct path and a sequence of image/multipath scatterers described by
%  delay times in vector DLAYS and corresponding attenuation factors in vector SCLS.
%  In addition, the generated impulse response accounts for the frequency
%  dependent attenuation in the air, based on vectors in data structure
%  SFAC.
%
%     [impres, t] = roomimpres(dlays, scls, c, fs, sfac)
%
%  The following parameters must be supplied at the input:
%  DLAYS  => Vector with sequence of delays in seconds representing multipath
%            interference
%  SCLS   => Vector with scale factors (reflection coefficients) correspoding to signals
%                    delayed by each element in DLAYS
%  C      => the speed of sound in meters per second
%  FS     => the sampling frequency for the impulse response
%  SFAC   => (optional) data structure with the frequency dependent attenuation
%                     vectors:
%                     sfac.freq = vector of frequencies where attenuation values are provided
%                     sfac.atten = vector of corresponding attenuation values in dB per meter (must be positive)
%                      see function atmAtten() for generating these points
%                      based on environmental conditions.
%                      The default value is set according the the ISO 9613-1 for temperature = 22 C, pressure
%                      29.92 inHg, Relative Humidity 38%.
%                      
% The outputs are:
% IMPRES => Vector containing the impulse response
% T      => optional parameter, a vector of time axis points corresponding
%           to the impuse response IMPRES
%
%     written by Kevin D. Donohue (donohue@engr.uky.edu) September 2005
%     (updated June 2009)
%

%   If attenuation factor not given
if nargin < 5
    sfac.freq = (fs/2)*[0:200]/200;
    %  Create attenuation values with nominal temperature, humidity, and
    %  pressure
    temp = 22; % Temperature centigrade
    press = 29.92; % pressure inHg
    hum = 38;  % humidity in percent
    dis = 1;  %  Distance in meters
    sfac.atten =  atmAtten(temp, press, hum, dis, sfac.freq);
end
if ~isfield(sfac,'freq')
    sfac.freq = (fs/2)*[0:200]/200;
end

%  The frequeny range must extend from 0 to fs/2 to be compatable with FIR2
if sfac.freq(1) ~= 0   %  Extend to 0 with nearest neighbor value
    sfac.freq = [0, sfac.freq];
    sfac.atten = [sfac.atten(1), sfac.atten];
end  %  Extend to fs/2 with nearest neighbor value
if sfac.freq(end) ~= fs/2
    sfac.freq = [sfac.freq, fs/2];
    sfac.atten = [sfac.atten, sfac.atten(end)];
end


len = 50;  %   FIR filter length minimum limit
nfax = sfac.freq/(fs/2);  %  Normalized frequency axis for FIR filter specs
%  Find maximum delay in samples (results in longest frequency response)
mxdly = ceil(max(dlays)*fs);

% Find filter spectral spread to determine time support for
freqsp = sum(sfac.freq.*10.^(-sfac.atten/20))/sum(10.^(-sfac.atten/20));
mxsprd = 1/freqsp;  %  Reciprocal of frequency at spread point (in seconds)

impres = zeros(1,ceil(5*fs*mxsprd+mxdly)); %  Set up output array for max delay and 5 time constants of the filter

%  Loop through each delay and assign impulse response corresponding to
%  image/multipath scatterer at that delay
for k=1:length(dlays)
    d = dlays(k)*c;  %  Compute distance traveled by signal
    attn = d*sfac.atten; %  Compute log slope of spectrum
    dt = abs(5*mxsprd);  %  Spread in time domain (5 time constants)
    filn = round(dt*fs);
    %  Time domain spread less than LEN samples, set equal to LEN (length of
    %  sample spectra).
    if filn < len;     %  limit filter order to greater than 50 points (spectral definition)
        filn = len;
    elseif filn > fs/2  %  Limit filter order to less than a half second in worse case
        filn = fs/2;
    end
    %  if filter order (same as spread) not even, increase by one sample to
    %  make order even
    if filn/2 ~= fix(filn/2)
        filn=filn+1;
    end
    %  Compute frequency spectrum for attenuation
    alph = 10.^(-attn/20);
    b = fir2(filn,nfax,alph);        %  FIR filter to implement attenuation effects
    addely = dlays(k)-(length(b)-1)/(2*fs);  % assign corresponding to delay to middle of filter coefficient
    %  Filter coefficients after centering extend beyond the beginning of time reference, truncate filter response
    %  trim filter to start at time = 0
    if addely < 0
        %  Trim filter to fit the beginning of time sequence
        trm = round(-addely*fs)+1;
        b = b(trm:end);
        addely = 0;
    end
    dum = scls(k)*delayt(b,fs,addely);  %  Delay by travel time
    %  If delay plus filter length exceeds beyond end of accumulation array, extend
    %  the accumulation array
    if length(dum) > length(impres)
        impres = [impres, zeros(1, length(dum)-length(impres))];
    end
    % Assign impulse response to proper delay in accumulation array.
    impres(1:length(dum)) = impres(1:length(dum)) + dum;
end
%  Trim off trailing zeros
greal = find(impres ~= 0);
trmval = max(greal);
impres = impres(1:trmval);
%  If requested, create a time axis for impulse response
if nargout == 2
    t = [0:length(impres)-1]/fs;
end