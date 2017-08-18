function sigout = tfmask(ina,soi,dt,dffac,fs)
% 
%  This function inputs an array of beamformed signals, creates a time
%  frequency mask for each one and applies it.
%
%     sigout = tfmask(ina, soi, dt, dffac, fs)
%
% Inputs: 
%   INA   => Columnwise array of beamformed signals
%   SOI   => Index for column in INA to identify speaker of interest
%   DT    => Delta time for time-frequency window (seconds)
%   DFFAC => Integer denoting the increase Delta frequency
%                   over the maximum possible based on DT. 
%                   DFFAC = 2 will result in padding with zeros 
%                   to double the grid point resolution.
%   FS    => Sampling Rate
%          
% Required Functions:
%               getMsk.m
%
% Written by Kevin D. Donohue (donohue@engr.uky.edu)    October 14 2011
% Updated August 2014 to include a different windowing for masker creation

wlen = dt;          % Window length for STFT (in seconds)
windinc = .25;       % Increment between overlapping windows in percent of window length 
zpadfac = dffac;    % Scale on length to zero pad for STFT                
gamma =1;           % Threshold for Binarizing the mask

%  Beamformed SOI
bsig = ina(:,soi);
siglen = length(bsig); 
%  Beamformed Interferers
bnos = [ina(:,1:(soi-1)) ina(:,(soi+1):end)];

%  Compute increment and tappering window
len = round(fs*wlen);
inc = (len-1)*windinc;                  % Increment in time
%  Symmetric flat window with cosine taper
tapwin = tukeywin(len,2*windinc);  
tapwin = (tapwin)/((1/windinc)-1);      % Scale window so overlapping window always add to 1
nfft = 2^(nextpow2(zpadfac*len));       % Round zero padding length upto next power of 2.
sigout = zeros(length(bsig)+nfft,1);    % Initialize output signal
tapwins = kaiser(len,12);  %  Stonger taper for masker
%  Increment through signals for STFT and masking
stp = 1;        % Starting Sample, integer
stpa = stp;     % Initialize starting sample without rounding
edp = len;      % Ending sample, integer
hlen =  floor(nfft/2);      %  Half the FFT length for mask development
while edp <= siglen
    % Window and take FFT of SOI signal at half FFT length
    sigspec = fft(bsig(stp:edp,1).*tapwins,hlen);
    Mask = ones(hlen,1);        % Initialize mask
    % Compute mask from interferers 
    for ii = 1:size(bnos,2) 
        % Window and take FFT of SOI signal at half FFT length
        sp = fft(bnos(stp:edp,ii).*tapwins, hlen);
        tMsk = getMsk(sigspec,sp,gamma);
        Mask = Mask .* tMsk;    % And binary mask values
    end
   
    %  Take mask into time domain for zero padding 
    timeMask = ifft(Mask);     % Make mask complex, magnitude 1
    %  Take mask back into frequency domain with zero padding at length NFFT
    Mask2 = fft(timeMask, nfft);
   
    %  Compute SOI frequency component with zero padding
    sigspec = fft(bsig(stp:edp,1),nfft);
    specout = sigspec.*Mask2;   % Apply mask
    %  Get time domain masked signal
    sigadd = real(ifft(specout,'symmetric'));
    %  Accumulate with overlap and zero padding
    sigout(stp:stp+len-1,1) = sigout(stp:stp+len-1,1) + sigadd(1:len,1).*tapwin;
    stpa = stpa+inc;       % Increment (not necessarily integer)
    stp = round(stpa);     % Round off to compute starting point of next segment 
    edp = stp+len-1;       % Compute end of next segment in STFT analysis
end
%  When finished, write masked signal
sigout = real(sigout);





 