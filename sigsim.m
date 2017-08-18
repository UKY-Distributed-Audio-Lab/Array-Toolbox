function [sigout, tax, td] = sigsim(f1,f2,winlen,fs,pcttap, sigpos, mpos, flims, c)
% This function will create an array of impulse-like signals associated
% with a sound eminating from spatial locations and received at an array of
% microphones.  Each microphone corresponds to a column in the output
% matrix "sigout."  ALso in the output is the corresponding time axis
% "tax" and time delays "td"
% 
%  [sigout, tax, td] = sigsim(f1,f2,winlen,fs,pcttap, sigpos, mpos, flims, c)
%
% The input arguments are sound source descriptors:
% f1 => Starting frequent of sound source
% f2 => Ending Frequency of sound source
% winlen => length of time window (should be at least 1/(f2-f1))
% pcttap => percentage of signal bandwitdth that tapers between 0 and 1 at
%           frequency edges
% fs =>  Signal sampling frequency
% sigpos => a matrix where each column represents the coordinates of a
%           signal source and each row is the dimension x,y, & z(if
%           present)
% mpos => a matrix where each column represents the coordinates of the mic
%         and each row is the dimension x,y, & z (if present)
% flims => a 2 column matrix with opposite corners over the field consider
%          for the presence of the sound source (note that the coordinates
%          of flims, mpos, and sigpos must be in the same frame of
%          reference).
%    c  => Speed of sound (FYI, c = 331.4+0.6*Temperature in Centigrade)


[fr, fc] = size(flims);  %  Get rows "fr" to determine how many dimensions are specified
fd = abs(flims(:,1)-flims(:,2));  % Compute Field distances for each dimension
% Obtain mic array information
[mr, mc] = size(mpos);  % Determine number of mics = mc
[sgr, sgc] = size(sigpos);  % Determine number of sound sources = sgc
rez = c/(2*fs); % Pick a superesolution value for field grid
%  Find max and min delays for all points in test field
%  Create Grid for all axes
for k=1:fr
    alen = ceil(fd(k)/rez);  %  Find number of point in axis
    st = min(flims(k,:));    %  Find starting point of axis
    ax{k} = st + (max(flims(k,:)-st))*[0:alen]/alen;  % Create axis points
end
% Find the minimum delay for between field point and mic positions
for k=1:mc
    for n=1:fr
       dum = min(abs([ax{n}- mpos(n,k)]));
       minax(n) = dum(1);
    end
    mdk(k) = norm(minax,2); % Find minimum distances for every mic 
end
mndf = min(mdk)/c;  %  Find minimum distance over all mics
% Find the maximum delay for between field point and mic positions
[gr, gc] = size(ax);
for k=1:mc
    for n=1:fr
       dum = max(abs([ax{n}- mpos(n,k)]));
       maxax(n) = dum(1);
    end
    mdk(k) = norm(maxax,2);  %  Find maximum distance for every mic
end
mxdf = max(mdk)/c;  % Find maximum distance over all mics

tdur = mxdf-mndf;  %  Find time difference for sound to travel between
                   %   minimum and maximum distnace in field to mic
                   %   array
% Compute time axis to include field travel time and time windows for every signal 
tax = [floor((mndf-tdur)*fs):ceil((mxdf+tdur+winlen)*fs)]/fs;
%  Create signal from frequency domain parameters
rng = pcttap*(f2-f1);  %  Find taper distance on frequency axis
tdur = max(tax)-min(tax);  %  Compute time duration for full axis
%  Compute window length for FFT processing
nwin = round(tdur*fs);  
%  Compute power of 2 bigger than twice the window length with zero padding 
nfft = 2;
while 2*nwin > nfft
    nfft=nfft*2;
end
fax = fs*[-nfft/2:nfft/2-1]/nfft;  %  Compute symeteric frequency axis
% Initialize signal vector
ft = zeros(size(fax));
indp = find((fax< f2) & (fax >= f1));  %  Find frequency range with non-zero values
indn = find((fax <= -f1) & (fax > -f2));
ft(indp) = 1;                          %  Assign bandwith a nonzero value
ft(indn) = 1;
%  Find taper ranges for smoothing square edges
tapl = find(f1 <= fax &  fax <= f1+rng);
taph = find(f2-rng <= fax & fax <= f2); 
tapwin = hanning(2*length(tapl));
ft(tapl) = tapwin(1:length(tapl));
tapwin = hanning(2*length(taph));
ft(taph) = tapwin(length(taph)+1:end);

taph = find(-f1 >= fax &  fax >= -f1-rng);
tapl = find(-f2+rng >= fax & fax >= -f2); 
tapwin = hanning(2*length(taph));
ft(taph) = tapwin(length(taph)+1:end);
tapwin = hanning(2*length(tapl));
ft(tapl) = tapwin(1:length(tapl));

sigout = zeros(length(fax),mc);
%  Loop for each sound source
for ks = 1:sgc
        %  Find the distance of the sound source from every mic
        for k=1:mc
           md(k) = norm((sigpos(:,ks) - mpos(:,k)),2);
        end
        td = md/c;  %  Find corresponding time delays

        td = td - min(tax);
        for k=1:mc
            phs = exp(j*2*pi*fax*td(k));
            ftdum = 10^(-6*log2(td(k))/20)*ft.*phs;
            sigout(:,k) = sigout(:,k) + real(ifft(ifftshift(ftdum)'));
        end
end
sigout = sigout(1:length(tax),:);
td = td + min(tax);