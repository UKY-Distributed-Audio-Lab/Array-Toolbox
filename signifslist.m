function strms1 = signifslist(fovp)
%  This function inputs a multichannel recording from a distributed
%  microphone system and parameters associated with the microphone
%  positions and the FOV.  These parameters are located in the data
%  structure FOVP.
%
%     strms1 = signifslist(fovp)
%
% Inputs:
% Field of FOVP:
%    mp     =>  Microphone positions, each column corresponds to a microphone [x,y,z]' in meters
%    sgrid  =>  2-dimensional cell array where each row contains x,y,z grid points of each FOV block
%    c      =>  Speed of sound
%    fs     =>  Sampling frequency for processing signal
%    pn     =>  string with path to microphone recording wave file
%    fn     =>  file name of microphone recroding wave file
%    fc     =>  High-pass filter cut-off for processin signal
%
% Optional fields:
%    trez   =>  Time window length (resolution) in seconds for detecting each block. Default: 20 ms
%    tinc   =>  Time increment in seconds.  Default: trez/2 (1/2 of time window length)
%    beta   =>  Partial whitening parameter (between 0 and 1).  Default: 0.85
%    pcfar  =>  False positive detection probability for threshold
%               computation.            Default: 1e-5
%    shape  =>  Shape parameter for Weibull distribution used in CFAR
%               threshold computation.  Default: 1.2
%    threshneigh =>  Vector for neighborhood used in threshold computation, plus/minus
%                    test point in meters.
%    graphicout_flag =>  if present and set to 1, acoustic images collapsed
%                        along the z-axis for each sub region will appear
%
% Output:
%  STRMS1   =>  Array of real numbers indicating all detected sources.  Each
%               row contains [time assoicated with midpoint of window (seconds),
%               detection statistic (unitless), x, y, z (coordinates in
%               meters)]
%
%  Written by Kevin Donohue (donohue@engr.uky.edu)  April 2012
%  Modified by Kirstin Brangers                     July 2012


strms1 = [];                    %  Initalize output array

[sizz, fso]= wavread([fovp.pn fovp.fn],'size');       %  Size and Sample Rate of sound file
samps =  sizz(1);                          %  Signal length
nm = sizz(2);                               %  Number of mic channels

%  Create signal conditioning filter coefficients
[b,a] = butter(4,fovp.fc/(fovp.fs/2),'high');

%  Check on time window parameter
if ~isfield(fovp,'trez')
    fovp.trez = .02;           % Default: 20 ms
end

%  Check on time increment parameter
if ~isfield(fovp,'tinc')
    fovp.tinc = fovp.trez/2;    % Default: trez/2 = 10 ms
end

%  Check on whitening parameter must be in [0, 1]
if ~isfield(fovp,'beta')
    fovp.beta = .85;          % Default: 0.85
end
if fovp.beta < 0 || fovp.beta > 1
    disp('Partial whitening factor out of range. Setting to default value!')
    fovp.beta = .85;            % Default: 0.85
end

%  Check on probability of false alarm
if ~isfield(fovp,'pcfar')
    fovp.pcfar = 1e-5;         % Default:  1e-5
end

%  Check on graphic output flag
if ~isfield(fovp,'graphicout_flag')
    fovp.graphicout_flag = 0;         % Default:  0, no images
end


%  Check on shape parameter for probability of false alarm computation
if ~isfield(fovp,'shape')
    fovp.shape = 1.2;          % Default: 1.2
end

%  Check on neighborhood limit around test for computing CFAR threshold
%  plus/minus in meters around test point
if ~isfield(fovp,'threshneigh')
    fovp.threshneigh = [7, 7, 1];  % Default: [7, 7, 1]
end



%  Step through each FOV region to determine locations of significant sound sources.
%  Create supporting data and vectors for performing SSL on all locations.

%  Compute window length in samples for segmenting time signal
winlen = ceil(fovp.fs*fovp.trez);
winleno = ceil(fso*fovp.trez);
%  Determine additional signal added to window length to ensure delays
%  do not shift zeros into the field
prs = mposanaly(fovp.mp,2);

%  Maximum delay in seconds needed to synchronize in Delay and Sum Beamforming
maxmicd = max(prs(:,3));

%  Extra delay time for padding window of data
textra = ceil(fovp.fs*maxmicd(1)/fovp.c);
textrao = ceil(fso*maxmicd(1)/fovp.c);
nwlen = winlen+textra;          % New window length in samples
nwleno = winleno+textrao;
%  Windowing function to taper edge effects in the whitening process
tapwin = flattap(nwlen,20);     % This creates a tapering window to reduce edge effects
winrec = tapwin*ones(1,nm);     % nm is number of mics
[dum, fnum] = size(fovp.sgrid); % Determine number of regions selected

%  Compute delays and weights for all points of interest
for kg=1:fnum
    [rstref{kg}, at{kg}] = voldelwts(fovp.sgrid{:,kg}, fovp.mp, fovp.c);
end

tst = 0;                        % Beginning of first window in seconds
ssta = [tst:fovp.tinc:(samps-nwleno-1)/fso];  % Starting points of sliding window (analog)
rad = fovp.threshneigh;         % Convert neighborhood in meters to samples

for cnt=1:length(ssta)
    ssti = 1 + round(fso*ssta(cnt));                    % Convert to index values
    sed = ssti+nwleno; % Window end
    [y,fso] = wavread([fovp.pn, fovp.fn], [ssti sed]);                % Read entire sound file
    y = resample(y,fovp.fs,fso);                % Resample to new sampling rate
    y = filtfilt(b,a,y);
    sigout = whiten(y(1:nwlen,:).*winrec, fovp.beta);   % Beta parameter (~.8)
    timestmp = ssta(cnt)+fovp.tinc;                         % Timestamp
    
    %  Create SRP Image from processed perimeter array signal
    for kg=1:fnum
        im = srpfast(sigout, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg});
        [mags,pos] = dent3dwei(im, rad ,fovp.pcfar);        %  Detect significant sources
        for so=1:length(mags)
            strms1 = [strms1; timestmp, mags(so), fovp.sgrid{kg}{1}(pos(so,1)), fovp.sgrid{kg}{2}(pos(so,2)), fovp.sgrid{kg}{3}(pos(so,3))];
        end
        
        %  Create Acoustic SRP Image Plot
        if fovp.graphicout_flag == 1
            figure(1+kg)
            imagesc(fovp.sgrid{kg}{1},fovp.sgrid{kg}{2},sum(im,3))
            colormap(jet); colorbar; axis('xy')
            title(['Number ' int2str(cnt) ' of ' int2str(length(ssta))])
            pause(.1)
        end
    end
end

end

