function strms1 = signifslist_par(fovp, numWorkers, workerID)
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
%    sgtwv  =>  The file name of the multichannel wav file recording
%    rflag  =>  A flag that tells whether or not the region being searched should
%               be a rectangular prism or a sphere (sphere=1, rectangular prism~=1)
%    sources => A matrix of 3 columns where each row specifies a point
%               representing a detected sound source (used in srpfast.m
%               when searching a spherical region) [x,y,z] (meters)
%
% Output:
%  STRMS1   =>  Array of real numbers indicating all detected sources.  Each
%               row contains [time assoicated with midpoint of window (seconds),
%               detection statistic (unitless), x, y, z (coordinates in
%               meters)]
%             
%  Written by Kevin Donohue (donohue@engr.uky.edu)  April 2012
%  Modified by Kirstin Brangers                     July 2012
%  Modified by Paul Griffioen (pmg6@students.calvin.edu)  August 2013


strms1 = [];                    %  Initalize output array
sgtwv = [fovp.pn, fovp.fn];
[sizz, fso]= wavread(sgtwv,'size');       %  Size and Sample Rate of sound file
samps =  sizz(1);                          %  Signal length
nm = sizz(2);                               %  Number of mic channels
rflag = fovp.rflag;
sources = fovp.sources;

%  Create signal conditioning filter coefficients
[b,a] = butter(4,fovp.fc/(fso/2));  

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
tapwin = flattap(nwleno,20);     % This creates a tapering window to reduce edge effects (default 20%)
winrec = tapwin*ones(1,nm);     % nm is number of mics
[~, fnum] = size(fovp.sgrid); % Determine number of regions selected

%  Compute delays and weights for all points of interest
rstref = cell(fnum);
at = cell(fnum);
for kg=1:fnum
  [rstref{kg}, at{kg}] = voldelwts(fovp.sgrid{:,kg}, fovp.mp, fovp.c);
end

tst = 2.6;                        % Beginning of first window in seconds
ssta = (tst:fovp.tinc:(samps-nwleno-1)/fso);  % Starting points of sliding window (analog)
rad = fovp.threshneigh;         % Convert neighborhood in meters to samples
steprem= mod(length(ssta),numWorkers);
disp(workerID);
tic
for cnt=workerID:numWorkers:length(ssta)-steprem
    
%     disp([num2str(cnt) ' of ' num2str(length(ssta))]);
    ssti =  1 + round(fso*ssta(cnt));                    % Convert to index values
    sed = ssti+nwleno-1; % Window end
    [y,fso] = wavread(sgtwv, [ssti sed]);                % Read entire sound file
    
    sigout = filtfilt(b,a,y);
    sigout= sigout.*winrec;
    
    
    
    
    sigout= denoisefast(sigout, .0004);
    ye = abs(hilbert(sigout));
    ye = ye - ones(length(ye),1)*mean(ye,1);                %subtract avg
    
    %sigout_test= resample(sigout,3,4);
    
    %ye_resamp = abs(hilbert(sigout_resamp));
    %ye_resamp = ye_resamp - ones(length(ye_resamp),1)*mean(ye_resamp,1); %subtract avg
    %sigout_resamp = resample(ye,fovp.fs,fso);                % Resample to new sampling rate
    
    %compute gradient 1 column at a time
    sig_grad= zeros( size(ye) );
    for i=1:size(ye,2)
        sig_grad(:,i)= gradient(ye(:,i));
    end
    
    sig_adjust = sig_grad.*sigout;
%    sigout_white = whiten(y(1:nwlen,:).*winrec, fovp.beta);   % Beta parameter (~.8)
    %determine envelope of gradient-adjusted signal
    %ye_adjust = abs(hilbert(sig_adjust));
    %ye_adjust = ye_adjust - ones(length(ye_adjust),1)*mean(ye_adjust,1);
    
    timestmp = ssta(cnt)+fovp.tinc;                         % Timestamp
    
    %  Create SRP Image from processed perimeter array signal
    %im={};
    imadj={};
    %imye={};
    %imyeadj={};
    
    %sig_collect = { sigout; sig_adjust; sig_grad }; %ye can be added before sig_grad
    %im_collect = { im; imadj; imyeadj }; %imye can be added before imyeadj
    
    %Change first argument of srpfast to whichever signal you want to
    %process (NOTE: the envelopes seem to be useless. use sig_grad or
    %sig_adjust)
    if (fnum==1) && (rflag~=1)
        kg=1;
        imadj= srpfast(sig_adjust, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg}, rflag, kg, sources, fnum);
    end
    
    %This is the slowest part by far
%     if (fnum==1) && (rflag~=1)      % If only the overall region is being searched, search the number of regions using a normal for loop
%         kg=1;
%             im{kg} = srpfast(sigout, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg}, rflag, kg, sources, fnum);
%             %         im = srpframenn(sigout, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez);
%             imadj{kg} = srpfast(sig_adjust, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg}, rflag, kg, sources, fnum);
%             %determine images for the envelopes
%             %imye{kg} = srpfast(ye, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg}, rflag, kg, sources, fnum);
%             imyeadj{kg} = srpfast(sig_grad, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg}, rflag, kg, sources, fnum);
%     else                            % If many regions are being searched, parallelize the regions being searched
% %         parfor kg=1:fnum
%         for kg=1:fnum
%             im{kg} = srpfast(sigout, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg}, rflag, kg, sources, fnum);
%             %         im = srpframenn(sigout, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez);
%         end
%     end
    
%     mags={};
%     pos={};
     magsadj={};
     posadj={};
    %magsye={};
    %posye={};
    %magsyeadj={};
    %posyeadj={};
    
    for kg=1:fnum
        %[mags{kg},pos{kg}] = dent3dweiyx(im_collect{1}{kg}, rad ,fovp.pcfar);        %  Detect significant sources
        %[magsg{kg},posg{kg}] = dent3dweiyx(im_collect{2}{kg}, rad, fovp.pcfar);
        %[magsye{kg},posye{kg}] = dent3dweiyx(imye{kg}, rad ,fovp.pcfar);
        [magsadj{kg},posadj{kg}] = dent3dweiyx(imadj, rad ,fovp.pcfar);
    end
    for kg=1:fnum
    for so=1:length(magsadj{kg})
            strms1 = [strms1; timestmp, magsadj{kg}(so), fovp.sgrid{kg}{1}(posadj{kg}(so,1)), fovp.sgrid{kg}{2}(posadj{kg}(so,2)), .0254];
       
            if size(strms1,1)/100 ==  round(size(strms1,1)/100)
                str= sprintf( 'dummy_workers_%i', numWorkers);
                str2= sprintf( '_ID_%i', workerID);
                str= [str,str2];
                save ( [str '.mat'], 'strms1');   %cannot be called in a parfor loop
                str= sprintf('Detection number %f', size(strms1,1));
                disp(str);
                toc;
            end
    end
    end
        
    %  Create Acoustic SRP Image Plot
    if fovp.graphicout_flag == 1
%         figure(1+kg)
%         imagesc(fovp.sgrid{kg}{1},fovp.sgrid{kg}{2},sum(im_collect{1},3))
%         colormap(jet); colorbar; axis('xy')
%         title(['Original Signal: Number ' int2str(cnt) ' of ' int2str(length(ssta))])
%         % axis([froom(1,1)-.25, froom(1,2)+.25, froom(2,1)-.25, froom(2,2)+.25])
%         pause(.1)
% 
%         figure(2+kg)
%         imagesc(fovp.sgrid{kg}{1},fovp.sgrid{kg}{2},sum(im_collect{2},3))
%         colormap(jet); colorbar; axis('xy')
%         title(['Gradient-Adjusted Signal: Number ' int2str(cnt) ' of ' int2str(length(ssta))])
%         pause(.1);

%         figure(3+kg)
%         imagesc(fovp.sgrid{kg}{1},fovp.sgrid{kg}{2},sum(imye{kg},3))
%         colormap(jet); colorbar; axis('xy')
%         title(['Original Signal Envelope: Number ' int2str(cnt) ' of ' int2str(length(ssta))])
%         % axis([froom(1,1)-.25, froom(1,2)+.25, froom(2,1)-.25, froom(2,2)+.25])
%         pause(.1) 

        figure(4+kg+workerID)
        imagesc(fovp.sgrid{kg}{1},fovp.sgrid{kg}{2},sum(imadj,3))
        colormap(jet); colorbar; axis('xy')
        title(['Signal Gradient: Number ' int2str(cnt) ' of ' int2str(length(ssta))])
        % axis([froom(1,1)-.25, froom(1,2)+.25, froom(2,1)-.25, froom(2,2)+.25])
        pause(.1) 
        %set breakpoint at blah for strms1 analysis
%         if ~isempty((magsadj{kg}))
%             blah=0;
%         end
    end
    
end
  
end

% [ytest,fso] = wavread([fovp.pn, fovp.fn]);                % Read entire sound file
% [aligny, delays, weights] = cuealign(ytest,fovp);
% save(['CueAlignResults.mat'],'aligny','delays','weights');

