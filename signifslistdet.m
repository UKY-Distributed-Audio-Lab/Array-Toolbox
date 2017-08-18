function strms1= signifslistdet(det, filename, nwlen,fovp, fsb, fsresample)

% plots the windows corresponding to predetermined detections DET in audio file
% FILENAME. NWLEN is the window length in samples

               %  Initalize output array
sgtwv = [fovp.fn];
[sizz, fso]= wavread(sgtwv,'size');       %  Size and Sample Rate of sound file
%samps =  sizz(1);                          %  Signal length
nm = sizz(2);                               %  Number of mic channels
rflag = fovp.rflag;
sources = fovp.sources;

%  Create signal conditioning filter coefficients
[b,a] = butter(2,fovp.fc/(fso/2)); 
if nargin <6
    fsresample= 44.1e3;
    if nargin <5
        fsb= 5e3;
    end
end
[bb,aa]= butter(6,fsb/fso/2);
fso= fovp.fs; %just to get the frequency fs
%  Compute window length in samples for segmenting time signal 
%winleno = ceil(fso*fovp.trez);
%  Determine additional signal added to window length to ensure delays
%  do not shift zeros into the field
%prs = mposanaly(fovp.mp,2);
fovp.fs=fso;
%  Maximum delay in seconds needed to synchronize in Delay and Sum Beamforming
%maxmicd = max(prs(:,3));

%  Extra delay time for padding window of data
%textrao = ceil(fso*maxmicd(1)/fovp.c);
%nwleno = winleno+textrao;
%  Windowing function to taper edge effects in the whitening process
tapwin = flattap(nwlen*fso,20);     % This creates a tapering window to reduce edge effects (default 20%)
winrec = tapwin*ones(1,nm);     % nm is number of mics
[~, fnum] = size(fovp.sgrid); % Determine number of regions selected

%  Compute delays and weights for all points of interest
rstref = cell(fnum);
at = cell(fnum);
for kg=1:fnum
  [rstref{kg}, at{kg}] = voldelwts(fovp.sgrid{:,kg}, fovp.mp, fovp.c);
end

tst = 1;                      % Beginning of first window in terms of index of det
rad = fovp.threshneigh;         % Convert neighborhood in meters to samples


for i=tst:size(det,1)
    bseg=floor(det(i,1)*fso);
    eseg=bseg+nwlen*fso-1;
    [y,fso]= wavread(filename, [bseg,eseg]);
    if fovp.graphicout_flag
        figure(1); plot(y);
    end
    %pause(.1);
    
    %now perform robust srp calculations
    sigout= y;
    %sigoutpera = sigoutper + nosamp*nosoutper + nos*sclnos*perpkpr;
    %sigout = whiten(sigout.*winrec, fovp.beta);
    sigfilt= [45e3 75e3];
    [sigout, ye, sig_grad, sig_adjust] = preprocess(sigout, winrec, fovp, fso, fsresample, sigfilt);
    toprocess= sig_adjust; %change this to change the signal used in srp
    if fovp.graphicout_flag
        figure(2); plot(toprocess);
    end
    %sig_adjust= whiten(sig_adjust, fovp.beta);
    timestmp = det(i,1);                         % Timestamp
    
    %  Create SRP Image from processed perimeter array signal
    imadj={};
    %Change first argument of srpfast to whichever signal you want to
    %process (NOTE: the envelopes alone seem to be useless. use sig_grad or
    %sig_adjust)
    if (fnum==1) && (rflag~=1)
        kg=1;
        %imadj= srpfast(sig_adjust, fovp.sgrid{:,kg}, fovp.mp, fovp.fs, fovp.c, fovp.trez, rstref{kg}, at{kg}, rflag, kg, sources, fnum);
        imadj= srpframenn(toprocess, fovp.sgrid{:,kg}, fovp.mp,fsresample, fovp.c, fovp.trez);
        %imadj= srpframenn(ye, fovp.sgrid{:,kg}, fovp.mp,fsresample, fovp.c, fovp.trez);
    end

     magsadj={};
     posadj={};
    
    for kg=1:fnum
        [magsadj{kg},posadj{kg}] = dent3dweiyx(imadj, rad ,fovp.pcfar);
    end
    strms1 = [];
    for kg=1:fnum
    for so=1:length(magsadj{kg})
            strms1 = [strms1; timestmp, magsadj{kg}(so), fovp.sgrid{kg}{1}(posadj{kg}(so,1)), fovp.sgrid{kg}{2}(posadj{kg}(so,2)), .0254];
       
            if size(strms1,1)/100 ==  round(size(strms1,1)/100)
                save dum2_sig_adjust.mat strms1   %cannot be called in a parfor loop
                str= sprintf('Detection number %f', size(strms1,1));
                disp(str);
                toc;
            end
    end
    end
        
    %  Create Acoustic SRP Image Plot
    if fovp.graphicout_flag == 1

        figure(5+kg)
        surf(imadj);
        figure(6+kg); imagesc(fovp.sgrid{kg}{1},fovp.sgrid{kg}{2},sum(imadj,3))
        
        colormap(jet); colorbar; axis('xy')
        title(['Signal Gradient: Number ' int2str(i) ' of ' int2str(size(det,1)) ' timestamp (s): ' num2str(timestmp,5)])
        % imaging for the video
        hold on;
        draw_circle(0,0, .2096, 'black');
        scatter( fovp.mp(1,:), fovp.mp(2,:), 60, 'blue');
        %now plot some blobs over the detections
        if ~isempty(strms1)
            holder= find(strms1(:,1)==timestmp,1);
            if ~isempty(holder)
                scatter( strms1(holder:end, 3), strms1(holder:end, 4), 80, 'white', 'filled');
                scatter( strms1(holder:end, 3), strms1(holder:end, 4), 20, 'red', 'filled');
            end
        end
        hold off;
        %frame= getframe;
        %writeVideo(writerObj, frame);
        %end imaging
        
        %pause(.1); 
        %set breakpoint at blah for strms1 analysis
%         if ~isempty((magsadj{kg}))
%             blah=0;
%         end
    end
%     pause
end