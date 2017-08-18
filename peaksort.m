function [dettruepks, locerr, tarmag, nospks, nosmag, tnstats] = peaksort(im,ax,ay,targetloc, pltflag)
%  This function inputs a reconstructed acoustic image IM, along with the corresponding axes
%  points AX and AY for the columns and rows, respectively, the known postions of the
%  target locations, TARGETLOC, and an option figure plot flag PLTFLAG
%
%     [dettruepks, locerr, tarmag, nospks, nosmag, tnstats]  = peaksort(im,ax,ay,targetloc,pltflag)
%  
%  The input PLTFLAG can be any integer value. If it is present a figure of the raw acoustic
%  image will be generated with superimposed target locations (red x's), detected target
%  locations (green o's), and other peaks not associated with the target location (yellow circles)
%  representing false alarms peaks for a threhsold set at 0.
%  Outputs:
%  DETTRUEPKS - matrix of the detected locations of the target
%  LOCERR     - the distance between the actual target and detected target
%  TARMAG     - magnitude at detected location
%  NOSPKS     - matrix of detected locations of noise peaks representing
%               false alarms
%  NOSMAG     - matrix of noise peaks magnitudes
%  TNSTATS    - vector of [target peaks mean, noise mean, std target peaks, std noise] 
%
%
%
% Sorting out a target peak from a false alarm peak is somewhat subjective.  The rules that this
% program uses is that a peak is considered a target if it is the nearest peak to the actual target
% location and any of the following hold:
%  1. It is in an adjacent pixel location to the true target, 
%  2. For a line between the target position the detected peak, the ratio between the smallest peak
%     height along the line and the height of the considered peak is
%     greater than 6 dB (assuming IM are coherent power values, not
%     amplitudes), and the smallest peak height along the line is positive.
% 
% False alarm peaks are determined as those peaks not associated directly with the target location.
% This excludes peaks that fall within the resolution area of the detected target.  The rules are
% peak is considered a false alarm peak if any of the following hold
%   1.  For a line between the closest detected target position (determined through the rules previously
%       stated) and the detected peak, the ratio between the smallest peak height along the line and the
%       maximum height of either the peak at the target position or the considered peak is
%       greater than 6 dB, or the smallest peak height along the line is negative.
%
%   Kevin D. Donohue (donohue@engr.uky.edu)  September 2005
%
% This function needs the functions "peakfind2d.m", "nearestpks.m", and "putpeaks.m" to run.



dettruepks = [];
locerr = [];
tarmag = [];
nospks = [];
nosmag =[];

%  Determine diagonal pixel width for testing closness
dpixwidth = sqrt((ax(2)-ax(1))^2+(ay(2)-ay(1))^2);
   
%  Find all peaks in image
[pxl,pyl,pml] = peakfind2d(im,ax,ay);
%  Trim the set of detected peaks to only those with positive amplitudes
trm = find(pml > 0);
pxl = pxl(trm);
pyl = pyl(trm);
pml = pml(trm);
peaks = [pxl; pyl];
%  Search for targets at each peak location
[dims, tarnum] = size(targetloc);  %  how many actual targets are present
[dims, pksnum] = size(peaks);      %  how many peaks must be searched through
tarcnt = 0;  %  Initalize Target Counter
for k=1:tarnum
    %  Zero out image at target locations for noise statistics:
    [dummin, xtl(k)] = min(abs(targetloc(1,k)-ax));
    [dummin, ytl(k)] = min(abs(targetloc(2,k)-ay));
    
    %  Find nearest peak to target location and characterize distance
    [ip, perr, valie, tval, pval]= nearestpk(peaks, targetloc(:,k), im, ax, ay);
    %  Determine if closest peak can be considered a target.
    %   If the peak in an adjacent pixel it is a target classify it as a
    %   target, or if the pixels values between the target and detected 
    %   peak do not drop below 6 dB
    if perr < 2*dpixwidth   % In adjacent pixel, consider it a target
        %  Build target array
        tarcnt=tarcnt+1;
        dettruepks(:,tarcnt) = peaks(:,ip);
        tarmag(tarcnt) = pml(ip);
        locerr(tarcnt) = perr;
        %  Exclude peak from further analysis
        if pksnum > 0
            peaks = peaks(:,[1:(ip-1),(ip+1):pksnum]);
            pml = pml(:,[1:(ip-1),(ip+1):pksnum]);
            [dims, pksnum] = size(peaks);    
        end           
    elseif  (valie >= eps) & (10*log10(pval/(valie+eps)) <= 6)  % if line is above zero
                                                       %and does not dip greater than 6 db
        %  Build target array
        tarcnt=tarcnt+1;
        dettruepks(:,tarcnt) = peaks(:,ip);
        tarmag(tarcnt) = pml(ip);
        locerr(tarcnt) = perr;
        %  Exclude peak from further analysis
        if pksnum > 0
            peaks = peaks(:,[1:(ip-1),(ip+1):pksnum]);
            pml = pml(:,[1:(ip-1),(ip+1):pksnum]);
            [dims, pksnum] = size(peaks);  
        end
   
    else 
        tarcnt=tarcnt+1;
        dettruepks(:,tarcnt) = [NaN; NaN];
        tarmag(tarcnt) = 0;
        locerr(tarcnt) = NaN;
    end
        
end
mtar = nanmean(tarmag);
sstdtar = nanstd(tarmag);

imnos = zeroout(im,targetloc, ax, ay);
mnos = nanmean(imnos);
sstdnos = nanstd(imnos);
tnstats = [mtar; mnos; sstdtar; sstdnos];


%  Examine peaks for false alarms, need to trim off secondary peaks
%  still on target rolloff
facnt = 0;  %  Initalize false alarm peak counter
dettar = find(tarmag > 0); %  Find detected targets
%  If undetected targets present, trim target array to only detected
%  targets
if ~isempty(dettar) & ~isempty(peaks)
    detnz = dettruepks(:,dettar);
    for k=1:pksnum
        %  Find nearest peak to detected noise location near the target and characterize
        %  distance
        [ip, perr, valie, tval, pval]= nearestpk(detnz, peaks(:,k), im, ax, ay);
        %   Determine if peak is close enought to a target to be considered part of it.
        %   If the peak pixels values between the target and detected 
        %   peak do not drop below 6 dB it is considered an independent peak
        if ~isempty(valie)
        if valie <= eps %  if line becomes negative, peak is sepearate from target
            %  Build false alarm array
            facnt=facnt+1;
            nospks(:,facnt) = peaks(:,k);
            nosmag(facnt) = pml(k);
        elseif 10*log10(max([pval,tval])/(valie+eps)) >= 6  % If line deviation 6db from either
                                                  %  end point consider them
                                                  %  separate
            %  build false alarm array
            facnt=facnt+1;
            nospks(:,facnt) = peaks(:,k);
            nosmag(facnt) = pml(k);
        end
        end
    end
elseif  ~isempty(peaks) %  if no target were detected then everything is a false alarm
    nospks = peaks;
    nosmag = pml;
else
    nospks = NaN;
    nosmag = 0;
end

%  If figure flag is provided, plot image with target and false alarm peaks 
if nargin == 5
    figure(pltflag)
    %  Create SRP image
    imagesc(ax,ay, im, [min(min(im)), max(max(im))]); colormap(gray); axis('xy')
    %  plot actual target positions
    if ~isempty(targetloc)
       putpeaks(pltflag,targetloc(1,:),targetloc(2,:),'ow');
    end
   % putpeaks(pltflag,dettruepks(1,:),dettruepks(2,:),'og');
   %putpeaks(pltflag,nospks(1,:),nospks(2,:),'oy');
end