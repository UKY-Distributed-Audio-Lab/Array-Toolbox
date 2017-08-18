function mkmovwdn(mldum,xax,yax,pos,frat,endrat,fname, rngout)
%  This function will take a 3-D array where the 3rd dimension is
%  time and it will make it into a movie.  It will perform various
%  filtering, image processing, and rate conversion according to
%  the input parameters
%
%     mkmovwdn(mldum,xax,yax,pos,frat,endrat,fname, rngout)
%
%  MLDUM =  3-D input matrix of movie frames, the 3rd dimension is time
%  XAX   =  axis points for the lateral positions
%  YAX  = axis point for the depth positions
%  POS = 4 element row vector for the postion and size of
%        the figure window in terms of screen pixels used to
%        make the movie [start_column, start_row, width, height]
%        where the start_ coordinates refers to the lower right
%        corner of the figure window.
%  FRAT = Frame rate in frames per second for frames in MLDUM
%  ENDRAT = Frame rate in frames per secone for the final image
%  FNAME = File name for stored AVI movie file
%  RNGOUT = (Optional) 2 element array with dynamic range for image
%           intensities [minimum, maximum].  Any values above or below
%           these values will be clipped.
%
%    Written by Kevin D. Donohue August 2005
%



%  Determine Size of 3-D acoustic movie
[r, c, n] = size(mldum);
%  Initalize AVI object
mv = avifile(fname, 'FPS', endrat, 'Compression', 'none');    

%  Create averaging kernal to smooth images over time
%   before subsampling down to ending frame rate
crat = frat/endrat;   %  framerate conversion factor
if crat > 1   %  Need to downsample to decrease rate
    winavg = ceil(3*crat);  %  Round window length up to next integer
    hw = kaiser(winavg,2);  %  Compute tapered window for averaging
    bwin = zeros(1,1,winavg); 
    % create averaging kernal in 3-dimensions
    for k =1:winavg
        bwin(1,1,k) = hw(k);  %  Average only along the time axis
    end
    bwin = bwin/(sum(abs(hw)));  % Normalize averaging kernal energy
else  %  if final rate is the same or higher than input, no smoothing needed.
    bwin(1,1,1) = 1;
end

%  Apply smoothing kernal before subsampling
mldum = convn(mldum,bwin,'same');

%  If range parameters for grayscale not given, use max and min of data 
if nargin < 8
    rngout(1) = min(min(min(mldum)));
    rngout(2) = max(max(max(mldum)));
end

% Initalize figure object
fig=figure;
%  Set figure properties for creating movie frames
set(fig,'DoubleBuffer','off','Position',pos);
set(gca,'xlim',[xax(1) xax(end)],'ylim',[yax(1) yax(end)],...
   'NextPlot','replace','Visible','off')
%  Assign color map and axes labels
colormap(gray)
axis('xy')
xlabel('meters')
ylabel('meters')
mcnt = 0;  %  initalize frame counter for ending fps image
%  Loop through each frame at ending rate 
for k=1:crat:n  % (note k is not necessarily an integer)
      gdum = squeeze(mldum(:,:,round(k)));  % remove singleton dimension
      imagesc(xax, yax, gdum, [rngout(1) rngout(2)])  %  Create picture
      mcnt = mcnt+1;  %  Increment final frame counter
      mmv = getframe(gcf);  % get movie frame from figure
      mv = addframe(mv,mmv);  %  add to movie file
end
mv = close(mv);  %  close and save movie object