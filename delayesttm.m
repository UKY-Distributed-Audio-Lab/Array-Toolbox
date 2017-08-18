function [d, cv, pro] = delayesttm(sig1,sig2,fs, pro)
%  This function estimates the time delay between 2 similar signal segments
%  (SIG1 and SIG2) of the same length and sampled on the same time reference
%  axis with rate (FS):
%
%       [d, cv, pro] = delayesttm(sig1,sig2,fs,pro)
%
%  The function works best when delays are small with respect to the segment
%  length (i.e. up to about 20% of segment length).  This version of the delay
%  estimator uses the unbiased cross correlation between the segments and
%  provides signal processing options such as detrending and clipping to reduce
%  influence of noise and other artifacts on the estimate.
%  The function inputs are:
%  SIG1  => Signal segment with assumed delay
%  SIG2  => Reference signal segment.  SIG1 and SIG2 must be the same size
%  FS    => Sampling frequency
%  PRO   => (optional) Data Structure with Signal processing options
%       Fields:  cliptype   =>  String indicating the type of clipping to be
%                               applied 'cnsc' = center soft clip
%                                       'cnhc' = center hard clip
%                                       'cnhl' = center hard limit
%                                       'none' = no clipping (default)
%                cliplevel  => number between 0 and 1 for computing the clip
%                              level a percent of the maximum absolute signal
%                              value.  Default is .0316.
%                detrendord => integer denoting the detrend order.
%                              -1 = no detrending
%                               0 = subtracts mean
%                               1 = subtracts best line fit (default)
%                interpord  => sinc function interpolation is performed
%                              to locate peak corelation position (delay)
%                              on a finer grid axis.  This is the number of
%                              neighborhood points used in the estimation.
%                              Usually 4 points is sufficient for a good
%                              interploation, much beyond 20 typically does not
%                              yield signficant accuracy even in a low noise 
%                              envornment.  Default is 10 (order 0 performs
%                              no interpolation).
%                interpfac  => Factor of increased resolution for estimating
%                              the peak. A factor of 10 implies the grid
%                              axis is sampled 10 times finer than the
%                              original.  Default is 10.
%
%  The outputs of the function are:
%  D     => relative time delay of SIG1 with respect to SIG2 (a negative
%           delay implies activity on SIG1 occurs before similar activity on
%           SIG2).
%  CV    => normalized value of cross correlation peak (correlation coefficient)
%  PRO   => Signal processing option data structure.  If not used for input,
%           this structure shows the default options that were used.
% 
%    written by Kevin D. Donohue, April 26, 2006 (donohue@engr.uky.edu)


%  if struct with processing parameters not provided create one with
%  the defaults
if nargin < 4
    %  Create Structure with defaults
    pro.cliptype = 'none';
    pro.cliplevel = 0;
    pro.detrendord = 1;
    pro.interpord = 4;
    pro.interpfac = 10;
else  %  if structure is provided determine the missing fields and provide defaults
    if ~isstruct(pro)  %  if not a structure, end program with error message
        error(['Error!  Argument' inputname(4) ' must be a data structure.'])
    else
        if ~isfield(pro,'cliptype')
            pro.cliptype = 'none';
        end
        if ~isfield(pro,'cliplevel')
            pro.cliplevel = .0316;
        else
            if pro.cliplevel >=1 | pro.cliplevel < 0
                error('Clip level must be beween 0 and 1')
            end
        end
        if ~isfield(pro,'detrendord')
            pro.detrendord = 1;
        else
            if ~(pro.detrendord == -1 |pro.detrendord == 0 | pro.detrendord == 1)
                error('Error!  Detrend order must be either -1, 0, or 1')
            end
        end
        if ~isfield(pro,'interpord')
            pro.interpord = 4;
        else
            if pro.interpord < 0  |  round(pro.interpord) ~= pro.interpord
                error(['Error! Order of interpolator must be a positive integer or 0.'])
            end
        end
        if ~isfield(pro,'interpfac')
            pro.interpfac = 10;
        else
            if pro.interpfac <= 0 |  round(pro.interpfac) ~= pro.interpfac
                 error(['Error! Order of interpolation factor an integer greater than 0.'])
            end
        end
    end
end


clipper = pro.cliplevel;  %  Set clipping threshold
upsamp = ceil(pro.interpfac);  %  Increase in resolution for time lag peak
interpint = pro.interpord;   % Number of samples to either side of the max point
                             %  in low resolution cross correlation
                  

%  Detrend signals to remove linear and DC offsets
if pro.detrendord == 1
    sig1 = detrend(sig1); 
    sig2 = detrend(sig2);
elseif pro.detrendord == 0
    sig1 = detrend(sig1,'constant'); 
    sig2 = detrend(sig2,'constant');
end

%  Clip signal according to specs in structure pro
if sum(pro.cliptype == 'none') ~=4
    sig1 = cclip(sig1,pro.cliptype,clipper);
    sig2 = cclip(sig2,pro.cliptype,clipper);
end

%  Compute cross correlation between signals
[xc, lags] = xcorr(sig1,sig2);
[mv, iv] = max(xc);  %  find max value in low resolution cross-correlation
mxpt = iv(1);        %  find position of max value
%  Extract points around maximum for upsampling and getting a higher
%  resolution delay estimate
st = max([1,(mxpt-interpint)]);
ed = min([length(xc),(mxpt+interpint)]);
segest = xc(st:ed);
% Apply a tapering window to limit edge effects
[r,c] = size(segest);
taperwin = hamming(length(segest));
if r>1   %  if column vector
    segest = segest.*taperwin;
else    % if row vector
    segest = segest.*taperwin';
end
%  Resample at higher rate
hrseg = resample(segest,upsamp,1);
[cv,iv] = max(hrseg);  %  find max on finer grid
cv = cv/(eps+(sqrt(sum(sig1.^2)*sum(sig2.^2))));  % Compute correlation coefficient
% add to beginning of extracted segment for actual delay
d = (lags(st) + (iv(1)-1)/upsamp)/fs;
%  Uncomment last set of lines to plot processed signals and show delay estimate
%figure(1);
%plot([0:length(sig1)-1], sig1, 'r', [0:length(sig1)-1], sig2, 'b');
%ds=d*fs*1000;
%title(['Delay = ' num2str(ds) ' ms'])
