function [d, stdds, pro] = delayestfr(sig1,sig2,fs, pro)
%  This function estimates the time delay between 2 similar signal segments
%  (SIG1 and SIG2) of the same length and sampled on the same time reference
%  axis with rate (FS):
%
%       [d, cv, pro] = delayesttm(sig1,sig2,fs,pro)
%
%  The function estimates the delay directly in the frequency domain based on
%  the slope of the group delay of the spectral products (this multiplication
%  where one signal is conjugated corresponds to correlation in time domain). 
%  This function will typically work best (at least better than in time
%  domain) when the signal is narowband and noisy since this program gives the option
%  of restricting the spectral region over which the group delay is estimated.
%  
%  The function inputs are:
%  SIG1  => Signal segment with assumed delay
%  SIG2  => Reference signal segment.  SIG1 and SIG2 must be the same size
%  FS    => Sampling frequency
%  PRO   => (optional) Data Structure with Signal processing options
%       Fields:  flow   =>  lower frequency limit in Hz of signal bandwidth
%                            (default is 100Hz)
%                fhigh =>   upper frequency limit in Hz of signal bandwidth
%                           (default is 15kHz)
%                detrendord => integer denoting the detrend order.
%                              -1 = no detrending
%                               0 = subtracts mean
%                               1 = subtracts best line fit (default)
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
    pro.flow = 100;
    pro.fhigh = 15000;
    pro.detrendord = 1;
    
else  %  if structure is provided determine the missing fields and provide defaults
    if ~isstruct(pro)  %  if not a structure, end program with error message
        error(['Error!  Argument' inputname(4) ' must be a data structure.'])
    else
        if ~isfield(pro,'flow')
            pro.flow = 100;
        end
        if ~isfield(pro,'fhigh')
            pro.fhigh = 15000;
        else
            if pro.flow < 0  | pro.flow >= pro.fhigh
                error('Low frequency limit must be positive and less than upper frequency limit')
            end
        end
        if ~isfield(pro,'detrendord')
            pro.detrendord = 1;
        else
            if ~(pro.detrendord == -1 |pro.detrendord == 0 | pro.detrendord == 1)
                error('Error!  Detrend order must be either -1, 0, or 1')
            end
        end
    end
end

       

%  Detrend signals to remove linear and DC offsets
if pro.detrendord == 1
    sig1 = detrend(sig1); 
    sig2 = detrend(sig2);
elseif pro.detrendord == 0
    sig1 = detrend(sig1,'constant'); 
    sig2 = detrend(sig2,'constant');
end


%plot(fax,phs);%  Compute FFTs
len1 = length(sig1);
len2 = length(sig2);
len = max([len1, len2]);  %  use longer signal length for zero pad
nfft = 2^nextpow2(2*len); %  pad with zero to next power of 2
fsig1 = fft(sig1,nfft);  %  FT of signal 1
fsig2 = fft(sig2,nfft);  %  FT of signal 2
%  Take product
pd = conj(fsig1).*(fsig2);
%  Compute frequency axis
fax = fs*[0:nfft-1]/nfft;
phs = unwrap(angle(pd));
%pause
gpd = gradient(phs,fs/nfft);
wrng = find(fax >= pro.flow & fax <= pro.fhigh);
if ~isempty(wrng)
ds = (abs(pd(wrng)).*gpd(wrng))/(mean(abs(pd(wrng)))*2*pi);
d = mean(ds);
stdds = std(ds);
else
    disp(['No points were found in the frequency range from ' num2str(pro.flow), ' to ' num2str(pro.fhigh)])
    d = nan;
    stdds = nan;
end
% Uncomment last set of lines to plot processed signals and show delay estimate
% figure(1);
% plot([0:length(sig1)-1], sig1, 'r', [0:length(sig1)-1], sig2, 'b');
% ds=d*fs/1000;
% title(['Delay = ' num2str(ds) ' ms'])
