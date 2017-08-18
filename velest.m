function [vest, cst]  = velest(ba,fs,d,vr);
%  This function estimates velocity of sound, where a propagating soundwave
%  is measured with an endfire linear array and its samples are stored in matrix BA.
%  The signal can be extended over time (a white noise burst is usually best) and
%  this program automatically steps through the sequence using a sliding window.
%  The values for sequential windows are in output vector VEST.  These values can be
%  averaged (or weighted and averaged) for a final velocity estimate.  The values from
%  each window can be use to obtain a standard deviation and provide an estimate of the
%  variability of the estimate.  The program assumes that the sonifying signal was broadband
%  (a white noise burst or a sequence of impulse like noises) and the speed of sound is close 
%  to 340 m/s (this can changed with optional parameter VR).  The program uses this information 
%  to decide on window size and filtering paramters. It also assumes that the signal after
%  filtering is clean enough to generate a .4 correlation coefficient.  For
%  any window with correlation coefficeint value less than this in the sequence, the
%  velocity estimate is set to NaN.
%
%              [vest,cst] = velest(ba,fs,d, vr);
%
%  Inputs:  BA is a matrix of timesamples where each column represents a signal
%           from one of the array microphones.  Adjacent microphones must
%           be adjacent columns in the matrix, where the microphones with
%           lower index values are insonified first
%           FS is the sampling frequency.
%           D is a vector denoting the distances in meters between each
%           adjacent pair of microphones (number of elements in D is one
%           less than number of microphones in array)
%           (Optional)  VR  initial velocity estimate so window sizes can
%           be appropiately chosen for the correlation operations.  Default
%           is 340 m/s.
% Outputs:  VEST is the actual velocity measurment matrix corresponding to each pair
%           of adjacent microphones (rows) and window position along the time 
%           signal (columns).
%           If cross-correlation did not exceed a particular threshold (.5) the output is
%           set to NaN.
%           CST is the peak correlation coefficient value associated with each velocity
%           estimate in VEST.
%
%   Written by Kevin D.Donohue (donohue@engr.uky.edu) June 14, 2006

[smps, ch] = size(ba);  % Get the samples and number of mic channels.
thr = .4;   %  Set threshold for rejecting estimates based on correlation coefficients values
td = sum(d);  %  Compute total distance spanned by mic array
if nargin == 3,
    vr = 340;   %  A rough estimate of the sound speed

end
%  Estimate a reasonable frame size that is big enough to include the furthest
%  pair of adjacent Mics

frm1 = round(100*fs*td/vr);   %  Smaller frame for first window
frms = round(200*fs*td/vr); %  Get larger frame for second window
olap = floor(frms/4);     %  Set overlap for steping through correlation signals

md = min(d);
fc = vr*2/md(1);  %  Find a cutoff frequency so that about 2 wavelenths or more will
                   %  fit in the smallest distance between microphones.
            
[b,a] = butter(6,2*fc/fs,'high');  %  High-pass filter coefficient for processing raw signal.

%  Perform high-pass filter on each channel
ba = filter(b,a,ba);


kw = 1;  %  Start at first sample:
ke1 = kw+frm1-1;  %  Extend to one window length for small
ke = kw+frms-1;  %  Extend to one window length for large

vest = zeros(ch-1,floor((smps-frms)/olap));
cst = zeros(ch-1,floor((smps-frms)/olap));

%  Step through entire signal 1 frame at a time
kstep=1;   %  Sample denoting begining of frame
while ke < smps;     %  Keep looping as long as end of frame
                     %   is less than total signal length
    % take correlation for all adjacent pairs
    for k=2:ch
        sigsm = detrend(ba(kw:ke1,k-1));  %   Extract Small Segment
        siglg = detrend(ba(kw:ke,k));     %   Extract Large Segment
        [dly, cv] = delayesttm(siglg,sigsm,fs);  %  Compute Delay
        %  Test to see if correlation value is sufficiently high for
        %  a reliable estimate
        if cv > thr
            vest(k-1,kstep) = d(k-1)/(dly+eps);   %  compute velocity for this pair
        else 
            vest(k-1,kstep) = NaN;   %  If correlation value too low set estimate to 0
        end
        cst(k-1,kstep) = cv;  % Normalize to create correlation coefficient
    end
    kstep = kstep+1;  %  Increment array index
    kw  = kw+olap;    %  Update frame starte sample
    ke1 = kw+frm1-1;  %  Extend to one window length for small
    ke  = kw+frms-1;  %  Extend to one window length for large
end

