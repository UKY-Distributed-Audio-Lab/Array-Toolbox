%  This scipt is used to demonstrate the use and performance of the
%  different delay functions in the Microphone Array Toolbox.  The
%  functions needed to run this script are:
%
%  delayint.m, delaytab.m delayt.m delayf.m, and subsamplefir.m
%
%  The script will graphically compare the original nondelayed
%  version of the signal with the delayed outputs from the functions,
%  delayt (time domain subsample delay) and delayf (frequency domain
%  subsample delay).  A comparison will also be presented to show
%  the error (rms difference) between the frequency subsample delay and
%  all other implementations.  Graphs of the mean execution times for each
%  of the algorithms will also be generated along with the normalized rms
%  error between the frequency delay and the other time domain delay
%  functions.
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu)  July 2005

fs = 16E3;  %  Sampling frequency in Hz
ns = round(fs/8);   %  Number of points  in signal.
siglen = .02; %length of signal in seconds
sigstart = .005;
[g, tax] = simimp(1000,2500,500,3000,sigstart,fs,siglen);
%g = randn(1,ns);  %  Generate white noise signal

% If the above variables are changed, be sure to uncomment the axis lines
% in the for loop


%  If an impulse signal is desired uncomment the next 2 lines
%g = zeros(1,ns);   %  Initialize array for impulse signal 
%g(23) = 1;         %  place impulse at some time index point

%  Generate FIR coefficient table for the "delaytab" function
subinc = 4;  % subsample increment
ord = 4;     % number of filter coefficients
tab = subsamplefir(subinc,ord,'cos');   %  Selected cosine squared-weighted sinc 

%  Loop to step through a series of delays while plotting comparisons 
inct = 1/(3*fs);  %  Delay increment seconds
inc = inct*fs; % Delay increment in samples
numd = 60;  %  Number of delays

% Position figure for simultaneously viewing 2 delay function outputs
 figure(1)
 set(1,'Position', [37   420   466   254])
 
 figure(2)
 set(2,'Position', [543   420   466   254])
 
 figure(3)
 set(3,'Position', [37   80   466   254])
 
 figure(4)
 set(4,'Position', [543   80   466   254])


for k=1:numd
    %  Frequency domain shift
    figure(1)
    tic;  %  Reference time
    [sdf, tn] = delayf(g,fs,inct*k, siglen);  %  Implement Delay 
    tf(k) = toc;  %  Record time  
    plot(tn,sdf,'r',[0:length(g)-1]/fs,g,'b')  % plot shifted version with original
    xlabel('seconds')
    axis([0 .02 -.2 .3])
    legend('Delayed','Original')
    title('Frequency Domain Subsample Delay with Original')
 
    %  Time domain shift with FIR fractional sample shift
    figure(2)
    tic;    %  Reference time
    [sdt, tn] = delayt(g,fs,inct*k, siglen);  %  Implement Delay
    tt(k) = toc;    %  Record time  
    plot(tn,sdt,'r',[0:length(g)-1]/fs,g,'b')  % plot shifted version with original
    xlabel('seconds')
    axis([0 .02 -.2 .3])
    legend('Delayed','Original')
    title('Time Domain Delay Subsample Filtering with Original')
 
    %  Time domain shift with Table of FIR fractional sample shifts (round
    %  off to nearest subincrement)
   figure(3)
    tic;    %  Reference time
    [sdb, tn] = delaytab(g,tab,fs,inct*k, siglen);  %  Implement Delay
    tb(k) = toc;   %  Record time  
   plot(tn,sdb,'r',[0:length(g)-1]/fs,g,'b')  % plot shifted version with original
   xlabel('seconds')
    legend('Delayed','Original')
   title('Time Domain Delay Subsample Roundoff Table with Original')

    %  Vector offset shift (round off to nearest integer sample)
   figure(4)
    tic;    %  Reference time
    [sdi, tn] = delayint(g,fs,inct*k, siglen);  %  Implement Delay
    ti(k) = toc;    %  Record time  
   plot(tn,sdi,'r',[0:length(g)-1]/fs,g,'b')  % plot shifted version with original
   xlabel('seconds')
    legend('Delayed','Original')
   title('Time Domain Delay Integer Roundoff with Original')
 
 %  Compute SNR from shift error (relative to frequency shift (in dB)
    errt(k) = 10*log10( sum((sdf-sdt).^2) /sum((sdf).^2)) ;  %  Time domain delay
    errb(k) = 10*log10( sum((sdf-sdb).^2)/ sum((sdf).^2));  %  Time domain table roundoff delay
    erri(k) = 10*log10(sum((sdf-sdi).^2) / sum((sdf).^2));  %  Time domain integer delay
    %pause(.4)  %  pause to look at shifted outputs
    pause(.10)
end

%  vector of mean time elapsed for each function
figure(5)
%set(gcf, 'Position', [100   398   763   305])
timelps = [mean(tf), mean(tt), mean(tb), mean(ti)]*1000;  %  mean elapsed time in milliseconds
stem(timelps);
set(gca,'XLim',[0,5],'XTick',[1:4])
set(gca,'XTickLabel',{'frequency', 'time subsample', 'table subsample', 'integer offset'})
ylabel('Milliseconds');
title('Comparison of elapsed time for each delay function execution')
% vector of errors
figure(6)
%set(gcf, 'Position', [273   151   623   305])
er = [mean(errt), mean(errb), mean(erri)];
plot([1:3],er,'x-');
set(gca,'XLim',[0,4],'XTick',[1:3])
set(gca,'XTickLabel',{'time subsample', 'table subsample', 'integer offset'})
ylabel('dB');
title('Comparison of normalized error with respect to frequency domain delays')



