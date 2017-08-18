%  This script will creat graphics to assess the potnential information
%  that can be extracted from the pop of a synovial joint.
clear
close all
figpos = [452   161   635   270];
poploc = [24680; 222500; 540300; 762900;  987000];  %  Location of pops in samples
winlen = 2048+1024;
pre = 300;
post = winlen-pre-1;
%  Path and filename of wave file with pops
fn = 'Dr Fryer Cracks.wav';
pn = 'C:\Users\Donohue\Documents\';
%  Read in file
[y,fs] = wavread([pn,fn]);
[bb,aa] = butter(4,2*[10000]/fs,'High');


c = 350;
spacerez = .002;  %  Desired Resolution in meters
timerez = spacerez/c;

%  Segment
tax = [0:winlen-1]/fs;
seg = zeros(winlen,length(poploc));
for k=1:length(poploc)
    seg(:,k) = y(poploc(k)-pre:poploc(k)+post);
    seg(:,k) = filtfilt(bb,aa,seg(:,k));
    %seg(:,k) = whiten(seg(:,k),.9);
    figure
    plot(tax*1000,seg(:,k))
    xlabel('ms')
    ylabel('Volts');
    title(['Time waveform for pop ' int2str(k)])
    set(gcf,'Position', figpos);
    xlim([0 20])
    print -dmeta
    pause
end
segr = resample(seg,10*fs,fs);
acs = zeros(10*winlen,length(poploc));
for k =1:length(poploc)
    [ac,lgs] = xcorr(segr(:,k),'coef');
    acs(:,k) = ac(10*winlen:end);
    lg = lgs(10*winlen:end);
    figure
    plot(1000*lg/(fs*10),acs(:,k))
    xlabel('ms')
    ylabel('Correlation Coefficient')
    title(['Autocorrelation for pop ' int2str(k) ' with resolutions limits'])
    hold on
    plot(1000*[timerez, timerez], [-.4 1], 'g--', 'Linewidth', 2)
    gg = find(acs(:,k) > .5);
    xlim([-.03, 2])
    plot(1000*[lg(gg(end)) lg(gg(end))]/(fs*10), [-.4 1], 'r--', 'Linewidth', 2)
    hold off
    legend('Autocorrelation','Desired Resolution', 'Possible Resolution')
    set(gcf,'Position', figpos);
    print -dmeta
    pause
end
%  Compute specta
fax = fs*[0:winlen-1]/(2*winlen);
for k =1:length(poploc)
    spec = abs(fft(seg(:,k),2*winlen));
    figure
    semilogx(fax,20*log10(spec(1:winlen)));
    title(['Spectrum for pop ' int2str(k)])
    xlabel('Hz');
    ylabel('dB')
    ylim([-60 40])
    set(gcf,'Position', figpos);
    print -dmeta
    grid
    pause
end
    
    

