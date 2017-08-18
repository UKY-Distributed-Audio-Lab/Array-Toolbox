%  This script simulates an impulse-like sound in a field of view (FOV) with a perimeter
%  array around the walls (offset by .25 meters toward the center) of the room of a
%  rectangular room.  In addition 2 white noise sources are simulated outside the FOV
%  on the actual wall of the room (representing fan noise or noise from
%  window).  These sources represent coherent noise.  The strongest signal on the mic
%  array is used for adjusting the noise power to achieve 10 dB SNR.  Therefore all
%  other mic signals have a less than 10 dB SNR.  In addition, a -30 dB white noise
%  signal is added to every mic with respect to the strongest signal to
%  represent low-level system or diffuse noise.
%
%  The simulated signal is then used to create a steered response power
%  (SRP) image using the function SRPFRAMENN.
%  Details and adjustments for the simulation are explained in the comments below
%
%   Written by Kevin D. Donohue (donohue@engr.uky.edu) October 2005.
clear
fno = 1;  %  Figure number for plot
fs= 16000;  %  Sample frequency in Hz
sigtot = 1;   %  Number of targets in FOV
numnos = 2;   %  Number of coherent targets on wall perimeter
              %  generate target bandwidths
snr = -25;  %  coherent noise sources SNR to be added relative to strongest target peaks
batar = .6; %  Beta values for PHAT processing
micnum = 8;  %  Number of mics in array to be tested
%  Target signal parameters
f12p = 3000;  %  Corresponding upper frequency limit
f11p = 100;  %  Lower frequency limit
%  White noise snr
wgnsnr = -30;
sclnos = 10^(wgnsnr/20);
%  Frequency dependent Attenuation
temp = 28; % Temperature centigrade
press = 29.92; % pressure inHg
hum = 80;  % humidity in percent
dis = 1;  %  Distance in meters (normalized to 1 meter)
prop.freq = fs/2*[0:200]/200;  %  Create 100 point frequency axis
prop.atten =  atmAtten(temp, press, hum, dis, prop.freq);  %  Attenuation vector
prop.c = SpeedOfSound(temp,hum,press);

%  Generate room geometry
%  Opposite corner points for room, also walls are locations for noise source placement
froom = [-3.5 -4 0; 3.5 4 3.5]';  % [x1, y1, z1; x2 y2 z2]'
% Opposite Corner Points of Perimeter for mic placement
fmics = [-3.25 -3.75 1.5; 3.25 3.75 1.5]';
%  Room reflection coefficients (walls, floor, ceiling)
bs = [.8 .8 .8 .8 .7 .7];
%  Field of view for reconstructing SRP image (opposite corner points)
fov = [-2.5 -2.5 1.5; 2.5 2.5 1.5]';

%  Time window for frequency domain block processing
trez = 20e-3;  %  In seconds
%  Room Resolution: Step through cartesion grid for mic and sound source
%  plane
rez = .04;  %  In meters

%  All vertcies in image plane
v = [fmics(1:2,1), [fmics(1,1); fmics(2,2)], fmics(1:2,2), [fmics(1,2); fmics(2,1)]];  
v = [v; ones(1,4)*1.5];
vn = [froom(1:2,1), [froom(1,1); froom(2,2)], froom(1:2,2), [froom(1,2); froom(2,1)]];  
vn = [vn; ones(1,4)*1.5];
%  Compute window length in samples for segmenting time signal 
winlen = ceil(fs*trez);
wininc = round(fs*trez/2);  %  Compute increment in sample for sliding time window along
%  Compute grid axis for pixel of the SRP image
gridax = {[fov(1,1):rez:fov(1,2)], [fov(2,1):rez:fov(2,2)], [fov(3,1):rez:fov(3,2)]}; 

%  Compute spacing for equal spacing of perimeter array
spm = (norm(v(:,2)-v(:,1),2)+norm(v(:,3)-v(:,2),2))/(micnum/2 - 2 +4/sqrt(2));
%  Compute starting point on perimeter from the first corner
stp = (spm/sqrt(2))*(norm(v(:,2)-v(:,1),2)/norm(v(:,3)-v(:,2),2));
%  Generate location of perimeter array microphones
mposperim = regmicsperim(v(:,[1,3]),spm,stp);
%  Find max distance (delay) over all mic pairs; this represents an upper bound
%  on all required relative delays when scanning over the FOV
[rm, nm] = size(mposperim);
prs = mposanaly(mposperim,2);
%  Maximum delay in seconds needed to synchronize in Delay and Sum beamforming
maxmicd = max(prs(:,3));
%  Extra delay time for padding window of data
textra = ceil(fs*maxmicd(1)/prop.c); 
%  Windowing function to taper edge effects in the whitening process
tapwin = flattap(winlen+textra,20);
winrec = tapwin*ones(1,micnum);

%  Simulate target signals
f11s = f11p-0.2*f11p;  %  Compute Stop bands from passbands
f12s = f12p+0.2*f12p;  %  Compute Stop bands from passbands
%  Ensure signal occurs at a late enough time to be included in
%  first window for processing
td = textra/(fs)+15*trez/2;
simsiglen = td+2*(4/(f12p-f11p) + 4/f11p);
%  Generate target waveforms
target1 = simimp(f11p,f12p,f11s,f12s,td,fs,simsiglen);
%  Expand to multiple targets if sigtot greater than 1
target1 = target1*ones(1,sigtot);

%  Random generation of signal position within FOV
sigpos = ((fov(:,2)-fov(:,1))*ones(1,sigtot)).*rand(3,sigtot) + fov(:,1)*ones(1,sigtot);
%  Compute array signals from target
[sigoutper, taxper] = simarraysigim(target1, fs, sigpos, mposperim, froom, bs, prop);
%  Random generation of coherent noise source positions on wall 
for knn=1:numnos
    randv = ceil(rand(1,1)*4);
    %  Noise source positions
    sigposn(:,knn) = vn(:,randv) + rand(1)*(vn(:,mod(randv,4)+1)-vn(:,randv));
end
% Create coherent white noise source with sample lengths as target signal
[rt,ct] = size(target1);
%  generate white noise 
onos = randn(rt,numnos);
%  place white noise target randomly on wall
[nosoutper, taxnosper] = simarraysigim(onos,fs, sigposn, mposperim, froom, bs, prop);
[mxp,cp] = max(max(abs(sigoutper)));  % Max point over all channels
envper = abs(hilbert(sigoutper(:,cp(1))));  % Compute envelope of strongest channel
%  Compute maximum envelope point for reference in SNRs
%  Also location of max point will be used to ensure time window processed includes
%  the target
[perpkpr, rpper] = max(envper);
%  Trim room signals to same length
[siglenper, mc] = size(sigoutper);
[noslenper, mc] = size(nosoutper);
siglen = min([siglenper, noslenper]);
sigoutper = sigoutper(1:siglen,:);
nosoutper = nosoutper(1:siglen,:);
%  Normalize noise power
nosoutper = nosoutper/sqrt(mean(mean(nosoutper.^2)));
%  Add coherent noise to target signals
nos = randn(siglen,mc);
asnr = 10^((snr/20));
nosamp = asnr*perpkpr;
sigoutpera = sigoutper + nosamp*nosoutper + nos*sclnos*perpkpr;
% Initialize signal window index to beginning index, offset to ensure it includes target
% signal
sst = 1+rpper(1)-fix(.9*winlen); 
sed = sst+min([winlen+textra, siglen]);   %  and end window end
%  create tapering window
tapwin = flattap(sed-sst+1,20);  %  One dimensional
wintap = tapwin*ones(1,micnum);  %  Extend to matrix covering all channels
%  Whiten signal (apply PHAT, with beta factor given at the begining)
sigout = whiten(sigoutpera(sst:sed,:).*wintap, batar);
%  Create SRP Image from processed perimeter array signals
im = srpframenn(sigout, gridax, mposperim, fs, prop.c, trez);
%  Set up figure for plotting
figure(fno)
%  Plot SRP image
imagesc(gridax{1},gridax{2}, im, [0, max(max(im))]);
colormap(jet); colorbar; axis('xy')
axis([froom(1,1)-.25, froom(1,2)+.25, froom(2,1)-.25, froom(2,2)+.25])
hold on
%  Mark actual target positions 
plot(sigpos(1,:),sigpos(2,:),'ok', 'MarkerSize', 18,'LineWidth', 2);
%  Mark coherenet noise positions
plot(sigposn(1,:),sigposn(2,:),'xb','MarkerSize', 18,'LineWidth', 2);  %  Coherent noise
%  Mark microphone positions
plot(mposperim(1,:),mposperim(2,:),'sr','MarkerSize', 12);
axis('tight')
%  Number them
for kn=1:length(mposperim(1,:))
    text(mposperim(1,kn),mposperim(2,kn), int2str(kn), 'HorizontalAlignment', 'center')
end
%  Draw Room walls
plot([vn(1,:), vn(1,1)],[vn(2,:), vn(2,1)],'k--')
% Label Plot
xlabel('Meters')
ylabel('Meters')
title(['SRP image (Mics at squares, Target in circle, Noise sources at Xs'] )
hold off
%  Plot signal array
figure(fno+1)
offset = zeros(1,micnum); % Initialize offset vector
for km=1:micnum
    %plot offset
    offset(km) = max(abs(sigoutpera([sst:sed],km))) + .1*std(sigoutpera([sst:sed],km));
end
fixoff = max(offset);
offt = 0;
for km=1:micnum
    offt = fixoff +offt;
    plot((sst:sed)/fs,sigoutpera(sst:sed,km)+offt)
    hold on
end
hold off
set(gca,'ytick',[])
xlabel('Seconds')
title('Array Signals, Mic 1 is on the bottom')
figure(fno)