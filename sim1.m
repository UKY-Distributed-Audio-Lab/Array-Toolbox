%  Simulation 1 - Array Uncorrelated White Noise Response 
%  mean and rms of field of view for white noise excitation
%  of microphone channels.  
%  Microphone geometries include even spacings of microphones over
%  a rectangular perimeter for 4, 8, and 16 microphones:
%  all on one side
%  split over opposite sides
%  spread over entire perimeter
%  No multipath is considred
%  In this first attempt Monte Carlo there will be 50 runs for each mic
%  configureation
%  variables computed will be:
%  mean pixel values
%  rms pixesl values
%  histograms of all pixel values
%  histograms of negative pixel values
%  histograms of positive pixel values
mruns = 4;

f11p = 1500;
f12p = 2600;
f11s = f11p-.1*f11p;
f12s = f12p+.1*f12p;

f21p = 2200;
f22p = 2600;
f21s = f21p-.1*f21p;
f22s = f22p+.1*f22p;

td = 0;
fs = 16000;
siglen = 2*(td+20/(f12s-f11s) + 5/f11p);
[target1, tax] = simimp(f11p,f12p,f11s,f12s,td,fs,siglen);
target2 = simimp(f21p,f22p,f21s,f22s,td,fs,siglen);
% load up signal array
sigarray = [target1, target2];

% SNR for uncorrelated noise
snr = 12;

%  Generate room geometry
c = 348;  %  Speed of sound in room
froom = [-3 0; 3 6]';
fov = [-1.5 1.5; 1.5 4.5]';
rez = c/(2*f12p);  % Grid resolution
trez = 20e-3;  %  Time resolution over which to compute power

gridax = {[fov(1,1):rez:fov(1,2)], [fov(2,1):rez:fov(2,2)]}; 
%  Reflection coefficient of multipath scatterers
rcoefinc = .2;
for mpk=2:2
rcoef = 2*(mpk-1)*rcoefinc;
%  Scatterer lines
fmscatbw1 = [-3 -.2; 3 -.2]';
fmscatfw1 = [-3 6.2; 3 6.2]';

fmscatbw2 = [-3 -6.2; 3 -6.2]';
fmscatfw2 = [-3 12.2; 3 12.2]';
%fmscatbw3 = [-3 -12; 3 -12]';
%fmscatfw3 = [-3 18; 3 18]';
%  Scatter spacings
spsc = 2;
mpscatbw1= regmicsline(fmscatbw1, spsc);
[r,s] = size(mpscatbw1);
mpscatbw1 = [ones(1,s)*rcoef; mpscatbw1];

mpscatbw2= regmicsline(fmscatbw2, spsc);
[r,s] = size(mpscatbw2);
mpscatbw2 = [ones(1,s)*rcoef^2; mpscatbw2];

%mpscatbw3= regmicsline(fmscatbw3, spsc);
%[r,s] = size(mpscatbw3);
%mpscatbw3 = [ones(1,s)*rcoef^3; mpscatbw3];

mpscatfw1= regmicsline(fmscatfw1, spsc);
[r,s] = size(mpscatfw1);
mpscatfw1 = [ones(1,s)*rcoef; mpscatfw1];

mpscatfw2= regmicsline(fmscatfw2, spsc);
[r,s] = size(mpscatfw2);
mpscatfw2 = [ones(1,s)*rcoef^2; mpscatfw2];

%mpscatfw3= regmicsline(fmscatfw3, spsc);
%[r,s] = size(mpscatfw3);
%mpscatfw3 = [ones(1,s)*rcoef^3; mpscatfw3];

mpscat = [];
% mpscat = [mpscatbw1,mpscatbw2,mpscatbw3, mpscatfw1,mpscatfw2,mpscatfw3];
mpscat = [mpscatbw1,mpscatbw2, mpscatfw1,mpscatfw2 ];

%  Set up microphone geometry
spm = 1;
fom = [-2 0; 2 0]';
% mpos1= regmicsline(fom, spm);
% fom2 = [-2 6;  2 6]';
% mpos2= regmicsline(fom2, spm);
% fom3 = [-2.5 0;  -2.5 6]';
% mpos3= regmicsline(fom3, spm);
% fom4 = [2.5 0;  2.5 6]';
% mpos4= regmicsline(fom4, spm);
% 
% 
% mpos = [mpos1, mpos2, mpos3, mpos4];
mpos = regmicsperim(froom,spm);
%  Position scatterers
st= fov(:,1);  %  Smallest X-Y position in FOV
delt = fov(:,2)-fov(:,1);  %  Delta Range for scatter position
[mr, mc] = size(sigarray);  %  Get number of scatterer positions to generate (columns of sigarray)

%  Set window and signal lengths
wlen = ceil(trez*fs);  %  Number of samples for computing power
corwin = trez;
%corwin = hanning(wlen); %  Tapering window
% Additional window length to compensate for maximum delay
textra = ceil(fs*norm(fom(:,1)-fom(:,2),2)/c); 

[rm, nm] = size(mpos);
prs = mposanaly(mpos,2);
[g, msrt] = sort(prs(:,3));
out = prs(msrt,1:3);
hh = find(out(:,3) < 2.5);
prs = out(hh,:);
nrm = max(prs(:,3));
prs(:,3) = (prs(:,3)/(nrm+eps));
%prs(:,3) = 1;
for monte=1:mruns
    % Generate random target positions 
    for k=1:mc
        sigpos(:,k) = st + rand(2,1).*delt;
    end
    [sigout, tax] = simarraysig(sigarray, fs, sigpos, mpos, c, mpscat); %, mpscat);
    sigrms = max(max(abs(sigout)))/sqrt(2);
    nosamp =10^(-snr/20)*sigrms(1);
    nos = nosamp*randn(size(sigout));
    sigout = sigout+nos;
    [siglen,cs] = size(sigout);

    imcount = 1;
    mleim = [];
    sigposs = []

    for nalpha=1:1:1
            % figure(1)
            %
            % [rm,cm] = size(mpos);
            % for k=1:cm
            % plot(mpos(1,k),mpos(2,k),'rx')
            % hold on
            % end
            % [rm,cm] = size(sigpos);
            % for k=1:cm
            % plot(sigpos(1,k),sigpos(2,k),'ko')
            % hold on
            % end
            % [rm,cm] = size(mpscat);
            % for k=1:cmc
            % plot(mpscat(2,k),mpscat(3,k),'gd')
            % hold on
            % end
            % axis([-4 4 -13 19]);
            % hold off
            %
            % pause(.1)

            sigst = 1;
            sigend = sigst+wlen;
            wint = hanning(sigend);
            winrec = wint*ones(1,nm);
            while sigend <= siglen
                %g = gradient(abs(hilbert(sigout(sigst:sigend,:))));
                g = whiten(sigout(sigst:sigend,:).*winrec,(nalpha-1)/10);
                %mleim(:,:,imcount) = srpframepairs(g, gridax, mpos, prs, fs, c, corwin);
                imo = srpframe(g, gridax, mpos, fs, c, corwin);
                sigposs(:,imcount) = [sigpos(:,1); sigpos(:,2)];
                mleim(:,:,imcount) = localdnoise(imo,6,10);
                imagesc(gridax{1}, gridax{2}, squeeze(mleim(:,:,imcount))); axis('xy'); colorbar
                hold on
                plot(sigposs(1,1),sigposs(2,1),'ko')
                plot(sigposs(3,1),sigposs(4,1),'kx')
                
                hold off

                imcount = imcount+1
                pause

                sigst = sigst + ceil(trez*fs/2);
                sigend = sigst + wlen;
            end
    end
    fname = ['monte' int2str(monte) '.mat']
    save(fname, 'mleim', 'sigposs', 'gridax')
end
end

    %    fname = ['monal' int2str(nalpha) 'mp' int2str(mpk) '.mat']
    %    save(fname, 'mleim', 'sigposs', 'gridax') 
    %    [r,c,dd] = size(mleim);
    %     for kk=1:dd
    %     imagesc(gridax{1}, gridax{2}, squeeze(mleim(:,:,kk))); axis('xy'); colorbar
    %     hold on
    %     plot(sigposs(1,1),sigposs(2,1),'ko')
    %     plot(sigposs(3,1),sigposs(4,1),'ko')
    %     hold off
    %     pause
    %     end
