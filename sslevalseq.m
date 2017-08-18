filenames(:,1) = { ...
     'nos___c_-1_-1.mat';... 
     'nos___c_-1_0.mat';... 
     'nos___c_-1_1.mat';... 
     'nos___c_0_-1.mat';... 
     'nos___c_0_0.mat';... 
     'nos___c_0_1.mat';... 
     'nos___c_1_-1.mat';... 
     'nos___c_1_0.mat';... 
     'nos___c_1_1.mat'};  %...
filenames(:,2) = { ...
     'nos_n_c_-1_-1.mat';... 
     'nos_n_c_-1_0.mat';... 
     'nos_n_c_-1_1.mat';... 
     'nos_n_c_0_-1.mat';...
     'nos_n_c_0_0.mat';... 
     'nos_n_c_0_1.mat';... 
      'nos_n_c_1_-1.mat';...
      'nos_n_c_1_0.mat';...
     'nos_n_c_1_1.mat'};  %...
filenames(:,3) = { ...     
     'words___c_-1_-1.mat';... 
     'words___c_-1_0.mat';... 
     'words___c_-1_1.mat';... 
     'words___c_0_-1.mat';... 
     'words___c_0_0.mat';... 
     'words___c_0_1.mat';... 
     'words___c_1_-1.mat';... 
     'words___c_1_0.mat';... 
     'words___c_1_1.mat'};  %...
filenames(:,4) = { ...
     'words_n_c_-1_-1.mat';... 
     'words_n_c_-1_0.mat';... 
     'words_n_c_-1_1.mat';... 
     'words_n_c_0_-1.mat';...
     'words_n_c_0_0.mat';... 
     'words_n_c_0_1.mat';... 
     'words_n_c_1_-1.mat';...
     'words_n_c_1_0.mat';...
     'words_n_c_1_1.mat'};
 
 fout = {'seqwnb.mat'; ...
         'seqwnbnos.mat'; ...
         'seqwrd.mat'; ...
         'seqwrdnos.mat'};
 
envparametersexp5_15_2006
pflag =0;
pars.snrtype = 'p';
batar =  [0 .2  .4  .6  .8  1]; %  Beta values to be tested
trez = 30e-3;  %  Processing Window Length In seconds
%  Room Resolution: Step through cartesion grid for mic and sound source
%  plane
rez = .02;  %  In meters


froom = [-3.25 -3.25 0; 3.25 3.25 3]';  % Opposite Corner Points of room

%  All vertcies in image plane
v = [froom(1:2,1), [froom(1,1); froom(2,2)], froom(1:2,2), [froom(1,2); froom(2,1)]];  
v = [v; ones(1,4)*mean(fov(3,:))];

gridax = {[fov(1,1):rez:fov(1,2)], [fov(2,1):rez:fov(2,2)], [fov(3,1):rez:fov(3,2)]}; 
cfreq = [300, 7000];  %  Critical frequency range where most signal energy resides

mposperim = mpos;
    %  Find max distance (delay) over all mic pair, this represents
    %  and upper bound on all requried relative delays when scanning over
    %  the FOV

    [rm, nm] = size(mposperim);
    prs = mposanaly(mposperim,2);
    maxmicd = max(prs(:,3));
    %  Extra delay time in seconds for padding window of data
    textra = maxmicd(1)/c; 
    [posnum, typenum] = size(filenames);
    
    
  snrs = [1 4 6];  %  Selecte particualr SNR values from recorded data
 for ktype = 1:typenum
 for k = 1:posnum
     filenames{k,ktype}
     load(['FoamOnly/' filenames{k,ktype}]);
     flim = 2*cfreq/fs;  %  Set frequency limit on PHAT transform     
     %  Initalize performance metric arrays for new file type
     if k == 1
       tmpr = cell(length(snrs), length(batar));  %  Initalize target magnitude array
       lerpr =cell(length(snrs), length(batar));  %  Initalize position error array
       nmpr = cell(length(snrs), length(batar));  %  Initalize noise magnitude array for target plus noise frames
       stpn = cell(length(snrs), length(batar));  %  Initalize mean peak target per channel to rms noise in target frame ratio
       nompr = cell(length(snrs), length(batar));  %  Initalize noise magnitude array for noise only frames
       stpno = cell(length(snrs), length(batar));  %  Initalize mean peak target per channel to rms noise in noise only frame ratio
     end
     %  Adjust target locations
     if tloc(1) == -1 & tloc(2) == -1
         tloc(1) = tloc(1); 
         tloc(2) = tloc(2) + .15;
     elseif tloc(1) == 0 & tloc(2) == -1;
         tloc(1) = tloc(1) - .015;
         tloc(2) = tloc(2) - .02+.15;
     elseif tloc(1) == 1 & tloc(2) == -1
         tloc(1) = tloc(1) + 0.015;
         tloc(2) = tloc(2) + .15;
     elseif tloc(1) == -1 & tloc(2) == 0
         tloc(1) = tloc(1) -.01;
         tloc(2) = tloc(2) + .15 + .01;
     elseif tloc(1) == -1 & tloc(2) == 1
         tloc(1) = tloc(1);
         tloc(2) = tloc(2) + .15 - .01;
     elseif tloc(1) == 0 & tloc(2) == 0
         tloc(1) = tloc(1) -.01;
         tloc(2) = tloc(2) +.15 -.005;
     elseif tloc(1) == 0 & tloc(2) == 1
         tloc(1) = tloc(1)-.015;
         tloc(2) = tloc(2) + .15 - .01;
     elseif tloc(1) == 1 & tloc(2) == 0
         tloc(1) = tloc(1)+.02;
         tloc(2) = tloc(2)+.15-.04;
     elseif tloc(1) == 1 & tloc(2) == 1
         tloc(1) = tloc(1)+.015;
         tloc(2) = tloc(2)+.15-.015;
     end
     %  Create corresponding time axis
     tax = [0:length(y)-1]/fs;
     winlen = ceil(fs*trez);  %  WIndow lenght in samples
     wininc = round(fs*trez);  %  WIndow increment in samples
     prez = ceil(fs*rez*sqrt(2)/c);  %  Time duration corresponding to 1 REZ pixel
     prez = [300]*prez;   %  Time interval for power integration (up to TREZ limit)
     [rsam, cmics] = size(y);  %  Get size in samples and mic channels of total mic recordings
     seglen = ceil(fs*(textra)) +winlen; %  Get total segment lenght to account of max delay differences of potential sounds
     tapwin = flattap(seglen,20);  %  Create a flat window to taper edges
     winrec = tapwin*ones(1,cmics); % Extend window to cover array channels
     [b,a] = butter(4,cfreq(1)/(2*fs), 'High');  %  Filter out room noise
     y = filtfilt(b,a,y);

    
     %  SNR value loop
     for kpn=1:length(snrs)
         kp = snrs(kpn)  %  Assign index to SNR location
         %  Get indecies for segment corresponding to SNR signal
         stax = find(tax >= par.targetb(kp) & tax <= par.targete(kp));
         sst = stax(1) + 2*winlen;
         sed = sst + seglen;
         segcnt = 1;  %  Initalize segment counter
         if kp > length(par.noiseb)
             ntax = find(tax >= par.noiseb(end) & tax <= par.noisee(end));
         else
             ntax = find(tax >= par.noiseb(kp) & tax <= par.noisee(kp));
         end
         nst =  ntax(1);  % Initial signal window index 
         ned = nst + seglen; 
         sigpow = [];
         nospow = [];
              while sed < stax(end) - 2*winlen 
              sigpow(segcnt) = snrarray(y([sst:sed],:),[],pars);
                       for batak = 1:length(batar)
                        bata = batar(batak);
                        %  Create Images from RF Perimeter Array
                        sigout = whiten(y([sst:sed],:),bata, flim);
                        mleim4 = srpframenn((sigout), gridax, mposperim, fs, c, trez, prez);
                        %  Analyze peaks
                        if pflag == 0
                            
                            [dettruepkspr, locerrpr, tarmagpr, nospkspr, nosmagpr] = peaksort(mleim4,gridax{1},gridax{2},tloc');
                            tmpr{kpn,batak} = [tmpr{kpn,batak}, tarmagpr];  %  Append target magnitudes
                            lerpr{kpn,batak} =[lerpr{kpn,batak}, locerrpr]; %  Append position errors
                            bignos = sort(nosmagpr, 'descend'); 
                            nmpr{kpn,batak} = [nmpr{kpn,batak}, bignos(1:8)];  %  Append noise magnitudes
                           
                            if ~isempty(nosmagpr)
                                stpn{kpn,batak} = [stpn{kpn,batak}, mean(tarmagpr)/(mean(nosmagpr)+eps)]; %  Save mean target to peak noise
                            else
                                stpn{kpn,batak} = [stpn{kpn,batak}, Inf];
                            end
                        else
                            figure(1); plot(y(sst:sed,:)); title('signal segement')
                            [dettruepkspr, locerrpr, tarmagpr, nospkspr, nosmagpr] = peaksort(mleim4,gridax{1},gridax{2},tloc', 4);
                            %  Plot coherenet noise positions
                            axis([-2 2 -2 2]);
                            hold on
                            %for knn=1:numnos
                            %putpeaks(4,tloc(1),tloc(2),'xb');  %  Coherent noise
                            %end
                            %  Plot mic array
                           % hold on
                             plot(mposperim(1,:),mposperim(2,:),'<r');
                             xlabel('Meters')
                             ylabel('Meters')
                             title(['Perimeter Array Target beta = ' num2str(bata) ' SNR = '  num2str(sigpow(segcnt))])
                             climit = get(gca,'Clim');
                             set(gca,'Clim', [0 climit(2)])
                            hold off
                            pause
                        end
                       end
                       segcnt = segcnt+1;
                       sst = sst+wininc;
                       sed = sst+seglen;
              end
              segcnt = 1;         
              nospow = [];
              while ned < ntax(end) 
              nospow(segcnt) = snrarray([],y([nst:ned],:),pars);
                       for batak = 1:length(batar)
                        bata = batar(batak);
                        sigout = whiten(y([nst:ned],:),bata,flim);
                        mleim4 = srpframenn((sigout), gridax, mposperim, fs, c, trez, prez);
                        %  Analyze peaks
                        if pflag == 0
                            [dettruepkspr, locerrpr, tarmagpr0, nospkspr, nosmagpr] = peaksort(mleim4,gridax{1},gridax{2},[]);
                            bignos = sort(nosmagpr,'descend');
                            nompr{kpn,batak} = [nompr{kpn,batak}, bignos(1:8)];  %  Append noise magnitudes
                            if ~isempty(nosmagpr)
                                stpno{kpn,batak} = [stpno{kpn,batak}, mean(tarmagpr)/(mean(nosmagpr)+eps)]; %  Save mean target to peak noise
                            else
                                stpno{kpn,batak} = [stpno{kpn,batak}, Inf];
                            end
                        else
                            figure(2); plot(y(nst:ned,:)); title('noise segement')
                            [dettruepkspr, locerrpr, tarmagpr, nospkspr, nosmagpr] = peaksort(mleim4,gridax{1},gridax{2},tloc', 5);
                            %  Plot coherenet noise positions
                            axis([-2 2 -2 2]);
                           % set(gca,'Clim', climit)
                            hold on
                            %for knn=1:numnos
                            %putpeaks(4,tloc(1),tloc(2),'xb');  %  Coherent noise
                            %end
                            %  Plot mic array
                            %hold on
                             plot(mposperim(1,:),mposperim(2,:),'<r');
                             xlabel('Meters')
                             ylabel('Meters')
                             title(['Perimeter Array beta = ' num2str(bata) ' SNR = ' num2str(nospow(segcnt))  ])
                             %climit = get(gca,'Clim');
                             set(gca,'Clim', [0 climit(2)])
                            hold off
                            pause
                        end
                       end  %  beta for loop
                       nst = nst+wininc;
                       ned = nst + seglen;
                       segcnt = segcnt+1;
              end  %  Sequence while loop
                        snrdum =  mean(sigpow)/mean(nospow);
              snrexp(k,kpn) = 10*log10(snrdum);
                 
                end  %  SNR Loop  
 end   %  File load loop of same type
         if pflag == 0
            save(fout{ktype}, 'tmpr', 'lerpr', 'nmpr', 'stpn', 'nompr', 'stpno', 'snrexp', 'filenames','batar')
         end
         snrexp = [];
 end
 exit
      