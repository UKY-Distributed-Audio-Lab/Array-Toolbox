function [ outstruct ] = powercalib( files, fovp, nwlen, posthresh, posnum, det, sourcepos)
% This function expects N many input audio files, all of which contain data
% recorded from a sound source in the same position for all files.
% For each file, the script will roll through a set of possible microphone
% positions and record the power each time. At the end, the sound source
% positions corresponding to the top 5% coherent power are averaged and
% used as the sound source position for that file. At the end, the
% positions from each file are averaged together to yield CALIBPOS.

% FILES is a structure containing the string filenames of the audio files
%  used for calibration.
% NWLEN is the window length in seconds.
% EXPECTPOS is the user input expected positions of the microphones
% POSTHRESH is the maxiumum allowed distance from the expectpos that any
%  given mic should be at
% POSNUM is the number of possible positions that should be tested for each
% microphone.
% DET is a detection matrix
if nargin <7
    sourcepos=[0,0];
    if nargin <6
        det=[50e-3];
        if nargin < 5
            posnum=10;
            if nargin < 4
                posthresh=.01;
                if nargin < 3
                    nwlen=10e-3;
                end
            end
        end
    end
end
outstruct=[];
expectpos= fovp.mp;
nmics= size( expectpos, 2); % # of mics is equal to # of columns

% change count sizes for more precise calibration
xcnt= posnum;
ycnt= posnum;
xstep= (2*posthresh)/(xcnt-1);
ystep= (2*posthresh)/(ycnt-1);

% convert sourcepos to grid space

sourceposx= (sourcepos(1,1)-fovp.sgrid{1,1}{1,1}(1))/(fovp.sgrid{1,1}{1,1}(size(fovp.sgrid{1,1}{1,1},2))-fovp.sgrid{1,1}{1,1}(1))*size(fovp.sgrid{1,1}{1,1},2);
sourceposy= (sourcepos(1,2)-fovp.sgrid{1,1}{1,2}(1))/(fovp.sgrid{1,1}{1,2}(size(fovp.sgrid{1,1}{1,2},2))-fovp.sgrid{1,1}{1,2}(1))*size(fovp.sgrid{1,1}{1,2},2);
sourcepos=[sourceposx, sourceposy];
allcoords= [];
%for each microphone
for mic=1:nmics
% Monte Carlo Generation
    %first step generates integers that correspond to number of steps to
    %take, ranging from 0 to posnum-1. This generation allows for each x
    %and y coordinate to be chosen about 1000 times.
    xcoords= randi(posnum,posnum*1000,1)-1;
    ycoords= randi(posnum,posnum*1000,1)-1;
    xcoords= expectpos(1,mic)-posthresh+xcoords*xstep;
    ycoords= expectpos(2,mic)-posthresh+ycoords*ystep;

    allcoords= cat(2,allcoords, [xcoords ycoords]);
end

%for each file in FILES
for ifile=1:size(files,2)
    %for each possible combo of positions, compute the power of the
    %calibration image at 0,0(srp)
    % CURRENTLY ONLY FUNCTIONS FOR 4 MICS
    savestruct=[];
    iterations= size(files,2)*size(allcoords,1);
    tic;
    for i=1:size(allcoords,1)
        fovp.mp= [allcoords(i,1) allcoords(i,3) allcoords(i,5) allcoords(i,7);
                  allcoords(i,2) allcoords(i,4) allcoords(i,6) allcoords(i,8);
                  0              0              0              0             ];
        [~, pow]=signifslistdetpow(det, char(files(ifile)), nwlen, fovp, 5e3, 44.1e3, sourcepos);
        
        savestruct=[ savestruct; fovp.mp(1,:) fovp.mp(2,:) pow];
        
        if mod(i,100)==0
           save([char(files(ifile)) '_savestructbackup.mat']);
        end
        t=toc;
        disp('. . .');
        disp( [int2str(t/60/60) ' hours elapsed for file' int2str(ifile) '.']);
        tpi=t/i;
        ileft= iterations-i;
        tleft= ileft*tpi;
        disp( [int2str(tleft/60/60/4) ' hours remaining for file' int2str(ifile) '.']);
    end
    
    outstruct= [outstruct savestruct];
    
end


end

