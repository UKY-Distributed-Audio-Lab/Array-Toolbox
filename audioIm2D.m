function mleim = audioIm2D(filename, env, pro, sflag)

% This function is used to create a 2D image array with respect to
% the frame number based on the inputs.
%
% Required files in the path to use this function:
%   1) flattap.m        (Array Toolbox)
%   2) srpframenn.m     (Array Toolbox)
%   3) whiten.m         (Array Toolbox)
%   4) winCalc.m        (Array Toolbox)
%
%   Inputs:
%       1) filename - wav file array from which to create to acoustic image
%
%       2) env - data structure describing the environmental parameters
%          fields include:
%               A) c - speed of sound
%               B) mpos - microphone positions
%               C) fov - field of view
%
%       3) pro - signal processing parameters
%          fields include:
%               A) trez - processing time window in seconds
%               B) rez - spacial resolution window
%               C) beta - whitening parameter
%               D) fs - signal processing sampling rate
%		E) lang - language being used (octave(1)/matlab(0))
%
%       4) sflag - set to 1 to save after each processed frame (optional -
%          default is 0)
%
%   Outputs:
%       1) mleim - 3 dimensional matrix, where the first 2 dimesions
%          contain the image information, while the 3rd contains the frame
%          number.
%          This matrix gets stored in a .mat file with the same name
%          as the input wav file.
%
%   Written by Kevin Donohue (2006)
%   Edited by Asher Finkel, Harikrishnan Unnikrishnan (5/28/08)

% Extract data from the structures env and pro

c=env.c;
fov=env.fov';
mposperim = env.mpos;

fs=pro.fs;
rez=pro.rez;
trez=pro.trez;
beta=pro.beta;
lang=pro.lang;

%  Set sflag to 0 when these parameters are not passed in
if nargin < 4
    sflag=0;
end

savename=[filename(1:end-3),'mat'];
gridax = {[fov(1,1):rez:fov(1,2)], [fov(2,1):rez:fov(2,2)], [fov(3,1):rez:fov(3,2)]};

[y,fso]=wavread(filename);
tax = [0:length(y)-1]/fs;
[rsam, cmics] = size(y);

% Calculate wininc and seglen
[wininc,seglen]=winCalc(env,fs,trez);

% Create a tapering window for all mic channels
tapwin = flattap(seglen,20);
winrec = tapwin*ones(1,cmics);

% Create a high pass filter
[b,a] = butter(4,200/(fs/2), 'high');

for ii = 1: cmics
	y(ii,:) = filtfilt(b,a,y(ii,:));
end;

sst = 1;  % Initial signal window index
sed = sst + seglen-1;

% iterate for different values of beta
for batak = 1:length(beta)
    bata = beta(batak);
    %  Create Images from RF Perimeter Array
    count=1;
    while sed<=length(y)
        % Taper window
        sigout = y([sst:sed],:).*winrec;

        % Whiten signal in the window
        sigout = whiten(sigout,bata);

        % Calculate magnitude for each pixel
        mleim(:,:,count) = srpframenn((sigout), gridax, mposperim, fs, c, trez);

        % Compute the window range of the next iteration
        sst = sst+wininc;
        sed = sst+seglen-1;

        if sflag == 1
		if lang == 0
			save(savename,'mleim');
		else
        		save('-mat-binary',[savename,'mleim']);
		end
        end

        count=count+1;
    end  %  while frame loop
end  %  Beta Loop

if lang == 0
	save(savename,'mleim');
else
	save('-mat-binary',[savename,'mleim']);
end
