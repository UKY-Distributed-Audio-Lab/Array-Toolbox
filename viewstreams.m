% This script converts the detection file to an auditory scene analysis
% (ASA) file and then creates a movie file (AVI) to display each of the 
% streams (location and time). 
%
% Required Programs:
%   convtpos.m      - Converts my detection file with time stamps to one
%                     compatible with Hari's program to create the ASA file
%   imsegm.m        - Hari's program to identify auditory streams from clustered
%                     detections - in time and space proximity
%   bfsegment.m     - Performs beamforming on detected streams and produces
%                     wave files for each beamformed stream
%   masksegment.m   - Performs auditory masking on detected streams, after
%                     beamforming
%   mkmov3d.m       - Converts ASA file to a stem plot movie of detected sources  
%   ListenMult.m    - Provides user interface allowing users to listen to
%                     found sources
%
% Written by Kevin D. Donohue (donohue@engr.uky.edu) June 2012
% Modified by Kirstin Brangers                       August 2012
%               Added function call to GUI - ListenMult({filename})
%               Modified save statement to include 'fovp'



%  After running ASA
moviemake_flag = 0;   %  Set flag 1 for a graphical output of detected streams 
% Get strings for file name and path to directory containing source
% detection information
[basenamep, pn] = uigetfile('*.mat');
basenam = basenamep(1:end-4);
load([pn,basenam '.mat'])   

%  Convert detection file format to one used by Hari's program
x = convtpos(strms1,fovp);
fovp.tinc = fovp.trez/2;
%  Now apply auditory stream identification using the proximity thresholds
timethr = 8;        % Time in seconds of tolerable nondetection (quiet period in speech) to link streams (1.1)
spacethr = .5;      %  Distance in meters where one source is allow(smaller detections are ignored) (.9)
toosmall = 0.15;       % Time in seconds to drop any stream shorter than this

tr = round(timethr/fovp.tinc) +1;       % Convert time to frame numbers
tsfrm = round(toosmall/fovp.tinc) +1;   % Convert time to frame numbers
stseg = imsegm(x,spacethr,tr,tsfrm);     % Hari's program to create true ASA matrix

strmnum = max(stseg(:,1));
% Create beamformed streams and save in separate files
[scriptout len] = bfsegment([basenam '.wav'],'',[1:strmnum(1)],stseg,fovp);
% Prompt user to listen to beamformed streams and delete streams that are
% multiples and noisy. Remaining streams will be masked.
bffname = [basenam 'bfscript.mat'];             
save(bffname, 'basenam', '', 'moviemake_flag','fovp', 'stseg', 'scriptout')     % Save beamformed data as .mat
ListenBeamform({bffname});                      % Call Beamform GUI 



% Create masked streams and save in separate files
% scriptm = masksegment(basenam,pn,scriptout,len);

%  Create graphical output ?
if moviemake_flag == 1
    %  Create avi file to for movie of detected streams (must create a unique
    %  avi file name each time.)  Will not overwrite
    fratein = 1/(fovp.trez/2);                  % Frames per second on input
    frateout = 15;                              % Output frame rate
    fname = [basenam '.avi'];                   % File name of movie
    sfactor = 200;                              %  Scale detection stat value to marker size
    mv = mkmov3d(stseg,fratein,frateout,fname,sfactor);
end
% Save scripting files for streams
 scfname = [basenam 'scripts.mat'];
% save(scfname,'fovp','stseg','scriptm') % Save file as .mat

% Pass data to GUI to display and listen to sources
% ListenMult({scfname});


