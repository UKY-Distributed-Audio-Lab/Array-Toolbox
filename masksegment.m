function  scriptm = masksegment(fnbase,pn,scriptf,streams,numofstrm)
% 
%  This function opens a sequence of files containing beamformed signals with
%  base name FNBASE and integer numberings 01, 02, 03, .... at the end of
%  each base name all located in a directory identified by string PN.
%  The test file SCRIPTF contains the starting and stopping locations
%  of each beamformed signal file along with the time it became active in
%  in the from the start of the recording.  The output is another set of
%  files on which auditory masking has been applied.
%
%       scriptm = masksegment(fnbase,pn,scriptf,numofstrm)
%   
%  Inputs:
%   FNBASE      =>  String identifying basename of files containing beamformed
%                    signals from each detected source
%   PN          =>  String with path name to directory containing beamformed file
%                    and also the directory to where the masked signals will also
%                    be written.
%   SCRIPTF     =>  Cell array containing the following information for each
%                    source:
%                    [segment#, xs,ys,zs, xe, ye, ze, time start, filename]
%   STREAMS     =>  The streams chosen by the user to be masked
%   NUMOFSTRM   =>  Number of valid streams detected during beamforming
%                    segment
%
% Output:
%   SCRIPTM     =>  Cell array the same as SCRIPTF
%                    File will be written to directory PN with the same 
%                    basename FNBASE and sequence number followed by 
%                    an 'm' to denote a masked sequence.
% 
% Required Functions:
%               tfmask.m
% 
% Written by Kevin D. Donohue (donohue@engr.uky.edu)    July 2012
% Modified by Kirstin Brangers                          August 2012
%           Added NUMOFSTRM as input parameter to determine
%               how many valid streams were found during beamforming
%           Added STREAMS as input parameter to determine user chosen
%               streams
%           Added/Modified lines 53-107 to account for checking existence of all
%               beamformed files in scriptf
 

% Initialized fileflag to determine number of beamformed signal in directory

scriptm = scriptf;      % Assign input scriptf to output scriptm for modifying later

dt = .05;               % Window length for spectrogram (in milliseconds)
%gamma = 0;             % Threshold for binarizing the mask
%mmode = 'DIFF';
dffac = 2;              % Zero padding ratio

[row col] = size(scriptf);
% Delete streams from scriptm that were not chosen by user
temp= zeros(row,1);       % Temp vector to hold values of scriptf
for i=1:row
    % Convert cell array to numeric array
        temp(i,1) = cell2mat(scriptf(i,1));
end

tempm = cell(numofstrm,9);
% tempscript = cell2mat(scriptf(1,:));
tempscript = temp;
hwait = waitbar(0,['Masking streams!  On stream ' int2str(0) ' of ' int2str(numofstrm)]);
for st=1:numofstrm
    strm = streams(st);   % Chosen stream number
%     tempscript = cell2mat(scriptf(1,:));
    a = find(tempscript(:,1) == strm);     % Find stream number in scriptm equal to chosen stream
    if ~isempty(a)      % If variable a is not empty, copy stream details over
        tempm{st,1} = scriptm{a,1};
        tempm{st,2} = scriptm{a,2};
        tempm{st,3} = scriptm{a,3};
        tempm{st,4} = scriptm{a,4};
        tempm{st,5} = scriptm{a,5};
        tempm{st,6} = scriptm{a,6};
        tempm{st,7} = scriptm{a,7};
        tempm{st,8} = scriptm{a,8};
        tempm{st,9} = scriptm{a,9};
    end    
end

% Save temp as scriptm = scriptm now has user chosen streams
scriptm = tempm;

fileflag = 2;
fcnt = 1;                           % Initialize file name count
count = 1;                          % Initialize count for existing files
for i = 1:row                       % Cycle through all beamformed files
        % Create string for new file name
        if fcnt < 10
            fntest = [pn, fnbase, '0' int2str(fcnt), '.wav'];   % Stream file name
        else
            fntest = [pn, fnbase,  int2str(fcnt), '.wav'];      % Stream file name
        end

        fileflag = exist(fntest,'file');  % Test to see if it exists

        if fileflag == 2
            % If so, add it to list of files to open and test for next one
            wfiles{count} = fntest; 
            count = count+1;
        end

        % Test for next file
        fcnt = fcnt+1;       % Increment file name count
end
% count = count-1;          % Correct file count number for last failed test
[r count] = size(wfiles);

for k=1:count
        [ydum, fs] =  audioread(wfiles{k});      % Read in wave files
        bsig(:,k) = ydum;  %  assign to array
end
[siglen,spnum] = size(bsig);
for k=1:count
    waitbar(k/count,hwait,['Masking streams!  On stream ' int2str(1) ' of ' int2str(count)])
    foutname = [wfiles{k}(1:end-4), 'm.wav'];   % Append 'm' to output filename
    scriptm{k,9} = foutname;
    sigout = tfmask(bsig,k, dt,dffac,fs);
    audiowrite(foutname, sigout / (10*std(sigout)+eps),fs)
end

% Check streams: If columns 1-7 are zeros, then delete streams
valid = cell(1,9) ;     % Initialize output cell array
temp= zeros(1,8);       % Temp vector to hold values of scriptout
count = 1;              % Initialize stream number 
for i=1:numofstrm
    % Convert cell array to numeric array
    for n=1:7
        temp(1,n) = cell2mat(scriptm(i,n));
    end
    % If cols1-7 are not all zeros, copy over to output cell array VALID
    if sum(temp) ~= 0
        valid{count,1} = temp(1,1);
        valid{count,2} = temp(1,2);
        valid{count,3} = temp(1,3);
        valid{count,4} = temp(1,4);
        valid{count,5} = temp(1,5);
        valid{count,6} = temp(1,6);
        valid{count,7} = temp(1,7);
        valid{count,8} = temp(1,8);
        valid{count,9} = scriptm{i,9};
        
        count = count+1;        % Increment stream number
    end
end
close(hwait)

scriptm = valid;








