function output = srpvarstest( fsbin, nwlenin, fcin, fsresamplein, betain, frames, filename, fovp )
% This script is designed to test a large number of possibilities of
% variable combinations used in srp processing for a small number of
% frames. The start and end times of frames should be given in structure
% FRAMES, as well as the target x,y location(in meters) and tolerance(in 
% meters). FILENAME is a string which provides the name of the audiofile to
% be processed. FOVP is a base fovp structure required for srp processing.
% The structure for FRAMES is as follows:
%       frames= 1 [frame_start, frame_end, targ_x, targ_y, tolerance
%               2  "                                                "
%               n  "                                                "]
%       where N is the number of frames to look at (probably only a few),
%       FRAME_START and FRAME_END are the start and end times (in sec.) of
%       the frame, TARG_X and TARG_Y are the approximate x and y coords of
%       the source, and TOLERANCE is the maximum acceptable distance away
%       from the source that a detection can be.
% The structure for the test variables is as follows:
%       variables= [low_end, high_end, step]
%       where LOW_END and HIGH_END are the minimum and maximum values to be
%       tested, and step is the increment size between the two. 
%       NOTE: this function expects step to divide evenly into the
%       difference of LOW_END and HIGH_END.
% The output structure is formatted as follows:
%       output= [frame #1, fsb, nwlen, fc, fsresample, beta, det_flag
%                "                                                  "
%               ...
%               frame # n, "                                        "
%               ...                                                  ]
%       The format of the output structure is confusing to look at.
%       Basically, the function will create a set of variables to test for,
%       and run it on all the subsamples of the first frame before moving
%       to the second, and so on. DET_FLAG indicates whether or not a
%       source was detected in the tolerable region.
%
% When this script is finished, the output structure should be fed into
% another script to analyze its contents. The intention is to pull out the
% variable combinations that yield a detection in all frames examined
% (hence why it is important to only use a small number of frames where
% the source is known to make a sound). These can then be examined for
% trends in the effectiveness of certain values of a particular variable
% and for the effectiveness of certain variable combinations.

fovp.graphicout_flag=0; %uncomment to disable showing images (saves cpu time and allows to run
                        % silently in background.
% Set up nested for loops to examine all possibilities.
cnt=0;
output= [];
tru_detect=0;
fsbcombos= (fsbin(2)-fsbin(1))/fsbin(3)+1;
nwlencombos= (nwlenin(2)-nwlenin(1))/nwlenin(3)+1;
fccombos= (fcin(2)-fcin(1))/fcin(3)+1;
fsresamplecombos= (fsresamplein(2)-fsresamplein(1))/fsresamplein(3)+1;
betacombos= (betain(2)-betain(1))/betain(3)+1;
varcombos= floor(fsbcombos*nwlencombos*fccombos*fsresamplecombos*betacombos);
%elapsed_previous=0;
tic;
for nwlen= nwlenin(1):nwlenin(3):nwlenin(2)
    % with a new window length, set up the detection matrix. If last
    % subsample window is longer than the rest of frame, don't clip it
    det= [];
    for i=1:size(frames,1)
        start= frames(i,1);
        while start<frames(i,2)
            det= [det; start];
            start= start+nwlen;
        end
    end
    for fc= fcin(1):fcin(3):fcin(2)
        fovp.fc= [ fc, fovp.fc(2) ];
        for fsb= fsbin(1):fsbin(3):fsbin(2)
            for fsresample= fsresamplein(1):fsresamplein(3):fsresamplein(2)
                for beta= betain(1):betain(3):betain(2)
                    fovp.beta=beta;
                    strms= signifslistdet(det, filename, nwlen*fovp.fs, fovp, fsb, fsresample);
                    cnt= cnt+1;
                    % Now add to output structure for each frame
                    for i=1:size(frames,1)
                        % Determine value of detection flag for each
                        % subsample window. Similar loop to before.
                        start= frames(i,1);
                        while start<frames(i,2)
                            %first check to see if there is a corresponding
                            %entry in strms
                            det_flag=0;
                            strms_flag=0;
                            dist= frames(i,5)+.01;
                            if ~isempty(strms)
                                ind= find(start== strms(:,1));
                                if ~isempty( ind )
                                    %there is an entry, so now determine if distance is tolerable
                                    %note: it is possible for more than one
                                    %detection to occur in a window
                                    strms_flag=1;
                                    dist= sqrt( (strms(ind, 3)-frames(i,3)).^2+(strms(ind,4)-frames(i,4)).^2);
                                    det_flag = any(dist <= frames(i,5));
                                    tru_detect= tru_detect+ det_flag;
                                    dist=min(dist);
                                end
                            else
                            end
                            output= [output; i, fsb, nwlen, fc, fsresample, beta, det_flag, strms_flag, dist];
                            start= start+nwlen;
                        end
                    end
                    % save output for every 100 combinations tested
                    if mod(cnt, 100)==0
                        bakname= 'output_bak.mat';
                        save output_bak.mat output;
                        disp(['Saving output as ' bakname ' . . .']);
                        
                    end
                    elapsed=toc;
                    disp( ['Finished combination ' num2str(cnt) ' of ' num2str(varcombos) '.']);
                    disp( ['Time Elapsed: ' num2str(elapsed) ' seconds']);
                    disp( ['         or   ' num2str(elapsed/60) ' minutes']);
                    disp( ['         or   ' num2str(elapsed/60/60) ' hours']);
                    disp( ['         or   ' num2str(elapsed/60/60/24) ' days']);
                    %srp_time= (elapsed-elapsed_previous)/size(det,2);
                    eta= (elapsed/cnt)*(varcombos-cnt);
                    disp( ['Estimated Time Remaining: ' num2str(eta) ' seconds.']);
                    disp( ['                     or   ' num2str(eta/60) ' minutes.']);
                    disp( ['                     or   ' num2str(eta/60/60) ' hours.']);
                    disp( ['                     or   ' num2str(eta/60/60/24) ' days.']);
                    disp( 'NOTE: ETA is an overestimate early in the program');
                    disp( [num2str(tru_detect) ' true detections so far. . .']);
                end
            end
        end
    end
end


end

