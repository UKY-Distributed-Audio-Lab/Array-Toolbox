function track_video(distances, circx, circy, angles, timestamps, video)

% uses drawTrack and some data to track head movement. creates an avi video
% of the results.
%
% Written by Kyle Ihli    July 2014
tic
elapsedTime=0;
frames=[];
writerObj= VideoWriter('rat_track-test2.avi');
open(writerObj);
for i=1:size(distances,1)
   figure(1);
   if floor(timestamps(i)/1000*30)<= size(timestamps,1)
    image(video{floor(timestamps(i)/1000*30)});
   end
   hold on;
   draw_circle(353,255,363/2);
   drawTrack_video(circx(1,:), circy(1,:), distances(i,:)/1000, angles(i,:));
   hold off;
   frame= getframe;
   writeVideo(writerObj, frame);
   if i~=size(distances,1)
       elapsedTime=toc;
       %pause((timestamps(i+1)/1000-elapsedTime));
       %pause(.1);
   end  
end

close(writerObj);