function track(distances, circx, circy, angles, timestamps)

% uses drawTrack and some data to track head movement
%plot(circx(1,:), circy(1,:));
for i=1:size(distances,1)
    tic
   figure(1);
   draw_circle(0,0,.2096);
   hold on;
   drawTrack_new(circx(1,:), circy(1,:), distances(i,:)/1000, angles(i,:));
   hold off;
   if i~=size(distances,1)
       elapsedTime=toc;
       pause((timestamps(i+1)-timestamps(i))/1000-elapsedTime);
       %pause(.1);
   end  
end