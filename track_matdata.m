function track_matdata(strms1,micx,micy)

% uses drawTrack and .mat file data to track head movement.
%
% Written by Kyle Ihli    July 2014
maxintens= max(strms1(:,2));
minintens= min(strms1(:,2));
slope= 65/(maxintens-minintens); %for linear dot size weighting
for i=1:size(strms1,1)
    tic
   figure(1);
   draw_circle(0,0,.2096);
   hold on;
   for k=1:4
       draw_circle(micx(k),micy(k),.02);
   end
   j=i+1;
   while (strms1(j,1) == strms1(i,1))
       j=j+1;
   end
   drawTrack_matdata(strms1(i:j,3), strms1(i:j,4), strms1(i:j,2), maxintens, minintens, slope);
   %hold off;
   if i~=size(strms1,1)
       elapsedTime=toc;
       pause((strms1(i+1,1)-strms1(i,1))/1000-elapsedTime); %only for viewing in real time
       %pause(.1);
   end  
end