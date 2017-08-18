function drawTrack(x,y,r,angles)
% Takes in the center x,y coordinates of circles, the angles given by tracking software
% and their radii r, then
% determines for what points the circles 'intersect' (i.e. when points from
% all circles are within THRESH of each other. It then draws these points.
%
% Function based on data collected from rat vocalizations in a cage.
% Written by Kyle Ihli    July 2014
xtoplot=[];
ytoplot=[];
xpairs=[];
ypairs=[];
xbpairs=[];
ybpairs=[];
xbtoplot=[];
ybtoplot=[];

% find distance from mic to center of body
% b= distance from mic to body
% a= length from center of body to head (assumed to be 8.2 cm)
% theta= angle from body-head line to mic-body line
% d= distance from mic to head
% b= sqrt(a^2*cos^2(theta)-a^2+d^2) +a*cos(theta)
a= 8.2e-2;
rads= angles*pi/180;

%for each mic
b= abs(sqrt( a^2*cos(rads).^2-a^2+r.^2 ) +a*cos(rads));

% find all intersections
% between 1 and others
for i=2:4
   [xout,yout]= circcirc(x(1,1),y(1,1),r(1,1), x(1,i),y(1,i),r(1,i));
   [xbout, ybout]= circcirc(x(1,1),y(1,1),b(1,1), x(1,i),y(1,i),b(1,i));
   if xout ~= NaN
        xpairs= [xpairs; xout];
        ypairs= [ypairs; yout];
   end
   if xbout~=NaN
       xbpairs= [xbpairs; xbout];
       ybpairs= [ybpairs; ybout];
   end
end

%between 2 and 3-4
for i=3:4
    [xout, yout]= circcirc( x(1,2),y(1,2),r(1,2), x(1,i),y(1,i),r(1,i));
    [xbout, ybout]= circcirc(x(1,2),y(1,2),b(1,2), x(1,i),y(1,i),b(1,i));
    if xout ~= NaN
        xpairs= [xpairs; xout];
        ypairs= [ypairs; yout];
    end
    if xbout~=NaN
       xbpairs= [xbpairs; xbout];
       ybpairs= [ybpairs; ybout];
   end
end

%between 3 and 4
[xout, yout]= circcirc( x(1,3), y(1,3), r(1,3), x(1,4), y(1,4), r(1,4));
[xbout, ybout]= circcirc( x(1,3), y(1,3), b(1,3), x(1,4), y(1,4), b(1,4));
if xout ~=NaN
    xpairs= [xpairs; xout];
    ypairs= [ypairs; yout];
end
if xbout ~=NaN
    xbpairs= [xbpairs; xbout];
    ybpairs= [ybpairs; ybout];
end

xpairs= [xpairs(:,1);xpairs(:,2)];
ypairs= [ypairs(:,1);ypairs(:,2)];
xbpairs= [xbpairs(:,1);xbpairs(:,2)];
ybpairs= [ybpairs(:,1);ybpairs(:,2)];

% now find where x,y pairs are within threshold of others
thresh= 3e-2;

for i=1:size(xpairs, 1)
       if sqrt(xpairs(i)^2+ypairs(i)^2) <.2096
          if any(abs(xpairs(i)-xpairs(i+1:end)) < thresh) && any(abs(ypairs(i)-ypairs(i+1:end)) < thresh)
            xtoplot= [xtoplot; xpairs(i)]; 
            ytoplot= [ytoplot; ypairs(i)];
          end
       end
end

for i=1:size(xbpairs, 1)
       if sqrt(xbpairs(i)^2+ybpairs(i)^2) <.2096
          if any(abs(xbpairs(i)-xbpairs(i+1:end)) < thresh) && any(abs(ybpairs(i)-ybpairs(i+1:end)) < thresh)
            xbtoplot= [xbtoplot; xbpairs(i)]; 
            ybtoplot= [ybtoplot; ybpairs(i)];
          end
       end
end

%refine results by "boxing" the average and excluding all points outside
%box (using threshold)


% for i=1:size(xbtoplot,1)
%     
% end

%now, rotate data to match up with axes
% xhold= xtoplot;
% xholdb= xbtoplot;
% xtoplot= -ytoplot;
% xbtoplot= -ybtoplot;
% ytoplot= xhold;
% ybtoplot= xholdb;

avgx= mean(xtoplot);
avgy= mean(ytoplot);
avgxb= mean(xbtoplot);
avgyb= mean(ybtoplot);

scatter(xtoplot,ytoplot);
scatter(xbtoplot, ybtoplot, 50, 'red');

plot( [avgx, avgxb], [avgy, avgyb]);