function draw_circle(x,y,r, color)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step; bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
%
% Written by Kyle Ihli    July 2014
if nargin==3
   color= 'blue'; 
end
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp, color);
end