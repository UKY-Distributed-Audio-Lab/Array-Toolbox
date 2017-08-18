%  This script demonstrates the function that places microphones around a
%  perimeter of a rectangular area.  It requires the use of the function
%  "regmicsperim.m."  The function requires a spacing and a starting point with
%  respect to the initial corner.  The values are derived in the script to maintain
%  a symetery around the perimeter.  A special formula was derived for this script
%  so that as mics are placed around the rectangle at the corners, the distance is
%  maintained diagonally from the mic on one side to the next. 
%  
%
%  written by Kevin D. Donohue (donohue@engr.uky.edu) October 2005

%  Select number of microphones
micnum =8;
% Opposite Corner Points of room (Columns of x,y,z) 
froom = [-3.5 -4 0; 3.5 4 3]';
%  Convert opposite corners to vertcies around rectangle halfway up the walls
v = [froom(1:2,1), [froom(1,1); froom(2,2)], froom(1:2,2), [froom(1,2); froom(2,1)]]; % (x-y plane coordinates) 
v = [v; ones(1,4)*((froom(3,2)-froom(3,1))/2 + froom(3,1))];  %  Set all z coordinates to half way up wall

%  Compute spacing that corresponds to an equal spacing of perimeter array microphones
%   (Specially derived formula)
spm = 2*(norm(v(:,2)-v(:,1),2)+norm(v(:,3)-v(:,2),2))/(micnum-4 +4*sqrt(2));
%  Starting position for first element along wall (so that if mics are a multiple of
%  4 it will look symetric on each wall) (Another specially derived formula)
stp = (spm/sqrt(2))*(norm(v(:,2)-v(:,1),2)/norm(v(:,3)-v(:,2),2));
%  Use opposite corners of rectangular perimeter, spacing, and starting
%  point on a perimeter wall as input arguments
mposperim = regmicsperim(v(:,[1,3]),spm,stp);
% Plot results mic postions
figure()
plot(mposperim(1,:),mposperim(2,:),'<r');
xlabel('Meters')
ylabel('Meters')
title(['Perimeter Array with Spacing = ' num2str(spm)])
axis([-4.5 4.5 -4.5 4.5]);
