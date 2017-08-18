%  This script demonstrates the use of functions associated with microphone
%  array geometry.  The functions needed to run this script are:
%
%     regmicsline.m, regmicsplane.m, and mposanaly.m
% 
%   Written by Kevin D. Donohue (donohue@engr.uky.edu) July 2005
%
%  Generate and plot the coordinates of a linear mic array over a 2-D space

fom = [-2 2; 0 1];  %  Define mics over x dimension from -2 to 2 meters
                    %  and 0 to 1 over the y dimension
sp = .3;  % define the spacing increment

mpos = regmicsline(fom,sp);  % Generate coordinates of mic position
figure(1)
set(gcf,'Position', [170   315   560   420])
plot(mpos(1,:), mpos(2,:), 'x')
xlabel('Meters')
ylabel('Meters')
title(['Microphone Positions with Spacing ' num2str(sp) ' in 2-D space'])
disp('Hit any key to continue')
pause
fom = [-2 2 2; 0 0 4];  %  Define mic plane over x dimension from -2 to 2 meters
                        %  and 0 to 4 over the y dimension
sp = .3;  % define the spacing increment

mpos = regmicsplane(fom,sp,'Rect');  % Generate coordinates of mic positions
                                     % with rectilinear geometry
figure(2)
set(gcf,'Position', [198   285   560   420])
plot(mpos(1,:), mpos(2,:), 'x')
xlabel('Meters')
ylabel('Meters')
title(['Microphone Positions with Rectangular Spacing ' num2str(sp) ' in 2-D space'])
disp('Hit any key to continue')
pause
mpos = regmicsplane(fom,sp,'Hex');  % Generate coordinates of mic positions
                                     % with Hexagonal geometry           
figure(3)
set(gcf,'Position', [235   247   560   420])
plot(mpos(1,:), mpos(2,:), 'x')
xlabel('Meters')
ylabel('Meters')
title(['Microphone Positions with Hexagonal Spacing ' num2str(sp) ' in 2-D space'])
disp('Hit any key to continue')
pause
%  Now generate 10 mics randomly distribued over a 4x4 plane and rank every
%  pair of mics by thier distance from each other

mcnum = 8;      %  Number of random mics to generate
n = 2;  %  number of mics in subsets for analysis
rng = [0 4];      %  Square range over which to distribute
dst = (rng(2)-rng(1));   %  Compute Square distance
mpos = rng(1)+rand(2,mcnum)*dst;     % Generate random number over array

figure(4)
set(gcf,'Position', [272   209   560   420])
plot(mpos(1,:), mpos(2,:), 'x')  %  Plot mic position
axis([rng, rng])  % set plot axis to full range
% label mics with index from MPOS
hold on
offs = dst*.01;
for k=1:mcnum
    text(mpos(1,k)+offs, mpos(2,k),int2str(k));  % place mic number near mic  
end
hold off
xlabel('Meters')
ylabel('Meters')
title(['Random Positions of mics over a ' num2str(dst) ' by ' num2str(dst) ' area'])
disp('Hit any key to continue')
pause
%  Generate analysis table
pa = mposanaly(mpos,n);   %  find all possible combinations of mics taken "n" 
                          %   at a time and compute statistics
%  Sort based on minimum distance (first column after the mic indices)
[rnk, indx] = sort(pa(:,n+1));
pamind = pa(indx,1:n+1);   %   Generate a sorted matrix with minimum distance
tabout = num2str(pamind);  %   Output number is table fashion in a dialog box
[r,c] = size(tabout);
%  Create header labels
nms = find(tabout(1,:) ~= ' ');
stp =[1];
for k=2:length(nms)
    if nms(k) ~= (nms(k-1)+1)
        stp = [stp, nms(k)-1];
    end
end
header = [];
for k=1:c
    header = [header, ' '];
end
header(stp(1):stp(1)+5-1) = 'Mic 1';
header(stp(2)-2:stp(2)+5-3) = 'Mic 2';
header(stp(3)-1:stp(3)+8-2) = 'Distance';
%  Output message box
msgbox([header; tabout], 'Table of Ranked Mic Distance Pairs')
