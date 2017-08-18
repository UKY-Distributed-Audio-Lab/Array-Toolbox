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

mpos = regmicsplane(fom,sp,'Hex');  % Generate coordinates of mic positions
                                     % with Hexagonal geometry           
figure(3)
set(gcf,'Position', [235   247   560   420])
plot(mpos(1,:), mpos(2,:), 'x')
xlabel('Meters')
ylabel('Meters')
title(['Microphone Positions with Hexagonal Spacing ' num2str(sp) ' in 2-D space'])



