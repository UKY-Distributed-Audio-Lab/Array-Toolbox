function [cw] = SpeedOfSound(temp,RH,p)
% This function computes the speed of sound given room temperature,
% humidity, and air pressure
%
%   [cw] = SpeedOfSound(temp,RH,p)
%
% with inputs, TEMP the temperature in Centigrade, RH the relative humidity in percent
% and P the air pressure in mmHg.
% Forumla used:
% cw=(331.45*sqrt(1+(temp/273)))*(4.5513*sqrt(((7+((0.01*RH*vp)/p))/(5+((0.01*RH*vp)/p)))/(29-11*((0.01*RH*vp)/p)))
% simplified is:
% cw=91.30029124144035*sqrt(273+temp)*sqrt((p*(7.*p+0.01*RH*vp))/((29.*p-0.11*RH*vp)*(5.*p+0.01*RH*vp)))
% where p is the ambient pressure, t is the ambient temperature, RH is the relative humidity, and vp is
% the vapor pressure of water as a function of t. 
%
%  Written by Jon Soli 6/2010 and Kevin D. Donohue (donohue@engr.uky.edu)
%

cw=91.30029124144035*sqrt(273+temp)*sqrt((p*(7.*p+0.01*RH*vp(temp)))/((29.*p-0.11*RH*vp(temp))*(5.*p+0.01*RH*vp(temp))));

function [vapor] = vp(temp)
vapor=(10^(((-7.90298)*(373.16/(temp+273.16)-1))+(5.02808*log10(373.16/(temp+273.16)))-(1.3816*10^(-7)*(10^(11.334*(1-(temp+273.16)/373.16))-1))+(8.1328*10^(-3)*(10^(-3.49149*(373.16/(temp+273.16)-1))-1))))*(29.9213);


