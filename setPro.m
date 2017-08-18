function pro = setPro(trez,rez,beta,fs,rad,lang)
% function setPro returns data structure of signal processing parameters
% 	pro = setPro(trez,rez,beta,fs,rad)
% input paramaters are:
%	1) trez - processing time window in seconds
%	2) rez - spacial resolution window
%	3) beta - whitening parameter
%	4) fs - signal processing sampling rate
%	5) lang - language being used (octave(1)/matlab(0))
%	6) rad - rad value for dent3d (optional, only used for 3D)

if nargin == 5
	pro = struct('trez',[trez],'rez',[rez],'beta',[beta],'fs',[fs],'lang',[lang]);
else
	pro = struct('trez',[trez],'rez',[rez],'beta',[beta],'fs',[fs],'lang',[lang],'rad',[rad]);
end
