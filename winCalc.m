function [wininc,seglen] = winCalc(env,fs,trez)

% This function computes the values of wininc and seglen from the input
%   parameters.  
%
%   [wininc,seglen] = winCalc(env,fs,trez)
%
% Required files:
% 1) mposanaly.m

c=env.c;
mposperim = env.mpos;

%  Find max distance (delay) over all mic pair, this represents
%  and upper bound on all requried relative delays when scanning over
%  the FOV

prs = mposanaly(mposperim,2);
maxmicd = max(prs(:,3));

%  Extra delay time in seconds for padding window of data
textra = maxmicd(1)/c;
winlen = ceil(fs*trez);
wininc = round(winlen/2);
seglen = ceil(fs*(textra))+winlen;
