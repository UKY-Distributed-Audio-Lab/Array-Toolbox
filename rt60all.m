clear
close all
testrt60pf
xx = sort(vec2a60);
rt60pf = mean(xx(1:6));
xx = sort(vec2anpow);
nospf = mean(xx(1:6));
xx = sort(vec2aspow);
sigpf = mean(xx(1:6));

% close all
% testrt60p
% xx = sort(vec2a60);
% rt60p = mean(xx(1:6));
% xx = sort(vec2anpow);
% nosp = mean(xx(1:6));
% xx = sort(vec2aspow);
% sigp = mean(xx(1:6));

close all
testrt60f
xx = sort(vec2a60);
rt60f = mean(xx(1:6));
xx = sort(vec2anpow);
nosf = mean(xx(1:6));
xx = sort(vec2aspow);
sigf = mean(xx(1:6));
