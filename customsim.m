% Simulates a room recording with adjustable reverb levels and mic
% placement.

nmics=8;
reverb=0.7; %adjust between 0.7 and 0.9
forv= [0,0,0;
       5,5,6]; %opposite corners of room
   
%ceiling planar mic array
% fom= [0,0,6;
%       2.5,2.5,6;
%       5,5,6];
% sp=5/(1+sqrt(2));
% mpos = regmicsplane(fom, sp);
mpos= [0 0 5 5 1.4645 1.4645 2.929 2.929;
             0 5 5 0 1.4645 2.929 2.929 1.4645;
             6 6 6 6 6      6     6     6]

% perimeter array
%mpos = regmicsperim2(fom, sp, num, stp);

% audio files to use
% man1.wav
% man2.wav
% man3.wav
% woman1.wav
% woman2.wav
% woman3.wav

% all the audio files should be at the same fs
sigarray=[];
[man1,fs]=audioread('man1.wav');
[man2,fs]=audioread('man2.wav');
[man3,fs]=audioread('man3.wav');

[woman1,fs]=audioread('woman1.wav');
[woman2,fs]=audioread('woman2.wav');
[woman3,fs]=audioread('woman3.wav');

sigarray= man1;
sigarray= [sigarray,man2];
sigarray= [sigarray,man3];
sigarray= [sigarray,woman1];
sigarray= [sigarray,woman2];
sigarray= [sigarray,woman3];

% Randomly generate signal x,y,z positions within certain ranges
% x,y: bound at least 0.5 from walls
%   z: bound between 1.2 and 2.3 meters
xrange= [0.5,4.5];
yrange= [0.5,4.5];
zrange= [1.2,2.3];

nsigs= size(sigarray,2);
sigx= rand(1,nsigs)*(xrange(2)-xrange(1))+xrange(1);
sigy= rand(1,nsigs)*(yrange(2)-yrange(1))+yrange(1);
sigz= rand(1,nsigs)*(zrange(2)-zrange(1))+zrange(1);

sigpos= [sigx; sigy; sigz];
% reflection coefficients
bs= [.90, .95 .90 .95 .2 .4];

[sigout,tax]= simarraysigim(sigarray,fs,sigpos,mpos,forv,bs,prop);

