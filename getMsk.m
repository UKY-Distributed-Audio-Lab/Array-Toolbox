function Msk = getMsk(spectar,specinter,gamma)
% 
% Creates a Time-Frequency binary map by taking element wise 
% difference between a estimate of target and interference.
% The spectrograms are normalized for unit power.
%
% syntax : MASK = getMsk(SPECTAR,SPECINTER,GAMMA)
%
% Inputs:
%       SPECTAR     : Estimated Target Spectrogram (TF representation)
%       SPECINTER   : Estimated Interference Spectrogram
%       GAMMA       : (optional) Binarizing threshold for power ratio (default 1). 
%                    higher the value more TF elements will be dropped
% Output:
%       MASK        : Binary mask. 1  to let TF pass through 0  to drop it.
%
% Written by Harikrishnan Unnikrishnan      April 5 2010
% Modified by Kevin D. Donohue              August 15 2011
% Modified by Kevin D. Donohue   August 2014 to remove thresholding modes

% Initializing
Msk = ones(size(spectar));
[r c] = size(spectar);

if nargin == 2
    gamma = 1;  %  Default
end

% compute magnitude of specta 
tartemp  =  abs(spectar);     
intertemp  = abs(specinter);  

% Compare power in time frequency bins
%  KD modified to prevent divide by zero
dimage(intertemp ~= 0) = tartemp(intertemp ~= 0)./(intertemp(intertemp ~= 0));    
Msk(dimage<=gamma) = 0;

%  KD addded this to drop isolated (impulse) spectral regions since they
%  may create significant musical noise (i.e. tonal blips)
for k=2:length(Msk)-1
    if Msk(k-1) == 0 && Msk(k+1) == 0
        Msk(k) = 0;
    end
end