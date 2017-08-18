function S = sii(varargin)

% "Methods for calculation of the Speech Intelligibility Index" 
% (ANSI S3.5-1997)
%
% MATLAB implementation of Section 4.
%
% Note: The remaining sections of the standard, which provide means
%       to calculate input parameters required by the "core" SII
%       procedure of Section 4, are implemented in seperate scripts:
%
%       Section 5.1 in script Input_5p1.m "method based on the direct
%          measurement/ estimation of noise and speech spectrum levels
%          at the listener's position"
%       Section 5.2 in script Input_5p2.m "method based on MTFI/CSNSL
%          measurements at the listener's position"
%       Section 5.3 in script Input_5p3.m "method based on MTFI/CSNSL
%          measurements at the eardrum of the listener"
%
% Parameters are passed to the procedure through pairs of "identifier"
% and corresponding "argument" Identifiesrs are always strings. Possible
% identifiers are:
%
%     'E' Equivalent Speech Spectrum Level (Section 3.6 in the standard)
%     'N' Equivalent Noise Spectrum Level (Section 3.15 in the standard)
%     'T' Equivalent Hearing Threshold Level [dBHL] (Section 3.23 in the standard)
%     'I' Band Importance function (Section 3.1 in the standard)
%
% Except for 'E', which must be specified, all parameters are optional. If an
% identifier is not specified a default value will be used. Paires of identifier
% and argument can occur in any order. However, if an identifier is listed, it
% must be followed immediately by its argument.
%
% Possible arguments for the identifiers are:
%
%   Arguments for 'E': 
%           A row or column vector with 18 numbers stating the Equivalent Speech Spectrum Levels in dB in bands 1 through 18. 
%
%   Arguments for 'N': 
%           A row or column vector with 18 numbers stating the Equivalent Noise Spectrum Levels in dB in bands 1 through 18.
%           If this identifier is omitted, a default Equivalent Noise Spectrum Level of -50 dB is assumed in all 18 bands (see note in Section 4.2).
%
%   Arguments for 'T': 
%           A row or column vector with 18 numbers stating the Equivalent Hearing Threshold Levels in dBHL in bands 1 through 18.
%           If this identifier is omitted, a default Equivalent Hearing Threshold Level of 0 dBHL is assumed in all 18 bands .
%
%   Arguments for 'I': 
%           A scalar having a value of either 1, 2, 3, 4, 5, 6, or 7. The Band-importance functions associated with each scalar are
%                       1:	Average speech as specified in Table 3 (DEFAULT)
%		                2:	various nonsense syllable tests where most English phonemes occur equally often (as specified in Table B.2)
%		                3:	CID-22 (as specified in Table B.2)
%		                4:	NU6 (as specified in Table B.2)
%		                5:	Diagnostic Rhyme test (as specified in Table B.2)
%		                6:	short passages of easy reading material (as specified in Table B.2)
%		                7:	SPIN (as specified in Table B.2)
%
%
% The function returns the SII of the specified listening condition, which is a value in the interval [0, 1].
%
%
% REMINDER OF DEFINITIONS & MEANINGS: 
% 
% Equivalent Speech Spectrum Level, E-prime
%           The SII calculation is based on free-field levels, even though the quantity relevant for perception and intelligibility 
%           is the level at the listener's eardrum. 
%           The Equivalent Speech Spectrum Level is the speech spectrum level at the center of the listener's
%           head (when the listener is temporarily absent) that produces in an average human with unoccluded ears an eardrum speech level equal
%           to the eardrum speech level actually present in the listening situation to be analyzed. 
%           Before the SII can be applied to a given listening situation, the corresponding Equivalent Speech 
%           Spectrum Level must be derived. For example, when speech is presented over insert earphones (earphones inside the earcanal), 
%           only the speech spectrum level at the eardrum is known. Using the inverse of the freefield-to-eardrum transfer function (Table 3 of the standard)
%           this eardrum level must be "projected" into the freefield, yielding the Equivalent Speech Spectrum Level. 
%           Sections 5.1, 5.2, and 5.3 of the standard give three examples of how to derive the Equivalent Speech Spectrum Level from physical measurements.
%           The standard allows the use of alternative transforms, such as the one illustrated above, where appropriate.
%
% Equivalent Noise Spectrum Level, N-prime
%           Similar to the Equivalent Speech Spectrum Level, the Equivalent Noise Spectrum Level is the noise spectrum level at the center of the listener's
%           head (when the listener is temporarily absent) that produces an eardrum noise level equal to the eardrum noise level actually present in the
%           listening situation to be analyzed. Sections 5.1, 5.2, and 5.3 give three examples of how to derive the Equivalent Speech Spectrum
%           Level from physical measurements.
%
% Hannes Muesch, 2003 - 2005


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERIFY INTEGRITY OF INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,Nvar] = size(varargin);
CharCount = 0;
Ident = [];
for k = 1:Nvar
    if ischar(varargin{k})&(length(varargin{k})==1)
        CharCount = CharCount + 1;
        Ident = [Ident; k];
    end
end
if Nvar/CharCount ~= 2
    error('Every input must be preceeded by an identifying string')
else
    for n = 1:length(Ident)
            if upper(varargin{Ident(n)}) == 'N'     % Equivalent Noise Spectrum Level (3.15)
                N = varargin{Ident(n)+1};
            elseif upper(varargin{Ident(n)}) == 'E' % Equivalent Speech Spectrum Level (3.6)
                E = varargin{Ident(n)+1};
            elseif upper(varargin{Ident(n)}) == 'T' % Equivalent Hearing Threshold Level [dBHL] (3.23)
                T = varargin{Ident(n)+1};
            elseif upper(varargin{Ident(n)}) == 'I' % Band Importance function (3.1)
                I = varargin{Ident(n)+1};
            else
                error('Only ''E'', ''I'', ''N'', and ''T'' are valid identifiers');
            end;
    end;
end;

if isempty(who('E')),  error('The Equivalent Speech Spectrum Level, ''E'', must be specified');     end
if isempty(who('N')),  N = -50*ones(1,18);                                                          end;
if isempty(who('T')),  T = zeros(1,18);                                                             end;
if isempty(who('I')),  I = 1;                                                                       end;


N = N(:)'; T = T(:)'; E = E(:)';
if length(N) ~= 18, error('Equivalent Noise Spectrum Level: Vector size incorrect');                end;
if length(T) ~= 18, error('Equivalent Hearing Threshold Level: Vector size incorrect');             end;
if length(E) ~= 18, error('Equivalent Speech Spectrum Level: Vector size incorrect');               end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPLEMENTATION OF SPEECH INTELLIGIBILITY INDEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =======================    THE NUMBERS IN PARENTHESIS REFER TO THE SECTIONS IN THE ANSI STANDARD  =======================  

% Define band center frequencies for 1/3rd octave procedure (Table 3)
f = [160 200 250 315 400 500 630 800 1000 1250 1600 2000, ...
     2500 3150 4000 5000 6300 8000];

% Define Internal Noise Spectrum Level (Table 3) 
X = [0.6 -1.7 -3.9 -6.1 -8.2 -9.7 -10.8 -11.9 -12.5 -13.5 -15.4 -17.7, ...
	-21.2 -24.2 -25.9 -23.6 -15.8 -7.1];

% Self-Speech Masking Spectrum (4.3.2.1 Eq. 5)
V = E - 24;

% 4.3.2.2	
B = max(V,N);
	
% Calculate slope parameter Ci (4.3.2.3 Eq. 7)
C = 0.6.*(B+10*log10(f)-6.353) - 80;

% Initialize Equivalent Masking Spectrum Level (4.3.2.4)
Z = [];
Z(1) = B(1);

% Calculate Equivalent Masking Spectrum Level (4.3.2.5 Eq. 9)
for i = 2:18
	Z(i) = 10*log10(10.^(0.1*N(i)) + ...
	sum(10.^(0.1*(B(1:(i-1))+3.32.*C(1:(i-1)).*log10(0.89*f(i)./f(1:(i-1)))))));
end;	

% Equivalent Internal Noise Spectrum Level (4.4 Eq. 10)
X = X + T;
	
% Disturbance Spectrum Level (4.5)
D = max(Z,X);
	
% Level Distortion Factor (4.6 Eq. 11)
L = 1 - (E - SpeechSptr('normal') - 10)./160;
L = min(1,L);

% 4.7.1 Eq. 12
K = (E-D+15)/30;
K = min(1,max(0,K));

% Band Audibility Function (7.7.2 Eq. 13)
A = L.*K;

% Speech Intelligibility Index (4.8 Eq. 14)
S = sum(BndImp(I).*A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           PRIVATE FUNCTIONS                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = BndImp(tst);
% Band importance functions:
% tst = 1:	Average speech as specified in Table 3
%		2:	various nonsense syllable tests where most English
%			phonemes occur equally often
%		3:	CID-22
%		4:	NU6
%		5:	Diagnostic Rhyme test
%		6:	short passages of easy reading material
%		7:	SPIN

if (nargin ~= 1),
	error('Incorrect # of input args to BndImp');
end;

if ~((tst==1)|(tst==2)|(tst==3)|(tst==4)|(tst==5)|(tst==6)|(tst==7)),
	error('Band Importance function must be integer between 1 and 7');
end;
	
BIArr= [0.0083	0		0.0365	0.0168	0		0.0114	0
		0.0095	0		0.0279	0.013	0.024	0.0153	0.0255
		0.015	0.0153	0.0405	0.0211	0.033	0.0179	0.0256
		0.0289	0.0284	0.05	0.0344	0.039	0.0558	0.036
		0.044	0.0363	0.053	0.0517	0.0571	0.0898	0.0362
		0.0578	0.0422	0.0518	0.0737	0.0691	0.0944	0.0514
		0.0653	0.0509	0.0514	0.0658	0.0781	0.0709	0.0616
		0.0711	0.0584	0.0575	0.0644	0.0751	0.066	0.077
		0.0818	0.0667	0.0717	0.0664	0.0781	0.0628	0.0718
		0.0844	0.0774	0.0873	0.0802	0.0811	0.0672	0.0718
		0.0882	0.0893	0.0902	0.0987	0.0961	0.0747	0.1075
		0.0898	0.1104	0.0938	0.1171	0.0901	0.0755	0.0921
		0.0868	0.112	0.0928	0.0932	0.0781	0.082	0.1026
		0.0844	0.0981	0.0678	0.0783	0.0691	0.0808	0.0922
		0.0771	0.0867	0.0498	0.0562	0.048	0.0483	0.0719
		0.0527	0.0728	0.0312	0.0337	0.033	0.0453	0.0461
		0.0364	0.0551	0.0215	0.0177	0.027	0.0274	0.0306
		0.0185	0		0.0253	0.0176	0.024	0.0145	0];

I = BIArr(:,tst)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = SpeechSptr(VclEfrt);
% This function returns the standard speech spectrum level from Table 3

Ei=[32.41	33.81	35.29	30.77;
    34.48	33.92	37.76	36.65;
    34.75	38.98	41.55	42.5;
    33.98	38.57	43.78	46.51;
    34.59	39.11	43.3	47.4;
    34.27	40.15	44.85	49.24;
    32.06	38.78	45.55	51.21;
    28.3	36.37	44.05	51.44;
    25.01	33.86	42.16	51.31;
    23		31.89	40.53	49.63;
    20.15	28.58	37.7	47.65;
    17.32	25.32	34.39	44.32;
    13.18	22.35	30.98	40.8;
    11.55	20.15	28.21	38.13;
    9.33	16.78	25.41	34.41;
    5.31	11.47	18.35	28.24;
    2.59	7.67	13.87	23.45;
    1.13	5.07	11.39	20.72];

switch lower(VclEfrt)
case 'normal', E = Ei(:,1)';
case 'raised', E = Ei(:,2)';
case 'loud',   E = Ei(:,3)';
case 'shout',  E = Ei(:,4)';
otherwise, error('Identifyer string to ''E'' not recognized')
end;


% EOF