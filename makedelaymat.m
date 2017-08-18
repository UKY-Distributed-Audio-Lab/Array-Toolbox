function [ delaymat ] = makedelaymat(a,mpos,fs,c)

% Creates a delaymatrix DELAYMAT for use with TDOA script.
% Requires user input a matrix with microphone signal start times,
% and their coordinates.
% 
%     EX:    __________________
%            start1 | x1 | y1 |
%            start2 | x2 | y2 |
%            start3 | x3 | y3 |
%            start4 | x4 | y4 |
% number of rows is (n choose 2) where n is number of rows in inmat
%
% Written by Kyle Ihli    July 2014
numrows= nchoosek(size(mpos,2), 2);
% number of dimensions is the number of columns minus 1
dim= size(mpos,1)-1;

[nwlen, channels]=size(a);
tapwin = flattap(nwlen,20);
fsresample= 44.1e3;
winrec = flattap(nwlen*(fsresample/300e3)+1,20)*ones(1,4);
[b,aa] = butter(4,[45e3,70e3]/(fs/2));  
a = filtfilt(b,aa,a);
ha= hilbert(a);
tax= [0:length(a)-1]/fs;
dmod= exp(j*2*pi*45e3*tax)';
x=[];
for i=1:channels
    dmodsig= tapwin.*real(dmod.*ha(:,i));
    x= [x, resample(dmodsig,fsresample,300e3)];
end
x = whiten(x.*winrec, .5, [1e3 14e3]/(fs/2));
ye = abs(hilbert(x));
ye = ye - ones(length(ye),1)*mean(ye,1);
sig_grad= zeros( size(ye) );
for q=1:size(ye,2)
    sig_grad(:,q)= gradient(ye(:,q));
end
sig_adjust= x.*sig_grad;
%only use one of the following
%[y,ind]=max(sig_grad.*a);
%[y,ind]=max(a);

% use xcorr
%change following line for the signal to be processed
proc= x;
corr=cell( [1, numrows]);
lags=corr;
corcolumn=1;
mpose=mpos;
while size(proc,2) >1
    for i=1:size(proc,2)-1
        maxlags= floor(sqrt( (mpose(1,1)-mpose(1,i+1))^2+(mpose(2,1)-mpose(2,i+1))^2 )/c*fsresample);
        [corr{corcolumn} lags{corcolumn}]= xcorr( proc(:,1),proc(:,i+1),maxlags,'coeff');
        corcolumn=corcolumn+1;
    end
    proc(:,1)=[];
    mpose(:,1)=[];
end
coeff=zeros(numrows,1);
ind=coeff;
for i=1:size(corr,2);
    [coeff(i),ind(i)]= max(corr{i});
    ind(i)= lags{i}(ind(i));
    ind(i)=ind(i)/fsresample; %convert to seconds
end
tst= ind/fsresample; %vector of peak times relative to start
inmat = zeros(size(mpos,2), 3); %xy only
for i=1:size(mpos,2)
   inmat(i,:)= [tst(i), mpos(1,i), mpos(2,i)]; 
end
% initialize delaymat to numrows x (2*dim+2) matrix
delaymat= zeros(numrows, 2*dim+2);
dlyrow=1;
while size(inmat,1) >1
    for i=1:size(inmat,1)-1
        delaymat(dlyrow,:)= [ind(dlyrow), inmat(1,2), inmat(1,3), inmat(i+1,2), inmat(i+1,3),coeff(dlyrow)];
        dlyrow= dlyrow+1;
    end
    inmat(1,:)=[];
end

end

