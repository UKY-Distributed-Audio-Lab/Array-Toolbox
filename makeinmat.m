function [ inmat ] = makeinmat( a, mpos, fs )
% makes inmat for use in script makedelaymat. Reads in audio sample A and
% determines the max values for each channel and when they occur. Compiles
% these values into inmat along with microphone positions MPOS for each channel.

%filter and whiten
nwlen=size(a,1);
tapwin = flattap(nwlen,20);
winrec = tapwin*ones(1,4);
[b,aa] = butter(4,[40e3,130e3]/(fs/2));  
a = filtfilt(b,aa,a);
a = whiten(a.*winrec, .7, [45e3 60e3]/(fs/2));

ye = abs(hilbert(a));
ye = ye - ones(length(ye),1)*mean(ye,1);

sig_grad= zeros( size(ye) );
for q=1:size(ye,2)
    sig_grad(:,q)= gradient(ye(:,q));
end
%only use one of the following
%[y,ind]=max(sig_grad.*a);
[y,ind]=max(a);

% use xcorr
corr=[];


[c12,lags]= xcorr(a(:,1),a(:,2),'coef');
[c13,lags]= xcorr(a(:,1),a(:,3),'coef');
[c14,lags]= xcorr(a(:,1),a(:,4),'coef');
[c23,lags]= xcorr(a(:,2),a(:,3),'coef');
[c24,lags]= xcorr(a(:,2),a(:,4),'coef');
[c34,lags]= xcorr(a(:,3),a(:,4),'coef');

tst= ind/fs; %vector of peak times relative to start
inmat = zeros(size(mpos,2), 3); %xy only
for i=1:size(mpos,2)
   inmat(i,:)= [tst(i), mpos(1,i), mpos(2,i)]; 
end
end

