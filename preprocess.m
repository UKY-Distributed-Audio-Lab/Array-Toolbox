function [ sigout, ye, sig_grad, sig_adjust ] = preprocess( sigout, winrec, fovp, fso, fsresample, sigfilt)
% This function applies pre-processing techniques to a signal before it is
% used for srp processing. Since these techniques are used in many scripts,
% it is useful to use a single function for them all.

% 1. -whiten original signal
%    -butterworth filter on resulting signal
% 2. -determine envelope
%    -butterworth filter envelope (low pass)
% 3. -trim samples from edges of signal and envelope
% 4. -compute the gradient of the envelope
% 5. -compute product of gradient and signal
%
%   Return the signal, envelope, gradient, and adjusted signal

[nwlen, nm]=size(sigout);
tapwin = flattap(nwlen,20);
%winrec = flattap(ceil(nwlen*(fsresample/fso)),20)*ones(1,size(sigout,2));
[b,aa] = butter(4,sigfilt/(fso/2));
sigout = filtfilt(b,aa,sigout);
sigout = whiten(sigout.*winrec, fovp.beta);
% ha= hilbert(sigout);
% tax= [0:length(sigout)-1]/fso;
% dmod= exp(j*2*pi*45e3*tax)';
%x=[];
% for i=1:channels
%     dmodsig= tapwin.*real(dmod.*ha(:,i));
%     x= [x, resample(dmodsig,fsresample,fso)];
% end
% x = whiten(x.*winrec, fovp.beta, [1e3 14e3]/(fso/2));
%x = whiten(x.*winrec, fovp.beta);
%ye = abs(hilbert(x));
%ye = ye - ones(length(ye),1)*mean(ye,1);
ye = abs(hilbert(sigout));
sig_grad= zeros( size(ye) );
for q=1:size(ye,2)
    sig_grad(:,q)= gradient(ye(:,q));
end
%sig_adjust= x.*sig_grad; %.*sig_grad;
sig_adjust= sigout.*sig_grad;

% sigout = whiten(sigout.*winrec, fovp.beta, [45e3 130e3]/(fso/2));
% sigout = filtfilt(b,a,sigout);
% 
% ye = abs(hilbert(sigout));
% ye = ye - ones(length(ye),1)*mean(ye,1);                %subtract avg
% ye = filtfilt(bb,aa,ye);
% sigout= sigout(50:end-50,:);
% ye= ye(50:end-50,:);
% 
% %compute gradient 1 column at a time
% sig_grad= zeros( size(ye) );
% for q=1:size(ye,2)
%     sig_grad(:,q)= gradient(ye(:,q));
% end
% sig_adjust = sig_grad.*sigout;

end

