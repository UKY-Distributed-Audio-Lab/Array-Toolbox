function denoised = denoisefast(sig,thres)
% Denoises the signal by checking if the absolute value of each element is less than the
% threshold THRES. If so, it brings the value to zero. else, it either adds or
% subtracts the threshold from the value to bring it closer to zero
% 
% Inputs:
%     sig - signal to be denoised. Each column represents a mic, rows represent samples
%     thres - the desired noise threshold
% Outputs:
%     denoised - the resulting denoised signal
%
% Written by Kyle Ihli    July 2014

denoised= zeros(size(sig,1),size(sig,2));

for row=1:size(sig,1)
    for column=1:size(sig,2)
       if abs(sig(row,column))<=thres
           denoised(row,column)=0;
       elseif sig(row,column)>0
           denoised(row,column)=sig(row,column)-thres;
       else
           denoised(row,column)= sig(row,column)+thres;
       end
    end
end