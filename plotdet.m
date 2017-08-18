function plotdet(det, filename, nwlen)

% plots the windows corresponding to predetermined detections DET in audio file
% FILENAME. NWLEN is the window length in samples
[~,fs]= audioread(filename, [1 2]);
for i=1:size(det,1)
    bseg=floor(det(i,1)*fs);
    eseg=bseg+nwlen-1;
    [y,fs]= audioread(filename, [bseg,eseg]);
    figure(1); plot(y);
    pause(.1);
end