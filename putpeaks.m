function putpeaks(fno,posx,posy,smk)
figure(fno)
hold on
for k=1:length(posx)
    plot(posx(k),posy(k),smk)
end
hold off
