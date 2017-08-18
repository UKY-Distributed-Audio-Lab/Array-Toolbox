ns = 140;
g = randn(1,ns)';
g = hilbert(g);
%g = zeros(1,ns);
%g(23) = 1;
inc = .005;
for k=1:40
    figure(1)
    tic
    [sd, tn] = delayt(g,10,inc*k);
    toc
    plot(tn,abs(sd),'r',[0:length(g)-1]/10,abs(g),'b')
    figure(2)
    tic
    [sd, tn] = delayf(g,10,inc*k);
    toc
    plot(tn,abs(sd),'r',[0:length(g)-1]/10,abs(g),'b')
    pause(.4)
end
