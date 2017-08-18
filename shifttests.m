%  This script tests the shift/delay fucntions

ns = 60;  %  Number of points in array
g = randn(1,ns)';


%g = hilbert(g);
%g = zeros(1,ns);
%g(23) = 1;
tab = subsamplefir(5,4,'c')
inc = .005;
for k=1:40
    figure(1)
    tic
    [sd, tn] = delayt(g,10,inc*k,50/10);
    toc
    plot(tn,sd,'r',[0:length(g)-1]/10,g,'b')
    figure(2)
    tic
    [sd, tn] = delayf(g,10,inc*k,50/10);
    toc
    plot(tn,sd,'r',[0:length(g)-1]/10,g,'b')
    figure(3)
    tic
    [sd, tn] = delaytab(g,tab,10,inc*k,50/10);
    toc
    plot(tn,sd,'r',[0:length(g)-1]/10,g,'b')
    
    pause(.4)
end
