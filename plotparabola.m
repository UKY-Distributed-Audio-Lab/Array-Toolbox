function plotparabola(delaymat, c, mpos,a,fs,xcart,ycart)
% Plot the parabolas around which a source may be located. This script uses
% theory behind Time Delay of Arrival (TDOA) to determine possible
% locations of a source based solely on the time delay between microphones.
% 
% Written by Kyle Ihli    July 2014
if nargin>3
   delaymat=makedelaymat(a,mpos,fs,c); 
end

coeffs= delaymat(:,6);
figure;
draw_circle(0,0,.2096);
hold on;
scatter(mpos(1,:),mpos(2,:));
drez= 50;
d2=(.4/drez):(.4/drez):.4;
dlyrows= size(delaymat,1);
pairs=[];

for i=1:dlyrows
    td=delaymat(i,1);
    d1=abs(td)*c+d2;
    ref= delaymat(i,2:3);
    nref= delaymat(i,4:5);
    if td<0 % then it is closer to ref mic & should be paired with d2
        for j=1:size(d1,2)
            [xout,yout]=circcirc( ref(1), ref(2), d2(1,j), nref(1), nref(2), d1(1,j));
            pairs= [pairs; xout(1),yout(1),delaymat(i,6);xout(2),yout(2),delaymat(i,6)];
        end
    else % pair non-ref mic with d2
        for j=1:size(d1,2)
            [xout,yout]=circcirc( ref(1), ref(2), d1(1,j), nref(1), nref(2), d2(1,j));
            pairs= [pairs; xout(1),yout(1),delaymat(i,6);xout(2),yout(2),delaymat(i,6)];
        end
    end
end

%now plot the parabolas in a weighted fashion
%uncomment to see full parabolas
colors={'red' 'magenta' 'yellow' 'cyan' 'blue' 'green'};
coeffrange= max(coeffs)-min(coeffs);
coeffstep= (max(coeffs)-min(coeffs))/6;
%6 different colors, so divide range into 6 segments
pairstart=1;
pairend= size(d1,2)*2;
for i=1:size(pairs,1)/pairend
%     coeffseg= ceil((coeffs(i)-min(coeffs))/coeffstep);
%     if coeffseg==0
%         coeffseg=1;
%     end
    scatter(pairs(pairstart:pairend,1),pairs(pairstart:pairend,2), 30, colors{1});
    pairstart=pairend+1;
    pairend=pairend+size(d1,2)*2;
end

thres=.007;
for i=1:size(pairs,1)-1
    if ~isnan(pairs(i,1))
        %check to see if x and y are within THRES of another xy
        flag= 0;
        for j=i+1:size(pairs,1)
            if all(abs(pairs(i,:)-pairs(j,:))<thres)
                flag=1;
                j=size(pairs,1);
            end
        end
        if ~flag
            %scatter(pairs(i,1),pairs(i,2),30,'white');
            pairs(i,:)=NaN;
        end
    end
end
%scatter(pairs(:,1),pairs(:,2), 30, 'red');
axis([-.21 .21 -.21 .21]);
scatter(xcart,ycart,70,'blue');
hold off;