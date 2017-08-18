function[source]=tdoafunction(delaymat,c,fov)
%This function inputs the matrix 'delaymat',speed of sound 'c' and field of
%view 'fov' and outputs the estimated source location using the TDOA method.

%The input matrix 'delaymat' will have the microphone pairs and their corresponding delays
%arranged along each row.The number of rows correspond to the number of tdoa estimates.
% The first column corresponds to the time delays between the paired off mics, the next with 
%the reference microphones positions followed by the non-reference mic positions of the pair.
%The size of the 'delaymat' will be given by n x (2m+1) where n is the
%number of time delay estimates and m is the dimension of the microphones.
%The fov gives the corner points of the field of view around which the microphones are placed.
%The output will be the coordinates of the estimated sound source location.

%[source]=tdoafunction(delaymat,c,fov)

 %Written by Arulkumaran Muthukumarasamy May 2007
 
delay1=delaymat(:,1);
mpos1=delaymat(:,2:5);
[ma1 mb1]=size(mpos1);%  Get the number of dimensions and mics
figure;
for i=1:ma1
    dist1 = ((mpos1(i,1)-mpos1(i,3))^2+(mpos1(i,2)-mpos1(i,4))^2)^0.5;%finding the distance between the microphone pairs
    center=[(mpos1(i,1)+mpos1(i,3))/2 (mpos1(i,2)+ mpos1(i,4))/2]';%finding the center of the position of  microphones
    %finding the direction of arrival
    w=(c*delay1(i))/dist1;
    theta1(i)=(asind(w));
    %adding an extra 90 degrees to slopes of the the microphones placed vertically
    if (mpos1(i,1)>=(fov(1,2)))|(mpos1(i,1)<=(fov(1,1)))
        angle=90;
    else
        angle=0;
    end
    %finding the slopes of the line connecting the potential source
    %locations through the center of the microphone pairs
    if ((mpos1(i,2))>=(fov(2,1))|mpos1(i,1)>=(fov(1,2)))
        slope(i)=tand(abs(90+theta1(i)+angle));    
    else  
        slope(i)=tand(abs(90-theta1(i)+angle));    
    end
       b(i)=center(2,:)-(slope(i)*center(1,:));%estimating the y-intercept
% uncomment these set of lines to plot the locus of the potential source locations for every pair 
   %figure;
   x=-10:1:10;
   y=slope(i)*x+b(i);
   hold on
   axis([-.21 .21 -.21 .21]);
   axis ij;
   plot(x,y);
   
end
 
 % finding the intersection of the lines by using the equation of the
 % line(slope-intercept form y=mx+b)where m is the slope obtained.
 D=[slope' -1*ones(ma1,1)];
 B=-1*b';
 m1=inv(D'*D);
 if (isinf(m1)==1)|(isnan(m1)==1)
     msgbox('ERROR:Matrix inverse does not exist')
 end
 m2=D'*B;
 source=m1*m2;%% returns the estimated source location
 