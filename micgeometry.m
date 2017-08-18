%  This is the set of the experiments done with real targets in the room.

micnum=8;  %  Number of microphones
fmic = [-3.62/2+0.225 -3.62/2+0.225 1.57; 3.62/2-0.225 3.62/2-0.225 2.22]';  % Opposite Corner Points of room
%  All vertcies im image plane
v = [fmic(1:2,1), [fmic(1,1); fmic(2,2)], fmic(1:2,2), [fmic(1,2); fmic(2,1)]];  
v = [v; ones(1,4)*1.57];
%  Generate wall geometry
fwall = [-3.62/2 -3.62/2 0; 3.62/2 3.62/2 2.22]';  % Opposite Corner Points of room
%  All vertcies im image plane
vw = [fwall(1:2,1), [fwall(1,1); fwall(2,2)], fwall(1:2,2), [fwall(1,2); fwall(2,1)]];  
vw = [vw; ones(1,4)*1.57];

fov = [-1.5 -1.5 1.57; 1.5 1.5 1.57]';
%  Compute spacing for equal spacing of perimeter array on wall
spm = (norm(v(:,2)-v(:,1),2)+norm(v(:,3)-v(:,2),2))/(micnum/2-2 +4/sqrt(2));
stp = (spm/sqrt(2));
mposdum = regmicsperim(v(:,[1,3]),spm,stp);
for msk=1:micnum
    mposperim(:,msk) = mposdum(:,mod(msk+micnum-2,micnum)+1);
end
