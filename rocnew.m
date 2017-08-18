function [pfa, pd, t, thresh, dist, roca,maxpd,minpfa] = rocnew(c1,c2,npt)
%  This function sweeps a threshold over the range of data in 'c1' and 
%  'c2', where 'c1' is the statistic related to the phenomenon of
%  interest (target), and 'c2' represents the absence of the phenomenon
%  (noise) to compute a receiver operating characteristic (ROC) curve.
%  It is assumed that the phenomenon of interest is greater than the 
%  absence of the phenomenon. The input parameter 'npt' is the
%  number of uniform increments for the applied thresholds between
%  the maximum and minimum values of the data in 'c1' and 'c2'.
%
%  [pfa, pd, t, dist, roca] = roc(c1,c2,npt)
%
%  'pfa' and 'pd' are a vectors containing the probability of false alarm and
%  probability of detection, respectively, for a series of thresholds 
%  corresponding to values stored in vector 't'.  The ROC curve can be
%  plotted as detection probability vs. specificity as:
%  >>plot(1-pfa,pd)
%  'roca' is the area under the ROC curve and 'dist' is the minimum
%  probability of error, overall thresholds assuming equal probability for
%  each class occuring.
%
%   Updated by Kevin D. Donohue (donohue@engr.uky.edu) May 2009.

%  Compute thresholds
bt = min(min(c1), min(c2));  %  Find minimum threshold value 
et = max(max(c1), max(c2));  %  Find maximum threshold value
t = et + (bt-et)*[0:npt-1]/(npt-1);  %  Compute thresholds
nc1 = length(c1);  %  Number of target samples
nc2 = length(c2);  %  Numver of noise samples
%  Loop applies sequecent of thresholds and computer corresponding
%  detection and false-alarm probabilities
for k=1:npt
   pd(k) = length(find(c1 > t(k)))/nc1;  % Probability of detection
   pfa(k) = length(find(c2 > t(k)))/nc2;  %  Probability of false alarm
end
% Compute ROC area
roca = sum(diff(pfa).*(pd(1:npt-1)+pd(2:npt)))/2;
%  Compute distances of curve to the 100%  Detection Probability
%  0% False Alarm probability
adist=(1-pd)+(pfa);
%  Fin the closest point on curve to the 100% / 0% point 
n=find(adist==min(adist));
dist=adist(n(1));
thresh=t(n(1));  %  Threshold that results in the minimum distance point
maxpd=pd(n(1));  %  Best probability of detection
minpfa=pfa(n(1));  %  Best false alarm rate
