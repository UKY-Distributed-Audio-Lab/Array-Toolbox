function tfa = gamthresh(pfa,kval)
%  This function computes the false-alarm threshold for a gamma distributed
%  function with k=.5. and mean 1.
%
%    tfa = gamthresh(pfa)
%
%  It iterates with the incomplete gamma until a treshold is found that
%  is accurate to 3 significant digits
%
%  Written by Kevin D. Donohue (donohue@engr.uky.edu) December 17, 2007

%  Check to make sure probabity values are in range
if pfa <=0 || pfa >=1
    disp('The False alarm probability must be between 0 and 1, Come on!')
    tfa = Nan; % If not output not a number
else
    %  If probabilities are in range set up tolerance and initals values
tol = pfa/1e5;  %  Compute tolerance for stopping rule
tfa = 1;  %  Initial threshold guess
epfa = 1 - gammainc(tfa,kval);  %  Initial error from guess
errr = epfa - pfa;
%  Set initial direction flag for incrementing threshold in initial
%  direction
if errr > 0
    dirflag = 1;   %  Increase treshold
else
    dirflag = 0;  %  Decrease treshold
end
inc = .5;  %  Start off with a threshold increment of .5
%  Loop increment threshold guess until tolerance is met
while abs(errr) > tol
    %  If positive error increase threshold value
    if errr > 0
        if dirflag == 0  %  If just changed directions ...
            inc = inc/2;   %  Cut threshold in half
            dirflag = 1;  % Set direction flag to increase
        end
        tfa = tfa + inc;  %  Increment threshold
        
    else  % if negative error decrease threshold value
        if dirflag == 1  % if just changed direction
            inc = inc/2;  % cut treshold in half
            dirflag = 0;  %  Set direction flag to decrease
        end
        tfa = tfa - inc;  %  Decrement threshold
    end
    epfa = 1 - gammainc(tfa,kval);  %  Compute new error
    errr = epfa - pfa;
end
tfa = tfa/kval;  %  Scale threshold to account for mean actually being 1/(.5)
end
       
