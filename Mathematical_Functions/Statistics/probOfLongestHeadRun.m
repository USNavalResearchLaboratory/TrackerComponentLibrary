function prob=probOfLongestHeadRun(n,p,r)
%%PROBLONGESTHEADRUN In a series of n independent trials, where each trial
%       has a probability p of success, find the probability of there being
%       a consecutive sequence of at least r successes.
%
%INPUTS: n The total number of trials.
%        p The probability of a success in each trial.
%        r The minimum number of CONSECUTIVE successes out of all trials.
%
%OUTPUTS: prob The probability of there being at least r successes out of n
%              trials.
%
%If r>n/2, then the simple formula of Equation 5 in [1] is used. Otherwise
%the formula of Equation 2 in [1] is used. If r>n, then the probability is
%just 0.
%
%EXAMPLE:
%This is the example given in the paper, which is originally from de
%Moivre. The solution is 65/128.
% prob=probOfLongestHeadRun(10,1/2,3)
%
%REFERENCES:
%[1] Y. Malinovsky, "A note on the closed-form solution for the longest head
%    run problem of Abraham de Moivre," The Mathematical Intelligencer,
%    vol. 44, no. 3, pp. 267-268, Sep. 2022.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

q=1-p;
if(r>n)
    prob=0;
    return;
elseif(r<0)
    error('r is invalid.')
end

if(r>n/2)
    %Equation 5.
    prob=p^r+(n-r)*p^r*q;
else
    pr=p^r;
    qpr=q*pr;
    betaValn=1;
    betaValnmr=1;
    for l=1:floor(n/(r+1))
        betaValn=betaValn+(-1)^l*binomial(n-l*r,l)*qpr^l;
        betaValnmr=betaValnmr+(-1)^l*binomial(n-r-l*r,l)*qpr^l;
    end
    %Equation 2
    prob=1-betaValn+pr*betaValnmr;
end
end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
