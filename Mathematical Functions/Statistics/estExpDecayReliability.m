function R=estExpDecayReliability(t,T,method)
%%ESTEXPDECAYRELIABILITY Determine the probability that a particular
%       example of something will be functional after time T given
%       independent the failure times of n other examples under the
%       assumption that the time of failure given something is functional
%       at time zero is exponentially distributed, so the probability of
%       success (it lasting longer than T) given that it works at time zero
%       is Pr{T>t}=exp(-t/theta) for some parameter theta. 
%
%INPUTS: t A 1Xn or nX1 vector of the times when something failed. All
%          elements should be >=0.
%        T The desired time after which a failure is acceptable. T>=0.
%   method The estimation approach. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            expected value estimate.
%          1 Use the maximum likelihood estimate. This estimate is biased.
%
%OUTPUTS: R The probability that  the system will still be operational at
%           some point after time T.
%
%This implements the methods of [1]. The basic reliability model is
%R=p*P
%where p is the probability of a failure at time 0 and P is the probability
%of a failure between time 0 and T. THe exponential decay model applies to
%P.
%
%EXAMPLE:
%We demonstrate here that the value of R at time T is more reliable when
%computing using method 0 than method 1 given n=30 samples of the data.
% p=0.9;
% T=2;
% theta=3;
% lambda=1/theta;
% RTrue=p*(1-ExponentialD.CDF(T,lambda))
% 
% numRuns=1e4;
% R0=0;
% R1=0;
% for curRun=1:numRuns
%     n=30;
%     t=zeros(n,1);
%     for k=1:n
%         if(rand()<p)%If working at time 0.
%            t(k)=ExponentialD.rand(1,lambda);
%         else
%            t(k)=0;
%         end
%     end
%     R0=R0+estExpDecayReliability(t,T,0);
%     R1=R1+estExpDecayReliability(t,T,1);
% end
% R0=R0/numRuns
% R1=R1/numRuns
%Typically, R0 will be closet to RTrue than R1.
%
%REFERENCES:
%[1] E. L. Pugh, "The best estimate of reliability in the exponential
%    case," Operations Research, vol. 11, no. 1, pp. 57-61, Jan.-Feb. 1963.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(method))
    method=0; 
end

n=length(t);
sel=(t~=0);
p=sum(sel)/n;

%Keep only failure times that are >0.
t=t(sel);
n=length(t);

if(n==0)
   R=0;
   return;
end

theta=sum(t)/n;%Equation 7.
switch(method)
    case 0%The expected value solution.
        R=p*(1-T/(n*theta))^(n-1);%Equation 17 (with p inserted).
    case 1%The ML solution.
        R=p*exp(-T/theta);%Equation 6 (with p inserted).
    otherwise
        error('Unknown method specified.')
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
