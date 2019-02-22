function [xb,Rb,d,minIdx]=meanCovBACON(x,m,version,gammaVal,maxIter)
%%MEANCOVBACON Given a set of samples of a distribution with outliers, 
%              compute the mean and covariance of the set of points using
%              the blocked adaptive computationally effcient outlier
%              nominators (BACON) algorithm of [1]. This algorithm tries to
%              efficient detect and eliminate the outliers.
%
%INPUTS: x A pXn set of n p-dimensional samples of the distribution. The
%          algorithm of [1] requires that n>1+3*p for the approximate
%          correction term cnp to be used. If n<=1+3*p, then the
%          approximate value of cnp given in the appendix of [2] is used.
%          It is still required that n>p.
%        m Optionally, one can specify the number of points to use in the
%          initial set of points used in the algorithm prior to iterating.
%          If omitted or an empty matrix is passed, the default is
%          m=max(p+1,min(fix(n/4),4*p));
%  version Two algorithms for obtaining an initial estimate are given in
%          [1]. This parameter chooses between them. Possible values are:
%          1 (The default if omitted or an empty matrix is passed) Use
%            algorithm 2, version V1 of [1].
%          2 Use algorithm 2, version V2 of [1].
% gammaVal The threshold for determining outliers. Due to the use of a
%          heuristic correction term (cnp) in the algorithm, this value
%          will only approximately be met. This is a value of the chi-
%          squared distribution with p degrees of freedom that is the PG
%          percentile for PG close to 1. If omitted or an empty matrix
%          is passed, then the value for the 99.97% value of the
%          distribution with p degrees of freedom is used. The heuristic
%          correction is inaccurate when PG is not close to 1.
%  maxIter The maximum number of iterations of the algorithm to perform.
%          The default if omitted or an empty matrix is passed is n+p+1.
%          and should usually be more than enough.
%
%OUTPUTS: xb The pX1 mean vector after approximately omitting outliers.
%         Rb The pXp covariance matrix, exceluding approximate outliers.
%          d The Mahalanobis distances values of all of the observations
%            using xb as the mean and Rb as the covariance matrix.
%     minIdx The indices of the observations that were used to compute xb
%            and Rb.
%
%This function implements the algorithm of [1]. The gammaVal threshold is
%only approximate. The cnp term in [1] is used to approximately make it
%hold. However, it cannot be used if n<=1+3*p. Thus, for n in that range,
%the heuristic cnp term of [2] is used.
%
%EXAMPLE 1:
%Here, we have a bivariate Gaussian distribution polluted with 10% of its
%samples from another Gaussian distribution.
% xBarTrue=[4;3];
% RTrue=[12,-7;
%        -7,14];
% SRTrue=chol(RTrue,'lower');
% xPoll=[-120;6];
% RPoll=[1,   0.05;
%        0.05,0.1];
% SRPoll=chol(RPoll,'lower');
% 
% numTrue=10000;
% numPoll=1000;
% N=numTrue+numPoll;
% x=zeros(2,numTrue+numPoll);
% for k=1:numTrue
%     x(:,k)=xBarTrue+SRTrue*randn(2,1);
% end
% for k=(numTrue+1):N
%     x(:,k)=xPoll+SRPoll*randn(2,1);
% end
% [xb,Rb]=meanCovBACON(x)
% [xAll, RAll]=calcMixtureMoments(x)
%One will see that xb and Rb have a number of digits in common with
%xBarTrue and RTrue, whereas RAll is massively larger than RTrue and xTrue
%is farther from xBarTrue.
%
%EXAMPLE 2:
%Here, we demonstrate the accuracy of the heuristic correction that makes
%the gammaVal term hold. In this instance, the distribution in question is
%Gaussian and is NOT pullted by anything else. Thus, if out initial PG
%value from which gammaVal is derived is accurately modeled in
%meanCovBACON, then the observed ratio of points omitted should be the
%same.
% N=100;
% xBar=[4;3];
% R=[12,-7;
%        -7,14];
% SR=chol(R,'lower');
% numRuns=1000;
% 
% PG=0.9997;
% gammaVal=ChiSquareD.invCDF(PG,2);
% x=zeros(2,N);%Allocate space.
% PGObs=0;
% for curRun=1:numRuns
%     for k=1:N
%         x(:,k)=xBar+SR*randn(2,1);
%     end
%     [~,~,~,minIdx]=meanCovBACON(x,[],[],gammaVal);
% 
%     PGObs=PGObs+length(minIdx);
% end
% PGObs=PGObs/(N*numRuns)
%One will see that PGObs is close to PG. If one were to reduce PG to 0.8
%one would see that PGObs is much smaller that it should be (around 0.3).
%On the other hand, if N=6, and PG unchanged, one would see a lesser
%underestimate with PGObs around 0.94 rather than 0.9997. The approximation
%for gammaVal is more accurate for higher numbers of samples and PG values
%closer to 1.
%
%REFERENCES:
%[1] N. Billor, A. S. Hadi, and P. F. Velleman, "BACON: Blocked adaptive
%    computationally efficient outlier nominators," Computational
%    Statistics and Data Analysis, vol. 34, no. 3, pp. 279-298, 28 Sep.
%    2000.
%[2] P. J. Rousseeuw and B. C. van Zomeren, "Unmasking outliers and
%    leverage points," Journal of the American Statistical Association,
%    vol. 85, no. 411, pp. 633-639, Sep. 1990.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

p=size(x,1);
n=size(x,2);

if(nargin<5||isempty(maxIter))
    maxIter=n+p+1;
end

if(n<2)
    error('There must be at least two points in x (n>=2).') 
end

if(n<=p)
    error('The number of columns of x must be > the number of rows.')
end

if(nargin<3||isempty(version))
    version=1;
end

if(nargin<4||isempty(gammaVal))
    PG=0.9997;
    gammaVal=ChiSquareD.invCDF(PG,p);
end

if(nargin<2||isempty(m))
    m=max(p+1,min(fix(n/4),4*p));
end

if(m>n)
    error('m must be <=n.');
end

chiVal=sqrt(gammaVal);

if(n<=1+3*p)
    %Use the value from Appendix 2 of [2].
    cnp=1+15/(n-p);
else
    %Equation 4 of [1].
    cnp=1+(p+1)/(n-p)+2/(n-1-3*p);
end

%Defined after Equation 4 of [1].
h=fix((n+p+1)/2);

%To hold the cost values (distances).
d=zeros(n,1);
switch(version)
    case 1%Algorithm 2, version 1
        [xBar, S]=calcMixtureMoments(x);
        SInv=inv(S);
        for i=1:n
            diff=x(:,i)-xBar;
            %Mahalanobis distances
            d(i)=sqrt(diff'*SInv*diff);
        end

        %Find the m smallest Mahalanobis distances. These are the potential
        %basic subset.
        [~,minIdx]=mink(d,m);
    case 2%Algorithm 2, version 2
        %Coordinate-wise median.
        mVec=median(x,2);
        
        for i=1:n
            d(i)=norm(x(:,i)-mVec);
        end
        %Find the m smallest distances.
        [~,minIdx]=mink(d,m);
    otherwise
        error('Unknown initial basic subset selection technique chosen.')
end

numMin=m;
numMinOld=0;

curIter=0;
while(numMin~=numMinOld&&curIter<maxIter)
    numMinOld=numMin;
    
    %Algorithm 3, step 2.
    [xb, Rb]=calcMixtureMoments(x(:,minIdx));
    RbInv=inv(Rb);
    for i=1:n
        diff=x(:,i)-xb;
        d(i)=sqrt(diff'*RbInv*diff);  
    end
    %Algorithm 3, Step 3.
    %Defined before Equation 4.
    chr=max(0,(h-numMin)/(h+numMin));
    cnpr=cnp+chr;

    thresh=chiVal*cnpr;

    minIdx=find(d<thresh);
    numMin=length(minIdx);
    curIter=curIter+1;
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
