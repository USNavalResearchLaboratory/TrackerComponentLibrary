function Qn=QnEst(x,scaleType,k,algorithm)
%%QNEST Compute the Q estimator of a real vector of points. This is a
%      robust estimation statistic, which is discussed as an alternative to
%      the median absolute deviation about the median (MAD) in [1]. As
%      described in [1] and [2], the Q-estimator is the kth order statistic
%      of abs(x(i)-x(j)) with i<j.
%
%INPUT: x A real nX1 vector of points. If an empty matrix is passed, the
%         NaN is returned.
% scaleType This parameter specifies how the Q estimate should be scaled.
%         Possible values are:
%         0 (The default if omitted or an empty matrix is passed) Do not
%            scale the estimate.
%         1 Scale the estimate to be asymptotically unbiased as an
%           estimator for the standard deviation of a Gaussian
%           distribution. The estimator is biased for small n.
%         2 Scale the estimate to be approximately unbiased for the given
%           n. This uses the low-precision empirically derived
%           multiplication factors given in [2].
%       k The order statistic to take. The default if omitted or an empty
%         matrix is passed is binomial(h-1,2), which is the value where the
%         scale values that one can choose above are valid. The maximum is
%         n*(n+1)/2-n. For binomial(h-1,2)+1<=k<=binomial(h,2), when no two
%         samples coincide, the breakdown point of the estimator is
%         maximized, as discussed in [1] and [2].
% algorithm An optional parameter specifying the algorithm used. Possible
%         values are:
%         0 (The default if omitted or an empty matrix is passed) Use the
%           efficient algorithm of [2].
%         1 Use a brute-force approach, implied by the definition of the Q
%          estimator.
%
%OUTPUTS: Qn The Q-estimator.
%
%As derived in [1], if one assumes that the underlying distribution is
%Gaussian and one takes the limit of n->Inf, then the S estimator can be
%scaled as c*Sn to be an estimator of the standard deviation of a normal
%distribution. Specifically, the scale value is a value d after Equation
%3.7 in [1]:
%d=1/(sqrt(2)*GaussianD.invCDF(5/8))=2.2191444659850757931851380420763
%Note that the numeric value given in the paper appears to be slightly in
%error. Other contants can be found for estimating the standard deviation
%of other estimators. For the Gaussian distribution, it was empirically
%demonstrated that the estimates produced with small values of n are
%biased. Thus, additional empirical constants to debias the estimates were
%found in [2].
%
%EXAMPLE:
%Here, we sample a normal distribution and then use the scaled Q estimator
%to estimate the standard deviation.
% n=1e4;
% sigma=20;
% mu=-8;
% x=mu+sigma*randn(n,1);
% Qn=QnEst(x,1)
%
%REFERENCES:
%[1] P. J. Rousseeuw and C. Croux, "Alternatives to the median absolute
%    deviation," Journal of the American Statistical Association, vol. 88,
%    no. 424, pp. 1273-1283, Dec. 1993.
%[2] C. Croux and P. J. Rousseeuw, "Time-efficient algorithms for two
%    highly robust estimators of scale," in Computational Statistics.
%    Berlin: Springer-Verlag, Aug. 1992, vol. 1: Proceedings of the 10th
%    Symposium on Computational Statistics, pp. 411-428, conference
%    location: Neuchâtel, Switzerland.
%
%August 2018 David Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(x);

if(n==0)
    Qn=NaN;
    return
elseif(n<2)
    Qn=0;
    return;
end

if(nargin<4||isempty(algorithm))
    algorithm=0;
end

if(nargin<2||isempty(scaleType))
    scaleType=0;%No scaling.
end

h=fix(n/2)+1;
if(nargin<3||isempty(k))
    k=binomial(h,2);
end

switch(algorithm)
    case 0
        work=zeros(n,1);
        weight=zeros(n,1);
        P=zeros(n,1);
        Q=zeros(n,1);

        y=sort(x,'ascend');

        left=n-(1:n).'+2;
        right=n*ones(n,1);

        jHelp=fix(n*(n+1)/2);
        kNew=k+jHelp;
        nL=jHelp;
        nR=n^2;

        found=false;
        iter=0;
        while(1)
            iter=iter+1;
            if(nR-nL>n)
                j=1;
                for i=2:n
                    if(left(i)<=right(i))
                        weight(j)=right(i)-left(i)+1;
                        jHelp=left(i)+fix(weight(j)/2);
                        work(j)=y(i)-y(n+1-jHelp);
                        j=j+1;
                    end
                end
                trial=weightedMedianHi(work(1:(j-1)),weight(1:(j-1)));

                j=0;
                for i=n:-1:1
                    while((j<n)&&((y(i)-y(n-j))<trial))
                        j=j+1;
                    end

                    P(i)=j;
                end

                j=n+1;
                for i=1:n
                    while((y(i)-y(n-j+2))>trial)
                        j=j-1;
                    end
                    Q(i)=j; 
                end

                sumP=sum(P);
                sumQ=sum(Q)-n;

                if(kNew<=sumP)
                    right=P;
                    nR=sumP;
                else
                    if(kNew>sumQ)
                        left=Q;
                        nL=sumQ;
                    else
                        Qn=trial;
                        found=true;
                        break;
                    end
                end
            else
                break;
            end
        end

        if(found==false)
            j=1;
            for i=2:n
                if(left(i)<=right(i))
                    for jj=left(i):right(i)
                        work(j)=y(i)-y(n-jj+1);
                        j=j+1;
                    end
                end
            end
            Qn=kthOrderStat(work(1:(j-1)),kNew-nL);
        end
    case 1%Use brute-force distances.
        costs=abs(repmat(x(:),1,n)-repmat(x(:)',n,1));
        %We only want half of the matrix without the main diagonal.
        theVec=vech(costs,1,true);
        Qn=kthOrderStat(theVec,k);
        test=1
    otherwise
        error('Unknown algorithm specified.')
end
        
%The asymptotic scale value when using a Gaussian distribution as derived
%in [1].
c=2.2191444659850757931851380420763;

switch(scaleType)
    case 0%Do not scale
    case 1%Use only the asymptotic scaling.
        Qn=c*Qn;
    case 2%Use the heuristic scaling given in [2].
        switch(n)
            case 2
                dn=0.399;
            case 3
                dn=0.994;
            case 4
                dn=0.512;
            case 5
                dn=0.844;
            case 6
                dn=0.611;
            case 7
                dn=0.857;
            case 8
                dn=0.669;
            case 9
                dn=0.872;
            otherwise
                if(mod(n,2))
                    dn=n/(n+1.4);
                else
                    dn=n/(n+3.8);
                end
        end
        
        Qn=dn*c*Qn;
    otherwise
        error('Unknown scale type specified.')
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
