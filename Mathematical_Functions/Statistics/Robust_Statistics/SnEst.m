function Sn=SnEst(x,scaleType,algorithm)
%%SNEST Compute the S estimator of a real vector of points. This is a
%      robust estimation statistic, which is discussed as an alternative to
%      the median absolute deviation about the median (MAD) in [1]. As
%      described in [1] and [2], the S-estimator is
%      Sn=lomedian_{i=1:n}(himedian_{j=1:n}(abs(x(i)-x(j))));
%      Where low median is the fix((n+1)/2)-th order statistic and the high
%      median is the fix((n+1)/2)+1-th order statistic.
%
%INPUT: x A real nX1 vector of points. If an empty matrix is passed, the
%         NaN is returned.
% scaleType This parameter specifies how the S estimate should be scaled.
%         Possible values are:
%         0 (The default if omitted or an empty matrix is passed) Do not
%            scale the estimate.
%         1 Scale the estimate to be asymptotically unbiased as an
%           estimator for the standard deviation of a Gaussian
%           distribution. The estimator is biased for small n.
%         2 Scale the estimate to be approximately unbiased for the given
%           n. This uses the low-precision empirically derived
%           multiplication factors given in [2]. If n is even and n>9, then
%           this gives the same result as option 1.
% algorithm An optional parameter specifying the algorithm used. Possible
%         values are:
%         0 Use the naïve algorithm mentioned in [2].
%         1 (The default if omitted or an empty matrix is passed) Use the
%           efficient algorithm of [2].
%
%OUTPUTS: Sn The S-estimator.
%
%As derived in [1], if one assumes that the underlying distribution is
%Gaussian and one takes the limit of n->Inf, then the S estimator can be
%scaled as c*Sn to be an estimator of the standard deviation of a normal
%distribution. Specifically, the scale value is a value c such that 
%GaussianD.CDF(GaussianD.invCDF(3/4)+inv(c))-GaussianD.CDF(GaussianD.invCDF(3/4)-inv(c))==1/2
%An approximate solution is 
%c=1.192598553123208483257612563682716957731
%Other contants can be found for estimating the standard deviation of other
%estimators. For the Gaussian distribution, it was empirically demonstrated
%that the estimates produced with small values of n are biased. Thus,
%additional empirical constants to debias the estimates were found in [2].
%
%EXAMPLE:
%Here, we sample a normal distribution and then use the scaled S estimator
%to estimate the standard deviation.
% numSamples=1e4;
% sigma=20;
% mu=-8;
% x=mu+sigma*randn(numSamples,1);
% sigmaEst=SnEst(x,1)
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

if(nargin<2||isempty(scaleType))
    scaleType=0;
end

if(nargin<3||isempty(algorithm))
   algorithm=1; 
end

n=length(x);

if(n==0)%if(isempty(n))
    Sn=NaN;
    return;
end

n2=fix(n/2);
n12=fix((n+1)/2);

switch(algorithm)
    case 0%The naïve brute-force approach in [1].
        lowmed=n12;
        himed=n2+1;

        a1=zeros(n,1);
        a2=zeros(n,1);

        for i=1:n
            for j=1:n
                a1(j)=abs(x(i)-x(j));
            end

            a1=sort(a1,'ascend');
            a2(i)=a1(himed);
        end
        a2=sort(a2,'ascend');

        Sn=a2(lowmed);
    case 1%The efficient algorithm of [1].
        x=sort(x,'ascend');
        a2=zeros(n,1);

        a2(1)=x(n2+1)-x(1);
        for i=2:n12
            nA=i-1;
            nB=n-i;

            diff=nB-nA;
            leftA=1;
            leftB=1;
            rightA=nB;

            Amin=fix(diff/2)+1;
            Amax=fix(diff/2)+nA;
            while(leftA<rightA)
                span=rightA-leftA+1;
                even=1-mod(span,2);
                half=fix((span-1)/2);
                tryA=leftA+half;
                tryB=leftB+half;
                if(tryA<Amin)
                    leftA=tryA+even;
                else
                    if(tryA>Amax)
                        rightA=tryA;
                        leftB=tryB+even;
                    else
                        medA=x(i)-x(i-tryA+Amin-1);
                        medB=x(tryB+i)-x(i);
                        if(medA>=medB)
                            rightA=tryA;
                            leftB=tryB+even;
                        else
                            leftA=tryA+even;
                        end
                    end   
                end
            end

            if(leftA>Amax)
                a2(i)=x(leftB+i)-x(i);
            else
                medA=x(i)-x(i-leftA+Amin-1);
                medB=x(leftB+i)-x(i);
                a2(i)=min(medA,medB);
            end
        end

        for i=(n12+1):(n-1)
            nA=n-i;
            nB=i-1;
            diff=nB-nA;
            leftA=1;
            leftB=1;
            rightA=nB;
            Amin=fix(diff/2)+1;
            Amax=fix(diff/2)+nA;

            while(leftA<rightA)
                span=rightA-leftA+1;
                even=1-mod(span,2);
                half=fix((span-1)/2);
                tryA=leftA+half;
                tryB=leftB+half;
                if(tryA<Amin)
                    leftA=tryA+even;
                else
                    if(tryA>Amax)
                        rightA=tryA;
                        leftB=tryB+even;
                    else
                        medA=x(i+tryA-Amin+1)-x(i);
                        medB=x(i)-x(i-tryB);
                        if(medA>=medB)
                            rightA=tryA;
                            leftB=tryB+even;
                        else
                            leftA=tryA+even;
                        end
                    end
                end
            end

            if(leftA>Amax)
                a2(i)=x(i)-x(i-leftB);
            else
                medA=x(i+leftA-Amin+1)-x(i);
                medB=x(i)-x(i-leftB);
                a2(i)=min(medA,medB);
            end
        end

        a2(n)=x(n)-x(n12);

        Sn=kthOrderStat(a2,n12);
    otherwise
        error('Unknown algorithm specified.')
end

c=1.192598553123208483257612563682716957731;

switch(scaleType)
    case 0%Do not scale
    case 1%Use only the asymptotic scaling.
        Sn=c*Sn;
    case 2%Use the heuristic scaling given in [2].
        switch(n)
            case 2
                cn=0.743;
            case 3
                cn=1.851;
            case 4
                cn=0.954;
            case 5
                cn=1.351;
            case 6
                cn=0.993;
            case 7
                cn=1.198;
            case 8
                cn=1.005;
            case 9
                cn=1.131;
            otherwise
                if(mod(n,2))
                    cn=n/(n-0.9);
                else
                    cn=1;
                end
        end
        
        Sn=cn*c*Sn;
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
