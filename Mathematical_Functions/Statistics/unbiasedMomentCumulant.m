function [val,retData]=unbiasedMomentCumulant(x,r,isMoment,algorithm,retData)
%%UNBIASEDMOMENTCUMULANT Given a set of samples of a scalar distribution,
%           return an unbiased estimate of the rth order central moment (h
%           statistic), which is E{(x-mu)^r}, where E is the expected value
%           operation, when given only samples of x and not knowing the
%           mean mu. Alternatively, one could compute an unbiased estimate
%           of the rth order cumulant (k statistic). Note that for r=0, if
%           a moment is desired, then the mean is returned, rather than
%           returning 0. If r is a vector, one can estimate generalized h
%           statistics, which are an estimate of
%           E{(x-mu)^r(1)}*E{(x-mu)^r(2)}*...*E{(x-mu)^r(v)} However,
%           algorithm 2 is used for generalized h statistics and its
%           complexity increases rapidly with the number of samples.
%
%INPUTS: x A 1XN or NX1 array of scalar samples of the distribution.
%        r The scalar order of the cumulant or of the central moment.
%          Alternatively, if this is a vector, then one will get a
%          generalized h statistic, and in that instance, algorithm 2 is
%          required.
% isMoment A parameter selecting whether a moment or a cumulant is desired.
%          The default if omitted or an empty matrix is passed is true,
%          indicating that a central moment is desired.
% algorithm An optional parameter selecting which algorithm to use.
%          possible values are:
%          0 Use an explicit solution with precomputed coefficients. These
%            are only available for scalar r=1,...,5. The explicit
%            solutions use the parameters in Table I of [1]. This is the
%            default if r<=5 is scalar.
%          1 Use Dwyer's method in [1] to determine the coefficients of the
%            expression of the moment in terms of the sample power sums
%            sum(x^k) for k=1 to r.
%          2 Use the symmetric function based algorithm of [2]. This
%            algorithm can only estimate central moments, not cumulants. 
%            The complexity of this algorithm increases rapidly with N.
%            This can be used with a vector r and is thus the default when
%            r is a vector.
% retData If algorithm=1 is used, then the retData value is returned. If
%         one calls this function passing back retData and keeps all other
%         inputs the same except replacing x with a different set of
%         samples, then the function will execute faster. One should not
%         pass back retData if any inputs other than x have changed.
%
%OUTPUTS: val The scalar value of the moment.
%     retData If algorithm=1, then this is a structure that can be passed
%     back so that if the function is called a second time with the same
%     inputs, changing only x to different samples, then the algorithm will
%     be faster.
%
%When mu is not available and must be found from the samples, as is the
%case here, then the average mean((x-mean(x)).^p) to approximate
%E{(x-mu)^p} is a biased estimator. Hence, h statistics are more
%complicated.
%
%EXAMPLE 1:
%Here, we estimate the sixth central moment of the scalar Laplace
%distribution from 1000 samples. Since the distribution is symmetric and
%unbounded, we can assume that if the estimator is unbiased, the averaging
%a large number of these estimates will approach the true moment value.
%That is done here and compared to the true moment value. One can also see
%how retData is used to speed up the algorithm.
% n=1000;
% lambda=0.1;
% mu=2;
% Gamma=1;
% momentNum=6;
% numRuns=1e4;
% 
% hAvg=0;
% retData=[];
% for curRun=1:numRuns
%     x=LaplaceD.rand(n,lambda,mu,Gamma);
%     [val,retData]=unbiasedMomentCumulant(x,momentNum,true,[],retData);
%     hAvg=hAvg+val;
% end
% avgMomentEst=hAvg/numRuns
% muList=zeros(momentNum,1);
% for k=1:momentNum
%     muList(k)=LaplaceD.momentGenFun(lambda,mu,Gamma,k);
% end
% centralMoment=raw2CentralMoments(muList);
% trueMoment=centralMoment(end)
%One should see that the average moment and the true moment should agree to
%a few decimal places.
%
%EXAMPLE 2:
%Now, we estimate the product of the second central moment and the fourth
%central moment. This necessitates the use of algorith 2, which can be
%rather slow. This example will typically take a few seconds to run.
% rng(0);%To make the results reproducible
% n=10;
% lambda=0.25;
% mu=2;
% Gamma=1;
% p=[2,4];
% numRuns=1000;
% 
% hAvg=0;
% for curRun=1:numRuns
%     x=LaplaceD.rand(n,lambda,mu,Gamma);
%     hAvg=hAvg+unbiasedMomentCumulant(x,p);
% end
% avgMomentEst=hAvg/numRuns
% muList=zeros(max(p),1);
% for k=1:max(p)
%     muList(k)=LaplaceD.momentGenFun(lambda,mu,Gamma,k);
% end
% centralMoment=raw2CentralMoments(muList);
% trueMoment=prod(centralMoment(p))
%The average moment estimate will be about 0.0998 and the true moment will
%be about 0.0938.
%
%REFERENCES:
%[1] P. S. Dwyer, "Moments of any rational integral isobaric sample moment
%    function," The Annals of Mathematical Statistics, vol. 8, no. 1, pp.
%    21-65, Mar. 1937.
%[2] D. S. Tracy and B. C. Gupta, "Generalized h-statistics and other
%    symmetric functions," The Annals of Statistics, vol. 2, no. 4, pp. 
%    837-844, 1974.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<3||isempty(isMoment))
    isMoment=true;
end

if(nargin<4||isempty(algorithm))
    if(isscalar(r))
        if(r<=5)
            algorithm=0;%Use the explicit solution.
        else
            algorithm=1;%Use the power sum solution of Dwyer.
        end
    else
        algorithm=2;
        
        if(isMoment==false)
            error('Vector values of r are not allowed for algorithm~=2.')
        end
    end
end

if(nargin<5)
    retData=[];
end

if(isMoment)%h-statistics
    switch(algorithm)
        case 0
            val=hStatExplicit(x,r);
            retData=[];
        case 1
            [val,retData]=momCumDwyer(x,r,false,retData);
        case 2
            val=genHStat(x,r);
            
            retData=[];
        otherwise
            error('Unknown Algorithm specified.')
    end
else%k-statistics
    switch(algorithm)
        case 0
            val=kStatExplicit(x,r);
            retData=[];
        case 1
            [val,retData]=momCumDwyer(x,r,true,retData);
        otherwise
            error('Unknown Algorithm specified.')
    end
end
end

function [momentVal,retData]=momCumDwyer(x,r1,isCumulant,retData)
%%MOMCUMDWYER Compute unbiased central moments (h statistics) or cumulants
%             (k statistics) using Dwyer's algorithm, described in [1]. The
%             definition of the symmetric powers in the algorithm can be
%             bettwer understood when considering [2].
%
%The algorithm uses a weird definition of a product, which involves
%appending the indices of various "b" variables on each other. Thus, an
%initial phase of the algorithm involves just determinign the indices of
%the necessary coefficients. Sections 17 and 18 of [1] give the b variables
%for k statistics and h statistics. Section 16 discusses how to turn the b
%variables into a variables. The expression for turning the a variables
%into moments, F_r, is in Section 1.
%
%REFERENCES:
%[1] P. S. Dwyer, "Moments of any rational integral isobaric sample moment
%    function," The Annals of Mathematical Statistics, vol. 8, no. 1, pp.
%    21-65, Mar. 1937.
%[2] P. S. Dwyer, "Combined expansions of products of symmetric power sums
%    and of sums of symmetric power products with application to sampling,"
%    The Annals of Mathematical Statistics, vol. 9, no. 1, pp. 1-47, Mar.
%    1938.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=length(x);

numParts=numberOfPartitions(r1,true);

if(nargin<4||isempty(retData))
    %If we have to compute the coefficients.
    %AParts holds the indices of the b that go into the sum to form the a as
    %well as their coefficients.
    aParts=cell(r1,2);
    for r=1:r1
        %Allocate space for the indices of the b values needed and their
        %coefficients.
        bIdx=cell(numParts(r+1),1);
        bCoeff=zeros(numParts(r+1),1);
        [thePartition,theData]=getNextPartition(r);
        curB=1;
        while(~isempty(thePartition))
            s=theData.d;
            p=theData.r(1:s);
            piVals=theData.m(1:s);
            rho=sum(piVals);

            [p,idx]=sort(p,'descend');
            piVals=piVals(idx);

            bCoeff(curB)=partComboNum(r,p,piVals)*(-1)^(rho-1)*factorial(rho-1);
            idx=zeros(rho,1);
            startIdx=1;
            for k=1:s
                numRep=piVals(k);
                sel=startIdx:(startIdx+numRep-1);
                idx(sel)=p(k);
                startIdx=startIdx+numRep;
            end
            bIdx{curB}=idx;

            [thePartition,theData]=getNextPartition(r,theData);
            curB=curB+1;
        end
        aParts{r,1}=bIdx;
        aParts{r,2}=bCoeff;
    end

    %Now, determine the b indices needed to construct the sum in terms of the
    %as. (The Fr equation).
    aPartNeeded=cell(r1,3);
    %This stores the a indices, the coefficient of the a term and the powers of
    %the P sums for each a term.
    [thePartition,theData]=getNextPartition(r1);
    curA=1;
    while(~isempty(thePartition))
        s=theData.d;
        p=theData.r(1:s);
        piVals=theData.m(1:s);
        rho=sum(piVals);

        [p,idx]=sort(p,'descend');
        piVals=piVals(idx);
        coeffVal=partComboNum(r1,p,piVals);

        idx=zeros(rho,1);
        startIdx=1;
        for k=1:s
            numRep=piVals(k);
            sel=startIdx:(startIdx+numRep-1);
            idx(sel)=p(k);
            startIdx=startIdx+numRep;
        end

        aPartNeeded{curA,1}=idx;
        aPartNeeded{curA,2}=coeffVal;
        aPartNeeded{curA,3}=piVals(:);
        aPartNeeded{curA,4}=p(:);

        [thePartition,theData]=getNextPartition(r1,theData);
        curA=curA+1;
    end

    %To build each of the required "a" terms with the parts needed from the
    %single index a terms.
    numCoeffs=curA-1;
    %Allocate space for the a coefficients. These must be built from the suffix
    %multiplication rule of Section 16.
    aCoeffs=zeros(numCoeffs,1);
    for curPart=1:numCoeffs
        idx=aPartNeeded{curPart,1};
        r=length(idx);

        bIdxList=aParts{idx(1),1};
        bCoeffList=aParts{idx(1),2};
        for k=2:r
            bIdxCur=aParts{idx(k),1};
            bCoeffCur=aParts{idx(k),2};
            numB=length(bIdxList);
            numCur=length(bIdxCur);

            bIdxListNew=cell(numB*numCur,1);
            bCoeffListNew=zeros(numB*numCur,1);
            newCur=1;
            for i1=1:numB
                for i2=1:numCur
                    bIdxListNew{newCur}=[bIdxList{i1};bIdxCur{i2}];
                    bCoeffListNew(newCur)=bCoeffList(i1)*bCoeffCur(i2);

                    newCur=newCur+1;
                end
            end
            bIdxList=bIdxListNew;
            bCoeffList=bCoeffListNew;
        end

        %We now have the b indices and coefficients to construct the a
        %coefficients once we compute the b values.
        numB=length(bIdxList);
        aSum=0;
        for curB=1:numB
            bIdx=bIdxList{curB};
            bCoeff=bCoeffList(curB);

            numBIdx=length(bIdx);
            s=numBIdx;

            bIdx=sort(bIdx,'descend');
            p=zeros(s,1);
            piVals=zeros(s,1);

            p(1)=bIdx(1);
            numUnique=1;
            piVals(1)=1;
            for curIdx=2:numBIdx
                if(bIdx(curIdx)==bIdx(curIdx-1))
                    s=s-1;
                    piVals(numUnique)=piVals(numUnique)+1;
                else
                    numUnique=numUnique+1;
                    p(numUnique)=bIdx(curIdx);
                    piVals(numUnique)=1;
                end
            end
            %Note that s==numUnique.
            piVals=piVals(1:s);
            p=p(1:s);
            rho=sum(piVals);

            if(isCumulant)%For Fisher k statistics
                b=(-1)^(rho-1)*factorial(rho-1)/fallingFactorial(n,rho);
            else%For h statistics
                if(s==1&&piVals(1)==1&&p(1)==r1)
                    AVal=1;
                elseif(p(1)>1&&piVals(1)==1&&s==2&&p(2)==1)
                    AVal=(-1)^piVals(2);
                elseif(p(1)==1&&s==1&&piVals(1)==r1)
                    AVal=(-1)^(r1-1)*(r1-1);
                else
                    AVal=0;
                end
                b=AVal/fallingFactorial(n,rho);
            end

            aSum=aSum+bCoeff*b;
        end

        aCoeffs(curPart)=aSum;
    end

    retData.aPartNeeded=aPartNeeded;
    retData.aCoeffs=aCoeffs;
else
    aPartNeeded=retData.aPartNeeded;
    aCoeffs=retData.aCoeffs;
    numCoeffs=length(aCoeffs);
end

%We have all of the coefficients. To get this moment with a generic length
%n set of data, we just need the power sums up to degree r.
P=zeros(r1,1);
for k=1:r1
    P(k)=sum(x.^k);
end

momentVal=0;
for curA=1:numCoeffs
    aCur=aCoeffs(curA);
    coeff=aPartNeeded{curA,2};
    piVals=aPartNeeded{curA,3};
    idx=aPartNeeded{curA,4};

    momentVal=momentVal+coeff*aCur*prod(P(idx).^piVals);
end

end

function val=partComboNum(r,p,piVals)
%%PARTCOMBONUM The expression for this multinomial value is given on page
%              23 of [1].
%
%REFERENCES:
%[1] P. S. Dwyer, "Moments of any rational integral isobaric sample moment
%    function," The Annals of Mathematical Statistics, vol. 8, no. 1, pp.
%    21-65, Mar. 1937.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %For large r, use this form to help reduce overflow.
    val=exp(gammaln(r+1)-sum(piVals.*gammaln(p+1))-sum(gammaln(piVals+1)));
end

function val=hStatExplicit(x,r)
%%HSTATEXPLICIT This implements the explicit h statistics of orders 1 to 5
%               that are given in the table I of [1] using the formula for
%               Fr given in Section 1 of [1]. h statistics are unbiased
%               estimates of the moment about the mean.
%
%INPUTS: x A 1XN or NX1 array of scalar samples of the distribution.
%        r The order of the moment. 1<=r<=5.
%
%OUTPUT: val The rth order h statistic. This is a scalar value.
%
%REFERENCES:
%[1] P. S. Dwyer, "Moments of any rational integral isobaric sample moment
%    function," The Annals of Mathematical Statistics, vol. 8, no. 1, pp.
%    21-65, Mar. 1937.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=length(x);
switch(r)
    case 1
        a1=1/n;
        
        P1=sum(x);
        val=a1*P1;
    case 2
        a2=1/(n-1);
        a11=-1/(n*(n-1));
        
        P1=sum(x);
        P2=sum(x.^2);
        
        val=a2*P2+a11*P1^2;
    case 3
        a3=n/((n-1)*(n-2));
        a21=-1/((n-1)*(n-2));
        c21=3;
        a111=2/(n*(n-1)*(n-2));
        
        P1=sum(x);
        P2=sum(x.^2);
        P3=sum(x.^3);
        
        val=a3*P3+c21*a21*P1*P2+a111*P1^3;
    case 4
        n13=(n-1)*(n-2)*(n-3);
        n4=n*n13;
        
        a4=(n^2-2*n+3)/n13;
        a31=-(n^2-2*n+3)/n4;
        c31=4;
        a22=-(2*n-3)/n4;
        c22=3;
        a211=1/n13;
        c211=6;
        a1111=-3/n4;
        
        P1=sum(x);
        P2=sum(x.^2);
        P3=sum(x.^3);
        P4=sum(x.^4);
        
        val=a4*P4+c31*a31*P3*P1+c22*a22*P2^2+c211*a211*P2*P1^2+a1111*P1^4;
    case 5
        n14=(n-1)*(n-2)*(n-3)*(n-4);
        n5=n*n14;
        
        a5=n*(n^2-5*n+10)/n14;
        a41=-(n^2-5*n+10)/n14;
        c41=5;
        a32=-(n-2)/n14;
        c32=10;
        a311=(n^2-4*n+8)/n5;
        c311=10;
        a221=(2*n-4)/n5;
        c221=15;
        a2111=-1/n14;
        c2111=10;
        a11111=4/n5;
        
        P1=sum(x);
        P2=sum(x.^2);
        P3=sum(x.^3);
        P4=sum(x.^4);
        P5=sum(x.^5);
        val=P5*a5+c41*a41*P4*P1+c32*a32*P3*P2+c311*a311*P3*P1^2+c221*a221*P2^2*P1+c2111*a2111*P2*P1^3+a11111*P1^5;         
    otherwise
        error('The explicit solutions are only available up to order 5.');
end
end

function val=kStatExplicit(x,r)
%%KSTATEXPLICIT This implements the explicit Fisher k statistics of orders
%               1 to 5 that are given in the table I of [1] using the
%               formula for Fr given in Section 1 of [1]. h statistics are
%               unbiased estimates of the cumulants of the distribution. Up
%               to order 3, they are the same as the h statistics.
%
%INPUTS: x A 1XN or NX1 array of scalar samples of the distribution.
%        r The order of the moment. 1<=r<=5.
%
%OUTPUT: val The rth order h statistic. This is a scalar value.
%
%REFERENCES:
%[1] P. S. Dwyer, "Moments of any rational integral isobaric sample moment
%    function," The Annals of Mathematical Statistics, vol. 8, no. 1, pp.
%    21-65, Mar. 1937.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=length(x);
switch(r)
    case 1
        a1=1/n;
        
        P1=sum(x);
        val=a1*P1;
    case 2
        a2=1/(n-1);
        a11=-1/(n*(n-1));
        
        P1=sum(x);
        P2=sum(x.^2);
        
        val=a2*P2+a11*P1^2;
    case 3
        a3=n/((n-1)*(n-2));
        a21=-1/((n-1)*(n-2));
        c21=3;
        a111=2/(n*(n-1)*(n-2));
        
        P1=sum(x);
        P2=sum(x.^2);
        P3=sum(x.^3);
        
        val=a3*P3+c21*a21*P1*P2+a111*P1^3;
    case 4
        n22=(n-2)*(n-3);
        n13=(n-1)*n22;
        n4=n*n13;
        
        a4=n*(n+1)/n13;
        a31=-(n+1)/n13;
        c31=4;
        a22=-1/n22;
        c22=3;
        a211=2/n13;
        c211=6;
        a1111=-6/n4;
        
        P1=sum(x);
        P2=sum(x.^2);
        P3=sum(x.^3);
        P4=sum(x.^4);
        
        val=a4*P4+c31*a31*P3*P1+c22*a22*P2^2+c211*a211*P2*P1^2+a1111*P1^4;
    case 5
        n14=(n-1)*(n-2)*(n-3)*(n-4);
        n5=n*n14;
        
        a5=n^2*(n+5)/n14;
        a41=-n*(n+5)/n14;
        c41=5;
        a32=-n*(n-1)/n14;
        c32=10;
        a311=2*(n+2)/n14;
        c311=10;
        a221=2*(n-1)/n14;
        c221=15;
        a2111=-6/n14;
        c2111=10;
        a11111=24/n5;

        P1=sum(x);
        P2=sum(x.^2);
        P3=sum(x.^3);
        P4=sum(x.^4);
        P5=sum(x.^5);
        val=P5*a5+c41*a41*P4*P1+c32*a32*P3*P2+c311*a311*P3*P1^2+c221*a221*P2^2*P1+c2111*a2111*P2*P1^3+a11111*P1^5;         
    otherwise
        error('The explicit solutions are only available up to order 5.');
end
end

function h=genHStat(x,p)
%%GENHSTAT Given a set of scalar samples, obtain an unbiased estimate of
%          the pth central sample moment. That is, estimate E{(x-mu)^p},
%          where E is the expected value operation and mu is the mean of
%          the distribution. Alternatively, p can be a length-v vector, and
%          an unbiased estimate of
%          E{(x-mu)^p(1)}*E{(x-mu)^p(2)}*...*E{(x-mu)^p(v)} is obtained.
%          When mu is not available and must be found from the samples, as
%          is the case here, then the average  mean((x-mean(x)).^p) is a
%          biased estimator. Additionally multiplying unbiased estimators
%          can lead to a biased estimator. The unbiased estimator for this
%          type of problem when p is a scalar is called an h-statistic and
%          when p is a vector is a generalized statistic.
%
%INPUTS: x An nX1 or 1Xn vector of real samples. Note that this function
%          can be slow for values of p>4 when a large number of samples is
%          used. 
%        p The integer moment desired. p>=1. p can be a vector of all
%          unique elements if the product of central moments should be
%          estimated.
%
%OUTPUTS: h The unbiased estimate of the pth central moment given the
%           samples.
%
%The estimator of [1], which is an h-statistic, is rather complicated. An
%explicit expression for h-statistics in terms of symmetric means is given
%in [2], and an expression for generalized h statististics in terms of
%symmetric means is also provided. These are how this function is
%implemented. Equation 2.2 is used when p is a scalar (a central moment)
%and Equation 3.3 is used when p is a vector (a product of central
%moments). Note that the evaluation of symmetric means involves the
%evaluation of matrix permanents and thus the complexity scales
%exponentialy with the number of samples in x. Only small numbers of
%samples can be used, because the computational complexity increases too
%rapidly.
%
%REFERENCES:
%[1] P. R. Halmos, "The theory of unbiased estimation," The Annals of
%    Mathematical Statistics, vol. 17, no. 1, pp. 34-43, Mar. 1946.
%[2] D. S. Tracy and B. C. Gupta, "Generalized h-statistics and other
%    symmetric functions," The Annals of Statistics, vol. 2, no. 4, pp. 
%    837-844, 1974.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(x);

if(any(p>=n))
    error('p must be >=n for an unbiased estimator to exist.')
end

if(isscalar(p))
    if(p==1)%Return the mean instead of 0.
        h=mean(x);
        return
    end
    
    h=0;
    for r=0:p
        h=h+(-1)^r*binomial(p,r)*symmetricMean(x,[p-r,ones(1,r)]);
    end
else
    v=length(p);
    
    p=p(:)';
    h=0;
    
    r=zeros(1,v);
    while(~isempty(r))
        sumR=sum(r);
        
        prodVal=1;
        for k=1:v
            prodVal=prodVal*binomial(p(k),r(k)); 
        end

        %We want p-r without the zero elements.
        pV=p-r;
        pV(pV==0)=[];

        h=h+(-1)^sumR*prodVal*symmetricMean(x,[pV,ones(1,sumR)]);

        r=getNextTuple(r,p);
    end
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
