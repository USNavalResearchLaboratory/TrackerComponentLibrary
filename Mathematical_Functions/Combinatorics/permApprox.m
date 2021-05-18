function [permApprox,bound]=permApprox(Z,algorithm,numIter)
%%PERMAPPROX Approximate the permanent of a real matrix using one of a
%            number of algorithms. Target-measurement association
%            probabilities in the JPDAF can be computed using matrix
%            permanents.
%
%INPUTS: Z An mXn real matrix. If m<=n, then the standard matrix permanent
%          bound is found. If m>n, then the permanent bound of A' is found
%          to be consistent with the permanents of A and A' being equal in
%          square matrices. Empty matrices have a permanent of one by
%          definition.
% algorithm An optional parameter specifying the approximation algorithm to
%          use. Possible values are:
%          0 (The default if omitted or an empty matrix is passed and Z is
%            positive) Use the crude approximation using the product of the
%            arithmetic mean of the columns of Z (or Z' if Z has has more
%            columns than rows). This algorithm assumes that Z has all
%            positive elements. This algorithm returns an error bound as a
%            second output. In some special instances, permApprox+bound is
%            a better approximation than permApprox.
%          1 This is the same as 0, but the bound computed is done assuming
%            that Z is a binary matrix.
%          2 (The default if Z has negative elements) Use the stochastic
%            approximation algorithm of [2]. Due to the use of random
%            numbers, one will not get the same result with consecutive
%            runs of the algorithm.
%  numIter The number of iterations of the stochastic permanent
%          approximation algorithm to use. If this parameter is omitted or
%          an empty matrix is passed, then the default value of n*n*m is
%          used.
%
%OUTPUTS: permApprox The approximate of the matrix permanent.
%              bound The bound on the accuracy of the matrix permanent.
%                    This is such that abs(perm(Z)-permApprox)<=bound. This
%                    is only returned for algorithm=0,1 or is n=m=1.
%                    Otherwise, an empty matrix is returned.
%
%EXAMPLE:
%A special case with algorithms 0 and 1 is that permApprox+bound is the
%exact solution when Z is a matrix of all the same value.
% Z=5*ones(7,10);
% [val1,bound]=permApprox(Z,0);
% val1=val1+bound;
% val2=perm(Z);
% relErr=(val1-val2)/val2
%One will see that the relative error is on the order of the finite
%precision error (a few times eps()).
%
%REFERENCES:
%[1] B. Roos, "New permanent approximation inequalities via identities,"
%    arXiv, 21 Feb. 2018. [Online]. Available:
%    https://arxiv.org/abs/1612.03702
%[2] A. Barvinok, "Polynomial time algorithms to approximate permanents
%    and mixed discriminants within a simply exponential factor," Random
%    Structures and Algorithms, vol. 14, no. 1, pp. 29-61, Jan. 1999.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    if(any(Z(:)<0))
        algorithm=2;
    else
        algorithm=0;
    end
end

if(nargin<3||isempty(numIter))
    m=size(Z,1);
    n=size(Z,2);
    
    numIter=n*n*m;
end

%Empty matrices have a permanent of 1 by definition.
if(isempty(Z))
    permApprox=1;
    bound=0;
    return; 
elseif(numel(Z)==1) 
    permApprox=Z;
    bound=0;
    return; 
end 

switch(algorithm)
    case 0
        [permApprox,bound]=permApproxRoo(Z,false);
    case 1
        [permApprox,bound]=permApproxRoo(Z,true);
    case 2
        permApprox=permApproxS(Z,numIter);
        bound=[];
    otherwise
        error('Unknown algorithm specified.')
end
end

function [permApprox,bound]=permApproxRoo(Z,useBinAlg)
%%PERMAPPROXROO Approximate the matrix permanent using the algorithm of [1]
%               using the product of the arithmetic mean of the columns of
%               Z (or Z' if Z has has more columns than rows). This simple,
%               crude approximation has a decent accuracy bound associated
%               with it in Theorem 4.1 (and Theorem 4.2 specifically for
%               binary matrices) in [1].
%
%INPUTS: Z An NXn real matrix with all-positive elements.
%  useBinAlg An optional input that indicates whether Z is all binary and
%          the algorithm of Theorem 4.2 for binary matrices should be used.
%          The default if omitted or an empty matrix is passed is false.
%
%OUTPUTS: permApprox The approximate of the matrix permanent.
%              bound The bound on the accuracy of the matrix permanent.
%                    This is such that abs(perm(Z)-permApprox)<=bound.
%
%Often, a better approximation than permApprox is permApprox+bound.
%
%REFERENCES:
%[1] B. Roos, "New permanent approximation inequalities via identities,"
%    arXiv, 21 Feb. 2018. [Online]. Available:
%    https://arxiv.org/abs/1612.03702
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<2||isempty(useBinAlg))
    useBinAlg=false; 
end

N=size(Z,1);
n=size(Z,2);

%If it is not a thin matrix, then use the transpose.
if(n>N)
    Z=Z';
    temp=N;
    N=n;
    n=temp;
end

%Defined in Section 2.
zTilde=(1/N)*sum(Z,1);

%Defined in Section 2.
pnTilde=prod(zTilde);

permApprox=exp(gammaln(N+1)-gammaln(N-n+1))*pnTilde;

if(nargout>1)
    %If the accuracy bound should be computed.
    
    %Compute the row differences
    y=zeros(N,N,n);
    for r=1:n
        for k=1:N
            y(:,k,r)=Z(:,r)-Z(k,r);
        end
    end
    y=abs(y);%We only need the absolute value.
    
    if(useBinAlg)%If using Theorem 4.2.
        if(n==2)
            kappa=1;
        else
            eta=zeros(N,N,n);
            ZRowSums=sum(Z,1);
            for u=1:N
                for v=1:N
                    eta(u,v,:)=reshape(ZRowSums-Z(u,:)-Z(v,:),[1,1,n]);
                end
            end
            zetaln=(1./eta).*gammaln(eta+1);

            kappaTildeuvl=exp(2*zetaln-(2/(N-2))*gammaln(N-2+1));
            KappaTildeLSum=sum(kappaTildeuvl,3);
            maxVal=-Inf;
            for u=1:N
                for v=1:N
                    if(u==v)
                        continue;
                    end

                    for r=1:n
                        for s=1:n
                            if(r==s)
                                continue;
                            end
                            sumVal=KappaTildeLSum(u,v)-kappaTildeuvl(u,v,r)-kappaTildeuvl(u,v,s);

                            maxVal=max(sumVal,maxVal);
                        end
                    end
                end
            end
            kappa=(1/(n-2))*maxVal;
        end
    else%We are using Theorem 4.1.
        if(n==2)
            kappa=1;
        else
            Z2=abs(Z).^2;
            Z2Sum=sum(Z2(:));

            ZColSums=sum(Z2,2);
            ZRowSums=sum(Z2,1);

            maxVal=-Inf;
            for j1=1:N
                Z2SumMinJ1=Z2Sum-ZColSums(j1);
                for j2=(j1+1):N
                    Z2SumModJ12=Z2SumMinJ1-ZColSums(j2);
                    for r1=1:n
                        Z2SumModR=Z2SumModJ12-ZRowSums(r1)+Z2(j1,r1)+Z2(j2,r1);
                        for r2=(r1+1):n
                            ZSum=Z2SumModR-ZRowSums(r2)+Z2(j1,r2)+Z2(j2,r2);

                            maxVal=max(ZSum,maxVal);
                        end
                    end
                end
            end

            kappa=(1/((n-2)*(N-2)))*maxVal;
        end
    end

    diagIdx=diagElIdx(N);
    outerSum=0;
    for r=1:n
        for s=1:n
            if(r==s)
                continue;
            end

            Y=y(:,:,r).*y(:,:,s);
            Y(diagIdx)=0;

            innerSum=sum(Y(:));
            outerSum=outerSum+innerSum^2;
        end
    end
    theta=(1/(N*(N-1)*sqrt(n*(n-1))))*sqrt(outerSum);
    
    %Defined in Theorem 4.1
    beta=(1/n)*sum(abs(zTilde).^2);
    
    coeff=exp(gammaln(N+1)-gammaln(N-n+1));
    bound=coeff*(theta/(2*N))*fn(n,sqrt(beta),sqrt(kappa));
end
end

function val=fn(n,x1,x2)
%%FN This is the fn function defined in Theorem 4.1 of [1]. The sum was
%    directly expanded.
%
%REFERENCES:
%[1] B. Roos, "New permanent approximation inequalities via identities,"
%    arXiv, 21 Feb. 2018. [Online]. Available:
%    https://arxiv.org/abs/1612.03702
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=(x1^n-x2^(n-1)*(n*(x1-x2)+x2))/(x1-x2)^2;
if(~isfinite(val))
    val=1;
end
end

function retVal=permApproxS(A,numIter)
%PERMAPPROXS Compute a stochastic approximation to the permanent of a real
%            matrix.
%
%INPUTS: A An mXn matrix. If m<=n, then the standard matrix permanent bound
%          is found. If m>n, then the permanent bound of A' is found to be
%          consistent with the permanents of A and A' being equal in square
%          matrices. Empty matrices have a permanent of one by definition.
%  numIter The number of iterations of the stochastic permanent
%          approximation algorithm to use. If this parameter is omitted or
%          an empty matrix is passed, then the default value of n*n*m is
%          used. 
%
%OUTPUTS: retVal An approxmation of the matrix permanent of A.
%
%The algorithm is the real algorithm described in [1]. The algorithm
%utilizes random sampling, so the results will not be equal when the
%algorithm is run multiple times.
%
%REFERENCES:
%[1] A. Barvinok, "Polynomial time algorithms to approximate permanents
%    and mixed discriminants within a simply exponential factor," Random
%    Structures and Algorithms, vol. 14, no. 1, pp. 29-61, Jan. 1999.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

%Empty matrices have a permanent of 1 by definition.
    if(isempty(A))
        retVal=1;
        return; 
    elseif(numel(A)==1) 
        retVal=A; 
        return; 
    end 

    m=size(A,1);
    n=size(A,2);
    
    if(m>n)
        A=A';
        temp=m;
        m=n;
        n=temp;
    end
    
    if(nargin<2||isempty(numIter))
        numIter=n*n*m;
    end
    
    curComb=0:(m-1);
    
    retVal=0;
    while(~isempty(curComb))
        retVal=retVal+permApproxStochSquare(A(:,curComb+1),numIter);
        curComb=getNextCombo(curComb,n);
    end
end

function val=permApproxStochSquare(A,numIter)
    N=size(A,1);

    alpha=zeros(numIter,1);
    for curIter=1:numIter
        U=randn(N,N);
        B=U.*sqrt(A);
        alpha(curIter)=det(B)^2;
    end

    val=mean(alpha);
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
