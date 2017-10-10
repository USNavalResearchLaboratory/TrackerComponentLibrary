function val=HermitePoly(x,orders,V, type)
%%HERMITEPOLY Evaluate univariate or multivariate Hermite polynomials of a
%             given set of multivariate orders. The definition of Hermite
%             polynomials more common the statistics (probabilists) is used
%             rather than the definition that is more common in 
%
%INPUTS: x The xDimXnumPoints matrix of numPoints points at which the
%          multivariate Hermite polynomials should be evaluated.
%   orders A numDimX1 or 1XnumDim vector of the order of the Hermite
%          polynomial, as applied to each dimension. The elements should
%          be integers >=0.
%        V An xDimXxDim symmetric positive definite matrix (a covariance
%          matrix) with respect to which the multivariate Hermite
%          polynomial is taken. If this parameter is omitted or an empty
%          matrix is passed, the identity matrix will be used.
%     type An optional value indicating the type of Hermite polynomials
%          used (there are two different definitions). Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            definition of Hermite polynomials that is common in the
%            statistics literature. This definition is in terms of
%            derivatives of the multivariate normal PDF divided by the
%            normal PDF. Such Hermite polynomials are often denoted by He_n
%            in the literature.
%          1 Use the definition that is common in the physics literature.
%            This is similar to the definition in statistics, except the
%            division by 2 in the exponent within the normal PDF is
%            omitted. Such Hermite polynomials are often denoted by H_n in
%            the literature.
%
%OUTPUTS: val A 1XnumPoints vector of the Hermite polynomial values for
%             each of the points.
%
%Multivariate Hermite polynomials find application in multivariate
%Gram-Charlier series A expansions, among other places.
%
%The algorithm implements Equation 1.7 of [1]. Note that when comparing to
%the explicitly given formulae in [1], the solution for H_{22} differs.
%The 2*z_1*z_2*mu_{11} term in the expression for H_{22} should actually be
%4*z_1*z_2*mu_{11} as the product in the sum of Equation 1.7 involves two
%binomial terms that equal one. Thus, it appears that this implementation
%is correct, despite the disagrrement.

%The definite of Hermite polynomials common to statistics is used in [1].
%The relation to the Hermite polynomials commonly used in physics H^p_n(x)
%is H^p_n(x)=2^(sum(n)/2)*HermitePoly(sqrt(2)*x,n)
%where n is the orders vector.
%
%REFERENCES:
%[1] C. S. Withers, "A simple expression for the multivariate Hermite
%    polynomials," Statistics and Probability Letters, vol. 47, no. 2, pp.
%    165-169, Apr. 2000.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(type))
    type=0; 
end

if(nargin<3||isempty(V))
    xDim=size(x,1);
    V=eye(xDim,xDim);
end

%If the physics definition is being used.
if(type==1)
    x=sqrt(2)*x;
end

%Definition of z from Equation 1.7
z=V\x;

numX=size(x,2);

maxK=fix(sum(orders)/2);

%We will precompute all of the moments of a zero-mean multivariate normal
%random variable with covariance matrix inv(V) up to order 2*maxK.
momentMat=genNormalMoments(2*maxK,inv(V));

%This sum implements Equation 1.7
sumVal=zeros(1,numX);
signVal=1;
for k=0:maxK
    sumVal=sumVal+signVal*h(orders,k,z,momentMat);

    signVal=-signVal;
end

val=sumVal;

%If the physics definition is being used.
if(type==1)
     val=2^(sum(orders)/2)*val;
end

end

function momentMat=genNormalMoments(maxDeg,V)
%Generate all central moments from degree 0 to maxDeg of a multivariate
%normal distribution having covariance matrix V.

numDim=size(V,1);
mu=zeros(numDim,1);%Use mu=0 for central moments.

numMaxDims=repmat(maxDeg+1,[1,numDim]);
if(size(numMaxDims,2)==1)
    %Deal with having a singleton dimension.
    numMaxDims=[numMaxDims,1];
end
momentMat=zeros(numMaxDims);

%Go through all of the partitions of the elements for each order and then
%all of the permutations per partition.
for curDeg=0:maxDeg
    pRecur=getNextMPartition(numDim+curDeg,numDim);

    while(~isempty(pRecur))
        %Get all permutations of pRecur;
        permMat=genAllMultisetPermutations(pRecur);
        numPerm=size(permMat,2);
        for curPerm=1:numPerm
            idxVal=nDim2Index(numMaxDims,permMat(:,curPerm));

            momentMat(idxVal)=GaussianD.momentGenFun(mu,V,permMat(:,curPerm)-1);
        end

        pRecur=getNextMPartition(pRecur);
    end
end

end


function sumVal=h(v,k,z,momentMat)
%This function implements the second half of Equation 1.7, the function
%h_{v,k}(z). momMat must hold moments of a normal (0,V) random variable.

numZ=size(z,2);

numMomDim=size(momentMat,1);
p=length(v);

sumVal=zeros(1,numZ);
%We have to go through all partitions of p values such that they sum to
%2*k. We can have zero values.
curPar=getNextMPartition(2*k+p,p);
while(~isempty(curPar))
    %Now, we must go through all multiset permutations of the partition.
        permVals=genAllMultisetPermutations(curPar);
        numN=size(permVals,2);

        for curN=1:numN
            curPerm=permVals(:,curN);
            n=curPerm-1;
            %Evaluate the product.
            prodVal=ones(1,numZ);
            for j=1:p
                %The added test deals with when the binomial value is zero.
                %Testing for it eliminates a 0*Inf problem if the z term is
                %zero.
                if(v(j)>=n(j))
                    prodVal=prodVal.*binomial(v(j),n(j)).*z(j,:).^(v(j)-n(j));
                else
                    prodVal=0;
                    break;
                end
            end
            
            momIdx=nDim2Index(numMomDim,curPerm(:));
            sumVal=sumVal+momentMat(momIdx)*prodVal;
        end
    curPar=getNextMPartition(curPar);
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
