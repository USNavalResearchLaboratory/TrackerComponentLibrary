function PDFVals=GramCharlierAPDFApprox(xVals,quasiMomentMat,muN,SigmaN)
%%GRAMCHARLIERAPDFAPPROX Obtain approximate values for a multivariate PDF
%               based on its quasi-moments using a Gram-Charlier series of
%               type A approximation. This approximation does not always
%               asymptotically converge to the true PDF producing the
%               quasi-moments as the number of quasi-moments used
%               increases. For example, for a scalar PDF, the asymptotic
%               decrease in the distribution must be faster than
%               exp(-x^2/4) to assure convergence.
%
%INPUTS: xVals The xDimXN set of N points for which the PDF values are
%              desired.
% quasiMomentMat A matrix taking n indices, where 
%              quasiMomentMat(a1,a2,a3...an) corresponds to the quasi-
%              moment whose multivariate order is given by a1-1,a2-1, etc.
%              It is assumed that quasi moments only up to an order equal
%              to the size of the first dimension of the matrix are
%              available (other entries in quasiMomentMat are ignored).
%              All other dimensions must be at least the same size as the
%              first. Also, the zero-th order moment must be 1 (as the PDF
%              integrates to 1). If a univariate distribution is used,
%              then this is a column vector.
%          muN The nX1 mean vector with respect to with the quasi-moments
%              are computed. This should generally be the mean of the
%              distribution being approxiamted. If this parameter is
%              omitted or an empty matrix is passed, a zero mean vector is
%              used.
%       SigmaN The nXn covariance matrix with respect to which the
%              quasi-moments are computed. This should generally be the
%              covariance matrix of the distribution being approximated.
%              If this parameter is omitted or an empty matrix is passed,
%              the identity matrix is used.
%
%OUTPUTS: PDFVals The PDF approximated using the quasi-moments evaluated at
%               the points in xVals. If multiple points are given, the
%               result is a row vector. The values are not guaranteed to
%               always be positive.
%
%The Gram-Charlier A expansion is used in numerous papers on state
%estimation, but under various names, such as in [2,3,4,5]. A development
%of the expansion is given in [1]. The expressions here differ from [4] in
%that in [4] an incorrect normalization for multivariate distributions is
%used. Here, the multivariate Hermite polynomials of [6] are used. These
%polynomials specifically take a covariance matrix for a centered
%distribution.
%
%REFERENCES:
%[1] P. I. Kuznetsov, R. L. Stratonovich, and V. I. Tikhonov, "Quasi-moment
%    functions in the theory of random processes," Theory of Probability 
%    and its Applications, vol. V, no. 1, pp. 80-97, 1960.
%[2] H. Singer, "Generalized Gauss-Hermite filtering for multivariate
%    diffusion processes," FernUniversität Hagen, Tech. Rep., 2006.
%    [Online]. Available: http://deposit.fernuni-hagen.de/87/
%[3] H. Singer, "Generalized Gauss-Hermite filtering," AStA Advances in
%    Statistical Analysis, vol. 92, no. 2, pp. 179-195, May 2008.
%[4] S. Challa, Y. Bar-Shalom, and V. Krishnamurthy, "Nonlinear filtering
%    via generalized Edgeworth series and Gauss-Hermite quadrature," IEEE
%    Transactions on Signal Processing, vol. 48, no. 6, pp. 1816-1820, Jun.
%    2000.
%[5] W. W. Willman, "Edgeworth expansions in state perturbation estimation,"
%    IEEE Transactions on Automatic Control, vol. AC-26, no. 2, pp.
%    493-498, 1981.
%[6] C. S. Withers, "A simple expression for the multivariate Hermite
%    polynomials," Statistics and Probability Letters, vol. 47, no. 2, pp.
%    165-169, Apr. 2000.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDimList=size(quasiMomentMat);
maxDeg=numDimList(1)-1;
numIdx=length(numDimList);

%If the last dimension is just a singleton dimension, then shrink
%everything by one.
if(numIdx==2&&numDimList(2)==1&&(nargin<4||size(SigmaN,2)==1))
    numIdx=1;
end

%Do the conversion with respect to the standard normal distribution if no
%further details are given.
if(nargin<2||isempty(muN))
    muN=zeros(numIdx,1);
end
if(nargin<3||isempty(SigmaN))
    SigmaN=eye(numIdx,numIdx);
end

numX=size(xVals,2);

sumVal=ones(1,numX);%The zeroth-order term is always zero.
zVals=bsxfun(@minus,xVals,muN);%The argument to the Hermite polynomials

%We will precompute all of the inverse factorials needed in the products in
%the sum.
invFactList=1./factorial(0:maxDeg);

for curOrder=1:maxDeg
    %Now, we go through all compositions of curOrder items into  numIndex
    %slots.
    numCompositions=binomial(curOrder-1,numIdx-1);
    
    for k=0:(numCompositions-1)
        curComp=unrankTComposition(k,numIdx,curOrder);
        
        curEl=nDim2Index(numDimList,curComp);
        sumVal=sumVal+prod(invFactList(curComp))*quasiMomentMat(curEl)*HermitePoly(zVals,curComp-1,SigmaN);
    end
end

PDFVals=GaussianD.PDF(xVals,muN,SigmaN).*sumVal;

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
