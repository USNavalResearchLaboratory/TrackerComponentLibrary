function y=polyValVec(P,x)
%%POLYVALVEC Evaluate a set of polynomials at a given scalar point or set
%            of points. Evaluating a collection of polynomials produces a
%            vector estimate. This could, for example, be used to evaluate
%            polynomials to interpolate each component of a vector quantity
%            to a desired point given the ointerpolating polynomials for
%            each dimension.
%            
%INPUTS: P A numDimX(numDeg+1) matrix of numDim polynomial coefficients for
%          each of the numDim components. numDeg is the degree of the
%          polynomials. The coefficients for each row are arranged for
%          exponents in the order [x^numDeg, x^(numDeg-1),...,1,0].
%        x A scalar value of x or a linear vector of numX points at which
%          the polynomials should be evaluated.
%
%OUTPUTS: y The numDim X numX set of polynomials evaluated at each of the
%           numX points.
%
%The algorithm just calls Matlab's polyval function multiple times in a
%loop.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(P,1);
numX=length(x);

%Make sure that it is a row vector.
x=reshape(x,1,numX);

y=zeros(numDim,numX);
for curDim=1:numDim
    y(curDim,:)=polyval(P(curDim,:),x);
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
