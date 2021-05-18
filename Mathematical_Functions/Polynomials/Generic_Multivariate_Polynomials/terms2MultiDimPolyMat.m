function coeffs=terms2MultiDimPolyMat(termMat)
%%TERMS2MULTIDIMPOLYMAT Given a multidimensional polynomial represented as
%                    a 2D matrix of terms (monomials), convert it into a
%                    hypermatrix of coefficients suitable for use with
%                    functions such as polyValMultiDim and where convn can
%                    be used to multiply two such polynomials.
%
%INPUTS: termMat An (n+1)XnumTerms matrix such that
%                termMat(:,i)=[c,a1,a2,...,an] is a monomial term where c
%                is the value of of the monomial coefficient and the
%                monomial is
%                x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1). The ordering
%                of the terms in termMat does not matter.
%
%OUTPUTS: coeffs A hypermatrix of the coefficients for the multivariate
%               polynomial. These are arranged such that
%               coeffs(a1,a2,a3...an) corresponds to the coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in coeffs is
%               the maximum order of that dimension +1.
%
%This function is the opposite of multiDimHyperMat2Terms. See
%multiDimHyperMat2Terms for usage examples.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTerms=size(termMat,2);
numIdx=size(termMat,1)-1;

numDims=zeros(1,numIdx);

if(isscalar(numDims))
    numDims=[numDims,1];
end

for idx=1:numIdx
    numDims(idx)=max(termMat(idx+1,:))+1; 
end

coeffs=zeros(numDims);
for curTerm=1:numTerms 
    idx=nDim2Index(numDims,termMat(2:end,curTerm)+1);
    coeffs(idx)=coeffs(idx)+termMat(1,curTerm);
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
