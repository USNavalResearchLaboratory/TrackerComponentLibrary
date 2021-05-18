function [index,numBeforeDeg]=monomial2DegNegLexOrderIdx(monomial,numBeforeDeg)
%%MONOMIAL2DEGNEGLEXORDERIDX Given the nonegative integer exponents of a
%          monomial term (e.g. x^5*y^3*z^6 is [5;3;6]), determine the index
%          of that monomial term in a list of all monomial terms in degree
%          negative lexicographic ordering. The degree of a monomial is the
%          sum of the exponents of the terms.
%
%INPUTS: monomial A numDimXnumMonomials set of numMonomials monomial terms
%                 consisting of nonnegative integer vectors of exponents of
%                 the elements Alternatively, if one just wants the
%                 numBeforeDegVec vector, then this should be an empty
%                 matrix and numBeforeDegVec should contain the
%                 dimensionality and the maximum degree desired.
%    numBeforeDeg An optional parameter. Providing this parameter will
%                 speed up repeated calls to this function by avoiding the
%                 recalculation of certain values. This vector is such that
%                 numBeforeDeg(i+1) is the number of monomial terms
%                 before degree i starting with degree 0, where
%                 numBeforeDeg(0+1)=0. This parameter is returned by the
%                 function. Note that numBeforeDeg depends on numDim, so
%                 vectors cannot be reused between cals where numDim
%                 changes. If one just wants this parameter (for example
%                 in advance to speed up many repeated calls), then set
%                 monomial=[] and numBeforeDeg=[numDim;maxDeg], where
%                 maxDeg the maximum degree desired, >=0.
%
%OUTPUTS: index A numMonomialsX1 vector of the position of the monomials in
%              the ordering of monomials in degree negative lexicographic
%              order, starting from 1. If monomial=[] and numBeforeDeg
%              is set to the dimensionality and maximum degree, then this
%              is an empty matrix.
% numBeforeDeg The same vector as the optional input. If the input is
%              not provided, then numBeforeDeg will be filled with
%              values up to the degree of the monomial provided having the
%              highest degree. If monomials=[] and numBeforeDeg on the
%              input is set to the maximum degree desired, then
%              numBeforeDeg will be filled all the way up to that
%              degree.
%
%Various monomial orderings are used in algorithms that manipulate
%polynomials, such as methods for solving or minimizing multivariate
%polynomials as in [1] and [2].
%
%The degree negative lexicographic ordering, for a fixed degree, coincides
%with the ordering of the terms provided by unrankTComposition with
%firstElMostSig=true. However, that rank is only an offset for a given
%degree. The number of terms before that degree is the sum of the number of
%compositions of lower degrees. Thus, one sees the need for the
%numBeforeDeg vector. The number of compositions of a given degree into
%n parts, (each part starting from 0, not 1) is
%binomial(degree+numDim-1,numDim-1).
%
%EXAMPLE:
%Consider monomials like x^1*y^0*z^0, x^12*y^4*z^8;, x*y*z, x^0*y^0*z^0(=1)
% monomial=[1,12,1,0;
%           0, 4,1,0;
%           0, 8,1,0];
% [index,numBeforeDeg]=monomial2DegNegLexOrderIdx(monomial)
% %One will get the indices index=[2;2687;15;1].
% %The value numBeforeDeg can be reused with monomials having the same
% %number of variables, such as
% monomial=[0;1;1];
% index=monomial2DegNegLexOrderIdx(monomial,numBeforeDeg)
% %On the other hand, if one only wants the numBeforeDeg vector for up to
% %degree 24 for three variables, then one
% [~,numBeforeDeg]=monomial2DegNegLexOrderIdx([],[3,24])
%
%REFERENCES:
%[1] P. Dreesen, "Back to the roots: Polynomial system solving using linear
%    algebra," Ph.D. dissertation, Katholieke Universiteit Leuven, Leuven,
%    Flanders, Belgium, Sep. 2013.
%[2] A. Papachristodoulou and S. Prajna, "A tutorial on sum of squares
%    techniques for systems analysis," in Proceedings of the American
%    Control Conference, Portland, OR, 8-10 Jun. 2005, pp. 2686-2700.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(monomial))
    numDim=numBeforeDeg(1);
    maxDeg=numBeforeDeg(2);
else
    numDim=size(monomial,1);
    degs=sum(monomial,1);
    maxDeg=max(degs);
end

%Fill in the numBeforeDegVec vector.
if(isempty(monomial)||nargin<2||isempty(numBeforeDeg))
    numBeforeDeg=zeros(maxDeg+1,1);
    
    if(maxDeg>0)
        numBeforeDeg(1+1)=1;
    end
    
    for curDeg=2:maxDeg
        numBeforeDeg(curDeg+1)=numBeforeDeg(curDeg-1+1)+binomial(curDeg-1+numDim-1,numDim-1);
    end
    
    %If just the list numBeforeDegVec is desired.
    if(isempty(monomial))
        index=[];
        return;
    end
end

numMonomials=size(monomial,2);
%Allocate space.
index=zeros(1,numMonomials);

for curMonomial=1:numMonomials
    offset=numBeforeDeg(degs(curMonomial)+1);
    index(curMonomial)=rankTComposition(monomial(:,curMonomial)+1,true)+offset+1;
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
