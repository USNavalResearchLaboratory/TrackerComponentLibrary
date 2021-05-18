function coeffsRed=polyValMultiDim1Comp(coeffs,evalIdx,xEval)
%%POLYVALMULTIDIM1COMP Get the multidimensional polynomial (or scalar
%                 value) obtained by substituting in a value for one
%                 component of a multidimensional polynomial. Unlike
%                 polyValMultiDim, which substitutes all of the variables
%                 at once, this only substitutes for a single component and
%                 can thus be useful as an intermediate step in evaluating
%                 definite integrals over multivariate polynomials.
%
%INPUTS: coeffs A hypermatrix of the coefficients for the multivariate
%               polynomial. These are arranged such that
%               coeffs(a1,a2,a3...an) corresponds to the coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in coeffs is
%               the maximum order of that dimension +1.
%       evalIdx The index of the dimension of coeffs corresponding to the
%               variable that is to be evaluated. If evalIdx is greater
%               than the number of indices in coeffs, it is assumed to
%               correspond to a trailing singleton dimension (and thus the
%               variable does not appear in the polynomial).
%         xEval The value of the variable at index evalIdx that is to be
%               evaluated.
%
%OUTPUTS: coeffsRed The multivariate polynomial obtained by substituting in
%               xEval for the variable at evalIdx. Note that the number of
%               indices of the multivariate polynomial is not reduced (to
%               avoid confusion about which dimension corresponds to which
%               variable). Rather, the dimension in question is reduced to
%               a singleton. Note, however, that Matlab automatically
%               eliminates trailing singleton dimensions.
%
%This function just implements the basic rules of polynomial evaluation.
%
%As an example, consider the multivariate polynomial
%x1^3+12*x1+24*x2-36*x3-3*x1*x2*x3+16*x3^3+6*x1^2*x3-18*x1*x2+64
%where we would like to reduce it by substituting x3=-4. This is done as
% coeffs=zeros(4,2,4);
% coeffs(3+1,0+1,0+1)=1;
% coeffs(1+1,0+1,0+1)=12;
% coeffs(0+1,1+1,0+1)=24;
% coeffs(0+1,0+1,1+1)=-36;
% coeffs(1+1,1+1,1+1)=-3;
% coeffs(0+1,0+1,3+1)=16;
% coeffs(1+2,0+1,1+1)=6;
% coeffs(1+1,1+1,1+0)=-18;
% coeffs(0+1,0+1,0+1)=64;
% coeffsRed=polyValMultiDim1Comp(coeffs,3,-4)
%One will get a result corresponding to the polynomial
%-816+12*x1-24*x1^2+x1^3+24*x2-6*x1*x2
%
%Note that the convn function can be used for multivariate polynomial
%multiplication.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The maximum degrees of the polynomials are size(coeffs)-1.
highestDegsP1=size(coeffs);
numIdx=length(highestDegsP1);

if(evalIdx>numIdx||highestDegsP1(evalIdx)==1)
   %If the index of the variable corresponds to a trailing singleton
   %dimension, or other singleton dimension.
    coeffsRed=coeffs;
    return;
end

%Allocate space for the returned coefficients.
highestDegsPRed=highestDegsP1;
highestDegsPRed(evalIdx)=1;
coeffsRed=zeros(highestDegsPRed);

%First, we will compute all of the necessary powers of xEval.
maxDegP1=highestDegsP1(evalIdx);
powList=zeros(maxDegP1,1);
powList(1)=1;
powList(2)=xEval;
for curPow=3:maxDegP1
    powList(curPow)=powList(curPow-1)*xEval;
end

%Now, go through all of the values in the original set of coefficients,
%multiply by the appropriate power of xEval, and put the result into the
%proper place in the final coefficient set.
totalEls=prod(highestDegsP1);
for curEl=1:totalEls
    idxListOrig=index2NDim(highestDegsP1,curEl);
    powIdx=idxListOrig(evalIdx);
    
    idxListOrig(evalIdx)=1;
    newEl=nDim2Index(highestDegsPRed,idxListOrig);

    coeffsRed(newEl)=coeffsRed(newEl)+powList(powIdx)*coeffs(curEl);
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
