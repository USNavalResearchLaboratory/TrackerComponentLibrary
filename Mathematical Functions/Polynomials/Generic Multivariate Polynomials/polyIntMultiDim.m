function coeffsRet=polyIntMultiDim(coeffs,intDim,intConst)
%%POLYINTMULTIDIM Compute the indefinite integral of a multivariate
%                 polynomial with respect to a particular dimension given a
%                 hypermatrix of its coefficients (including cross terms).            
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
%        intDim The dimension of x (index of coeffs) with respect to which
%               the integral is to be taken.
%      intConst The integration constant. If this parameter is omitted or
%               an empty matrix is passed, a default value of 0 is used.
%
%OUTPUTS: coeffsRet The indefinite integral of the coeffs hypermatrix.
%
%This function is just implemented using the basic rules of polynomial
%integration.
%
%As an example, consider the multivariate polynomial:
%14+3*x1^2-18*x2+12*x1*x3-3*x2*x3
%The integrals with respect to the first, second and third variables are
% coeffs=zeros(3,2,2);
% coeffs(0+1,0+1,0+1)=14;
% coeffs(2+1,0+1,0+1)=3;
% coeffs(0+1,1+1,0+1)=-18;
% coeffs(1+1,0+1,1+1)=12;
% coeffs(0+1,1+1,1+1)=-3;
% intCoeffs1=polyIntMultiDim(coeffs,1);
% intCoeffs2=polyIntMultiDim(coeffs,2);
% intCoeffs3=polyIntMultiDim(coeffs,3);
%One will see that intCoeffs1 corresponds to
%14*x1+x1^3-18*x1*x2+6*x1^2*x3-3*x1*x2*x3
%intCoeffs2 corresponds to
%14*x2*+3*x1^2*x2-9*x2^2+12*x1*x2*x3-(3/2)*x2^2*x3
%and intCoeffs3 corresponds to.
%14*x3+3*x1^2*x3-18*x2*x3+6*x1*x3^2-(3/2)*x2*x3^2
%
%Note that the fucntion can be used with variables that are not present in
%the original polynomial. In the above example,
% intCoeffs4=polyIntMultiDim(coeffs,4);
%provides coefficients corresponding to
%14*x4+3*x1^2*x4-18*x2*x4+12*x1*x3*x4-3*x2*x3*x4
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(intConst))
   intConst=0; 
end

dimSizeList=size(coeffs);
numIdx=length(dimSizeList);
totalEls=prod(dimSizeList);

%If the dimension of integration is a singleton dimension that has been
%suppressed from the end of dimSizeList, then we must add in all of the
%missing singleton dimensions.
dimRetSizeList=[dimSizeList,ones(1,intDim-numIdx)];
%Integration will increase the maximum degree of the term over which
%integration is performed by one, meaning that the matrix has to be made
%larger.
dimRetSizeList(intDim)=dimRetSizeList(intDim)+1;

%Allocate space for the return coefficients.
coeffsRet=zeros(dimRetSizeList);

%We will now loop through all of the values in the original coefficient
%set, adjusting them and inserting them into the proper location in the new
%modified coefficient set, after multiplication by the appropriate
%constant.

if(intDim>numIdx)
%This is used when integrating over a trailing singleton dimension.
    indicesFull=ones(numIdx,1);
    indicesFull(intDim)=2;
end
for curEl=1:totalEls
    curCoeff=coeffs(curEl);
    indices=index2NDim(dimSizeList,curEl);
    
    if(intDim<=numIdx)
        constVal=1/indices(intDim);
        indices(intDim)=indices(intDim)+1;
        coeffsRet(nDim2Index(dimRetSizeList,indices))=constVal*curCoeff;
    else%If we are integration over a trailing singleton dimension.
        indicesFull(1:numIdx)=indices;
        coeffsRet(nDim2Index(dimRetSizeList,indicesFull))=curCoeff;
    end
end

%Add in the integration constant. It is the zeroth-order term.
coeffsRet(1)=intConst;

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
