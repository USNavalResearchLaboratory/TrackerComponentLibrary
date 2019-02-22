function sumCoeffs=polySumMultiDim(coeffs1,coeffs2)
%%POLYSUMMULTIDIM Add together two multidimensional polynomials given their
%                 coefficients. The ordering of the coefficients is as in
%                 the polyValMultiDim and polyDerMultiDim functions and for
%                 1D polynomials, is the opposite that used in polyval and
%                 polySum.
%
%INPUTS: coeffs1,coeffs2 Hypermatrices of the coefficients of the two
%               multidimensional polynomials. The coefficients are
%               arranged such that coeffs(a1,a2,a3...an) corresponds to the
%               coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in coeffs is
%               the maximum order of that dimension +1. 
%
%OUTPUTS: sumCoeffs A hypermatrix representing the sum of the two
%                   multivariate polynomials. Note that dimensions are not
%                   shrunk when the sum of high order equals zero (terms
%                   cancel).
%
%Polynomials are summed by adding corresponding terms.
%
%As an example, consider adding 
%x1^3+12*x1+24*x2-36*x3-3*x1*x2*x3+16*x3^3+6*x1^2*x3-18*x1*x2+64
%to
%-36+6*x1^2-3*x1*x2+48*x3^2
%This is done as
% coeffs1=zeros(4,2,4);
% coeffs1(3+1,0+1,0+1)=1;
% coeffs1(1+1,0+1,0+1)=12;
% coeffs1(0+1,1+1,0+1)=24;
% coeffs1(0+1,0+1,1+1)=-36;
% coeffs1(1+1,1+1,1+1)=-3;
% coeffs1(0+1,0+1,3+1)=16;
% coeffs1(1+2,0+1,1+1)=6;
% coeffs1(1+1,1+1,1+0)=-18;
% coeffs1(0+1,0+1,0+1)=64;
% 
% coeffs2=zeros(2,1,2);
% coeffs2(0+1,0+1,0+1)=-36;
% coeffs2(2+1,0+1,0+1)=6;
% coeffs2(1+1,1+1,0+1)=-3;
% coeffs2(0+1,0+1,2+1)=48;
% sumCoeffs=polySumMultiDim(coeffs1,coeffs2)
%One will see that the results corresponds to the polynomial
%x1^3+12*x1+24*x2-36*x3-3*x1*x2*x3+16*x3^3+6*x1^2*x3+6*x1^2-21*x1*x2+28+48*x3^2
%
%Note that when adding 1D polynomials, the ordering of the elements is
%revserse that used in polySum. Thus, to add
%x^2-3*x+2
%to
%x^5+3*x^3-4*x^2+1
%one uses
% coeffs1=[2;-3;1];
% coeffs2=[1;0;-4;3;0;1];
% sumCoeffs=polySumMultiDim(coeffs1,coeffs2)
%
%Note that the convn function can be used for multivariate polynomial
%multiplication.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dim1SizeList=size(coeffs1);
dim2SizeList=size(coeffs2);

numDims1=length(dim1SizeList);
numDims2=length(dim2SizeList);

%Determine the dimensionality of the sum.
dimSumSizeList=max([dim2SizeList,zeros(1,numDims1-numDims2)],[dim1SizeList,zeros(1,numDims2-numDims1)]);

numEls1=prod(dim1SizeList);
numEls2=prod(dim2SizeList);
numElsSum=prod(dimSumSizeList);
sumCoeffs=reshape(zeros(numElsSum,1),dimSumSizeList(:)');%Allocate space

%Go through all of the elements and find the sum. The simplest way is just
%to add all of coeffs1 to sumCoeffs and then add all of coeffs2 to
%sumCoeffs. 
for curEl=1:numEls1
    idxSum=nDim2Index(dimSumSizeList,index2NDim(dim1SizeList,curEl));
    sumCoeffs(idxSum)=coeffs1(curEl);
end

for curEl=1:numEls2
    idxSum=nDim2Index(dimSumSizeList,index2NDim(dim2SizeList,curEl));
    sumCoeffs(idxSum)=sumCoeffs(idxSum)+coeffs2(curEl);
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
