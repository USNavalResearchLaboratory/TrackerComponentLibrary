function coeffsRet=polyDerMultiDim(coeffs,derivDim)
%%POLYDERMULTIDIM Compute the derivative of a multivariate polynomial with
%                 respect to a particular dimension given a hypermatrix of
%                 its coefficients (including cross terms). 
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
%      derivDim The dimension of x (index of coeffs) with respect to which
%               the derivative is to be taken.
%
%OUTPUTS: coeffsRet The derivative of the coeffs hypermatrix. The
%                   dimensionality is sized to fit the number of nonzero
%                   coefficients. Note that Matlab suppresses trailing
%                   singleton dimensions (for more than a 2D matrix) when
%                   using the size command and appends one singleton
%                   dimension when dealing with scalars. 
%
%The function is just implemented using the basic rules of polynomial
%differentiation.
%
%As an example, consider the multivariate polynomial
%x1^3+12*x1+24*x2-36*x3-3*x1*x2*x3+16*x3^3+6*x1^2*x3-18*x1*x2+64
%The derivatives with respect to the first, second, and third coordinates
%are
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
% dcoeffsd1=polyDerMultiDim(coeffs,1);
% dcoeffsd2=polyDerMultiDim(coeffs,2);
% dcoeffsd3=polyDerMultiDim(coeffs,3);
%One will see that dcoeffsd1 corresponds to
%12+3*x1^2-18*x2+12*x1*x3-3*x2*x3
%dcoeffsd2 corresponds to
%24-18*x1-3*x1*x3
%and dcoeffsd3 corresponds to
%-36+6*x1^2-3*x1*x2+48*x3^2
%
%%Note that when evaluating a 1D polynomial the order of the coefficients
%is reverse that of the polyder function in Matlab. For example, to
%evaluate ther derivative of 2*x^2-3*x+4 one uses
% coeffs=[4;-3;2];
% dcoeffs=polyDerMultiDim(coeffs,1);
%to get dcoeffs=[-3;4;0].
%
%Note that the convn function can be used for multivariate polynomial
%multiplication.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dimSizeList=size(coeffs);
numIdx=length(dimSizeList);

%If we are trying to differentiate with respect to a singleton index, then
%the result is zero as no instances of that term arise.
if(derivDim>numIdx||dimSizeList(derivDim)==1)
    coeffsRet=0;
    return;
end

%Rather than explicitly programming the loops, which would be slow in
%Matlab, we take advantage of being able to index Matlab variables using
%cell arrays. Here, we zero all terms that do not include derivDim.
idxVec=cell(numIdx,1);
idxVec(1:(derivDim-1))={':'};
idxVec{derivDim}=1;
idxVec((derivDim+1):end)={':'};
coeffs(idxVec{:})=0;

%Shift everything in the derivative dimension by one to perform
%differentiation of components involving those elements, not including the
%constants by which the terms must be multiplied.
coeffs=circshift(coeffs,-1,derivDim);

%Now, we have to multiply all of the shifted terms by the appropriate
%constants of differentiation. It is just the previous order of the
%derivDim terms.
for i=1:(dimSizeList(derivDim)-1)
    idxVec{derivDim}=i+1;
    coeffs(idxVec{:})=coeffs(idxVec{:})*(i+1);
end

%Finally, the dimensionality of the derivDim has been reduced by one. Also,
%all of the terms not multiplied by a power of derivDim are gone, meaning
%that the dimensionality of their components might also be reduced.
coeffsRet=shrinkMultiDimPoly2Fit(coeffs);

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
