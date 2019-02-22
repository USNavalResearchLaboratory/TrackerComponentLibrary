function val=polyDefIntMultiDim(coeffs,xLows,xUps)
%%POLYDEFINTMULTIDIM Compute the value of the definite integral of a
%               multidimensional polynomial given lower and upper bounds of
%               integration for each variable.
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
%         xLows An nX1 or 1Xn set of values of the lower integration
%               bounds for each of the dimensions of the polynomial.
%         xUps An nX1 or 1Xn set of values of the upper integration
%               bounds for each of the dimensions of the polynomial.
%
%OUTPUTS: val The value of the definite integral.
%
%This function is just implemented using the basic rules of polynomial
%integration. Indefinite itnegrals are taken and appropraite bounds are
%substituted.
%
%As an example, consider the multivariate polynomial:
%14+3*x1^2-18*x2+12*x1*x3-3*x2*x3
%The definite integral from -5->5, 0->4, and -12->-3 for each of the
%variables is
% coeffs=zeros(3,2,2);
% coeffs(0+1,0+1,0+1)=14;
% coeffs(2+1,0+1,0+1)=3;
% coeffs(0+1,1+1,0+1)=-18;
% coeffs(1+1,0+1,1+1)=12;
% coeffs(0+1,1+1,1+1)=-3;
% xLows=[-5;0;-12];
% xUps=[5;4;-3];
% val=polyDefIntMultiDim(coeffs,xLows,xUps)
%One will obtain the correct result of 17280.
%
%Note that the convn function can be used for multivariate polynomial
%multiplication.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDims=length(xLows);

%Integrate over each dimension.
for curIntDim=numDims:-1:1
    %First, find the indefinite integral
    coeffs=polyIntMultiDim(coeffs,curIntDim,0);
    
    %Next, evaluate it at the bounds and take the difference.
    coeffs=polyValMultiDim1Comp(coeffs,curIntDim,xUps(curIntDim))-polyValMultiDim1Comp(coeffs,curIntDim,xLows(curIntDim));
end

%In the end, coeffs has become the scalar value to return.
val=coeffs;

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
