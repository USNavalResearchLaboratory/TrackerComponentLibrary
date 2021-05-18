function curPoly=multiDimPolyMat2String(coeffs,xVar,maxDigits,complexFormat)
%%MULTIDIMPOLYMAT2STRING Given a multivariate polynomial represented as a
%         hypermatrix of its coefficients (including cross terms), obtain a
%         string that is suitable for display or for using as an entry in a
%         cell array in the input to the function solvePolySysWithExtProg.
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
%      xVar An nX1 cell array where xVar{j} corresponds to the variable for
%           coefficient aj in termMat. These are the names of the variables
%           given as character strings, which are used by PHCpack to
%           represent the unknown random variables. The names cannot being
%           with a number.
% maxDigits The maximum number of digits to display for each number. The
%           default if this parameter is omitted or an empty matrix is
%           passed is 16, which is sufficient for floating point numbers.
% complexFormat A string indicating the value to use to specify an
%           imaginary number. By default, this is i, which is suitable for
%           Bertini and PHCpack. However, one might want to use ii if
%           complex numbers are to be used in Macaulay2.
%
%OUTPUTS: curPoly A string representing the polynomial expressed in coeffs
%                 using xVar as the variables. Note that this function does
%                 not combine like terms.
%
%This function just calls multiDimPolyMat2Terms and then terms2String.
%
%EXAMPLE:
% coeffs=zeros(2,2,2,2);
% coeffs(1+1,0+1,1+1,0+1)=1;
% coeffs(0+1,1+1,1+1,0+1)=-2;
% coeffs(1+1,0+1,0+1,0+1)=-6;
% coeffs(0+1,1+1,0+1,0+1)=-1;
% coeffs(0+1,0+1,0+1,1+1)=-4;
% coeffs(0+1,0+1,0+1,0+1)=12;
% 
% xVar={'xa','xb','xc','xd'};
% curPoly=multiDimPolyMat2String(coeffs,xVar)
% %The result is
% curPoly='12-4*xd-xb-2*xb*xc-6*xa+xa*xc'
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    maxDigits=[];%Use default in terms2String
end

if(nargin<4)
    complexFormat=[];%Use default in terms2String
end

curPoly=terms2String(multiDimPolyMat2Terms(coeffs),xVar,maxDigits,complexFormat);

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
