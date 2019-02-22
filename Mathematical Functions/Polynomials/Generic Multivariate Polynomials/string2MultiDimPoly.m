function coeffs=string2MultiDimPoly(theString)
%%STRING2MULTIDIMPOLY Given a string representing a multidimensional
%           polynomial equation, convert it into a hyoermatrix
%           representation of the polynomial, as can be used in function
%           such as polyValMultiDim and polyDerMultiDim.
%
%INPUTS: theString The string representing the multivariate polynomial.
%           All coefficients come before variables. Variables are all x
%           followed by a number, such as x3. All numbers for variables
%           must be >0. The format of this string is the same as that
%           returned by multiDimPolyMat2String. For example, 
%           '60*x2^2-4*x2^3+72*x1*x2+12' is a valid string.
%
%OUTPUTS: coeffs A hypermatrix of the coefficients for the multivariate
%                polynomial. These are arranged such that
%                coeffs(a1,a2,a3...an) corresponds to the coefficient of an
%                x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%                number of indices coeffs takes is equal to the
%                dimensionality of x (not counting singleton dimensions at
%                the end of coeffs). Note that this ordering is the reverse
%                that used in the 1D polyval function that is built into
%                Matlab. The number of elements for each index in coeffs is
%                the maximum order of that dimension +1.
%
%This function just calls string2Terms and terms2MultiDimPolyMat.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

coeffs=terms2MultiDimPolyMat(string2Terms(theString));

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
