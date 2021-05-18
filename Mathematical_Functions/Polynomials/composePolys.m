function cPoly=composePolys(a,b)
%%COMPOSEPOLYS Compose two polynomials given as power series. This means,
%              for polynomials A(x) and B(x) in a scalar x, evaluate the
%              composition C(x)=A(B(x)). C(x) itself is a polynomial. Thus,
%              the output of the function is the set of coefficients of C
%              (subject to finite precision limitations in the recursion).
%
%INPUTS: a Coefficients for a power series of the form
%          y(z)=a(end)+a(end-1)*z+a(end-2)*z^2+...
%          into which the power series in b is to be substituted. The
%          ordering of the coefficients is the same as in Matlab's polyval
%          function.
%        b Coefficients for the power series that is substitued into a.
%
%OUTPUTS: c Coefficients for the power series that represents the
%           substitution of the polynomial b into a.
%
%The substitution is performed by explicitly evaluating power of b (via
%repeated convolutions), multiplying the powers by the appropriate value in
%a and summing the results.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

orderA=length(a)-1;
orderB=length(b)-1;

%The order of the output.
prodOrder=orderA*orderB;

%Zero-pad b so that it is the same length as the length of the output. This
%simplifies assignments in the loop.
b=[zeros(prodOrder-orderB,1);b(:)];

%bProd will hold powers of b. Right now, it is just initialized to 1 (b^0).
bProd=[zeros(prodOrder,1);1];

%Initialize the output with the zeroth-order term.
cPoly=zeros(prodOrder+1,1);
cPoly=cPoly+a(end)*bProd;
for curPow=1:orderA
    %Get the polynomial for the next power of b. The convolution operation
    %is equivalent to polynomial multiplication. The extra leading terms
    %(which are all zeros) are just thrown out.
    bProd=conv(b,bProd);
    bProd=bProd((prodOrder+1):end);
    
    %Add in the next power term.
    cPoly=cPoly+a(end-curPow)*bProd;
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
