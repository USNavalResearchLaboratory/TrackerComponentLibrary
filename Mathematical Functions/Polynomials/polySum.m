function pSum=polySum(p1,p2)
%%POLYSUM Add together two vectors representing 1D polynomials as power
%         series. The format of the polynomial vectors is the same as that
%         used in Matlab's polyval function: The coefficients are ordered
%         y(z)=p(end)+p(end-1)*z+p(end-2)*z^2+...
%
%INPUTS: p1,p2 Coefficients for the two power series. The number of terms
%              in the different series can be different. Thus p1 and p2 do
%              not have the same length. p1(1) is the coefficient of the
%              highest order term. p(end) is the additive constant term. p1
%              and p2 can be row or column vectors.
%
%OUTPUTS: pSum A column vector holding the sum of p1 and p2. The first
%              element corresponds to the highest term, as used in Matlab's
%              polyval function. Note that the vector is not shortened when
%              the sum of high order terms equals zero (terms cancel).
%
%Polynomials are summed by adding corresponding terms.
%
%As an example, consider adding 
%x^2-3*x+2
%to
%x^5+3*x^3-4*x^2+1
%This is done as
% p1=[1;-3;2];
% p2=[1;0;3;-4;0;1];
% pSum=polySum(p1,p2)
%One will see that the result is pSum'=[1,0,3,-3,-3,3], which corresponds
%to x^5+3*x^3-3*x^2-3*x+3.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

pSum=[zeros(length(p1)-length(p2),1);p2(:)]+[zeros(length(p2)-length(p1),1);p1(:)];

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
