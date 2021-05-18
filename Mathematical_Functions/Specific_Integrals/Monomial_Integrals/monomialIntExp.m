function intVal=monomialIntExp(alpha)
%%MONOMIALINTEXP Evaluate the integral of w(x)*prod_{i=1}^N x(i)^alpha(i)
%       taken over all N-dimensional space with the weighting function
%       w(x)=exp(-norm(x)). That is, this function gives the value of the
%       desired monomial integral taken with an exponential  weighting
%       function. The region is designated as E_n^{r} in [1].
%
%INPUTS: alpha An NX1 or 1XN vector of the integer exponents of the
%              monomial term. All elements must be >=0.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formula is taken from Chapter 7.10 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% [xi,w]=seventhOrderExpCubPoints(3);
% alpha=[0;6;0];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntExp(alpha)
%One will find that intVal and theMoment are both about 7.2382e+04.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(all(mod(alpha,2)==0))
    n=length(alpha);
    sumVal=n+sum(alpha);

    intVal=2*exp(gammaln(sumVal)+sum(gammaln((alpha+1)/2))-gammaln(sumVal/2));
else
    intVal=0;
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
