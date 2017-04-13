function intVal=monomialIntSimplex(alpha)
%%MONOMIALINTSIMPLEX Evaluate  the integral of prod_{i=1}^N x(i)^alpha(i)
%           over the region defined by sum(x)<=1 and x>=0. That is, this
%           function gives the value of the desired monomial integral taken
%           over the n-dimensional simplex. The region is designated as T_n
%           in [1].
%
%INPUTS: alpha An NX1 or 1XN vector of the integer exponents of the
%              monomial term. All elements must be >=0.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formula is taken from Chapter 7.8 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% [xi,w]=fifthOrderSimplexCubPoints(3);
% alpha=[0;2;2];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntSimplex(alpha)
%One will see that the results are both near 7.9365e-04.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(alpha);

intVal=exp(sum(gammaln(alpha+1))-gammaln(n+sum(alpha)+1));

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
