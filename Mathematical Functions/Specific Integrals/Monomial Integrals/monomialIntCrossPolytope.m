function intVal=monomialIntCrossPolytope(alpha)
%%MONOMIALINTCROSSPOLYTOPE Evaluate the integral of
%          prod_{i=1}^N x(i)^alpha(i) over the space defined by
%          sum(abs(x))=1. That is, this function gives the value of the
%          desired monomial integral taken over an N-dimensional cross
%          polytope. The region is designated as G_n in [1]. In two
%          dimensions, this shape is a square, in three an octohedron, in
%          N, a cross polytope (or N-dimensional hyperoctahedron).
%
%INPUTS: alpha An NX1 or 1XN vector of the integer exponents of the
%              monomial term. All elements must be >=0.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formula is taken from Chapter 7.7 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% [xi,w]=fifthOrderCrossPolyCubPoints(3);
% alpha=[0;2;2];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntCrossPolytope(alpha)
%One will find that theMoment and intVal are both near 0.0063.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(all(mod(alpha,2)==0))
    n=length(alpha);

    intVal=exp(n*log(2)+sum(gammaln(alpha+1))-gammaln(n+sum(alpha)+1));
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
