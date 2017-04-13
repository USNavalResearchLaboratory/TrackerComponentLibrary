function intVal=monomialIntCubicalShell(alpha,K1,K2)
%%MONOMIALINTCUBALSHELL Evaluate the integral of prod_{i=1}^N x(i)^alpha(i)
%              taken over the region K1<=x(i)<=K2 with 0<K1<K2<Inf. That
%              is, this function gives the value of the desired monomial
%              integral over an N-dimensional hypercubic shell. The region
%              is designated as C_N^{shell} in [1].
%
%INPUTS: alpha An NX1 or 1XN vector of the integer exponents of the
%              monomial term. All elements must be >=0.
%       K1, K2 Two positive, finite, real values such that 0<K1<K2.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formula is taken from Chapter 7.3 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% K1=0.5;
% K2=1.25;
% [xi,w]=thirdOrderNDimShellCubPoints(3,K1,K2);
% alpha=[0;2;0];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntCubicalShell(alpha,K1,K2)
%One will find that both theMoment and intVal are about 8.0547.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(all(mod(alpha,2)==0))
    alphaSum=sum(alpha);
    n=length(alpha);
    
    intVal=(2^n/prod(alpha+1))*(K2^(n+alphaSum)-K1^(n+alphaSum));
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
