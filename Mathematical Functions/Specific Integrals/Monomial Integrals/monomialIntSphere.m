function intVal=monomialIntSphere(alpha,wExp)
%%MONOMIALINTSPHERE Evaluate the integral of w(x)*prod_{i=1}^Nx(i)^alpha(i)
%           where w(x)=sum(x.^2)^(wExp/2) over the region defined by
%           sum(x.^2)<=1. That is, this function gives the value of the
%           desired monomial integral times w(x) taken over an n-
%           dimensional hypersphere. The region is designated as S_n in
%           [1].
%
%INPUTS: alpha An NX1 or 1XN vector of the integer exponents of the
%              monomial term. All elements must be >=0.
%         wExp The exponent term in the weighting function w(x) described
%              above. This must be an integer >-N. The default if this
%              parameter is omitted or an empty matrix is passed is 0,
%              which means that w(x)=1.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formula is taken from Chapter 7.4 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% [xi,w]=seventhOrderSpherCubPoints(3);
% alpha=[0;2;2];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntSphere(alpha)
%One will see that intVal and theMoment are both about 0.1197.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(wExp))
    wExp=0; 
end

if(all(mod(alpha,2)==0))
    alphaSum=sum(alpha);
    n=length(alpha);
    
    intVal=(2/(n+wExp+alphaSum))*exp(sum(gammaln((alpha+1)/2))-gammaln((n+alphaSum)/2));
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
