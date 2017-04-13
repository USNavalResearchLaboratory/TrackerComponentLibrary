function intVal=monomialIntHexagon(alpha)
%%MONOMIALINTHEXAGON Evaluate the integral of prod_{i=1}^2 x(i)^alpha(i),
%             which is a monomial term, taken over a two-dimensional
%             hexagon with vertices (+/-1,0) and
%             (+/-(1/2),+/-(1/2)*sqrt(3)). The region is designated as H_2
%             in [1].
%
%INPUTS: alpha An 3X1 or 1X3 vector of the integer exponents of the
%              monomial term. All elements must be >=0. THis function is
%              only implemented for sum(alpha)<=6.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formulae are taken from Chapter 7.11 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% [xi,w]=fifthOrderHexagonCubPoints();
% alpha=[0;4];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntHexagon(alpha)
%One will see that theMoment and intval are both about 0.2273.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

V=3*sqrt(3)/2;

alpha=alpha(:);

if(any(mod(alpha,2)~=0))
    intVal=0;
elseif(length(alpha)~=2)
    error('Unsupported alpha vector provided.')
else
    if(all(alpha==[0;0]))
        intVal=V;
    elseif(all(alpha==[2;0])||all(alpha==[0;2]))
        intVal=(5/24)*V;
    elseif(all(alpha==[4;0])||all(alpha==[0;4]))
        intVal=(7/80)*V;
    elseif(all(alpha==[2;2]))
        intVal=(7/240)*V;
    elseif(all(alpha==[6;0]))
        intVal=(85/1792)*V;
    elseif(all(alpha==[0;6]))
        intVal=(81/1792)*V;
    elseif(all(alpha==[4;2]))
        intVal=(73/8960)*V;
    elseif(all(alpha==[2;4]))
        intVal=(93/8960)*V;
    else
       error('Unsupported alpha vector provided.')
    end
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
