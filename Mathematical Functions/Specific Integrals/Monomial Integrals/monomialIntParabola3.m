function intVal=monomialIntParabola3(alpha,a,b)
%%MONOMIALINTPARABOLA3 Evaluate the integral of prod_{i=1}^2 x(i)^alpha(i)
%               which is a monomial term, evaluated between the two
%               parabolas, one with equation y=b-b*x^2/a^2 and the other
%               with equation y=b*x^2/a^2-b. The region is designated as
%               Par^3 in [1].
%
%INPUTS: alpha An 2X1 or 1X2 vector of the integer exponents of the
%              monomial term. All elements must be >=0. This function is
%              only implemented for sum(alpha)<=4.
%          a,b The two real parameters of the two parabolas, whose
%              equations are given above.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formulae are taken from Chapter 7.15 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% a=1;
% b=2;
% [xi,w]=fifthOrderCubPointsParabola3(a,b);
% alpha=[0;4];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntParabola3(alpha,a,b)
%One will see that intVal and the Moment are abouth about 9.4569.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(length(alpha)~=2)
    error('Unsupported alpha vector provided.')
end

V=8*a*b/3;

alpha=alpha(:);

if(any(mod(alpha,2)~=0))
    intVal=0;
elseif(all(alpha==[0;0]))
    intVal=V;
elseif(all(alpha==[2;0]))
    intVal=(a^2/5)*V;
elseif(all(alpha==[0;2]))
    intVal=(8*b^2/35)*V;
elseif(all(alpha==[4;0]))
    intVal=(3*a^4/35)*V;
elseif(all(alpha==[2;2]))
    intVal=(8*a^2*b^2/315)*V;
elseif(all(alpha==[0;4]))
    intVal=(128*b^4/1155)*V;
else
    error('Unsupported alpha vector provided.')
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
