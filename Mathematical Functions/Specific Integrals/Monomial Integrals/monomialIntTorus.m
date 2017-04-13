function intVal=monomialIntTorus(alpha,rho1,rho2)
%%MONOMIALINTTORUS Evaluate the integral of prod_{i=1}^N x(i)^alpha(i),
%             which is a monomial term, taken over a three-dimensional
%             torus with a circular cross section. The region is specified
%             as (sqrt(x^2+y^2)-rho1)^2+z^2<=rho2^2 with 0<rho2<=rho1<Inf.
%             The region is designated as Tor_3 in [1].
%
%INPUTS: alpha A 3X1 or 1X3 vector of the integer exponents of the monomial
%              term. All elements must be >=0. THis function is only
%              implemented for sum(alpha)<=4.
%    rho1,rho2 The two defining parameters of the torus as given above with
%              0<rho2<=rho1<Inf. rho2 is the radius of the filled-in
%              portion of the ring and rho1 is the distance of the center
%              of the ring (the hole) from the center of the filled-in
%              portion of the ring
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formulae are taken from Chapter 7.18 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% rho1=1;
% rho2=1/3;
% [xi,w]=thirdOrderTorusCubPoints(rho1,rho2);
% alpha=[0;2;0];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntTorus(alpha,rho1,rho2)
%One will see that theMoment and intVal both are about 1.1880.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(any(mod(alpha,2)~=0))
    intVal=0;
elseif(any(alpha>4)||length(alpha)~=3||any(alpha<0))
    error('Unsupported alpha vector provided.')
else
    alpha=alpha(:);
    
    if(all(alpha==[0;0;0]))
        intVal=2*pi^2*rho1*rho2^2;
    elseif(all(alpha==[2;0;0])||all(alpha==[0;2;0]))
        intVal=pi^2*rho1*rho2^2*(rho1^2+(3/4)*rho2^2);
    elseif(all(alpha==[0;0;2]))
        intVal=(1/2)*pi^2*rho1*rho2^4;
    elseif(all(alpha==[4;0;0])||all(alpha==[2;2;0])||all(alpha==[0;4;0]))
        intVal=(3/4)*pi^2*rho1*rho2^2*(rho1^4+(5/2)*rho1^2*rho2^2+(5/8)*rho2^4);
        if(alpha(1)==2)
            intVal=intVal/3;
        end
    elseif(all(alpha==[0;0;4]))
        intVal=(1/4)*pi^2*rho1*rho2^6;
    elseif(all(alpha==[2;0;2])||all(alpha==[0;2;2]))
        intVal=(1/4)*pi^2*rho1*rho2^4*(rho1^2+(1/2)*rho2^2);
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
