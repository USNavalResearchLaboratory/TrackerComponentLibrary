function intVal=monomialWeightedEllipse(alpha,a,b)
%%MONOMIALWEIGHTEDELLIPSE Evaluate the integral of
%        w(x)*prod_{i=1}^N x(i)^alpha(i) taken over a 2D ellipse with
%        equation x^2/a^2+y^2/b^2=1 with a>b. The weighting function is
%        w(x)=sqrt((x-c)^2+y^2)*sqrt((x+c)^2+y^2) where c=sqrt(a^2-b^2).
%        The region is designated as Elp in [1].
%
%INPUTS: alpha An 2X1 or 1X2 vector of the integer exponents of the
%              monomial term. All elements must be >=0. This is only
%              implemented for the case where sum(alpha)<=6.
%          a,b The parameters for the ellipse such that Inf>a>=b>0.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formulae are taken from Chapter 7.12 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% a=5;
% b=3;
% [xi,w]=seventhOrderWeightEllipsCubPoints(a,b);
% alpha=[4;2];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialWeightedEllipse(alpha,a,b)
%One will find that theMoment and intVal are both around 351.4870.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(any(mod(alpha,2)~=0))
    intVal=0;
elseif(any(alpha>6)||length(alpha)~=2||any(alpha<0))
    error('Unsupported alpha vector provided.')
else
    alpha=alpha(:);
    
    c=sqrt(a^2-b^2);
    
    V=2*pi*log((a+b)/c);
    
    if(all(alpha==[0;0]))
        intVal=V;
    elseif(all(alpha==[2;0]))
        intVal=(pi/2)*(a*b+(1/(2*pi))*c^2*V);
    elseif(all(alpha==[0;2]))
        intVal=(pi/2)*(a*b-(1/(2*pi))*c^2*V);
    elseif(all(alpha==[4;0]))
        intVal=(3*pi/8)*((1/2)*a*b^3+(5/4)*c^2*a*b+(3/(8*pi))*c^4*V);
    elseif(all(alpha==[2;2]))
        intVal=(pi/8)*((1/2)*a*b^3+(1/4)*c^2*a*b-(1/(8*pi))*c^4*V);
    elseif(all(alpha==[0;4]))
        intVal=(3*pi/8)*((1/2)*a*b^3-(3/4)*c^2*a*b+(3/(8*pi))*c^4*V);
    elseif(all(alpha==[6;0]))
        intVal=(15*pi/48)*((1/3)*a*b^5+(13/12)*c^2*a*b^3+(11/8)*c^4*a*b+(5/(16*pi))*c^6*V);
    elseif(all(alpha==[4;2]))
        intVal=(3*pi/48)*((1/3)*a*b^5+(7/12)*c^2*a*b^3+(1/8)*c^4*a*b-(1/(16*pi))*c^6*V);
    elseif(all(alpha==[2;4]))
        intVal=(3*pi/48)*((1/3)*a*b^5+(1/12)*c^2*a*b^3-(1/8)*c^4*a*b+(1/(16*pi))*c^6*V);
    elseif(all(alpha==[0;6]))
        intVal=(15*pi/48)*((1/3)*a*b^5-(5/12)*c^2*a*b^3+(5/8)*c^4*a*b-(5/(16*pi))*c^6*V);
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
