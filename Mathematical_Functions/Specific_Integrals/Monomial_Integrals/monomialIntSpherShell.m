function intVal=monomialIntSpherShell(alpha,rho)
%%MONOMIALINTSPHERSHELL Evaluate the integral of prod_{i=1}^N x(i)^alpha(i)
%           over the region defined by rho<=sum(x.^2)<=1. That is, this
%           function gives the value of the desired monomial integral taken
%           over an n-dimensional hyperspherical shell. The region is
%           designated as S_n^{shell} in 1.
%
%INPUTS: alpha An NX1 or 1XN vector of the integer exponents of the
%              monomial term. All elements must be >=0.
%          rho The lower bound of the shell region, as described above,
%              0<rho<1.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formula is taken from Chapter 7.5 of [1]. These types of explicit
%moment formulae are useful for testing cubature integration points.
%
%%EXAMPLE:
%Here we verify that the moment value produced here equals that obtained
%using cubature points of an appropriate order.
% rMin=0.75;
% [xi,w]=seventhOrderSpherShellCubPoints(3,rMin);
% alpha=[0;2;4];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntSpherShell(alpha,rMin)
%One will see that theMoment and intVal are both about 0.0369.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(alpha);
intVal=monomialIntSphere(alpha,0)*(1-rho^(N+sum(alpha)));
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
