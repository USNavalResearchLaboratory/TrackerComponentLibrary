function [xi,w]=seventhOrderSpherShellCubPoints(numDim,rMin,algorithm)
%%SEVENTHORDERSPHERSHELLCUBPOINTS Generate seventh-order cubature points
%               for integration over a shell of the unit (hyper-)sphere
%               (weighting function is just 1). The shell is defined as the
%               region scuh that rMin<=sum(x.^2)<=1.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. The limits for this depend on the algorithm
%               used.
%          rMin The radius of the inner hollow section of the shell.
%               0<=rMin<1.
%     algorithm A value indicating which algorithm should be used. This
%               function implements formula S_n^{shell} 7-1 in [1], pg.
%               293. This formula requires a set of cubature points for a
%               spherical surface. This algorithm parameter is thus the
%               same as that in seventhOrderSpherSurfCubPoints, which is
%               called for the necessary points.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%This implements formula S_n^{shell} 7-1 in [1], pg. 293, 2*numDim points.
%Chapter 2.8 of [1] also describes a method for generating arbitrary-order
%spherical shells.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=[]; 
end

n=numDim;
rho=rMin;

c0=(1-rho^(n+0))/(n+0);
c2=(1-rho^(n+2))/(n+2);
c4=(1-rho^(n+4))/(n+4);
c6=(1-rho^(n+6))/(n+6);

ab=[-c2,c0;
    -c4,c2]\[-c4;
             -c6];
a=ab(1);
b=ab(2);

%Find the r values.
r2=roots([1;-a;b]);
r2=real(r2);%Deal with any finite precision issues.

if(any(r2<0))
    error('Finite precision errors led to an incorrectly determined value')
end

r=sqrt(r2);

%Determine A1 and A2
A=[1,         1;
   r2(1), r2(2)]\[c0;
                  c2];

[u,B]=seventhOrderSpherSurfCubPoints(n,algorithm);

numB=length(B);

M=2*numB;
xi=zeros(numDim,M);
w=zeros(M,1);

w(1:numB)=B*A(1);
w((numB+1):M)=B*A(2);

xi(:,1:numB)=r(1)*u;
xi(:,(numB+1):M)=r(2)*u;

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
