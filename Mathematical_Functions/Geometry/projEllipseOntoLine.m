function [sMin,sMax,s0]=projEllipseOntoLine(z,A,gammaVal,x0,v)
%%PROJELLIPSEONTOLINE Determine the orthogonal projection of an ellipse or
%                     ellipsoid onto a line.
%
%INPUTS: z The numDimX1 center of the ellipsoid.
%        A A numDimXnumDim symmetric, positive definite matrix that
%          specifies the size and shape of the ellipse or ellipsoid, where
%          a point zp is on the ellipse/ellipsoid if
%          (zp-z)'*A*(zp-z)=gammaVal.
% gammaVal The threshold for declaring a point to be in the ellipsoid. If
%          an empty matrix is passed, the default value of 1 is used.
%     x0,v Two numDimX1 vectors that specify the line. A point y is on the
%          line if y=x0+s*v, where s is a scalar quantity.
%
%OUTPUTS: sMin,sMax The range of values of s on the line y=x0+s*v onto
%                   which the ellipsoid is orthogonally projected. These
%                   are scalar quantities.
%
%The mathematics for projecting an ellipsoid onto a line are described in
%Section 12 of [1].
%
%EXAMPLE:
% A=[1,0;
%    0,10];
% M=[0.413074198133900,  0.910697373904216;
%    -0.910697373904216,   0.413074198133900];%A rotation matrix
% A=M*A*M';
% z=[5;6];
% x0=[0;0];
% v=[1;2];
% [sMin,sMax,s0]=projEllipseOntoLine(z,A,[],x0,v);
% figure(1)
% clf
% hold on
% axis([-0.5, 7.5, 2, 10])%Same size in x and y
% axis square
% drawEllipse(z,A,1,'b','linewidth',2)
% %Draw a segment of the line.
% xStart=x0+v;
% xEnd=x0+5*v;
% plot([xStart(1);xEnd(1)],[xStart(2);xEnd(2)],'-k','linewidth',4)
% %Draw the segment of the orthogonal projection onto the line.
% xStart=x0+sMin*v;
% xEnd=x0+sMax*v;
% plot([xStart(1);xEnd(1)],[xStart(2);xEnd(2)],'-r','linewidth',2)
% %Draw the line from the center of the ellipse to the place on the line onto
% %which it is orthogonally projected.
% xS0=x0+s0*v;
% plot([z(1);xS0(1)],[z(2);xS0(2)],'--c')
%
%REFERENCES:
%[1] S. B. Pope, "Algorithms for ellipsoids," Cornell University, Tech.
%    Rep. FDA-08-01, Feb. 2008. [Online].
%    Available: https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(gammaVal))
   gammaVal=1; 
end

A=A/gammaVal;

L=chol(A,'lower');

s0=v'*(z-x0)/(v'*v);
w=norm(L\v/(v'*v));

sMin=s0+w;
sMax=s0-w;
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
