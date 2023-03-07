function coeffs=scalRotTransEllipsForm2Poly(A,B,h,k,phi)
%%SCALEROTTRANSEELIPSFORM2POLY Given the parameters of an ellipse in 2D
%   represented in parameteric form as
%   x=A*cos(phi)*cos(t)-B*sin(phi)*sin(t)+h;
%   y=A*sin(phi)*cos(t)+B*cos(phi)*sin(t)+k;
%   where t is the parameteric parameter, which can range from 0 to 2*pi,
%   obtain a set of coefficients for a bivariante quadratic polynomial such
%   that points on the surface of the ellipse are given by f(y,y)=0. The
%   coefficients are returned in a format that can be used by
%   polyValMultiDim.
%
%INPUTS: A, B The scalar positive scale factors.
%        h, k The scalar offsets in the x and y coordinates.
%         phi The scalar rotation angle in radians.
%
%OUTPUTS: coeffs A 3X3 matrix of the coefficients for the bivariate
%               quadratic polynomial. These are arranged such
%               that coeffs(a1,a2) corresponds to the coefficient
%               of an x1^(a1-1)*x2^(a2-1) term. This is the format used in
%               polyValMultiDim.
%
%This function implements the expressions given in Equation 5 of [1]. The
%
%EXAMPLE:
%This example demonstrates that this function is consistent with
%polyEllipsForm2ScalRotTrans: Given an ellipse parameterized by A,B,h,k and
%phi, we obtain polynomial coefficients. Those coefficients are then
%converted back into this form using polyEllipsForm2ScalRotTrans. The
%converted-back values are (within finite precision limits) the same as the
%original value... except A and B are switched, but phi is flipped by pi/2,
%so the results are actually the same. The original and converted back
%values are displayed. Then, the two ellipses are parametrically plotted
%so that one can verify that the results are the same. 
% A=2;
% B=1/10;
% h=-2;
% k=3;
% phi=0.1268;
% coeffs=scalRotTransEllipsForm2Poly(A,B,h,k,phi);
% [A1,B1,h1,k1,phi1]=polyEllipsForm2ScalRotTrans(coeffs);
% v=[A;B;h;k;phi];
% v1=[A1;B1;h1;k1;phi1];
% %Display the two sets of parameters side by side.
% [v,v1]
% %Now, obtain x and y parameterically for each set of parameters. These x
% %and y values are then plotted.
% numPoints=500;
% t=linspace(0,2*pi,numPoints);
% x=A*cos(phi)*cos(t)-B*sin(phi)*sin(t)+h;
% y=A*sin(phi)*cos(t)+B*cos(phi)*sin(t)+k;
% x1=A1*cos(phi1)*cos(t)-B1*sin(phi1)*sin(t)+h1;
% y1=A1*sin(phi1)*cos(t)+B1*cos(phi1)*sin(t)+k1;
% figure(1)
% clf
% hold on
% plot(x,y,'linewidth',6)
% plot(x1,y1,'linewidth',2)
%
%REFERENCES:
%[1] G. B. Hughes and M. Chraibi, "Calculating ellipse overlap areas,"
%    Computing and Visualization in Science, vol. 15, no. 5, pp. 291-301,
%    Oct. 2012.
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

sinPhi=sin(phi);
cosPhi=cos(phi);
A2=A^2;
B2=B^2;

%Equation 5a
AA=cosPhi^2/A2+sinPhi^2/B2;
%Equation 5b
BB=2*sinPhi*cosPhi*(1/A2-1/B2);
%Equation 5c
CC=sinPhi^2/A2+cosPhi^2/B2;
%Equation 5d
DD=-2*cosPhi*(h*cosPhi+k*sinPhi)/A2+2*sinPhi*(k*cosPhi-h*sinPhi)/B2;
%Equation 5e
EE=-2*sinPhi*(h*cosPhi+k*sinPhi)/A2+2*cosPhi*(h*sinPhi-k*cosPhi)/B2;
%Equation 5f, modified with c.
FF=(h*cosPhi+k*sinPhi)^2/A2+(h*sinPhi-k*cosPhi)^2/B2-1;

coeffs(3,1)=AA;
coeffs(2,2)=BB;
coeffs(1,3)=CC;
coeffs(2,1)=DD;
coeffs(1,2)=EE;
coeffs(1,1)=FF;

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
