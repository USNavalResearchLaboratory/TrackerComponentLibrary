function [A,B,h,k,phi]=polyEllipsForm2ScalRotTrans(coeffs)
%%POLYELLIPSFORM2SCALROTTRANS Given a 2D ellipse represented as a quadratic
%   polynomial such that f(x,y)=0 is points on the ellipse, this function
%   obtains the parameterization for a parameteric representation of the
%   points in terms of scale factors (A, and B), a rotation angle (phi) and
%   translation coordinates (h, and k). Thus, if t is a parameteric
%   parameter that runs from 0 to 2*pi, x and y coordinates of points on
%   the ellipse are given by
%   x=A*cos(phi)*cos(t)-B*sin(phi)*sin(t)+h;
%   y=A*sin(phi)*cos(t)+B*cos(phi)*sin(t)+k;
%
%INPUTS: coeffs A 3X3 matrix of the coefficients for the bivariate
%               quadratic polynomial. These are arranged such
%               that coeffs(a1,a2) corresponds to the coefficient
%               of an x1^(a1-1)*x2^(a2-1) term. This is the format used in
%               polyValMultiDim.
%
%OUTPUTS: A, B The scalar positive scale factors.
%         h, k The scalar offsets in the x and y coordinates.
%          phi The scalar rotation angle in radians.
%
%The relations here come from inverting Equation 5 in [1]. Using the
%notation of [1], we note that we can write BB=(1/A2-1/B2)*sin(2*phi)
%and also AA-CC==(1/A^2-1/B^2)*cos(2*phi), so we use an inverse transgent to
%find phi. After that, the other terms are found by substituting and
%solving each of the equations. In the final equation for FF, we repalce
%the 1 with a c that we solve for. That is a scale factor (since one can
%multiply coeffs by any constant without changing the result). Given the
%scale factor, we then replace A with A/sqrt(c) and B with B/sqrt(c) to
%undo the scaling.
%
%EXAMPLE:
%To show that this produces correct results, we plot an ellipse using the
%drawEllipse function, which uses a different parameterization. We then
%convert that parameterization into polynomial coefficients using
%quadEllipsForm2Poly. Those coefficients are then transformed into the
%parameters for the scaled, rotated, and translated parametric form of
%this function. Being in parametric form, we plot the ellipse over the
%other one. One can see that they are the same.
% A=inv([4, -1.5;
%       -1.5,1]);
% x0=[2;3]*10;
% c=2.5;
% %Plot the original ellipse
% figure(1)
% clf
% drawEllipse(x0,A,c,'-k','linewidth',4)
% hold on
% %Convert the ellipse into polynomial form.
% coeffs=quadEllipsForm2Poly(A,x0,c)/1000;
% %Convert the polynomial form into the scaled, rotated, translated form.
% [A,B,h,k,phi]=polyEllipsForm2ScalRotTrans(coeffs);
% %This is a parameteric form, so plot the points on the ellipse:
% numPoints=500;
% t=linspace(0,2*pi,numPoints);
% x=A*cos(phi)*cos(t)-B*sin(phi)*sin(t)+h;
% y=A*sin(phi)*cos(t)+B*cos(phi)*sin(t)+k;
% plot(x,y,'-y','linewidth',1)
% legend('Original Ellipse','Parametric Ellipse')
% xlabel('x')
% ylabel('y')
%
%REFERENCES:
%[1] G. B. Hughes and M. Chraibi, "Calculating ellipse overlap areas,"
%    Computing and Visualization in Science, vol. 15, no. 5, pp. 291-301,
%    Oct. 2012.
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

AA=coeffs(3,1);
BB=coeffs(2,2);
CC=coeffs(1,3);
DD=coeffs(2,1);
EE=coeffs(1,2);
FF=coeffs(1,1);

phi=atan2(BB,AA-CC)/2;
sinPhi=sin(phi);
cosPhi=cos(phi);
%Using a double angle identity:
cos2Phi=cosPhi^2-sinPhi^2;

A2=2*cos2Phi/(AA-CC+(AA+CC)*cos2Phi);
B2=2*cos2Phi/(CC-AA+(AA+CC)*cos2Phi);

h=(1/2)*(-A2*DD+(B2-A2)*cosPhi*EE*sinPhi+(A2-B2)*DD*sinPhi^2);
k=(1/2)*(-B2*EE+(B2-A2)*cosPhi*DD*sinPhi-(A2-B2)*EE*sinPhi^2);

c=(h*cosPhi+k*sinPhi)^2/A2+(h*sinPhi-k*cosPhi)^2/B2-FF;

%Adjust the scaling of everything.
A2=A2*c;
B2=B2*c;

A=sqrt(A2);
B=sqrt(B2);

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
