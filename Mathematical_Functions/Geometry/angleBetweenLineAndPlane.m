function angVal=angleBetweenLineAndPlane(a,n)
%%ANGLEBETWEENLINEANDPLANE Given a line specified in parameteric form as
%       x=x0+a*t, where t is the scalar parametric parameter, and a plane
%       given in parameteric form as x'*n=d, return the (acute) angle
%       between the line and the plane.
%
%INPUTS: a The xDimX1 parameteric slope vector of the line.
%        n An xDimX1 normal to the plane. If this is omitted or an empty
%          matrix is passed, then n is a vector that is all zero except for
%          the last element, which is 1. In 3D that is equivalent to using
%          the x-y plane.
%
%OUTPUTS: angVal The angle between the line and the plane given in radians.
%
%This just uses the identity that
%sin(angVal)=abs(dot(a,n)/(norm(a)*norm(n)))
%
%EXAMPLE:
%An easy test for the function is to use the x-y plane and just have a line
%in the x-z plane defined at a specified angle to the x axis. In that
%instance, the result of this function can be quickly compared ot the
%truth.
% n=[0;0;1];
% theta=pi/3;
% a=[cos(theta);0;sin(theta)];
% AbsErr=abs(angleBetweenLineAndPlane(a,n)-theta)
%In this instance the absolute error is 0.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

angVal=asin(abs(dot(a,n)/(norm(a)*norm(n))));

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
