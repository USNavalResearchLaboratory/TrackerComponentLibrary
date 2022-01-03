function interpVal=quadraticInterp3Pt(x,x1,x2,x3,y1,y2,y3)
%%QUADRATICINTERP3PT Given a function evaluated at three specified points,
%               interpolate the function at the points specified in xDes.
%
%INPUTS: xDes A matrix of points of the independent variable at which the
%             interpolation should be performed.
%    x1,x2,x3 The independent variables values to which y1, y2, and y3
%             correspond.
%  y1, y2, y3 The values of the function to be interpolated.
%
%OUTPUTS: interpVal The interpolated values. This matrix has the same
%                   dimensions as xDes.
%
%A quadratic equation is defined by three points. The interpolation formula
%herein is just derived by solving for the coefficients given three
%equations.
%
%EXAMPLE:
%In this instance, we show that the interpolated value returned by the
%function matches the true function when when the points provided are
%actually from a quadratic polynomial.
% a=4;
% b=-6;
% c=-3;
% f=@(x)a*x.^2+b*x+c;
% 
% x1=-2;
% x2=0;
% x3=4;
% f1=f(x1);
% f2=f(x2);
% f3=f(x3);
% xDes=1.5;
% %One will see that the interpolate value and the true value are the same.
% yDesInterp=quadraticInterp3Pt(xDes,x1,x2,x3,f1,f2,f3)
% yDesTrue=f(xDes)
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

delta12=x1-x2;
delta13=x1-x3;
delta23=x2-x3;

deltaX1=x-x1;
deltaX2=x-x2;
deltaX3=x-x3;

interpVal=y1.*deltaX2.*deltaX3./(delta12.*delta13)-y2.*deltaX3.*deltaX1./(delta23.*delta12)+y3.*deltaX1.*deltaX2./(delta13.*delta23);

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
