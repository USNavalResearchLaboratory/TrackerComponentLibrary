function interpVal=quadraticInterpUnif3Pt(xDes,x0,deltaX,y0,y1,y2)
%%QUADRATICINTERPUNIF3PT Given a function evaluated at three uniformly
%               spaced points, interpolate the function at the points
%               specified in xDes.
%
%INPUTS: xDes A matrix of points of the independent variable at which the
%             interpolation should be performed.
%          x0 The point at which y0 is evaluated.
%      deltaX The spacing between the points on the grid.
%  y0, y1, y2 The values of the function to be interpolated evaluated at
%             x0, x0+deltaX and x0+2*deltaX.
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
% deltaX=2;
% x0=-2;
% x1=x0+deltaX;
% x2=x1+deltaX;
% f0=f(x0);
% f1=f(x1);
% f2=f(x2);
% xDes=1.5;
% %One will see that the interpolate value and the true value are the same.
% yDesInterp=quadraticInterpUnif3Pt(xDes,x0,deltaX,f0,f1,f2)
% yDesTrue=f(xDes)
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

x=(xDes-x0)./deltaX;

deltaF0=y1-y0;
deltaF1=y2-y1;
d2f0=deltaF1-deltaF0;

interpVal=y0+x.*deltaF0+(1/2).*x.*(x-1).*d2f0;

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
