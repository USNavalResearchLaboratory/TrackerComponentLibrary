function [xMin,yMin]=quadraticInterpUnif3PtMinMax(x0,deltaX,y0,y1,y2)
%%QUADRATICINTERPUNIF3PTMIN Given a function evaluated a 3 uniformly spaced
%           points, approximate its stationary point (minimum or maximum
%           depending on the orientation), using quadratic interpolation.
%
%INPUTS: x0 The point at which y0 is evaluated.
%    deltaX The spacing between the points on the grid.
%  y0, y1, y2 The values of the function to be interpolated evaluated at
%           x0, x0+deltaX and x0+2*deltaX.
%
%OUTPUTS: xMin The independent variable value to get the stationary point
%              (which is the finite minimum or maximum).
%         yMin The interpolated value of the function at the stationary
%              point.
%
%A quadratic equation is defined by three points. The interpolation formula
%herein is just derived by solving for the coefficients given three
%equations. The stationary point was derived by setting the derivative
%equal to zero and solving for x. We note that if y1 is the stationary
%point, then y0 and y2 will be symmetric about it and a 0/0 can occur in
%the implemented formula. Thus if a non-finite number arises, the assumed
%correct value of x is inserted.
%
%EXAMPLE:
%A quadratic equation is plotted. Three points are selected and the minimum
%is found and plotted and one can see that it coincides with the minimum of
%the quadratic equation.
% a=4;
% b=-6;
% c=-3;
% f=@(x)a*x.^2+b*x+c;
% 
% numPlotPts=500;
% xPlot=linspace(-10,10,numPlotPts);
% figure(1)
% clf
% hold on
% plot(xPlot,f(xPlot),'-b','linewidth',2)
% deltaX=2;
% x0=-2;
% x1=x0+deltaX;
% x2=x1+deltaX;
% y0=f(x0);
% y1=f(x1);
% y2=f(x2);
% [xMin,yMin]=quadraticInterpUnif3PtMinMax(x0,deltaX,y0,y1,y2);
% scatter(xMin,yMin,400,'.k')
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

deltaF0=y1-y0;
deltaF1=y2-y1;
d2f0=deltaF1-deltaF0;

xMin=(1/2)-deltaF0/d2f0;
if(~isfinite(xMin))
    %If symmetric about the minimum or all the values are the same.
    xMin=1/2;
end

if(nargout>1)
    yMin=y0+xMin.*deltaF0+(1/2).*xMin.*(xMin-1).*d2f0;
end

%Transform to the global coordinate system.
xMin=xMin*deltaX+x0;

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
