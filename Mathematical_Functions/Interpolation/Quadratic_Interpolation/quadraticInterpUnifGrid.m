function yInterp=quadraticInterpUnifGrid(xDes,x0,deltaX,yPts)
%%QUADRATICINTERPUNIFGRID Given a set of function values evalauted at
%       uniformly spaced points, interpolate the function values at a set
%       of desired points using quadratic interpolation. This function
%       selects the nearest set of basis points for the interpolation of
%       each point given the values in xPts and yPts.
%
%INPUTS: xDes A vector or matrix of real values of the independent variable
%             at which interpolation should be performed.
%          x0 The value of the independent variable at which yPts(1) was
%             evaluated.
%      deltaX The spacing of the independent variable  between values of
%             the dependent variable in yPts. deltaX>0.
%        yPts A numBasisX1 or 1XnumBasis array of values of the real
%             function to be interpolated evaluated starting with argument
%             at x0 and continuing with increments of deltaX. numBasis>=3.
%
%OUTPUTS: yInterp The interpolated values. This matrix has the same
%                 dimensions as xDes.
%
%This function finds brackets the region .for each points using binSearch
%and then performs interpolation using quadraticInterpUnif3Pt.
%
%EXAMPLE:
%A sine wave is plotted in blue. Then, a small set of points for
%interpolation are chosen. The line is then interpolated between the points
%and plotted in red. One can see the discontinuities when the basis points
%for the quadratic interpolation change. Increasing the number of basis
%points can decrease the discontinuity.
% numPtsPlot=3000;
% numPtsBasis=15;
% xPlot=linspace(-4,4,numPtsPlot);
% yPlot=sin(xPlot);
% xBasis=linspace(-4,4,numPtsBasis);
% yBasis=sin(xBasis);
% 
% figure(1)
% clf
% hold on
% plot(xPlot,yPlot,'-b','linewidth',4)
% scatter(xBasis,yBasis,400,'.k')
% x0=xBasis(1);
% deltaX=xBasis(2)-xBasis(1);
% yInterp=quadraticInterpUnifGrid(xPlot,x0,deltaX,yBasis);
% plot(xPlot,yInterp,'-r','linewidth',2)
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPts2Interp=numel(xDes);
yInterp=zeros(size(xDes));

numY=length(yPts);

for curPt=1:numPts2Interp
    idx=round((xDes(curPt)-x0)/deltaX)+1;

    %The points to use as the basis of the interpolation are one before
    %and one after the index of the closes point, unless the closest point
    %is an endpoint.
    idx=idx+(idx==1)-(idx==numY);
    
    x0Cur=x0+deltaX*(idx-2);
    y0=yPts(idx-1);
    y1=yPts(idx);
    y2=yPts(idx+1);
    yInterp(curPt)=quadraticInterpUnif3Pt(xDes(curPt),x0Cur,deltaX,y0,y1,y2);
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
