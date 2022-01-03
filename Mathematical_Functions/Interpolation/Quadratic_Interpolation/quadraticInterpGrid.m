function yInterp=quadraticInterpGrid(xDes,xPts,yPts)
%%QUADRATICINTERPGRID Given a set of points and a set of function values
%       evaluated at those points, interpolate the function values at a set
%       of desired points using quadratic interpolation. This function
%       selects the nearest set of basis points for the interpolation of
%       each point given the values in xPts and yPts. The points in xPts do
%       not have to be uniformly spaced.
%
%INPUTS: xDes A vector or matrix of real values of the independent variable
%             at which interpolation should be performed.
%        xPts A numXX1 or 1XnumX vector holding points at which the
%             reference function values are evaluated. These must be given
%             in increasing order. numX>=3.
%        yPts A numXX1 or 1XnumX vector holding the values of the real
%             function to be interpolated evaluated at the points in xPts.
%
%OUTPUTS: yInterp The interpolated values. This matrix has the same
%                 dimensions as xDes.
%
%This function finds brackets the region .for each points using binSearch
%and then performs interpolation using quadraticInterp3Pt.
%
%EXAMPLE:
%A sine wave is plotted in blue. Then, a small set of non-uniformly spaced
%points for interpolation are chosen. The line is then interpolated between
%the points and plotted in red. One can see the discontinuities when the
%basis points for the quadratic interpolation change. Increasing the number
%of basis points can decrease the discontinuity.
% numPtsPlot=4000;
% xPlot=linspace(-4,4,numPtsPlot);
% yPlot=sin(xPlot);
% 
% %Non-uniformly spaced basis points.
% xBasis=[-3.6031, -3.1943, -3.1040, -2.7751, -1.4236, -1.2073, -0.3550, -0.1226, 0.0287, 1.2500, 1.8025, 1.9833, 2.9627, 3.1371, 4.0091];
% yBasis=sin(xBasis);
% 
% figure(1)
% clf
% hold on
% plot(xPlot,yPlot,'-b','linewidth',4)
% scatter(xBasis,yBasis,400,'.k')
% yInterp=quadraticInterpGrid(xPlot,xBasis,yBasis);
% plot(xPlot,yInterp,'-r','linewidth',2)
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPts2Interp=length(xDes);
yInterp=zeros(numPts2Interp,1);

numXGrid=length(xPts);
for curPt=1:numPts2Interp
    [~,idx]=binSearch(xPts,xDes(curPt),0);

    %The points to use as the basis of the interpolation are one before
    %and one after the index of the closes point, unless the closest point
    %is an endpoint.
    idx=idx+(idx==1)-(idx==numXGrid);
    
    x1=xPts(idx-1);
    x2=xPts(idx);
    x3=xPts(idx+1);
    y1=yPts(idx-1);
    y2=yPts(idx);
    y3=yPts(idx+1);

    yInterp(curPt)=quadraticInterp3Pt(xDes(curPt),x1,x2,x3,y1,y2,y3);
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
