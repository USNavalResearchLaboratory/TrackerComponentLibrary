function zInterp=biquadraticInterp(xyPts,x0,y0,deltaX,deltaY,zGrid)
%%BIQUADRATICINTERP Biquadratic interpolation is interpolation of a
%       function that depends on two independent variables, by first
%       performing quadratic interpolation across one dimension, and then
%       performing quadratic interpolation across the next dimension. This
%       function assumes that in x and y generating the values in zGrid is
%       uniform.
%
%INPUTS: xyPts A 2XnumPts set of [x;y] pairs at which the value of z should
%              be interpolated.
%       x0, y0 These are values of the independent variables corresponding
%              to the entry in zGrid(1,1).
% deltaX, deltaY The uniform spacings between the independent variables for
%              gridded elements in zGrid. Thus, zGrid(i,j) holds the z
%              value associated with x=x0+deltaX*(i-1) and
%              y=y0+deltaY*(i-1).
%        zGrid A numXXnumY matrix of values of the independent variable
%              for different values of the dependent variables. The first
%              index selects the x value and the second index the y value.
%              numX>=3 and numY>=3.
%
%OUTPUTS: zInterp A 1XnumPts vector of interpolated values.
%
%Biquadratic interpolation is described in [1]. However, this function does
%not use the implementation that is given at the end of [1].
%
%EXAMPLE:
%A 2D function is evaluated and plotted. Then, the function is evaluated at
%a much smaller number of points and those are used to interpolate and plot
%the function again. The basis points are ploted in black. The
%small discontinuities between interpolation regions can be seen.
% numXPts=400;
% numYPts=401;
% plotXPts=linspace(-1,1,numXPts);
% plotYPts=linspace(-1,1,numYPts);
% z=@(xy)sum(sin(1-(16/15)*xy).^2-(1/50)*sin(4-(64/15)*xy)-sin(1-(16/15)*xy),1);
% [Y,X]=meshgrid(plotYPts,plotXPts);
% xy=[X(:).';Y(:).'];
% Z=reshape(z(xy),[numXPts,numYPts]);
% 
% figure(1)
% clf
% hold on
% surface(X,Y,Z,'edgeColor','None')
% 
% %Get a small number of points to use for interpolation and mark the
% %points on the plot.
% numXBasis=10;
% numYBasis=11;
% x0=-1;
% y0=-1;
% xBasis=linspace(x0,1,numXBasis);
% yBasis=linspace(y0,1,numYBasis);
% deltaX=xBasis(2)-xBasis(1);
% deltaY=yBasis(2)-yBasis(1);
% [YB,XB]=meshgrid(yBasis,xBasis);
% xyB=[XB(:).';YB(:).'];
% ZGrid=reshape(z(xyB),[numXBasis,numYBasis]);
% scatter3(XB(:),YB(:),ZGrid(:),200,'.k')
% view(25,50)
% 
% %Now, interpolate the surface using the above basis points.
% ZInterp=biquadraticInterp(xy,x0,y0,deltaX,deltaY,ZGrid);
% ZInterp=reshape(ZInterp,[numXPts,numYPts]);
% 
% %Plot the biquadratically interpolated surface and the control points.
% figure(2)
% clf
% hold on
% surface(X,Y,ZInterp,'edgeColor','None')
% scatter3(XB(:),YB(:),ZGrid(:),200,'.k')
% view(25,50)
%
%REFERENCES:
%[1] D. Smith, "NOAA technical memorandum NOS NGS 84: Biquadratic
%    interpolation," National Oceanic and Atmospheric Administration,
%    National Geodetic Survey, Silver Spring, MD, Tech. Rep., Sep. 2020.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPts=size(xyPts,2);
zInterp=zeros(1,numPts);

numX=size(zGrid,1);
numY=size(zGrid,2);

for curPt=1:numPts
    xCur=xyPts(1,curPt);
    yCur=xyPts(2,curPt);
    
    idxX=round((xCur-x0)/deltaX)+1;
    idxY=round((yCur-y0)/deltaY)+1;
    
    %The points cannot be endpoints.
    idxX=idxX+(idxX==1)-(idxX==numX);
    idxY=idxY+(idxY==1)-(idxY==numY);
    
    %First, interpolate across X for each Y.
    x0Cur=x0+deltaX*(idxX-2);
    y0Cur=y0+deltaY*(idxY-2);
    
    zInterp1=quadraticInterpUnif3Pt(xCur,x0Cur,deltaX,zGrid(idxX-1,idxY-1),zGrid(idxX,idxY-1),zGrid(idxX+1,idxY-1));
    zInterp2=quadraticInterpUnif3Pt(xCur,x0Cur,deltaX,zGrid(idxX-1,idxY),zGrid(idxX,idxY),zGrid(idxX+1,idxY));
    zInterp3=quadraticInterpUnif3Pt(xCur,x0Cur,deltaX,zGrid(idxX-1,idxY+1),zGrid(idxX,idxY+1),zGrid(idxX+1,idxY+1));
    
    zInterp(curPt)=quadraticInterpUnif3Pt(yCur,y0Cur,deltaY,zInterp1,zInterp2,zInterp3);
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
