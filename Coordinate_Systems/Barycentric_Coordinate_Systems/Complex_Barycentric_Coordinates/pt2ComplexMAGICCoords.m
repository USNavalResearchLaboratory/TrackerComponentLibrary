function c=pt2ComplexMAGICCoords(x,v)
%%PT2COMPLEXMAGICCOORDS Convert a 2D point in or around a polygon to Made-
%       to-order Angle Guided Interpolating Coordinates (MAGIC). This is a
%       is a complex barycentric coordinate system in 2D, which is defined 
%       in [1]. These coordinates can be useful for Barycentric mapping
%       and for interpolation. These coordiantes are not valid on the
%       boundary of the defining polygon.
%
%INPUTS: x The 2XnumPts set of points to convert into MAGIC coordiantes.
%        v The 2XnumBases set of vertices of the polygon defining the
%          coordinate system. This does not need to be convex. The vertices
%          must be given in COUNTERCLOCKWISE order. It is assumed that the
%          first vertex is not repeated at the end. Whether the points are
%          in counterclockwise order can be deduced from the sign of
%          signedPolygonArea.
%
%OUTPUTS: c A numBasesXnumPts set of complex barycentric weights forming
%           the MAGIC coordiantes for each point.
%
%EXAMPLE 1:
%In this example, the coordinates of two points are found and the points
%are recreated from the c values. The residual error of the reverse
%conversion is seen to be on the order of finite precision errors.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% x=[[6;6],[1;2]];
% c=pt2ComplexMAGICCoords(x,v);
% xBack=complexBarycentricCoords2Pt(c,v);
% ResidualErr=max(max(abs(xBack-x)))
%
%EXAMPLE 2:
%This draws a polygon and a set of horizontal and vertical lines in the
%the polygon. Then, the points of the lines are converted to complex MAGIC
%coordinates and the vertices of the polygon are moved to make a non-convex
%polygon. After conversion back to Cartesian coordinates, one can see how
%the original grid in the polygon has been warped.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;9];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% 
% vertLines=-1:1:10;
% numVertLines=length(vertLines);
% horizLines=1:10;
% numHorizLines=length(horizLines);
% numLinPts=100;
% xVert=zeros(2,numLinPts,numVertLines);
% xHoriz=zeros(2,numLinPts,numHorizLines);
% for k=1:numVertLines
%     xVert(1,:,k)=vertLines(k);
%     xVert(2,:,k)=linspace(0,11,numLinPts);
%     sel=pointIsInPolygon(v,xVert(:,:,k),false)==0;
%     xVert(:,sel,k)=NaN;
% end
% for k=1:numHorizLines
%     xHoriz(1,:,k)=linspace(-2,11,numLinPts);
%     xHoriz(2,:,k)=horizLines(k);
%     sel=pointIsInPolygon(v,xHoriz(:,:,k),false)==0;
%     xHoriz(:,sel,k)=NaN;
% end
% 
% figure(1)
% clf
% hold on
% plot([v(1,:),v(1,1)],[v(2,:),v(2,1)],'-c','linewidth',2)
% axis([-2, 12, 0, 12])
% for k=1:numVertLines
%     plot(xVert(1,:,k),xVert(2,:,k),'-r')
% end
% for k=1:numHorizLines
%     plot(xHoriz(1,:,k),xHoriz(2,:,k),'-b')
% end
% cVert=pt2ComplexMAGICCoords(xVert(:,:),v);
% cHoriz=pt2ComplexMAGICCoords(xHoriz(:,:),v);
% 
% %Move the vertices, but keep the shape convex and the vertices still in
% %counterclockwise order.
% v(:,1)=[3;2];
% v(:,2)=[8;0];
% v(:,3)=[6;6];
% v(:,4)=[9;11];
% v(:,5)=[3;10];
% 
% figure(2)
% clf
% hold on
% plot([v(1,:),v(1,1)],[v(2,:),v(2,1)],'-c','linewidth',2)
% axis([-2, 12, 0, 12])
% %Synthesize the corresponding points in the transformed coordinate
% %system.
% xTransVert=reshape(complexBarycentricCoords2Pt(cVert,v),[2,numLinPts,numVertLines]);
% xTransHoriz=reshape(complexBarycentricCoords2Pt(cHoriz,v),[2,numLinPts,numHorizLines]);
% for k=1:numVertLines
%     plot(xTransVert(1,:,k),xTransVert(2,:,k),'-r')
% end
% for k=1:numHorizLines
%     plot(xTransHoriz(1,:,k),xTransHoriz(2,:,k),'-b')
% end
%
%REFERENCES:
%[1] O. Weber, M. Ben-Chen, C. Gotsman, and K. Hormann, "A complex view of
%    barycentric mappings," Computer Graphics Forum, vol. 30, no. 5, pp.
%    1533-1542, Aug. 2011.
%
%October 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVert=size(v,2);
numPts=size(x,2);

%Complexify the points.
zPts=x(1,:)+1j*x(2,:);

%Complexify the vertices.
zj=v(1,:)+1j*v(2,:);

e=zeros(numVert,1);
for k=1:(numVert-1)
    e(k)=zj(k+1)-zj(k);
end
e(numVert)=zj(1)-zj(numVert);

gammaj=zeros(numVert,1);
c=zeros(numVert,numPts);
for curZ=1:numPts
    z=zPts(curZ);
    r=zj-z;
    
    %Equation 18
    for k=1:(numVert-1)
        gammaj(k)=abs(e(k)/(r(k)*r(k+1)))*(1/(pi-imag(log(r(k+1)/r(k)))));
    end
    gammaj(numVert)=abs(e(numVert)/(r(numVert)*r(1)))*(1/(pi-imag(log(r(1)/r(numVert)))));
    
    %Equation 10
    c(1,curZ)=gammaj(1)*(r(2)/e(1))-gammaj(numVert)*(r(numVert)/e(numVert));
    for k=2:(numVert-1)
        c(k,curZ)=gammaj(k)*(r(k+1)/e(k))-gammaj(k-1)*(r(k-1)/e(k-1));
    end
    c(numVert,curZ)=gammaj(numVert)*(r(1)/e(numVert))-gammaj(numVert-1)*(r(numVert-1)/e(numVert-1));
    c(:,curZ)=c(:,curZ)/sum(c(:,curZ));
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
