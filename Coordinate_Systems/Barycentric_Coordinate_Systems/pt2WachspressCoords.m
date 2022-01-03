function [phi,dphi]=pt2WachspressCoords(x,v,algorithm)
%%PT2WACHSPRESSCOORDS Convert a point in a convex polygon to Wachspress
%       coordinates. This is a type of barycentric coordinate system in 2D
%       that applies to points in convex polygons and is described in
%       Section 4 of [1]. These coordinates can be useful for Barycentric 
%       mapping and for interpolation.
%
%INPUTS: x The 2XnumPts set of points to convert into Wachspress
%          coordinates.
%        v The 2XnumBases set of vertices of the convex polygon defining
%          the coordinate system. The vertices must be given in
%          COUNTERCLOCKWISE order. It is assumed that the first vertex is
%          not repeated at the end. Convexity can be checked using the
%          polygonIsConvex function and whether the points are in
%          counterclockwise order can be deduced from the sign of
%          signedPolygonArea. 
% algorithm An optional parameter specifying how the coefficients are
%          computed. Possible values are:
%          0 Use the formula of Equation 4.2 in [1].
%          1 Use the formula of Equation 4.7 in [1].
%          2 (The default if omitted or an empty matrix is passed) Use the
%            method based on normals the the edges given in Sections
%            4.2-4.3 of [1]. This algorithm is also given in [2].
%
%OUTPUTS: phi A numBasesXnumPts set of the Wachpress coordinates of the
%             points given in x. Note that sum(phi(:,k))=1 for any k and
%             for points in the polygon, it is guaranteed that
%             all(phi(:)>=0).
%        dphi A numBasesX2XnumPts set of derivatives of the elements of phi
%             (rows) taken with respect to the elements of x (columns) for
%             each input measurement (3rd dimension).
%
%EXAMPLE 1:
%In this example, the coordinates of two points are found and the points
%are recreated from the phi values. The residual error of the reverse
%conversion is seen to be on the order of finite precision errors.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% x=[[6;6],[1;2]];
% phi=pt2WachspressCoords(x,v);
% xBack=barycentricCoords2Pt(phi,v);
% ResidualErr=max(max(abs(xBack-x)))
%
%EXAMPLE 2:
%This draws a convex polygon and a set of horizontal and vertical lines in
%the polygon. Then, the points of the lines are converted to Wachspress
%coordinates and the vertices of the polygon are moved. After conversion
%back to Cartesian coordinates, one can see how the original grid in the
%polygon has been warped.
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
% phiVert=pt2WachspressCoords(xVert(:,:),v);
% phiHoriz=pt2WachspressCoords(xHoriz(:,:),v);
% 
% %Move the vertices, but keep the shape convex and the vertices still in
% %counterclockwise order.
% v(:,1)=[3;2];
% v(:,2)=[8;0];
% v(:,3)=[10;10];
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
% xTransVert=reshape(barycentricCoords2Pt(phiVert,v),[2,numLinPts,numVertLines]);
% xTransHoriz=reshape(barycentricCoords2Pt(phiHoriz,v),[2,numLinPts,numHorizLines]);
% for k=1:numVertLines
%     plot(xTransVert(1,:,k),xTransVert(2,:,k),'-r')
% end
% for k=1:numHorizLines
%     plot(xTransHoriz(1,:,k),xTransHoriz(2,:,k),'-b')
% end
%
%EXAMPLE 3:
%In this example, the derivative of phi is verified with respect to finite
%differencing. The relative error is on the order of what one would expect
%due to finite precision limitations.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% x=[6;5];
% [phi,dphi]=pt2WachspressCoords(x,v);
% f=@(y)pt2WachspressCoords(y,v);
% dPhiNumDiff=numDiff(x,f,numVertices);
% RelErr=max(max(abs((dphi-dPhiNumDiff)./dPhiNumDiff)))
%
%EXAMPLE 4:
%Wachspress coordinates can be used for interpolation. In this examples,
%the coordinates are fit based on an initial set of vertices and then used
%to interpolate other values. The gradients returnd by the Wachspress
%coordinates can be used to interpolate the partial derivatives of the
%interpolated value with respect to the original point, as is shown here
%with the agreement of the interpolated value and the name value found via
%numeric differentiation. The relative error is small.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% iVals=[v(1,:).^2+v(2,:).^2;
%        v(1,:).^3-v(2,:).^2];
% x=[6;5];
% [~,dphi]=pt2WachspressCoords(x,v);
% f=@(x)barycentricCoords2Pt(pt2WachspressCoords(x,v),iVals);
% dInterpdxNumDiff=numDiff(x,f,2);
% dInterpdx=barycentricCoords2Pt(dphi,iVals);
% RelErr=max(max(abs((dInterpdx-dInterpdxNumDiff)./dInterpdxNumDiff)))
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%[2] M. S. Floater, A. Gillette, and N. Sukumar, "Gradient bounds for
%    Wachspress coordinates on polytopes," SIAM Journal on Numerical
%    Analysis, vol. 52, no. 1, pp. 515-532, 2014.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=2;
end

numPts=size(x,2);
numBases=size(v,2);

if(numBases<3)
    error('A minimum of 3 basis vertices must be provided.')
end

if(algorithm<2)
    C=zeros(numBases,1);
    C(1)=triangleArea(v(:,numBases),v(:,1),v(:,2));
    for k=2:(numBases-1)
        C(k)=triangleArea(v(:,k-1),v(:,k),v(:,k+1));
    end
    C(numBases)=triangleArea(v(:,numBases-1),v(:,numBases),v(:,1));
else
    w=zeros(numBases,1);
    p=zeros(2,numBases);
    u=getPolygonNormals(v);
end

phi=zeros(numBases,numPts);
if(nargout>1)
    dphi=zeros(numBases,2,numPts);
    R=zeros(numBases,2);
end
for curPt=1:numPts
    if(algorithm<2&&algorithm>=0)
        A=zeros(numBases,1);
        for k=1:(numBases-1)
            A(k)=triangleArea(x(:,curPt),v(:,k),v(:,k+1));
        end
        A(numBases)=triangleArea(x(:,curPt),v(:,numBases),v(:,1));

        if(nargout>1)
            dA=zeros(numBases,2);
            for k=1:(numBases-1)
                dA(k,:)=[v(2,k)-v(2,k+1),v(1,k+1)-v(1,k)]/2;
            end
            dA(numBases,:)=[v(2,numBases)-v(2,1),v(1,1)-v(1,numBases)]/2;
        end
    
        if(algorithm==0)%Equation 4.2 in [1].
            phi(1,curPt)=C(1)/(A(numBases)*A(1));
            for k=2:numBases
                phi(k,curPt)=C(k)/(A(k-1)*A(k)); 
            end
            denom=sum(phi(:,curPt));

            if(nargout>1)
                %Derivatives of the non-normalized phi values.
                R(1,:)=-(C(1)./(A(numBases)*A(1)).^2).*(A(numBases).*dA(1,:)+dA(numBases,:).*A(1,:));
                for k=2:numBases
                    R(k,:)=-(C(k)./(A(k-1)*A(k)).^2).*(A(k-1).*dA(k,:)+dA(k-1,:).*A(k,:));
                end

                %Derivatives of the normalized phi values.
                dphi(:,:,curPt)=(denom.*R-bsxfun(@times,phi(:,curPt),sum(R,1)))./denom.^2;
            end

            phi(:,curPt)=phi(:,curPt)/denom;
        else%Equation 4.7 in [1].
            phi(1,curPt)=C(1)*prod(A(2:(numBases-1)));
            for k=2:numBases
                phi(k,curPt)=C(k)*prod(A([1:(k-2),(k+1):numBases]));
            end
            denom=sum(phi(:,curPt));

            if(nargout>1)
                %Derivatives of the non-normalized phi values.
                idx=2:(numBases-1);
                R(1,1)=C(1)*productDeriv(A(idx),dA(idx,1));
                R(1,2)=C(1)*productDeriv(A(idx),dA(idx,2));
                for k=2:numBases
                    idx=[1:(k-2),(k+1):numBases];
                    R(k,1)=C(k)*productDeriv(A(idx),dA(idx,1));
                    R(k,2)=C(k)*productDeriv(A(idx),dA(idx,2));
                end

                %Derivatives of the normalized phi values.
                dphi(:,:,curPt)=(denom.*R-bsxfun(@times,phi(:,curPt),sum(R,1)))./denom.^2;
            end

            phi(:,curPt)=phi(:,curPt)/denom;
        end
    elseif(algorithm==2)
        for k=1:numBases
            h=(v(:,k)-x(:,curPt))'*u(:,k);
            p(:,k)=u(:,k)/h;
        end
        w(1)=det([p(:,numBases),p(:,1)]);
        
        for k=2:numBases
            w(k)=det([p(:,k-1),p(:,k)]);
        end
        phi(:,curPt)=w/sum(w);
        
        if(nargout>1)
            R(1,:)=p(:,numBases)+p(:,1);
            for k=2:numBases
                R(k,:)=p(:,k-1)+p(:,k);
            end
            phiRSum=sum(bsxfun(@times,phi(:,curPt),R),1);
            dphi(:,:,curPt)=phi(:,curPt).*bsxfun(@minus,R,phiRSum);
        end
    else
        error('Unknown algorithm specified.')
    end
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
