function [phi,dphi]=pt2MeanValCoords(x,v,algorithm)
%%PT2MEANVALCOORDS Convert a point to mean value coordinates. This is a
%       type of barycentric coordinate system in 2D that applies to points
%       in convex and non-convex polygons and is described in Section 4 of
%       [1]. These coordinates can be useful for interpolation and for
%       barycentric mapping of convex polygons.
%
%INPUTS: x The 2XnumPts set of point to convert into mean value
%          coordinates.
%        v The 2XnumBases set of vertices of the convex polygon defining
%          the coordinate system in clockwise or counterclockwise order. It
%          is assumed that the first vertex is not repeated at the end.
% algorithm An optional parameter specifying how the coefficients are
%          computed. Possible values are:
%          0 Use Equation 5.1 of [1].
%          1 (The default if omitted or an empty matrix is passed) Use
%            the Equation in Section 5.2 of [1]. This is more robust at the
%            edges of the polygon.
%
%OUTPUTS: phi A numBasesXnumPts set of the man value coordinates of the
%             points given in x. Note that sum(phi(:,k))=1 for any k.
%        dphi A numBasesX2XnumPts set of derivatives of the elements of phi
%             (rows) taken with respect to the elements of x (columns) for
%             each input measurement (3rd dimension).
%
%EXAMPLE 1:
%In this example, the coordinates of two points are found and the points
%are recreated from the phi values. The residual error of the reverse
%conversion is seen to be on the order of finite precision errors. Note
%that the residual error if a point were given outside of the polygon would
%be large --the mapping is not valid for points outside the polygon. With
%algorithm=2, they get mapped inside the polygon.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% x=[[6;6],[7;5]];
% phi=pt2MeanValCoords(x,v);
% xBack=barycentricCoords2Pt(phi,v);
% ResidualErr=max(max(abs(xBack-x)))
%
%EXAMPLE 2:
%This creates a grid on a non-convex polygon. The vertices are moved to
%create a convex polygon and this plots how the points move. Unlike
%Wachspress coordinates, this can handle non-convexivity in the original
%polygon. However, if one transforms a convex polgon into a non-convex one,
%some points might end up out of the transformed non-convex polygon.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[5;5];
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
%    sel=pointIsInPolygon(v,xHoriz(:,:,k),false)==0;
%    xHoriz(:,sel,k)=NaN;
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
% phiVert=pt2MeanValCoords(xVert(:,:),v);
% phiHoriz=pt2MeanValCoords(xHoriz(:,:),v);
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
% x=[6;6];
% [phi,dphi]=pt2MeanValCoords(x,v,0);
% f=@(y)pt2MeanValCoords(y,v,0);
% dPhiNumDiff=numDiff(x,f,numVertices);
% RelErr=max(max(abs((dphi-dPhiNumDiff)./dPhiNumDiff)))
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=1;
end

numPts=size(x,2);
numBases=size(v,2);

phi=zeros(numBases,numPts);
dphi=zeros(numBases,2,numPts);
switch(algorithm)
    case 0%The algorithm of Section 5.1
        t=zeros(numBases,1);
        w=zeros(numBases,1);
        cPerp=zeros(2,numBases);
        R=zeros(2,numBases);
        for curPt=1:numPts
            e=bsxfun(@minus,v,x(:,curPt));
            r=sqrt(sum(e.*e,1));
            e=bsxfun(@rdivide,e,r);

            for k=1:(numBases-1)
                cosAlpha=dot(e(:,k),e(:,k+1));
                sinAlpha=det([e(:,k),e(:,k+1)]);
                t(k)=sinAlpha/(1+cosAlpha);

                cPerp(1,k)=(e(2,k+1)/r(k+1)-e(2,k)/r(k))/sinAlpha;
                cPerp(2,k)=(e(1,k)/r(k)-e(1,k+1)/r(k+1))/sinAlpha;
            end
            cosAlpha=dot(e(:,numBases),e(:,1));
            sinAlpha=det([e(:,numBases),e(:,1)]);
            t(numBases)=sinAlpha/(1+cosAlpha);
            cPerp(1,numBases)=(e(2,1)/r(1)-e(2,numBases)/r(numBases))/sinAlpha;
            cPerp(2,numBases)=(e(1,numBases)/r(numBases)-e(1,1)/r(1))/sinAlpha;

            w(1)=(t(numBases)+t(1))/r(1);
            for k=2:numBases
                w(k)=(t(k-1)+t(k))/r(k);
            end
            phi(:,curPt)=w/sum(w);

            if(nargout>1)
                R(:,1)=(t(numBases)/(t(numBases)+t(1)))*cPerp(:,numBases)+(t(1)/(t(numBases)+t(1)))*cPerp(:,1)+e(:,1)/r(1);
                for k=2:numBases
                    R(:,k)=(t(k-1)*cPerp(:,k-1)+t(k)*cPerp(:,k))/(t(k-1)+t(k))+e(:,k)/r(k);
                end

                phiRSum=sum(bsxfun(@times,phi(:,curPt),R'),1);
                dphi(:,:,curPt)=phi(:,curPt).*bsxfun(@minus,R',phiRSum);
            end
        end
    case 1%The algorithm of Section 5.2.
        w=zeros(numBases,1);
        prodTerms=zeros(numBases,1);
        
        prodCoeffs=zeros(numBases,1);
        prodVals=zeros(numBases,1);
        
        if(nargout>1)
           dProdArg=zeros(numBases,2); 
           dwArg=zeros(numBases,2);
        end

        for curPt=1:numPts
            d=bsxfun(@minus,v,x(:,curPt));
            r=sqrt(sum(d.*d,1)).';

            for k=1:(numBases-1)
                prodTerms(k)=r(k)*r(k+1)+dot(d(:,k),d(:,k+1));
            end
            prodTerms(numBases)=r(numBases)*r(1)+dot(d(:,numBases),d(:,1));
            
            prodVals(1)=prod(prodTerms(2:(numBases-1)));
            prodCoeffs(1)=r(numBases)*r(2)-dot(d(:,numBases),d(:,2));
            w(1)=sqrt(prodCoeffs(1)*prodVals(1));
            for k=2:(numBases-1)
                prodVals(k)=prod(prodTerms([1:(k-2),(k+1):numBases]));
                prodCoeffs(k)=(r(k-1)*r(k+1)-dot(d(:,k-1),d(:,k+1)));
                w(k)=sqrt(prodCoeffs(k)*prodVals(k));
            end
            prodVals(numBases)=prod(prodTerms(1:(numBases-2)));
            prodCoeffs(numBases)=(r(numBases-1)*r(1)-dot(d(:,numBases-1),d(:,1)));
            w(numBases)=sqrt(prodCoeffs(numBases)*prodVals(numBases));

            denom=sum(w);
            phi(:,curPt)=w/denom;

            if(nargout>1)
                dr=-[d(1,:)'./r,d(2,:)'./r];

                for k=1:(numBases-1)
                    dProdArg(k,:)=dr(k,:)*r(k+1)+r(k)*dr(k+1,:)-[d(1,k+1),d(2,k+1)]-[d(1,k),d(2,k)];
                end
                dProdArg(numBases,:)=dr(numBases,:)*r(1)+r(numBases)*dr(1,:)-[d(1,1),d(2,1)]-[d(1,numBases),d(2,numBases)];

                prodCoeffDeriv=dr(numBases,:)*r(2)+r(numBases)*dr(2,:)+[d(1,2),d(2,2)]+[d(1,numBases),d(2,numBases)];
                sel=2:(numBases-1);
                dwArg(1,:)=prodCoeffDeriv*prodVals(1)+prodCoeffs(1)*[productDeriv(prodTerms(sel),dProdArg(sel,1)),productDeriv(prodTerms(sel),dProdArg(sel,2))];
                for k=2:(numBases-1)
                    prodCoeffDeriv=dr(k-1,:)*r(k+1)+r(k-1)*dr(k+1,:)+[d(1,k+1),d(2,k+1)]+[d(1,k-1),d(2,k-1)];
                    sel=[1:(k-2),(k+1):numBases];
                    dwArg(k,:)=prodCoeffDeriv*prodVals(k)+prodCoeffs(k)*[productDeriv(prodTerms(sel),dProdArg(sel,1)),productDeriv(prodTerms(sel),dProdArg(sel,2))];
                end
                prodCoeffDeriv=dr(numBases-1,:)*r(1)+r(numBases-1)*dr(1,:)+[d(1,1),d(2,1)]+[d(1,numBases-1),d(2,numBases-1)];
                sel=1:(numBases-2);
                dwArg(numBases,:)=prodCoeffDeriv*prodVals(numBases)+prodCoeffs(numBases)*[productDeriv(prodTerms(sel),dProdArg(sel,1)),productDeriv(prodTerms(sel),dProdArg(sel,2))];

                dw=bsxfun(@rdivide,dwArg,2*w);
                dphi(:,:,curPt)=(denom.*dw-bsxfun(@times,w,sum(dw,1)))./denom.^2;
            end
        end
    otherwise
        error('Unknown Algorithm Specified.')
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
