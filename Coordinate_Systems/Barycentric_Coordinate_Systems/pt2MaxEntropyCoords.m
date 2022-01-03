function phi=pt2MaxEntropyCoords(x,v,algorithm,maxIter,epsScal)
%%PT2MAXENTROPYCOORDS Convert a point in a polygon to maximum entropy
%       coordinates. This is a type of barycentric coordinate system in 2D
%       that applies to points in polygons and is described in Section 7 of
%       [1]. These coordinates can be useful for Barycentric mapping and
%       for interpolation.
%
%INPUTS: x The 2XnumPts set of points to convert into maximum entropy
%          coordinates.
%        v The 2XnumBases set of vertices of the polygon defining the
%          coordinate system. The vertices can be given in clockwise or
%          counterclockwise order. It is assumed that the first vertex is
%          not repeated at the end.
% algorithm An optional parameter specifying the edge weight function to
%          use. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the 
%            edge weight of Hormann and Sukumar, which is the first one
%            given in [1].
%          1 Use the second edge weight given in [1].
%  maxIter The maximum number of iterations of Newton's method to use. The
%          default if omitted or an empty matrix is passed is 100 and
%          generally, far fewer iterations are actually needed.
%  epsScal Convergence is determined when
%          abs(stepSize)<=epsScal*eps(lambda), where lambda is the 2X1
%          vector parameter being estimated from which the weights phi are
%          derived in [1]. The default if omitted or an empty matrix is
%          passed is 1.
%
%OUTPUTS: phi A numBasesXnumPts set of the maximum entropy coordinates of
%             the points given in x. Note that sum(phi(:,k))=1 for any k.
%             Points on the edge and outside of the polygon will typically
%             return a vector of NaNs.
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
% x=[[6;5],[8;2]];
% phi=pt2MaxEntropyCoords(x,v);
% xBack=barycentricCoords2Pt(phi,v);
% ResidualErr=max(max(abs(xBack-x)))
%
%EXAMPLE 2:
%This draws a convex polygon and a set of horizontal and vertical lines in
%the polygon. Then, the points of the lines are converted to maximum
%entropy coordinates and the vertices of the polygon are moved. After
%conversion back to Cartesian coordinates, one can see how the original
%grid in the polygon has been warped.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% 
% vertLines=-1:1:10;
% numVertLines=length(vertLines);
% horizLines=1:10;
% numHorizLines=length(horizLines);
% numLinPts=99;
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
% phiVert=pt2MaxEntropyCoords(xVert(:,:),v);
% phiHoriz=pt2MaxEntropyCoords(xHoriz(:,:),v);
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
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(epsScal))
    epsScal=1;
end

if(nargin<4||isempty(maxIter))
    maxIter=100;
end

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

numPts=size(x,2);
numBases=size(v,2);
phi=zeros(numBases,numPts);
for curPt=1:numPts
    if(any(isnan(x(:,curPt))))
       phi(:,curPt)=NaN;
       continue;
    end

    diffs=bsxfun(@minus,x(:,curPt),v);
    mags=sqrt(sum(diffs.*diffs,1));

    rhoVals=zeros(numBases,1);
    if(algorithm==0)
        for k=1:(numBases-1)
            rhoVals(k)=mags(k)+mags(k+1)-norm(v(:,k)-v(:,k+1));
        end
        rhoVals(numBases)=mags(numBases)+mags(1)-norm(v(:,numBases)-v(:,1));
    elseif(algorithm==1)
        for k=1:(numBases-1)
            rhoVals(k)=mags(k)*mags(k+1)+dot(diffs(:,k),diffs(:,k+1));
        end
        rhoVals(numBases)=mags(numBases)*mags(1)+dot(diffs(:,numBases),diffs(:,1));
    else
       error('Unknown algorithm specified.') 
    end

    piVals=zeros(1,numBases);
    piVals(1)=prod(rhoVals(2:(numBases-1)));
    for k=2:(numBases-1)
        piVals(k)=prod(rhoVals([1:(k-2),(k+1):numBases]));
    end
    piVals(numBases)=prod(rhoVals(1:(numBases-2)));

    mi=piVals/sum(piVals);
    d=-diffs;

    %lambda can start uniform.   
    lambda=[1/2;1/2];
    %Use Newton's method.
    for curIter=1:maxIter
        [F,H]=getFAndH(lambda,d,mi);
        
        if(any(~isfinite(H(:))))
           lambda=[NaN;NaN];
           break;
        end
        
        stepVal=pinv(H)*F;
        
        if(all(abs(stepVal)<=epsScal*eps(lambda)))
            %Convergence to some multiple of finite precision limits.
            break; 
        end
        lambda=lambda-stepVal;
    end

    wi=mi.*exp(sum(bsxfun(@times,lambda,d),1));
    phi(:,curPt)=wi/sum(wi);
end
end

function [gradF,H]=getFAndH(lambda,d,mi)
%%GETFANDH Get the the gradient and the Hessian of the cost function.

mExpdTerms=bsxfun(@times,mi.*exp(sum(bsxfun(@times,lambda,d),1)),d);

%The gradient.
gradF=sum(mExpdTerms,2);
H=zeros(2,2);
H(1,1)=sum(mExpdTerms(1,:).*d(1,:));
H(2,1)=sum(mExpdTerms(1,:).*d(2,:));
H(1,2)=H(2,1);
H(2,2)=sum(mExpdTerms(2,:).*d(2,:));

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
