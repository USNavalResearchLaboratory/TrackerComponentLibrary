function phi=pt2ThreePtCoords(x,v,p)
%%PT2THREEPTCOORDS Convert a 2D point to coordinates in the power family of
%      barycentric conversions, which is given in Corollary 2 in [1].
%      Theses are often called three-point coordiantes. Within a convex
%      polygon, this covers Wachpress coordinates, mean value coordinates,
%      discrete harmonic coordiantes, and more.
%
%INPUTS: x The 2XnumPts set of points to convert into a type of three-point
%          (power) coordiantes.
%        v The 2XnumBases set of vertices of the polygon defining the
%          coordinate system. This does not need to be convex. The vertices
%          must be given in COUNTERCLOCKWISE order. It is assumed that the
%          first vertex is not repeated at the end. Whether the points are
%          in counterclockwise order can be deduced from the sign of
%          signedPolygonArea.
%        p The power in the three-point formula in Corollary 2 in [1]. The
%          default if omitted or an empty matrix is passed is 0.
%
%OUTPUTS: phi A numBasesXnumPts set of the coordinates of the
%             points given in x. Note that sum(phi(:,k))=1.
%
%EXAMPLE 1:
%This example shows that the coefficients returned with p=0 are the same as
%Wachspress coordinates (withing finite precision bounds) both for a point
%inside the polygon as well as for a point outside the polygon.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% x=[[6;6],[1;2]];
% phi0=pt2WachspressCoords(x,v);
% c0=pt2ThreePtCoords(x,v,0);
% CoordDiff0=max(abs(phi0(:)-c0(:)))
%
%EXAMPLE 2:
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
% phi=pt2ThreePtCoords(x,v,2);
% xBack=barycentricCoords2Pt(phi,v);
% ResidualErr=max(max(abs(xBack-x)))
%
%REFERENCES:
%[1] O. Weber, M. Ben-Chen, C. Gotsman, and K. Hormann, "A complex view of
%    barycentric mappings," Computer Graphics Forum, vol. 30, no. 5, pp.
%    1533-1542, Aug. 2011.
%
%October 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(p))
    p=0;
end

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
phi=zeros(numVert,numPts);

for curZ=1:numPts
    z=zPts(curZ);
    r=zj-z;
    
    %Equation 13
    for k=1:(numVert-1)
        gammaj(k)=e(k)/imag(r(k)'*r(k+1))*(abs(r(k+1))^p/r(k+1)-abs(r(k))^p/r(k));
    end
    gammaj(numVert)=e(numVert)/imag(r(numVert)'*r(1))*(abs(r(1))^p/r(1)-abs(r(numVert))^p/r(numVert));
    
    %Equation 10
    phi(1,curZ)=gammaj(1)*(r(2)/e(1))-gammaj(numVert)*(r(numVert)/e(numVert));
    for k=2:(numVert-1)
        phi(k,curZ)=gammaj(k)*(r(k+1)/e(k))-gammaj(k-1)*(r(k-1)/e(k-1));
    end
    phi(numVert,curZ)=gammaj(numVert)*(r(1)/e(numVert))-gammaj(numVert-1)*(r(numVert-1)/e(numVert-1));
    phi(:,curZ)=phi(:,curZ)/sum(phi(:,curZ));
end

%If it is complex, that is only due to finite precision errors.
phi=real(phi);
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
