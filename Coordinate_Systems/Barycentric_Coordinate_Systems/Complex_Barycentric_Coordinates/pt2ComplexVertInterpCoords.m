function c=pt2ComplexVertInterpCoords(x,v)
%%PT2COMPLEXVERTINTERPCOORDS Convert a 2D point to the holomorphic
%       vertex-interpolating mapping that is defined in Section 4.2 of [1].
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
%           the complex vertex interpolating coordiantes for each point.
%
%EXAMPLE:
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
% c=pt2ComplexVertInterpCoords(x,v);
% xBack=complexBarycentricCoords2Pt(c,v);
% ResidualErr=max(max(abs(xBack-x)))
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
    
    %If we are given a vertex.
    zeroPt=find(r==0,1);
    if(~isempty(zeroPt))
        c(zeroPt,curZ)=1;
        continue; 
    end
    
    %Equation 17
    for k=1:(numVert-1)
        gammaj(k)=r(k+1)/r(k)-r(k)/r(k+1);
    end
    gammaj(numVert)=r(1)/r(numVert)-r(numVert)/r(1);
    
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
