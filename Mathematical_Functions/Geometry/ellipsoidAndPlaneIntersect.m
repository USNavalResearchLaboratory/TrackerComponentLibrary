function [boolVal,qNew]=ellipsoidAndPlaneIntersect(z,Q,n,q)
%%ELLIPSOIDANDPLANEINTERSECT Given an ellipsoid and the parameters for a
%                   plane, determine whether the ellipsoid and the plane
%                   intersect and if they intersect, return a point on the
%                   plane inside of the ellipse.
%
%INPUTS: z The numDimX1 center of the ellipsoid. This is a real vector.
%        Q A numDimXnumDim real,symmetric, positive definite matrix that
%          specifies the size and shape of the ellipsoid, where a point zp
%          is inside of the ellipsoid if
%          (zp-z)'*A*(zp-z)<=1.
%        n A numDimX1 unit vector that is normal to the plane.
%        q A point on the plane. n and q uniquely define the plane.
%
%OUTPUTS: boolVal This is true if the plane intersects the ellipsoid and is
%                 false otherwise.
%            qNew This is the numDimX1 point that minimizes
%                 (qNew-z)'*A*(qNew-z) such that qNew=q+x'*null(n'), which
%                 is the equation for a point on the plane. Thus, if
%                 boolVal is true, then qNew is a point inside of the
%                 ellipsoid that is on the plane. If they do not intersect,
%                 then qNew is the smallest value of (zp-z)'*A*(zp-z) that
%                 intersects the plane.
%
%It is clear that if we can find any point zp on the plane such that
%(zp-z)'*Q*(zp-z)<=1, then the plane and the ellipsoid intersect. To find
%such a point, we minimize (zp-z)'*Q*(zp-z) subject to the constraint that
%the desired point z is on the plane. Here, we derive the solution in 2D.
%The extension to more dimensions is straightforward.
%
%For 3D, first we subtract z from the position and from q. Here, we assume
%that the subtraction has already been performed. Thus, the optimization
%problem can be written
%minimize x'*Q*x
%such that x=q+t*r+u*s
%where [T]=null(n'), r=T(:,1) and s=T(:,2) are two vectors in the plane
%(orthogonal to n) and t and u are two scalar values. The minimization is
%thus
%minimize (q+t*r+u*s)'*Q*(q+t*r+u*s)
%over t and u
%We do out the multiplication to get
%minimize q'*Q*q+t^2*r'*Q*r+u^2*s'*Q*s+2*t*q'*Q*r+2*u*q'*Q*s+2*t*u*r'*Q*s
%The derivatives with regard to t and u respectively are 
%d/dt=2*t*r'*Q*r+2*q'*Q*r+2*u*r'*Q*s
%d/du=2*u*s'*Q*s+2*q'*Q*s+2*t*r'*Q*s
%Setting the derivatives equal to zero, we get the system of linear
%equations
%-[q'*Q*r| = [r'*Q*r r'*Q*s|*[t|
% |q'*Q*s]   |r'*Q*s s'*Q*s] |u]
%Thus, the solution for the point minimzing the ellipsoid equation is
%simple. Find t and u and the point is x=q+t*r+u*s --Plus add back the
%value z that was subtracted at the beginning (to get qNew). If x'*Q*x<=1,
%then the point is within the ellipsoid and the plane and the ellipsoid
%intersect. Otherwise, they do not intersect.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The null space vectors will be orthonormal.
orthoVecs=null(n');

numDim=length(n)-1;

q=q-z;

b=zeros(numDim,1);
for curDim=1:numDim
    b(curDim)=-q'*Q*orthoVecs(:,curDim);
end

A=zeros(numDim,numDim);
for curDim1=1:numDim
    A(curDim1,curDim1)=orthoVecs(:,curDim1)'*Q*orthoVecs(:,curDim1);
    
    for curDim2=(curDim1+1):numDim
        A(curDim1,curDim2)=orthoVecs(:,curDim1)'*Q*orthoVecs(:,curDim2);
        A(curDim1,curDim2)=A(curDim1,curDim2);
    end
end

params=A\b;

qNew=q+sum(bsxfun(@times,params',orthoVecs),2);

boolVal=qNew'*Q*qNew<=1;
qNew=qNew+z;

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
