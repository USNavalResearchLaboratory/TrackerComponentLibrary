function [m,r,s]=findEllipsPlaneIntersect3D(z,Q,n,q)
%%FINDELLIPSPLANEINTERSECT3D The intersection of an ellipsoid and a plane
%               in three dimensions is an ellipse. This function determines
%               the parameters of the ellipse of intersection in 3D.
%
%INPUTS: z The 3X1 center of the ellipsoid. This is a real vector.
%        A A 3X3 real,symmetric, positive definite matrix that
%          specifies the size and shape of the ellipsoid, where a point zp
%          is inside of the ellipsoid if
%          (zp-z)'*A*(zp-z)<=1.
%        n A 3X1 unit vector that is normal to the plane.
%        q A point on the plane. n and q uniquely define the plane.
%
%OUTPUTS: m, r, s These are three 3X1 vectors that uniquely specify
%                 the ellipse of intersection. The equation of a point on
%                 the ellipse of intersection is
%                 x=m+cos(theta)*r+sin(theta)*s for theta from 0 to 2*pi.
%                 Thus, m is the center of the ellipse and r and s are
%                 the two axes of the ellipse in 3D. 
%
%This function implements the algorithm of [1]. The
%ellipsoidAndPlaneIntersect function is used to move q to a point within
%the ellipsoid, as is required by the algorithm. The algorithm of [1] is
%only for axis-aligned ellipsoids that are centered at zero. Here, the
%provided ellipsoid is moved to the origin (z=0) and the plane is similarly
%shifted with q=q-z. This means that the final results must be shifted
%back. Also, an eigenvalue decomposition is used to obtain a rotation
%matrix to make Q diagonal and to appropriately rotate n and q. This must
%also be undone in the end.
%
%EXAMPLE:
%This example draws the ellipsoid, the plane, and the ellipse of
%intersection.
% Q=[36, 4, 10;
%     4, 30, 16;
%     10, 16, 24];
% z=[0;1;0];
% n=[1;0;0];
% q=z+[0.05;0;0];
% 
% figure(1)
% clf
% hold on
% h=drawEllipse(z,Q,1);
% alpha(h{1},0.5);%Make transparent.
% 
% vecs=null(n');
% numPoints1=10;
% x=linspace(-1,1,numPoints1);
% [X,Y]=meshgrid(x,x);
% 
% numPoints=numPoints1^2;
% 
% xp=zeros(numPoints1,numPoints1);
% yp=zeros(numPoints1,numPoints1);
% zp=zeros(numPoints1,numPoints1);
% for curPoint=1:numPoints
%     v=q+vecs(:,1)*X(curPoint)+vecs(:,2)*Y(curPoint);
%     xp(curPoint)=v(1);
%     yp(curPoint)=v(2);
%     zp(curPoint)=v(3);
% end
% 
% view(-45,20)
% h=surf(xp,yp,zp);
% alpha(h,0.5);%Make transparent.
% 
% %Find the ellipse of intersection
% [m,r,s]=findEllipsPlaneIntersect3D(z,Q,n,q);
% 
% %Draw the ellipse.
% theta=linspace(0,2*pi);
% xyzPoints=bsxfun(@plus,m,bsxfun(@times,cos(theta),r)+bsxfun(@times,sin(theta),s));
% plot3(xyzPoints(1,:),xyzPoints(2,:),xyzPoints(3,:),'-r','linewidth',10)
% 
%REFERENCES:
%[1] P. P. Klein, "On the ellipsoid and plane intersection equation,"
%    Applied Mathematics, vol. 3, no. 11, pp. 1634-1640, Nov. 2012.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Make sure that the ellipsoid and the plane intersect and get a point on
%the plane within the ellipsoid. This point will replace the original q.
[boolVal,q]=ellipsoidAndPlaneIntersect(z,Q,n,q);

%If the ellipsoid and the plane do not intersect.
if(boolVal==false)
    m=[];
    r=[];
    s=[];
    return;
end

%Move the origin to the center of the ellipsoid.
q=q-z;

%Next, rotate the ellipsoid to be axis-aligned.
[V,D]=eig(Q);
d1=sqrt(diag(D));%D1 is defined before Equation 3.
D1=diag(d1);

%Rotate q and n into the axis-aligned, ellipsoid-centered coordinate
%system.
q=V'*q;
n=V'*n;

%The vectors (r and s) orthogonal to n can be obtained using the null
%function.
orthogVecs=null(n');
r=orthogVecs(:,1);
s=orthogVecs(:,2);

%Rotate the vectors in the plane such that dot(D1*r,D1*s)=0. This is from
%Section 2 of [1].
D1rD1s=dot(D1*r,(D1*s));
if(abs(D1rD1s)>eps())
    omega=1/2*atan(2*D1rD1s/(dot(D1*r,D1*r)-dot(D1*s,D1*s)));
    if(~isfinite(omega))
       omega=pi/4; 
    end
    
    sinOmega=sin(omega);
    cosOmega=cos(omega);
    
    r=cosOmega*r+sinOmega*s;
    s=-sinOmega*r+cosOmega*s;
end

%Equation 18
D1r=D1*r;
beta1=dot(D1r,D1r);
D1s=(D1*s);
beta2=dot(D1s,D1s);

%Equation 12
kappa=dot(q,n);

%Equation 25
d=kappa^2/dot((1./d1).^2,n.^2);

%Equation 10
A=sqrt((1-d)/beta1);
B=sqrt((1-d)/beta2);

%Equation 39 for the center of the ellipse
D1n=D1*n;
m=kappa*(n-(dot(D1n,D1r)/beta1)*r-(dot(D1n,D1s)/beta2)*s);

%Incorporate A and B into r and s
r=A*r;
s=B*s;

%Finally, undo the rotations
r=V*r;
s=V*s;
m=V*m;

%Undo the origin shift.
m=m+z;

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
