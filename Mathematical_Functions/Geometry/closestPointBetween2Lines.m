function [midpoint,minDist,t]=closestPointBetween2Lines(a1,b1,a2,b2)
%%CLOSESTPOINTBETWEEN2LINES Given two non-parallel lines in numDim
%              dimensions that do not necessarily intersect, find the point
%              closest to both lines. This is done by finding the point on
%              each line that is closest to the other line, and then
%              averaging the points (get the point halfway between the
%              lines).
%
%INPUTS: a1,b1 numDimXN matrix of parameters for the first line in N pairs
%              of lines, the closes point between each being desired. The
%              ith line is represented in parametic form as
%              r=a1(:,i)*t+b1(:,i), where t is a scalar parameter and p is
%              a point on the line.
%        a2,b2 The numDimXN set of parameters for the second lines.
%                           
%OUTPUTS: midpoint The numDimXN set of closest points between each of the N
%                  pairs of lines.
%          minDist The NX1 set of minimum distances between the pairs of
%                  lines.
%                t The 2XN set of values of t used in the linear equations
%                  to get the minimum points. These values can be positive
%                  or negative.
%
%Line 1 is r1=a1*t1+b1 and line 2 is r2=a2*t2+b2. We want to find the
%values of t1 and t2 such that norm(r1-r2)^2 is minimized.
%norm(r1-r2)^2=norm(a1*t1+b1-a2*t2-b2)^2
%             =norm(A*T+c)^2 where A=[a1,-a2] and c=b1-b2 and T=[t1;t2]
%             =T'*A'*A*T+2*c'*A*T+c'*c
%Setting the derivative with respect to the vector T equal to zero, we get
%             2*A'*A*T+2*A'*c=0
%The solution is T=-inv(A'*A)*A'*c, which (being the zero deivative point)
%is an extrema and being the only extrema (unless the lines are parallel)
%is the closest point between the lines. The pinv function is used to
%evaluate (A'*A)\A' as it is more stable.
%
%EXAMPLE 1:
%This is a simple 2D case, so the lines intersect. The starting points are
%drawn in blue and lines with positive t are drawn in blue. The
%intersection point is marked in red.
% x1=[4;3];
% u1=[3;2];
% u1=u1/norm(u1);
% x2=[6;7];
% u2=[6;-1];
% u2=u2/norm(u2);
% r=10;
% p1=[x1,x1+r*u1];
% p2=[x2,x2+r*u2];
% figure(1)
% clf
% hold on
% scatter(x1(1),x1(2),400,'.b')
% scatter(x2(1),x2(2),400,'.b')
% plot(p1(1,:),p1(2,:),'-b')
% plot(p2(1,:),p2(2,:),'-b')
% pt=closestPointBetween2Lines(u1,x1,u2,x2);
% scatter(pt(1),pt(2),300,'.r')
%
%EXAMPLE 2:
%This is the same as example 1, but is in 3D. This time, the lines do not
%intersect in 3D, so we use t to compute the two nearest points on either
%line and we draw a red line between them.
% x1=[4;3;8];
% u1=[3;2;0];
% u1=u1/norm(u1);
% x2=[6;7;-2];
% u2=[6;-1;4];
% u2=u2/norm(u2);
% 
% r=20;
% p1=[x1,x1+r*u1];
% p2=[x2,x2+r*u2];
% figure(1)
% clf
% hold on
% scatter3(x1(1),x1(2),x1(3),400,'.b')
% scatter3(x2(1),x2(2),x2(3),400,'.b')
% plot3(p1(1,:),p1(2,:),p1(3,:),'-b')
% plot3(p2(1,:),p2(2,:),p2(3,:),'-b')
% [p,~,t]=closestPointBetween2Lines(u1,x1,u2,x2);
% scatter3(p(1),p(2),p(3),200,'.r')
% pt1=x1+t(1)*u1;
% pt2=x2+t(2)*u2;
% plot3([pt1(1);pt2(1)],[pt1(2);pt2(2)],[pt1(3);pt2(3)],'-r')
% view(8,90)
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(a1,2);
numDim=size(a1,1);

minDist=zeros(N,1);
midpoint=zeros(numDim,N);
t=zeros(2,N);
for i=1:N
    A=[a1(:,i),-a2(:,i)];
    c=b1(:,i)-b2(:,i);

    %pinv(A)=(A'*A)\A', we are doing pinv(A)*c using lsqminnorm.
    TMin=-lsqminnorm(A,c);

    t(:,i)=TMin;
    
    r1Min=a1*TMin(1)+b1;
    r2Min=a2*TMin(2)+b2;
    
    minDist(i)=norm(r1Min-r2Min);
    midpoint(:,i)=(r1Min+r2Min)/2;
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
