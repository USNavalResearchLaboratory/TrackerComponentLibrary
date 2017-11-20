function [midpoint,minDist]=closestPointBetween2Lines(a1,b1,a2,b2)
%%CLOSESTPOINTBETWEEN2LINES Given two non-parallel lines in 3D that do not
%                           necessarily intersect, find the point closest
%                           to both lines. This is done by finding the
%                           point on each line that is closes to the other
%                           line, and then averaging the points (get the
%                           point halfway between the lines).
%
%INPUTS: a1,b1  3XN matrix of parameters for the first line in N pairs of
%               lines, the closes point between each being desired. The ith
%               line is represented in parametic form as
%               r=a1(:,i)*t+b1(:,i), where t is a scalar parameter and r is
%               a point on the line.
%        a2,b2  The 3XN set of parameters for the second lines.
%                           
%OUTPUTS: midpoint The 3XN set of closest points between each of the N
%                  pairs of lines.
%         minDist  NX1 set of minimum distances between the pairs of lines.
%
%Line 1 is r1=a1*t1+b1 and line 2 is r2=a2*t2+b2. We want to find the
%values of t1 and t2 such that norm(r1-r2)^2 is minimized.
%norm(r1-r2)^2=norm(a1*t1+b1-a2*t2-b2)^2
%             =norm(A*T+c)^2 where A=[a1,-a2] and c=b1-b2 and T=[t1;t2]
%             =T'*A'*A*T+2*c'*A*T+c'*c
%Setting the derivative with respect to the vector T equal to zero, we get
%             2*A'*A*T+2*A'*c=0
%The solution is T=-inv(A'*A)*A'*c, which (being the zero deivative point) is
%an extrema and being the only extrema (unless the lines are parallel) is
%the closest point between the lines. The pinv function is used to evaluate
%(A'*A)\A' as it is more stable.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(a1,2);
minDist=zeros(3,N);
midpoint=zeros(N,1);

for i=1:N
    A=[a1(:,i),-a2(:,i)];
    c=b1(:,i)-b2(:,i);

    %pinv(A)=(A'*A)\A', we are doing pinv(A)*c using lsqminnorm.
    TMin=-lsqminnorm(A,c);

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
