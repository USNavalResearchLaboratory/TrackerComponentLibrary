function hullIdx=findConvexHull2DJarvis(points)
%%FINDCONVEXHULL2DJARVIS Given a set of points, find the smallest convex
%          polygon enclosing the points in 2D. This is am implementing of
%          Jarvis' algorithm, also known as the gift wrapping algorithm, is
%          used. The complexity is O(n*h), where n is the number of points
%          and h is the number of points in the hull. The indices of the
%          vertices are returned in counterclockwise order.
%
%INPUTS: points A 2XN set of N finite points of the form [x;y], the convex
%               hull about which one desires. Repeated points do not appear
%               to be a problem.
%
%OUTPUTS: hullIdx A numVertX1 set of indices of the points forming the
%                 vertices of the convex hull. The first point is not
%                 repeated at the end even though it is a closed shape.
%                 The points are provided in counterclockwise order.
%
%The original algorithm is given in [1]. As in findConvexHull2D, this
%function relies on the accuracy of the turnOrientation function.
%
%EXAMPLE:
%We generate some random points and then find the hull. The results are
%plotted.
% numPts=1000;
% x=randn(2,numPts);
% 
% idx=findConvexHull2DJarvis(x);
% 
% figure(1)
% clf
% hold on
% scatter(x(1,:),x(2,:))
% numInHull=length(idx);
% for k=1:(numInHull-1)
%     plot([x(1,idx(k));x(1,idx(k+1))],[x(2,idx(k));x(2,idx(k+1))],'-b')
% end
% plot([x(1,idx(numInHull));x(1,idx(1))],[x(2,idx(numInHull));x(2,idx(1))],'-b')
%
%REFERENCES:
%[1] R. A. Jarvis, "On the identification of the convex hull of a finite
%    set of points in the plane," Information Processing Letters, vol. 2,
%    no. 1, pp. 18-21, Mar. 1973.
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPoints=size(points,2);

%First, find the lowest-most point (smallest first dimension).
lowestPoint=findLowestPoint(points);

hullIdx=zeros(numPoints,1);
%Add the leftmost point to the hull.
numInHull=0;
p=lowestPoint;
while(1)
    numInHull=numInHull+1;
    hullIdx(numInHull)=p;
    
    q=p+1;
    if(q>numPoints)
        q=1;
    end
    
    for k=1:numPoints
        if(k==p||k==q)
            continue;
        end
        if(turnOrientation(points(:,p),points(:,k),points(:,q))==1)
            %If point k is more counterclockwise than point q.
            q=k;
        end
    end
    
    %q is the most counterclockwise point from p, so select it as the
    %current point, which will get it added to the hull the next time
    %around the loop.
    p=q;
    
    %If we have gone all the way around the hull of points.
    if(p==lowestPoint)
        break;
    end
end

%Size to fit.
hullIdx=hullIdx(1:numInHull);

end

function lowestIdx=findLowestPoint(points)
%%Find the point with the smallest y value. In case of ties, choose the one
%%with the largest x value. Swap this point with the beginning of the list
%%of points.
    numPoints=size(points,2);

    lowestIdx=1;
    for i=2:numPoints
        if((points(2,i)<points(2,lowestIdx))||((points(2,i)==points(2,lowestIdx))&&(points(1,i)>points(1,lowestIdx))))
            lowestIdx=i;
        end
    end
end

