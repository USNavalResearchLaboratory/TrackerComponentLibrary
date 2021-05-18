function vertices=findConvexHull2D(points)
%%FINDCONVEXHULL2D Given a set of points, find the smallest convex polygon
%                  enclosing the points in 2D (The convex hull). A form of
%                  Graham's algorithm, which has O(nlog(n)) complexity, is
%                  used. The vertices of the polygon (which coincide with
%                  points) are returned in counterlockwise order.
%
%INPUTS: points A 2XN set of N finite points of the form [x;y], the convex
%               hull about which one desires. Repeated points are allowed.
%
%OUTPUTS: vertices A 2XnumVert set of numVert vertex points of the convex
%                  hull polygon given in order. The first point is not
%                  repeated at the end even though it is a closed shape.
%                  The points are provided in counterclockwise order
%
%Graham's algorithm is described in Chapter 33.3 of[1], though no method
%for sorting the points in the necessary counterclockwise order is given
%and the sorting step is the most difficult part. The sorting can be done
%without explicitly computing angles as discussed in Chapter 3.5 of [2].
%However, the author suggests that one only use integers to avoid horrible
%finite precision problems. Here, the method of determining whether points
%are oriented  left uses the function turnOrientation, which makes use of
%the error-free function exactSignOfSum. Thus, as long as no precision is
%lost in any intermediate step, then the results are exact. This means that
%if points are floating point doubles that have been truncated to single
%precision (points cannot be singles or integer types for this function),
%then the results will be exact. Given doubles, steps were taken in the
%turnOrientation function to reduce the likelihood of finite precision
%errors. However, such errors cannot be completely eliminated without
%using quadruple precision arithmetic in an intermediate step.
%
%The implementation discussed in Cormen's book is simpler than Graham's
%original implementation, which is given in [3], as it does not use complex
%numbers. The implementation here does not use complex numbers and none of
%the points can be complex.
%
%EXAMPLE:
% numPts=1000;
% x=randn(2,numPts);
% xHull=findConvexHull2D(x);
% 
% figure(1)
% clf
% hold on
% scatter(x(1,:),x(2,:))
% numInHull=size(xHull,2);
% for k=1:(numInHull-1)
%     plot([xHull(1,k);xHull(1,k+1)],[xHull(2,k);xHull(2,k+1)],'-b')
% end
% plot([xHull(1,numInHull);xHull(1,1)],[xHull(2,numInHull);xHull(2,1)],'-b')
%
%REFERENCES:
%[1]T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein, Introduction
%   to Algorithms, 2nd ed. Cambridge, MA: The MIT Press, 2001.
%[2] J. O'Rourke, Computational Geometry in C, 2nd ed. Cambridge, United
%    Kingdom: Cambridge University Press, 1998.
%[3] R. L. Graham, "An efficient algorithm for determining the convex hull
%    of a finite planar set," Information Processing Letters, vol. 1, no.
%    4, pp. 132-133, Jun. 1972.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(isempty(points))
        vertices=[];
        return;
    end

    %First, rearrange things so that the point with the smallest y value is
    %at the beginning of the array. In case of ties in the y-coordinate,
    %the point with the largest x value is chosen. Thus, the angle from any
    %given point through the first point with respect to a parallel to the
    %x-axis going through the first point wil be in the range 0 to pi.
    points=findLowestPoint(points);

    numPoints=size(points,2);
    %If only two points are provided, then the convex hull is either just
    %the two points or just one of the points if they are equal. In this
    %instance, the  findLowestPoint function would have put the points in
    %the proper order from right to left.
    if(numPoints<=2)
        if(all(p0==points(:,2)))
            vertices=points(:,1);
        else
            vertices=points;
        end
        return
    end
    
    %Now, sort everything going counterclockwise from the first point in
    %points.
    gtCompare=@(p1,p2)comparePoints(p1,p2,points(:,1));
    points(:,2:end)=heapSort(points(:,2:end),true,gtCompare);

    %Once the points are sorted, scan through looking for ones that are at
    %the same angle with respect to points(:,1) and the parallel to the
    %x-axis going through points(:,1). Only the most distant point is
    %kept when two points are at the same angle. When two points are equal,
    %then keep the one with the lowest index as a convention.
    selDelete=findDuplicatePointsToRemove(points);
    %Get rid of the points that should be removed.
    points=points(:,~selDelete);
    
    %We have to make sure that the problem was not degenerate after points
    %were removed. That is, that we still have enough points to enter the
    %loop.
    numPoints=size(points,2);
    if(numPoints<=2)
        vertices=points;
        return
    end
    
%Now the ordering of the points has completed and Graham's algorithm for
%finding the convex hull of the points sorted in a counterclockwise manner
%about the lowest point points(:,1) can be used.
    
    %Allocate space for the stack.
    stack=zeros(2,numPoints);

    %Push the first three things onto the stack.
    stack(:,1:3)=points(:,1:3);
    curStackIdx=3;
    for curIdx=3:numPoints
        pCur=points(:,curIdx);
        %While the orientation of the turn is counterclockwise
        while(turnOrientation(stack(:,curStackIdx-1),stack(:,curStackIdx),pCur)~=1)
            curStackIdx=curStackIdx-1;%Pop the point off of the stack.
            
            %The first two things in the stack should never be removed.
            %In the unlikely event that finite precision errors cause it to
            %want to remove one of them, throw an error. This should never
            %happen if the points all are truncated to single precision but
            %extended to doubles), because the turnOrientation function
            %should be exact and that is the deciding factor in this
            %algorithm.
            if(curStackIdx<2)
               error('Finite precision errors prevent the computation of the convex hull')
            end
        end
        curStackIdx=curStackIdx+1;
        %Push point pCur onto the stack.
        stack(:,curStackIdx)=pCur;
    end
    
    vertices=stack(:,1:curStackIdx);
end

function points=findLowestPoint(points)
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

    %Swap the first point in the list with the item at lowestIdx.
    temp=points(:,1);
    points(:,1)=points(:,lowestIdx);
    points(:,lowestIdx)=temp;
end

function compVal=comparePoints(p1,p2,p0)
%%Returns 1 if p1 is at a greater counterclockwise angle with respect to p0
%and the x-axis going through p0 than pj, and 0 if it is a smaller angle or
%if the points are equal.

%The orientation of the turn from p0->p1->p2. If negative, then p1 comes
%before p2 (less than). If positive, then p1 comes after p2 (greater than).
%If zero, then they are collinear or equal.
A=turnOrientation(p0,p1,p2);

if(A>0)
    compVal=0;%Less than
elseif(A<0)
    compVal=1;%Greater than
else%The points are collinear. The one that is closer is "less than".
    r1=norm(p1-p0);
    r2=norm(p2-p0);
    
    if(r1<r2)
        compVal=0;%Less than
    else%Greater than or equal to
        compVal=1; 
    end
end
end

function selDelete=findDuplicatePointsToRemove(points)
    numPoints=size(points,2);
    selDelete=false(numPoints,1);
    
    p0=points(:,1);
    
    %The first point after p0 will stay unless it is equal to p0.
    if(all(p0==points(:,2)))
        selDelete(2)=true;
    end
    
    for i=2:(numPoints-1)
        p1=points(:,i);
        p2=points(:,i+1);
        
        %If the points are collinear
        if(turnOrientation(p0,p1,p2)==0)
            r1=norm(p1-p0);
            r2=norm(p2-p0);
            
            %Mark the closer point for deletion.
            if(r1<r2)
                selDelete(i)=true;
            else%r2 is either closer or the points are equal. The
            %convention for equal points is to keep the one with the lowest
            %index, so either way, p2 is marked for deletion.
                selDelete(i+1)=true;
            end
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
