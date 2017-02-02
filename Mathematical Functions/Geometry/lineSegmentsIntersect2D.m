function boolVal=lineSegmentsIntersect2D(line1,line2)
%%LINESEGMENTSINTERSECT2D Determine whether two 2D line segments intersect.
%                    The intersection is determined without actually
%                    finding the point of intersection.
%
%INPUTS: line1  A 2X2 matrix of two points in the first line where
%               line1(1,:) are the x-coordinates and line1(2,:) are the
%               y-coordinates.
%        line2  A 2X2 matrix of two points in the second line defined the
%               same way as line1.
%
%The algorithm of Chapter 1.5 of [1] is used. It does not compute slopes
%and thus avoid issues with infinite sloped when lines are vertical.
%
%REFERENCES:
%[1] J. O'Rourke, Computational Geometry in C, 2nd ed. Cambridge, United
%    Kingdom: Cambridge University Press, 1998.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[boolVal,isImproper]=lineSegmentsIntersectProper(line1,line2);

if(boolVal)
   return; 
end

if(isImproper)
  if(cIsBetweenaAndb(line1(:,1),line1(:,2),line2(:,1))...
     ||cIsBetweenaAndb(line1(:,1),line1(:,2),line2(:,2))...
     ||cIsBetweenaAndb(line1(:,1),line2(:,1),line2(:,2))...
     ||cIsBetweenaAndb(line1(:,2),line2(:,1),line2(:,2)))
    boolVal=true;
  end
else
    boolVal=false;
end

end

function [boolVal,isImproper]=lineSegmentsIntersectProper(line1,line2)
%Determine whether the line segments have a proper intersections. That is,
%no segment has an endpoint on the other.

a=line1(:,1);
b=line1(:,2);
c=line2(:,1);
d=line2(:,2);


%Eliminate improper cases, which must be checked separately.
if(pointsAreCollinear(a,b,c)||pointsAreCollinear(a,b,d)||pointsAreCollinear(c,d,a)||pointsAreCollinear(c,d,b))
     boolVal=false;
     isImproper=true;
     return;
end

isImproper=false;
boolVal=xor(turnOrientation(a,b,c)==1,turnOrientation(a,b,d)==1)&&xor(turnOrientation(c,d,a)==1,turnOrientation(c,d,b)==1);
end

function boolVal=cIsBetweenaAndb(a,b,c)
%Determine whether point c is between points a and b.

%Assume that the points are collinear, so that the following commented-out
%statement would never execute.
%     if(~pointsAreCollinear(a,b,c))
%         boolVal=false;
%         return;
%     end

%If segment a-b is not vertical, check the betweenness on x; otherwise
%check it on y.
    if(a(1)~=b(1))
        boolVal=((a(1)<=c(1))&&(c(1)<=b(1)))||((a(1)>=c(1))&&(c(1)>=b(1)));
    else
        boolVal=((a(2)<=c(2))&&(c(2)<=b(2)))||((a(2)>=c(2))&&(c(2)>=b(2)));
    end
end

function boolVal=pointsAreCollinear(a,b,c)
%The following is the same as
%boolVal=(b(1)-a(1))*(c(2)-a(2))-(c(1)-a(1))*(b(2)-a(2))==0
%but should be less susceptible to finite preciion errors.

S=[-a(2)*b(1);
    a(1)*b(2);
   -a(1)*c(2);
    a(2)*c(1);
   -b(2)*c(1);
    b(1)*c(2)];

boolVal=exactSignOfSum(S)==0;
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
