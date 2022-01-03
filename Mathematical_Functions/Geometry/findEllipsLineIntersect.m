function [xInt,t]=findEllipsLineIntersect(xc,A,x0,a,aType)
%%FINDELLIPSLINEINTERSECT Given the parameters for an ellipsoid (in 2D an
%       ellipse) of the form (x-xc)'*A*(x-xc)=1 and a line either given
%       parameterically as x=x0+a*t or whether x0 and a are two points on
%       the line, find the points of intersection between the ellipsoid and
%       the line.
%
%INPUTS: xc The xDimX1 center of the ellipsoid.
%         A The xDimXxDim symmetric matrix defining the above equation for
%           an ellipsoid.
%     x0, a These xDimX1 values define a line and the type of line is given
%           by aType. If aType=0, then the line is parameteric as x=x0+a*t
%           where t is a scalar value. If aType=1, then x0 and a are both
%           points on the line. 
%     aType As mentioned above, this specifies how the line is
%           parameterized. The default if omitted or an empty matrix is
%           passed is 0.
%
%OUTPUTS: xInt Either xDimX2 (or rarely xDimX1) matrix holding the vector
%              points of intersection of the line and the ellipsoid or an
%              empty matrix if there are no real solutions (the line does
%              not intersect the ellipsoid).
%            t If solutions exist, this is the 2X1 (or rarely 1X1) set of
%              parametric parameters. If aType=0, then these are the t
%              values in the parameteric equation. If aType=1, then these
%              are the t values in the parameteric equation x=x0+(a-x0)*t,
%              in which case values between 0 and 1 indicate that the
%              intersection occured between the two specified points.
%
%Start with the equations x=x0+a*t and (x-xc)'*A*(x-xc)=1. Substitite the
%first into the second and one ends up with a quadratic equation
%c1*t^2+c2*t+c3=0 with c1=a'*A*a, c2=2*xTilde'*A*a, c3=xTilde'*A*xTile-1
%where xTilde=x0-xc.
%
%EXAMPLE:
%Find and plot the intersection points of a line an an ellipse. Also
%consider a line that does not intersect and one sees that nothing is
%plotted from that, because no solution exists.
% %The ellipse.
% xc=[1;2];
% A=[1,0.5;
%    0.5,2];
% %Two points defining a line that does intersect.
% x01=[0;3];
% x02=[2;1];
% %Two points defining a line that does not intersect.
% x11=[3;3];
% x12=[3;1];
% 
% aType=1;
% xInt0=findEllipsLineIntersect(xc,A,x01,x02,aType);
% xInt1=findEllipsLineIntersect(xc,A,x11,x12,aType);%This is empty.
% 
% figure(1)
% clf
% hold on
% drawEllipse(xc,A,1,'linewidth',2)
% plot([x01(1);x02(1)],[x01(2);x02(2)],'-k','linewidth',2)
% plot([x11(1);x12(1)],[x11(2);x12(2)],'-k','linewidth',2)
% 
% numSol0=size(xInt0,2);
% numSol1=size(xInt1,2);%Should be 0.
% for k=1:numSol0
%     scatter(xInt0(1,k),xInt0(2,k),400,'.r')
% end
% for k=1:numSol1
%     scatter(xInt1(1,k),xInt1(2,k),400,'.g')
% end
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(aType))
    aType=0;
end

if(aType==1)
    %Make it parameteric.
    a=a-x0; 
end

xTilde=x0-xc;

c1=a'*A*a;
c2=2*xTilde'*A*a;
c3=xTilde'*A*xTilde-1;

t=roots([c1;c2;c3]);

if(isempty(t)||any(~isreal(t)))
    t=[];
    xInt=[];
    return 
end

xDim=size(x0,1);
numSol=length(t);

xInt=zeros(xDim,numSol,1);
for curSol=1:numSol
    xInt(:,curSol)=x0+a*t(curSol);
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
