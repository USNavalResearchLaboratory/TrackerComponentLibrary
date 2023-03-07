function [xAL,xBL,doesClip]=clipLine2ConvexPolygon(xAL,xBL,v)
%%CLIPLINE2CONVEXPOLYGON Given the vertices of a convex polygon in 2D in
%       either clockwise or counterclockwise order, clip one or more lines
%       to the polygon. The lines are defined by two points on the lines,
%       though they extend beyond the two points.
%
%INPUTS: xAL, xBL These are 2 2XnumLines matrices where xAL(:,k) and
%                 xBL(:,k) are two points on the kth line.
%               v A 2XnumVertices set of vertices defining the perimeter of
%                 the convex polygon. The last vertex should not be a
%                 repeat of the first vertex.
%
%OUTPUTS: xAL, xBL These are 2XnumLines set of the starting and ending
%                  points of the lines clipped to the polygon. If the kth
%                  line is complete outside the polygon, then xAL(:,k) and
%                  xBL(:,k) will be NaNs.
%         doesClip A numLinesX1 boolean vector indicating whether or not
%                  each line is in or on the polygon at all.
%
%EXAMPLE:
%Draw a polygon and a few line segments in black. The line segments
%represent infinte lines. Clip the lines to the polygon and draw the
%clipped lines in red. One of the lines does not clip, so a black segment
%is drawn, but nothing is drawn in red.
% numVertices=5;
% v=zeros(2,numVertices);
% v(:,1)=[1;1];
% v(:,2)=[15;0];
% v(:,3)=[10;10];
% v(:,4)=[7;10];
% v(:,5)=[0;8];
% 
% %Two line segments.
% xA=[[-3;10],[0;0],[4;6],[-2;1]];
% xB=[[12;10],[16;8],[2;7],[-3;3]];
% [xAT,xBT]=clipLine2ConvexPolygon(xA,xB,v);
% 
% figure(1)
% clf
% hold on
% plot([v(1,:),v(1,1)],[v(2,:),v(2,1)],'linewidth',2)
% for k=1:size(xA,2)
%     plot([xA(1,k),xB(1,k)],[xA(2,k),xB(2,k)],'-k','linewidth',4)
%     
%     if(~isnan(xAT(1,k)))
%         plot([xAT(1,k),xBT(1,k)],[xAT(2,k),xBT(2,k)],'-r','linewidth',2)
%     end
% end
%
%REFERENCES:
%[1] V. Skala, "An efficient algorithm for line clipping by convex
%    polygon," Computer and Graphics, vol. 17, no. 4, pp. 417-421, 1993.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(v,2);%N is the number of edges and vertices.

%If many lines are clipped at once, precomputing this outside the loop can
%be fastest.
siHat=zeros(2,N);
for k=1:(N-1)
    siHat(:,k)=v(:,k+1)-v(:,k);
end
siHat(:,N)=v(:,1)-v(:,N);

numPts=size(xAL,2);
doesClip=false(numPts,1);

for curPt=1:numPts
    xA=xAL(:,curPt);
    xB=xBL(:,curPt);

    si=bsxfun(@minus,v,xA);
    index=zeros(2,1);
    prevSpecialIdx=-2;

    k=0;
    i=N-1;
    j=0;
    s=xB-xA;
    xi=(si(2,i+1)*s(1)-si(1,i+1)*s(2));
    while((j<N)&&(k<2))
        eta=(si(2,j+1)*s(1)-si(1,j+1)*s(2));
        if(xi*eta<0)%There is an intersection
            index(k+1)=i+1;%The edge having the intersection.
            k=k+1;
            xi=eta;
        elseif(xi*eta==0)%xi and/or eta=0
            if(~(prevSpecialIdx==N-1&&i==0)&&~(prevSpecialIdx==i-1))        
                if(eta==0)
                    if(xi~=0)%If not special case m.
                        index(k+1)=i+1;%The edge having the intersection.
                        k=k+1;
                        xi=eta;
                        prevSpecialIdx=i;
                    end
                else%xi==0 and eta is not zero.
                    index(k+1)=i+1;%The edge having the intersection.
                    k=k+1;
                    xi=eta;
                    prevSpecialIdx=i;
                end
            end
        end
        i=j;
        j=j+1;
    end

    if(k==0)
        doesClip(curPt)=false;%There is no intersection.
        xA=NaN;
        xB=NaN;
    else
        doesClip(curPt)=true;%There is an intersection.
        tMin=-Inf;
        tMax=Inf;
        i=index(1);
        t1=det([si(:,i),-siHat(:,i)])/det([s,-siHat(:,i)]);
        i=index(2);
        t2=det([si(:,i),-siHat(:,i)])/det([s,-siHat(:,i)]);

        if(t1<t2)
           temp=t2;
           t2=t1;
           t1=temp;
        end

        if(t2<tMax)
            xB=xA+s*t2;
        end

        if(t1>tMin)
            xA=xA+s*t1;
        end
    end

    xAL(:,curPt)=xA;
    xBL(:,curPt)=xB;
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
