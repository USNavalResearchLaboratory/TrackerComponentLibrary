function [isInPolygon,omegas]=pointIsInPolygon(vertices,point,boundaryIsImportant)
%%POINTISINPOLYGON Given a polygon specified by a number of vertices,
%                determine whether a point is in the polygon. A simple
%                polygon has lines that do not cross. If the polygon is not
%                simple (is self-intersecting), then the even-odd rule is
%                used to determine whether a point is in the polygon or
%                not. As the algorithm is a type of winding number
%                algorithm, the winding number omega can also be obtained.
%
%INPUTS: vertices A 2XN matrix of the N vertices making up the polygon in
%                 order. Edges are between neighboring vertices. The last
%                 vertex can be the same as the first vertex. If not, it is
%                 assumed that an edge exists between the last and first
%                 vertices.
%           point A 2XnumPoints set of points that will be determined to be
%                 inside or outside of the polygon given by vertices.
% boundaryIsImportant An optional boolean variable indicating whether the
%                 boundary of the polygon is important. If
%                 boundaryIsImportant=true, then an algorithm that will
%                 correctly indicate points on the boundary as being in the
%                 polygon will be used. If it is false, then the results
%                 for points on the boundary can be inconsistent, though
%                 the algorithm will be slightly faster. The default is
%                 true.
%
%OUTPUTS: isInPolygon An NX1 vector where the ith element is true if the
%                     ith is true if the ith point is in the polygon, false
%                     otherwise. If boundaryIsImportant=false, the results
%                     can be inconsistent on the boundary of the polygon.
%              omegas An NX1 vector of the integer winding number obtained
%                     for each point. This is not meainingful for points on
%                     the boundary of the polygon. 
%
%The implementation is Algorithms 6 and 7 from [1].
%
%REFERENCES:
%[1] K. Hormann and A. Agathos, "The point in polygon problem for arbitrary
%    polygons," Computational Geometry, vol. 20, no. 3, pp. 131-144, Nov.
%    2001.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    boundaryIsImportant=true;
end
P=vertices;
%Put the first point on the end to close the polygon. It does not matter if
%the user has already done this.
P(:,end+1)=P(:,1);
%The number of vertices (It ends at the point at which it starts).
n=size(P,2)-1;
numPoints=size(point,2);

%Allocate space.
isInPolygon=zeros(numPoints,1);
omegas=zeros(numPoints,1);
for curPoint=1:numPoints
    R=point(:,curPoint);
    
    %omega is the winding number that we are finding.
    omega=0;
    if(boundaryIsImportant==false)
        %This is Algorithm 6, which is faster, but which does not return
        %consistent results for points that are on the boundary.
        for i=1:n
            if((P(2,i)<R(2))~=(P(2,i+1)<R(2)))
            %If crossing
                if(P(1,i)>=R(1))
                    if(P(1,i+1)>R(1))
                        %Modify_omega
                        omega=omega+2*(P(2,i+1)>P(2,i))-1;
                    else
                        detVal=(P(1,i)-R(1))*(P(2,i+1)-R(2))-(P(1,i+1)-R(1))*(P(2,i)-R(2));
                        if((detVal>0)==(P(2,i+1)>P(2,i)))
                        %If right_crossing
                            %Modify_omega
                            omega=omega+2*(P(2,i+1)>P(2,i))-1;
                        end
                    end
                else
                    if(P(1,i+1)>R(1))
                        detVal=(P(1,i)-R(1))*(P(2,i+1)-R(2))-(P(1,i+1)-R(1))*(P(2,i)-R(2));
                        if((detVal>0)==(P(2,i+1)>P(2,i)))
                        %If right_crossing
                            %Modify_omega
                            omega=omega+2*(P(2,i+1)>P(2,i))-1;
                        end
                    end
                end
            end
        end
        omegas(curPoint)=omega;
        isInPolygon(curPoint)=(mod(omega,2)~=0);
    else
        %This is algorithm 7 which is slower, but which correctly identifies
        %points on the boundary as being inside the polygon.
        isEdgeOrVertex=false;
        if((P(1,1)==R(1))&&(P(2,1)==R(2)))
            %If it is on the first vertex
            isEdgeOrVertex=true;
        else
            for i=1:n
                if(P(2,i+1)==R(2))
                    if(P(1,i+1)==R(1))
                        %If it is on a vertex
                        isEdgeOrVertex=true;
                        break;
                    else
                        if((P(2,i)==R(2))&&((P(1,i+1)>R(1))==(P(1,i)<R(1))))
                            %If it is on an edge
                            isEdgeOrVertex=true;
                            break;
                        end
                    end
                end

                if((P(2,i)<R(2))~=(P(2,i+1)<R(2)))
                %If crossing
                    if(P(1,i)>=R(1))
                        if(P(1,i+1)>R(1))
                            %Modify_omega
                            omega=omega+2*(P(2,i+1)>P(2,i))-1;
                        else
                            detVal=(P(1,i)-R(1))*(P(2,i+1)-R(2))-(P(1,i+1)-R(1))*(P(2,i)-R(2));
                            if(detVal==0&&~all(P(:,i)==P(:,i+1)))
                                isEdgeOrVertex=true;
                                break;
                            end
                            if((detVal>0)==(P(2,i+1)>P(2,i)))
                            %If right_crossing
                                %Modify_omega
                                omega=omega+2*(P(2,i+1)>P(2,i))-1;
                            end
                        end
                    else
                        if(P(1,i+1)>R(1))
                            detVal=(P(1,i)-R(1))*(P(2,i+1)-R(2))-(P(1,i+1)-R(1))*(P(2,i)-R(2));
                            if(detVal==0&&~all(P(:,i)==P(:,i+1)))
                                isEdgeOrVertex=true;
                                break;
                            end
                            if((detVal>0)==(P(2,i+1)>P(2,i)))
                            %If right_crossing
                                %Modify_omega
                                omega=omega+2*(P(2,i+1)>P(2,i))-1;
                            end
                        end
                    end
                end  
            end
        end
        omegas(curPoint)=omega;
        isInPolygon(curPoint)=(isEdgeOrVertex==true||(mod(omega,2)~=0));
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
