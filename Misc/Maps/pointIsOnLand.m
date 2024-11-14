function [isOnLand,vertexIdx]=pointIsOnLand(latLonPointDeg,mapData)
%%POINTISONLAND Given a point in terms of WGS-84 latitude and longitude in
%    degrees, determine whether it is on land. This uses the Global Self-
%    consistent, Hierarchical, High-resolution Geography Database (GSHHG),
%    developed by the University of Hawai'i, Honolulu and the National
%    Oceanic and Atmospheric Administration (NOAA) as returned by the
%    getGSHHGBoundaryData function.
%
%INPUTS: latLonPointDeg A 2XnumPoints set of [latitude;longitude;] points
%                in degrees to test if they are on land. The points must be
%                within the bounding box used when calling
%                getGSHHGBoundaryData.
%        mapData This is the mapData structure returned by the
%                getGSHHGBoundaryData function; see the comments to that
%                function for more details on the format. This function
%                supports dataType=0 maps from that function for
%                determining if a point is on land, or assuming the
%                point is on land, one can pass dataType=1 data to
%                determine whether or not the point is on a navigable river
%                (so isOnLand=false) or is not on the river
%                (isOnLand=true).
%
%OUTPUTS: isOnLand An NX1 boolean vector indicating whether each point is
%                  on land or in a body of water given in the data passed.
%        vertexIdx An NX1 vector such that if the ith point is on something
%                  (lake, island, land), vertexIdx(i) is the index of the
%                  polygon in vertices on which it resides. Note that a
%                  point in a lake will have isOnLand=false. If the point
%                  is in not on land or in some body of water/ island
%                  therein/ pond therein on land, then the corresponding
%                  vertexIdx will be zero.
%
%The header information in contains bounding rectangles for each of the
%polygons in vertices. If the point is within the bounding rectangle, then
%the pointIsInPolygon function is used to see if it is in the polygon. Only
%top-level (e.g. land, not island or lake) polygons are searched prior to
%declaring a point definitely not on land with vertexIdx empty as all
%lower-level polygons reside within the top-level ones. If a point is found
%in a top level polygon, then the lower-levels ones are hierarchically
%searched to determine whether it is one land, in a lake, etc.
%
%EXAMPLE 1:
%This initial example plots Europe without lakes and tests a grid of points
%for being on land or not. Those points that are on land are displayed in
%magents and those in the water are cyan.
% latLonRecDeg=[30;60;-30;30];
% figure(1)
% clf
% hold on
% dataType=0;
% resolutionLevel=1;
% excludeTypes=[2;3;4;6];
% mapData=getGSHHGBoundaryData(latLonRecDeg,dataType,resolutionLevel,excludeTypes);
% plotGSHHGBoundaryMap(mapData)
% numLat=51;
% numLon=50;
% latPts=linspaceNoEnd(latLonRecDeg(1),latLonRecDeg(2),numLat,true);
% lonPts=linspaceNoEnd(latLonRecDeg(3),latLonRecDeg(4),numLon,true);
% [Lat,Lon]=meshgrid(latPts,lonPts);
% ll=[Lat(:).';Lon(:).'];
% isOnLand=pointIsOnLand(ll,mapData);
% llOnLand=ll(:,isOnLand);
% llNotOnLand=ll(:,~isOnLand);
% scatter(llOnLand(2,:),llOnLand(1,:),'.m')
% scatter(llNotOnLand(2,:),llNotOnLand(1,:),'.c')
%
%EXAMPLE 2:
%This example is near the great lakes, so all the lake data is loaded and
%a grid of points is tested. Those in magenta are on land and those in cyan
%are on water.
% latLonRecDeg=[42;46.5;-88;-78];
% figure(1)
% clf
% hold on
% dataType=0;
% resolutionLevel=2;
% excludeTypes=[];
% mapData=getGSHHGBoundaryData(latLonRecDeg,dataType,resolutionLevel,excludeTypes);
% plotGSHHGBoundaryMap(mapData)
% numLat=51;
% numLon=50;
% latPts=linspaceNoEnd(latLonRecDeg(1),latLonRecDeg(2),numLat,true);
% lonPts=linspaceNoEnd(latLonRecDeg(3),latLonRecDeg(4),numLon,true);
% [Lat,Lon]=meshgrid(latPts,lonPts);
% ll=[Lat(:).';Lon(:).'];
% isOnLand=pointIsOnLand(ll,mapData);
% llOnLand=ll(:,isOnLand);
% llNotOnLand=ll(:,~isOnLand);
% scatter(llOnLand(2,:),llOnLand(1,:),'.m')
% scatter(llNotOnLand(2,:),llNotOnLand(1,:),'.c')
%
%April 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

vertices=mapData.vertices;
headerInfo=mapData.headerInfo;
shapeInfo=mapData.shapeInfo;
childStructureInfo=mapData.childStructureInfo;

numPoints=size(latLonPointDeg,2);

if(numPoints==0)
    isOnLand=[];
    vertexIdx=[];
    return;
end

isOnLand=false(numPoints,1);
vertexIdx=zeros(numPoints,1);

for curPoint=1:numPoints
    [isOnLand(curPoint), vertexIdx(curPoint)]=onePointIsOnLand(latLonPointDeg(:,curPoint),vertices,headerInfo,shapeInfo,childStructureInfo);
end

end

function [isOnLand, vertexIdx]=onePointIsOnLand(latLonPoint,vertices,headerInfo,shapeInfo,childStructureInfo)
%We only need to search the parent polygons before declaring a point to not
%be on land, because it cannot be in a child polygon (e.g. a lake) unless
%it is in the parent polygon.
if(~isempty(childStructureInfo))
    %If we are given land polygons rather than navigable river polygons.
    numPoly=childStructureInfo.numRootParents;

    for curPoly=1:numPoly
        isInPolygon=pointIsInGeographicPolygon(latLonPoint,shapeInfo(1:4,curPoly),vertices{curPoly});
        
        %If the point is in a parent polygon, we have to search all of the
        %child polygons, to see if it is actually on land or if it is in a lake or
        %on an island in a lake, etc.
        if(isInPolygon)
            curLevel=1;
            curParentIdx=curPoly;
            while(curLevel<childStructureInfo.maxChildLevel)
                numChildren=childStructureInfo.numChildList(curParentIdx);
                if(numChildren==0)
                    break;
                end
                
                %If it has children, look at them. If one finds a child in
                %which the point resides, then move down one level to
                %investigate the children of that child.
                curChildIdx=childStructureInfo.firstChildIdx(curParentIdx);
                for curChild=1:numChildren
                    isInChild=pointIsInGeographicPolygon(latLonPoint,shapeInfo(1:4,curChildIdx),vertices{curChildIdx});
                    
                    if(isInChild)
                         break;
                    end
                    
                    curChildIdx=curChildIdx+1;
                end
                
                %If it is in the child polygon, then go down another level to
                %investigate children of the child.
                if(isInChild)
                    curParentIdx=curChildIdx;
                    curLevel=curLevel+1;
                else
                %Otherwise, we have reached the deepest level in which the
                %point resides.
                    if(curLevel==1||curLevel==5||curLevel==6)
                        isOnLand=true;%Root levels are land.
                    else%For other levels, it is just odd is water, even is land.
                        isOnLand=mod(curLevel,2);
                    end
                    vertexIdx=curParentIdx;
                    return;
                end
            end
            %If we get here, then the point is not in one or more sublevels, or
            %we have reched the bottom sublevel. We will determine whether it
            %is on land based on what the last parent is.
            curLevel=headerInfo(3,curParentIdx);
            if(curLevel==1||curLevel==5||curLevel==6)
                isOnLand=true;
            else%For other levels, it is just odd is water, even is land.
                isOnLand=mod(curLevel,2);
            end
            
            vertexIdx=curParentIdx;
            return;
        end
    end

    isOnLand=false;
    vertexIdx=0;
else
    %If we are assumed to have already tested the point and found it to be
    %on land and we are calling this function a second time with rivers.
    numPoly=length(vertices);
    for curPoly=1:numPoly
        if(iscell(vertices{curPoly}))
            %Not a big navigable river.
            continue;
        end
        isInPolygon=pointIsInGeographicPolygon(latLonPoint,shapeInfo(1:4,curPoly),vertices{curPoly});
        if(isInPolygon)
            isOnLand=false;
            vertexIdx=curPoly;
            return
        end
    end
    %It is not on any of the navigable rivers, so we shall assume it is on
    %land (assume before calling this function, being on land had already
    %been tested).
    isOnLand=true;
    vertexIdx=0;
end
end

function isInPolygon=pointIsInGeographicPolygon(latLonPoint,boundRect,vertexArray)

    %First, see if the point, or the point aliased by 360 degrees in either
    %direction, falls in the bounding rectangle.
    isInBoundingBox=false;
    for curOffset=[0,-360,360]
        if(~((latLonPoint(1)>boundRect(2)||latLonPoint(1)<boundRect(1)||latLonPoint(2)+curOffset>boundRect(4)||latLonPoint(2)+curOffset<boundRect(3))))
            isInBoundingBox=true;
            break;
        end
    end

    %If no aliasing interval is in the rectangle, then it cannot be in the
    %polygon.
    if(isInBoundingBox==false)
        isInPolygon=false;
        return
    end

    %Now, we can check if the point is in the polygon.
    isInPolygon=pointIsInPolygon(vertexArray,latLonPoint,true);
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
