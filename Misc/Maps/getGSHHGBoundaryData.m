function mapData=getGSHHGBoundaryData(latLonRecDeg,dataType,resolutionLevel,excludeTypes)
%%GETGSHHGBOUNDARYDATA Read in geographic boundary data from the Global
%   Self-consistent, Hierarchical, High-resolution Geography Database
%   (GSHHG), developed by the University of Hawai'i, Honolulu and the
%   National Oceanic and Atmospheric Administration (NOAA). The data
%   consists of shorelines, rivers, and political boundaries on land. Note
%   that the political boundaries (i.e. country/ state borders) do not
%   include shoreline information, so a country will not look right unless
%   the shore is also drawn and the boundaries do nor form closed polygons.
%   This does not include ice near the North pole, which resularly changes
%   shape. The data loaded by this function can be used with the
%   plotGSHHGBoundaryMap and pointIsOnLand functions.
%
%INPUTS: latLonRec A 4X1 vector of [minLat;maxLat;minLon;maxLon] of WGS-84
%                 North latitude and East longitude in degrees indicating a
%                 rectangle in which data should be loaded. The latitudes
%                 must be between -90 and 90 degrees and the longitudes
%                 must start or end at a point between -180 and 180 degrees
%                 and cannot span more than 360 degrees. The ending
%                 longitude can be over 180 degrees. If one wants a window
%                 spanning the poles, the data should be loaded twice (one
%                 rectangle on either side). The polygons for land get
%                 clipped to the edges of the rectangle, but not broken
%                 into multiple polygons, so tests at the boundary of the
%                 rectangle will not be accurate as lines from clipped
%                 polygons will collect there.
%        dataType A value indicating which of the three data types to
%                 load. Possible values are:
%                 0 (The default if omitted or an empty matrix is passed)
%                   Load shoreline information.
%                 1 Load river information.
%                 2 Load political boundary information.
% resolutionLevel A value indicating the level of resolution of the data to
%                 load. Each step down in resolution represents an
%                 approximate 80% reduction in size and quality. Thus, high
%                 resolution is 80% less than full resolution and
%                 intermediate resolution is 80% less than high resolution.
%                 Possible values are:
%                 0 Load crude resolution data.
%                 1 (The default if omitted) Load low resolution data.
%                 2 Load intermediate resolution data.
%                 3 Load high resolution data.
%                 4 Load full resolution data.
%    excludeTypes An optional array indicating which hierarchical level
%                 types should be omitted. For example, if one might not
%                 want to include small details like the boundaries between
%                 a lake on an island which itself is in a lake on land.
%                 This is  an array of values of the hierarchical_level
%                 parameter, which is described below. An empty matrix
%                 means that nothing should be omitted. If the excludeTypes
%                 parameter itself is omitted, then if dataType=1 or 2,
%                 then nothing is omitted and if dataType=0, then
%                 excludeTypes=[2;3;4;6], which means that lakes and
%                 anything in them are omitted as well as the boundary
%                 between the ground of Antarctica and the ice, since the
%                 boundary between Antarctica's ice and the ocean should
%                 suffice.
%
%OUTPUTS: mapData This holds a data structure that can then be passed to
%           the plotGSHHGBoundaryMap and pointIsOnLand functions. The
%           members are latLonRecDeg (the same as the input), vertices,
%           headerInfo, and shapeInfo. The members are:
%           -headerInfo A 10XN matrix of parameters where the ith column
%            describes the ith set of vertices.  The values are saved
%            as int32 types rather than doubles to make it more efficient
%            to use the levels for indexation within the childStructureInfo
%            class. The rows are:
%            1) ID: A unique integer number identifying the polygon.
%            Note that this can be repeated if the polygone ended up
%            gatting split across the beginning and end of the window in
%            longitude compared to the GSHHG source data.
%            2) nPoints: The number of points in the polygon not counting
%            the repeated starting point for closed shapes.
%            3) hierarchical_level: This indicates the type of object being
%            traced out, such as a shoreline in the ocean, or a lake on
%            land, and is described in more detail below.
%            4) version This is the release number of the GSHHG data being
%            loaded. This should be the same for all data entries.
%            5) greenwich 0 If the original (non-split) polygon does not
%            cross the 0 degree longitude or the international data line.
%            1 if 0 degrees longitude is crossed, 2 if the international
%            dateline is crossed and 3 if both are crossed.
%            6) 0=not set; 1=river-lake and GSHHG level=2.
%            7) source: This indicates the underlying dataset that provided
%            the coastline data for the shape in the GSHGG. The values are
%            characters converted into int32 values. Thus, if one wants to
%            test for the letter W, then one should compare to int32('W').
%            The values are normally uppercase. However, if a lowercase
%            value is used, it means that the body of water being bounded
%            is a river-lake-- a part of a river that is so large it is
%            best represented by a closed polygon. Possible values for the
%            data sources are int32 versions of the ASCII characters
%            W,w World Vector Shorelines
%            C,c CIA World Data Bank II
%            A,a Atlas of the Cryosphere 
%            8) container: If dataType=0. This is the ID of the polygon
%            containing this polygon or a -1 if this polygon is not
%            contained by anything. For example, if this is the shoreline of
%            a lake, then this would be the ID of the land mass on which the
%            lake resides.
%            9) ancestor: If dataType=0 and resolutionLevel<4. This is the
%            ID of the polygon in the full-resolution dataset to which this
%            reduced resolution shape corresponds.
%            10) Split number. This starts at 1 the first time a polygone
%            with a particular ID is added. If the polygone has to be split
%            for some reason, then this is incremented.
%           -shapeInfo A 5XN matrix holding information regarding the
%            bounds and area of each polygon. The order of the entries is
%            the same as in headerInfo. The entries are stores as doubles.
%            1-4) South, North, West, East,: The North latitudes (South,
%            North) and East longitudes (West, East) that rectangularly
%            bound the shape. These are the actual bounds of the clipped
%            polygon.
%            5) f_area: If dataType=0 this is the area of the corresponding
%            (non-split) polygon in the full dataset.
%           -vertices A length N cell array of vertices of the N different
%            objects that were found. For dataType=0 and for dataType=1
%            when the thing in question is a river-lake, each entry is a
%            2XnumVertex set of points in WGS-84 [latitude;longitude] North
%            and East in degrees. For rivers as lines and political
%            boundaries as lines, each entry holds a length-2 cell array.
%            The first element in the array is a list of vertices in the
%            lines, but since the lines can be clipped by the bounding box,
%            the points represent multiple segments. the second cell entry
%            holds a (numSegments+1) vector of starting indices of the
%            segments, with the last element being the total number of
%            points.
%
%This function loads the data in the Global Self-consistent, Hierarchical,
%High-resolution Geography Database (GSHHG), which was formerly called the
%Global, Self-consistent, Hierarchical, High-resolution Shoreline Database
%(GSHHS). The original version of the database is described in [1].
%The latest version can be obtained from
%http://www.soest.hawaii.edu/pwessel/gshhg/
%The data that is loaded here is from the binary release of the GSHHG and
%the data file is assumed to be located in a "data" folder in the same
%folder as this function and the file is assumed to be named gshhg-bin.zip.
%That is, the files have not been unzipped from the distributed version,
%though the zip file will have had to be renamed as it is distributed with
%the version number in the name. This function has been tested using
%version 2.3.7.
%
%The values of the hierarchical_level parameter are described in 
%http://www.soest.hawaii.edu/pwessel/gshhg/
%and depend on the type of dataset being loaded. The possible values depend
%on whether data for shorelines, rivers, or borders are being loaded and
%are:
%Shorelines are organized into 6 hierarchival levels:
%1. L1: boundary between land and ocean, except Antarctica.
%2. L2: boundary between lake and land.
%3. L3: boundary between island-in-lake (on land) and lake.
%4. L4: boundary between pond-in-island and island (in a lake on land).
%5. L5: boundary between Antarctica ice and ocean.
%6. L6: boundary between Antarctica grounding-line and ocean.
%Rivers are organized into 11 classification levels:
%1. L0: Double-lined rivers (river-lakes).
%2. L1: Permanent major rivers.
%3. L2: Additional major rivers.
%4. L3: Additional rivers.
%5. L4: Minor rivers.
%6. L5: Intermittent rivers - major.
%7. L6: Intermittent rivers - additional.
%8. L7: Intermittent rivers - minor.
%9. L8: Major canals.
%10. L9: Minor canals.
%11. L10: Irrigation canals.
%Borders are organized into three levels
%1. L1: National boundaries.
%2. L2: State boundaries within the Americas.
%3. L3: Marine boundaries.
%
%Sheline data can be passed to pointIsOnLand to determine if a point is on
%land. If a point has been determined to be on land, one can then pass the
%river data to pointIsOnLand to determine if the point is on a river-lake
%(other rivers are not considered).
%
%Only Antarctica and the ice sheet associated with Antarctica do not form
%closed-shapes in latitude and longitude coordinates, but when they are
%loaded by this function, the shapes are closed.
% 
%REFERENCES:
%[1] P. Wessel and W. H. F. Smith, "A global, self-consistent,
%    hierarchical, high-resolution shoreline database," Journal of
%    Geophysical Research, vol. 101, no. B4, pp. 8741-8743, 10 Apr. 1996.
%
%April 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<1||isempty(latLonRecDeg))
        %Take everything.
        latLonRecDeg=[-90;90;-180;180];
    end

    if(nargin<2||isempty(dataType))
        dataType=0;
    end

    if(nargin<3||isempty(resolutionLevel))
        resolutionLevel=1;
    end

    if(nargin<4)
       if(dataType~=0)
           %Don't exclude anything when getting political boundaries or
           %rivers.
           excludeTypes=[];
       else
           %Exclude islands in lakes and ponds on islands in lakes when
           %by default when getting coastlines. Also exclude the boundary
           %between Antarctica's ground and ocean as the ice boundary
           %sufficies.
           excludeTypes=[2;3;4;6];
       end
    end
    
    clipToView=true;
    
    %[latitude;longitude];
    boundRectMin=[latLonRecDeg(1);
                  latLonRecDeg(3)];
    boundRectMax=[latLonRecDeg(2);
                  latLonRecDeg(4)];

    %These are the four vertices of the bounding rectangle in latitude
    %and longitude.
    boundRectPointsLatLon=[[latLonRecDeg(1);latLonRecDeg(3)],[latLonRecDeg(2);latLonRecDeg(3)],[latLonRecDeg(2);latLonRecDeg(4)],[latLonRecDeg(1);latLonRecDeg(4)]];

    switch(dataType)
        case 0
            switch(resolutionLevel)
                case 0
                    file2Decompress='gshhs_c.b';
                case 1
                    file2Decompress='gshhs_l.b';
                case 2
                    file2Decompress='gshhs_i.b';
                case 3
                    file2Decompress='gshhs_h.b';
                case 4
                    file2Decompress='gshhs_f.b';
                otherwise
                    error('Invalid resolution level specified')
            end
        case 1
            switch(resolutionLevel)
                case 0
                    file2Decompress='wdb_rivers_c.b';
                case 1
                    file2Decompress='wdb_rivers_l.b';
                case 2
                    file2Decompress='wdb_rivers_i.b';
                case 3
                    file2Decompress='wdb_rivers_h.b';
                case 4
                    file2Decompress='wdb_rivers_f.b';
                otherwise
                    error('Invalid resolution level specified')
            end
        case 2
            switch(resolutionLevel)
                case 0
                    file2Decompress='wdb_borders_c.b';
                case 1
                    file2Decompress='wdb_borders_l.b';
                case 2
                    file2Decompress='wdb_borders_i.b';
                case 3
                    file2Decompress='wdb_borders_h.b';
                case 4
                    file2Decompress='wdb_borders_f.b';
                otherwise
                    error('Invalid resolution level specified')
            end
        otherwise
            error('Invalid data type specified')
    end
    
    %Load the file
    ScriptPath=mfilename('fullpath');
    ScriptFolder=fileparts(ScriptPath);
    unzippedFiles=readZipArchive([ScriptFolder,'/data/gshhg-bin.zip'],file2Decompress);
    binaryData=unzippedFiles{2};
    numBytes=length(binaryData);
    
    %First, scan through the file to see how many polygons there are and
    %what the total number of points is.
    numPolygons=0;
    curStartByte=1;
    while(curStartByte<=numBytes)
        numPolygons=numPolygons+1;
        pointsInPolygon=typecast(binaryData((curStartByte+7):-1:(curStartByte+4)),'uint32');
        curStartByte=curStartByte+44;%The beginning of the points.
        %The points are 8 bytes each, so this is what we need to add to get
        %to the start of the next header.
        curStartByte=curStartByte+pointsInPolygon*8;
    end
    %Allocate space for the vertices.
    vertices=cell(numPolygons,1);
    %The *2 is an upper bound to deal with the fact that some polygons
    %might split in two across an aliased boundary region. Only Antarctica
    %will have numGatedOffsets=3 below, and we will only be using 1 piece.
    headerInfo=zeros(10,numPolygons*2,'int32');
    shapeInfo=zeros(5,numPolygons*2,'double');

    curPolyTaken=1;
    curStartByte=1;
    for curPoly=1:numPolygons
        %First, check whether the polygon is a level that we even want to
        %extract.
        pointsInPolygon=typecast(binaryData((curStartByte+7):-1:(curStartByte+4)),'uint32');
        flagInt=typecast(binaryData((curStartByte+11):-1:(curStartByte+8)),'uint32'); 
        level=bitand(flagInt,255);

        if(any(level==excludeTypes))
            %Go to the next header. We do not want to take polygons form
            %this level.
            curStartByte=curStartByte+44;
            curStartByte=curStartByte+pointsInPolygon*8;
            continue;
        end

        %Extract the bounds of the polygon. This forms a box. We will see
        %if this box overlaps with the window that we are trying to extract
        %stuff from.
        west=double(typecast(binaryData((curStartByte+15):-1:(curStartByte+12)),'int32'))/1e6;
        east=double(typecast(binaryData((curStartByte+19):-1:(curStartByte+16)),'int32'))/1e6;
        south=double(typecast(binaryData((curStartByte+23):-1:(curStartByte+20)),'int32'))/1e6;
        north=double(typecast(binaryData((curStartByte+27):-1:(curStartByte+24)),'int32'))/1e6;

        %The bounds if we consider aliasing in longitude.
        rectMin1=[south;
                   west];
        rectMax1=[north;
                   east];

        %Record all aliased offsets where an intersection occurs with the
        %bounding box. A landmass that spans both sides of an aliased
        %region will gate with more than one thing. (Antarctica is a
        %special case). If we clip, we are going to have to break the
        %landmass into two pieces (one for each side of the region).
        offsetsThatGate=zeros(3,1);
        numGatedOffsets=0;
        for offset=[0,-360,360]
            doesIntersect=rectsIntersect(rectMin1+[0;offset],rectMax1+[0;offset],boundRectMin,boundRectMax,true);
            if(doesIntersect)
                numGatedOffsets=numGatedOffsets+1;
                offsetsThatGate(numGatedOffsets)=offset;
            end
        end

        if(numGatedOffsets==0)
            %Go to the next header. This polygon is not going to be
            %taken.
            curStartByte=curStartByte+44;
            curStartByte=curStartByte+pointsInPolygon*8;
            continue;
        end

        %Get all of the other header information for this polygon.
        polygonNumber=typecast(binaryData((curStartByte+3):-1:(curStartByte)),'uint32');

        version=bitand(bitshift(flagInt,-8),255);
        greenwich=bitand(bitshift(flagInt,-16),3);
        source=bitand(bitshift(flagInt,-24),1);
        river=bitand(bitshift(flagInt,-25),1);
        p=bitshift(flagInt,-26);%Area magnitude scale

        %area=double(typecast(binaryData((curStartByte+31):-1:(curStartByte+28)),'uint32'))/double(10^p);
        area_full=double(typecast(binaryData((curStartByte+35):-1:(curStartByte+32)),'uint32'))/double(10^p);
        container=typecast(binaryData((curStartByte+39):-1:(curStartByte+36)),'int32');
        ancestor=typecast(binaryData((curStartByte+43):-1:(curStartByte+40)),'int32');

        %Save the header and shape information.
        headerInfoCur=[polygonNumber;
                      pointsInPolygon;
                      level;
                      version;
                      greenwich;
                      river;
                      source;
                      container;
                      ancestor;
                      1];%The last one is the split number.

        %Extract all of the points.
        points=zeros(2,pointsInPolygon);
        %Go to the first byte of the points.
        curStartByte=curStartByte+44;
        offset=0;
        for curPt=1:pointsInPolygon
            longitude=double(typecast(binaryData((curStartByte+3):-1:curStartByte),'int32'))/1e6+offset;
            latitude=double(typecast(binaryData((curStartByte+7):-1:(curStartByte+4)),'int32'))/1e6;
            
            %Get rid of discontinuities.
            if(curPt>1)
                lonDiff=longitude-prevLon;
                if(lonDiff>180)
                    longitude=longitude-360;
                    offset=offset-360;
                elseif(lonDiff<-180)
                    longitude=longitude+360;
                    offset=offset+360;
                end
            end
            prevLon=longitude;
            points(:,curPt)=[latitude;longitude];

            curStartByte=curStartByte+8;
        end

        %Given all those points, in the polygon, clip it to the bounding
        %box if needed If clipped to the bounding box, then we might end up
        %splitting the polygon into TWO pieces, one that ends up on either
        %side of the bounding box.

        if(abs(points(2,1)-points(2,end))==360)
            %If Antarctica (ground or ice) is taken, then we will do
            %something special to be able to handle the fact that it spans
            %all longitudes and that it isn't closed in latitude: We
            %duplicate and alias the points so that there is no issue
            %with boxes going past the end of the observable region.
            %Then, we close the polygon, so other algorithms for
            %testing whether something is in a polygon will still work.
            points=[bsxfun(@plus,points,[0;360]),points,bsxfun(@plus,points,[0;-360])];
            points=[points,[-90;points(2,end)],[-90;points(2,1)],points(:,1)];
            numGatedOffsets=1;
            offsetsThatGate(1)=0;
        end

        if(clipToView)
            if(dataType==0||dataType==1&&level==0)
                for k=1:numGatedOffsets
                    pointsCur=clipPolygonSH2D(bsxfun(@plus,points,[0;offsetsThatGate(k)]),boundRectPointsLatLon,true);
                    if(isempty(pointsCur))
                        continue;
                    end
    
                    headerInfo(:,curPolyTaken)=headerInfoCur;
                    shapeInfo(:,curPolyTaken)=[min(pointsCur(1,:));%South
                                               max(pointsCur(1,:));%North
                                               min(pointsCur(2,:));%West
                                               max(pointsCur(2,:));%East
                                               area_full];
                    vertices{curPolyTaken}=pointsCur;
    
                    headerInfoCur(10)=headerInfoCur(10)+1;
                    curPolyTaken=curPolyTaken+1;
                end
            else
                %Determine which line segments are in the box. We mark when
                %line segments start.
                pointsInLine=zeros(2,pointsInPolygon);
                for k=1:numGatedOffsets
                    pointsCur=bsxfun(@plus,points,[0;offsetsThatGate(k)]);

                    numPtsInLine=0;
                    prevPt=pointsCur(:,1);
                    pointIsLineStart=false(pointsInPolygon,1);
                    prevWasInBox=false;
                    for curPt=2:pointsInPolygon
                        [p1Clip,p2Clip]=clipLineSegment2Rect(latLonRecDeg,prevPt,pointsCur(:,curPt));
                        if(isempty(p1Clip))
                            prevWasInBox=false;
                            prevPt=pointsCur(:,curPt);
                            continue;
                        end

                        %The line segment is in the box, or is crossing
                        %the box.
                        if(prevWasInBox==false)
                            %Add the full segment.
                            pointsInLine(:,numPtsInLine+1)=p1Clip;
                            pointsInLine(:,numPtsInLine+2)=p2Clip;
                            pointIsLineStart(numPtsInLine+1)=true;
                            numPtsInLine=numPtsInLine+2;
                            prevWasInBox=true;
                        else
                            %The previous point was already in the box, so
                            %we just have to add this point.
                            numPtsInLine=numPtsInLine+1;
                            pointsInLine(:,numPtsInLine)=p2Clip;
                        end
                        prevPt=pointsCur(:,curPt);
                    end

                    %Add all the line segments.
                    if(numPtsInLine>0)
                        lineStartIdx=[find(pointIsLineStart);numPtsInLine];
                        headerInfo(:,curPolyTaken)=headerInfoCur;
                        shapeInfo(:,curPolyTaken)=[min(pointsInLine(1,1:numPtsInLine));%South
                                                   max(pointsInLine(1,1:numPtsInLine));%North
                                                   min(pointsInLine(2,1:numPtsInLine));%West
                                                   max(pointsInLine(2,1:numPtsInLine));%East
                                                   area_full];

                        %Save the points in the line as well as the indices
                        %specifying where the lines begin and end.
                        vertices{curPolyTaken}={{pointsInLine(:,1:numPtsInLine)},{lineStartIdx}};
                        headerInfoCur(10)=headerInfoCur(10)+1;
                        curPolyTaken=curPolyTaken+1;
                    end
                end
            end
        else
            headerInfo(:,curPolyTaken)=headerInfoCur;
            shapeInfo(:,curPolyTaken)=[south;
                                       north;
                                       west;
                                       east;
                                       area_full];

            vertices{curPolyTaken}=points;
            curPolyTaken=curPolyTaken+1;
        end
    end
    numPolyTaken=curPolyTaken-1;

    %Shrink to fit.
    headerInfo=headerInfo(:,1:numPolyTaken);
    shapeInfo=shapeInfo(:,1:numPolyTaken);
    vertices=vertices(1:numPolyTaken);

    mapData.latLonRecDeg=latLonRecDeg;


    %Sort the vertices and First sort by level number. Levels 1, 5,and 6
    %stay at the top, then levels 2 and 3.
    if(dataType==0)
        fiveSel=(headerInfo(3,:)==5);
        sixSel=(headerInfo(3,:)==6);
        
        %Change them to 1 for sorting purposes.
        headerInfo(3,fiveSel|sixSel)=1;
        
        %The number of things of types 1, 5 and 6.
        numRootParents=sum(headerInfo(3,:)==1);
        
        %The number of polygons that could have possible children is all
        %those that are not of whatever the highest type is, not counting 5
        %and 6. For example, type 4 if all types are represented.
        for maxChildLevel=4:-1:2
            numBase=sum(headerInfo(3,:)==maxChildLevel);
            if(numBase~=0)
                break
            end
        end
        numParents=numPolyTaken-numBase;
        if(numBase==0)
            maxChildLevel=1;%There are no children;
        end
        
        %Sort so that all level 1's come first, all two's second, etc.
        %What were 5 and 6 are placed among the 1's. The sort is
        %hierarchical so that among values the same hierarchical level, it
        %is sorted by container. Thus, all ancestors of ID 1, then all
        %ancestors of ID2, etc.
        [~,idx]=sortrows(headerInfo',[3,8]);
        headerInfo=headerInfo(:,idx);
        vertices=vertices(idx);%Sort the vertices the same way.
        shapeInfo=shapeInfo(:,idx);%Sort the shape info the same way
        
        %Replace the fives and sixes in the proper places.
        headerInfo(3,fiveSel(idx))=5;
        headerInfo(3,sixSel(idx))=6;
        
        %If there are any children, then figure out what they are.
        %otherwise, numChildren and childIdx are empty as there are no
        %children.
        if(maxChildLevel~=1)
            %Allocate space to store information on the children of each
            %node.
            numChildList=zeros(numParents,1);
            firstChildIdx=zeros(numParents,1);
            
            %Store the indices of the first child of each node, or leave it
            %zero if the node has no children. To do this, we want the
            %index of when the next level starts, which will be level 2 in
            %this case.
            curChildIdx=find(headerInfo(3,:)==2,1);
            for curParentIdx=1:numParents
                parentID=headerInfo(1,curParentIdx);
                if(headerInfo(8,curChildIdx)==parentID)
                    %The parent has children. Record the index and count
                    %them.
                    firstChildIdx(curParentIdx)=curChildIdx;
                    numChildList(curParentIdx)=0;
                    while(headerInfo(8,curChildIdx)==parentID)
                        numChildList(curParentIdx)=numChildList(curParentIdx)+1;
                        curChildIdx=curChildIdx+1;
                        if(curChildIdx>numPolyTaken)
                            break;
                        end
                    end

                    if(curChildIdx>numPolyTaken)
                        break;
                    end
                end
                %If the parent has no children, then numChildren and
                %childIdx for that entry remain zero.
            end
        else
            numChildList=[];
            firstChildIdx=[];
        end
        
        mapData.vertices=vertices;
        mapData.headerInfo=headerInfo;
        mapData.shapeInfo=shapeInfo;

        %Fill the structure to return.
        childStructureInfo.numRootParents=int32(numRootParents);
        childStructureInfo.maxChildLevel=int32(maxChildLevel);
        childStructureInfo.numChildList=int32(numChildList);
        childStructureInfo.firstChildIdx=int32(firstChildIdx);
        mapData.childStructureInfo=childStructureInfo;
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
