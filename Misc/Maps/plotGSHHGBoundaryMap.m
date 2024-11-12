function plotGSHHGBoundaryMap(mapData,varargin)
%%PLOTGSHHGBOUNDARYMAP Plot the map information data returned by the
%    getGSHHGBoundaryData function on an equirectangular projection world
%    map using the current figure (without clearing it). The x-axis is
%    longitude in degrees. A call to this function without parameters will
%    plot a low resolution outline map of the continents and other major
%    landmasses. Multiple calls using the different data types from the
%    getGSHHGBoundaryData function as inputs can build up a detailed map
%    including rivers and poltical boundaries. The plot is drawn in the
%    active figure. Large highly-detailed plotting can be slow. The
%    outer boundary of the plot can have random overlapping clipped lines
%
%INPUTS: mapData This is the mapData structure returned by the
%                getGSHHGBoundaryData function; see the comments to that
%                function for more details ont he format. If omitted or an
%                empty matrix is passed, then the function is used to get
%                mapData using latLonRecDeg=[-90;90;-180;180], dataType=0,
%                resolutionLevel=1 and excludeTypes=[2;3;4;6];
%       varargin The same variable-length set of arguments that one would
%                pass to format the plot function. This is just passed on
%                to the plot function to format the lines. The default if
%                omitted is '-b'.
%
%OUTPUTS: None. The function just draws in the current plot.
%
%EXAMPLE:
%This example plots the land with country boundaries and the borders of
%major lakes (both at a low resolution) with the center of the map not
%being at 0 degrees longitude, but shifted.
% latLonRecDeg=[-90;90;-180-90;180-90];
% theLand=getGSHHGBoundaryData(latLonRecDeg,0,1,[3;4;6]);
% theCountries=getGSHHGBoundaryData(latLonRecDeg,2,1,[2;3]);
% figure(1)
% clf
% hold on
% plotGSHHGBoundaryMap(theLand)
% plotGSHHGBoundaryMap(theCountries)
%
%April 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<1||isempty(mapData))
        latLonRecDeg=[-90;90;-180;180];
        dataType=0;
        resolutionLevel=1;
        excludeTypes=[2;3;4;6];
        mapData=getGSHHGBoundaryData(latLonRecDeg,dataType,resolutionLevel,excludeTypes);
    end

    latLonRecDeg=mapData.latLonRecDeg;
    vertices=mapData.vertices;
    
    if(nargin<2)
        varargin{1}='-b';
    end

    holdVal=ishold();
    
    %Clear the axes if hold was not already on for the active figure.
    if(holdVal==false)
        clf
    end
    hold on

    numItems=length(vertices);

    for curItem=1:numItems
        if(iscell(vertices{curItem}))
            %If we are plotting a river or a state boundary as a line.
            pointsInLine=vertices{curItem}{1}{:};
            lineStartIdx=vertices{curItem}{2}{:};
            numSegments=length(lineStartIdx)-1;
            for k1=1:numSegments
                idx=lineStartIdx(k1):lineStartIdx(k1+1);
                plot(pointsInLine(2,idx),pointsInLine(1,idx),varargin{:})
            end
        else
            pointsCur=vertices{curItem};
            plot(pointsCur(2,:),pointsCur(1,:),varargin{:})
        end
    end

    axis(latLonRecDeg([3;4;1;2]));
    
    if(holdVal==false)
        hold off
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
