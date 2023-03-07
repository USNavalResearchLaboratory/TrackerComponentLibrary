function vals=linspaceNoEnd(minVal,maxVal,numPoints,bothEnds)
%LINSPACENOEND Generate numPoints values uniformly spaced starting at
%              minVal going up to but not including maxVal. If bothEnds is
%              true, then neither minVal nor maxVal will be included. The
%              distance from the last point in vals to maxVal equals the
%              spacing between all other points (plus finite precision
%              errors).
%
%INPUTS: minVal,maxVal The real scalar minimum value and maximum values
%              across which uniform points will be generated.
%    numPoints The number of uniformly spaced points to generate. if this
%              parameter is omitted or an empty matrix is passed, then a
%              default of 100 will be used.
%     bothEnds If false, minVal wil be included in the points, but maxVal
%              won't. If true, the both minVal and maxVal will be omitted
%              from the points. The default if omitted or an empty matrix
%              is passed is false.
%
%OUTPUTS: vals A 1XnumPoints vector of equispaced points.
%
%This function can be useful if one wishes to generates uniform angles
%since linspaceNoEnd(0,2*pi,numPoints) will generate numPoints uniformly
%spaced points without including 2*pi, which simply aliased back to 0.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(bothEnds))
    bothEnds=false;
end

if(nargin<3||isempty(numPoints))
    numPoints=100;
end

numPoints=floor(double(numPoints));

if(bothEnds)
    delta=(maxVal-minVal)/(numPoints+1);
    vals=minVal+delta+(0:(numPoints-1))*delta;
else
    delta=(maxVal-minVal)/numPoints;
    vals=minVal+(0:(numPoints-1))*delta;
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
