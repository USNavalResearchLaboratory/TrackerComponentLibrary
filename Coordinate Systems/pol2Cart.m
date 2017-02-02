function cartPoints=pol2Cart(polarPoints,systemType)
%%POL2CART  Convert points from polar coordinates to 2D Cartesian
%           coordinates. The angle can measured either ounterclockwise from
%           the x-axis, which is standard in mathematics, or clockwise from
%           the y axis, which is more common in navigation.
%
%INPUTS: polarPoints A 2XN matrix of points in polar coordinates. Each
%                    column of the matrix has the format [r;azimuth], with
%                    azimuth given in radians.
%      systemType   An optional parameter specifying the axes from which
%                   the angles are measured. Possible vaues are
%                   0 (The default if omitted) The azimuth angle is
%                     counterclockwise from the x axis.
%                   1 The azimuth angle is measured clockwise from the y
%                     axis.
%
%OUTPUTS: cartPoints A 2XN or matrix of the points transformed into
%                    Cartesian coordinates. Each columns of cartPoints is of
%                    the format [x;y].
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    systemType=0;
end

%Extract the coordinates
r=polarPoints(1,:);
azimuth=polarPoints(2,:);

switch(systemType)
    case 0
        x=r.*cos(azimuth);
        y=r.*sin(azimuth);
    case 1
        x=r.*sin(azimuth);
        y=r.*cos(azimuth);
    otherwise
        error('Invalid system type specified.')
end

cartPoints=[x;y];
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
