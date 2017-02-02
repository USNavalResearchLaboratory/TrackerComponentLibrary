function points=Cart2Pol(cartPoints,systemType)
%%CART2POL Convert a set of 2D Cartesian coordinates to 2D polar
%          coordinates. The angle can be measured either counterclockwise
%          from the x-axis, which is standard in mathematics, or clockwise
%          from the y axis, which is more common in navigation.
%
%INPUTS: cartPoints A 2XN matrix of Cartesian points that are to be
%                   transformed into polar coordinates. Each columns of
%                   cartPoints is of the format [x;y] or [x;y;z].
%      systemType   An optional parameter specifying the axes from which
%                   the angles are measured. Possible vaues are
%                   0 (The default if omitted) The azimuth angle is
%                     counterclockwise from the x axis.
%                   1 The azimuth angle is measured clockwise from the y
%                     axis.
%
%OUTPUTS: points    A 2XN matrix of the converted points. Each column of
%                   the matrix has the format [r;azimuth], with azimuth
%                   given in radians.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2)
        systemType=0;
    end

    %Extract the coordinates and find r.
    x=cartPoints(1,:);
    y=cartPoints(2,:);
    r=sqrt(x.^2+y.^2);
    
    switch(systemType)
        case 0
            azimuth=atan2(y,x);
        case 1
            azimuth=atan2(x,y);
        otherwise
            error('Invalid system type specified.')
    end

    points=[r;azimuth];
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
