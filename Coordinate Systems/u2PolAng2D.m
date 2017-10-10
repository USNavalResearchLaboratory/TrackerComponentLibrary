function azimuth=u2PolAng2D(uVals,systemType)
%%U2POLANG2D Convert the direction cosine value u in 2D into a polar
%            angle. The "direction cosine" u is just the x coordinate of a
%            unit vector from the receiver to the target in the coordinate
%            system at the receiver. Optionally, a full [u;v] unit vector
%            in 2D can be provided. The sign of the other component
%            eliminates ambiguity about whether the target is in front of
%            the sensor. If omitted, the missing v component is assumed to
%            be positive.
%
%INPUTS: uVals A 1XnumPoints (for only u) or a 2XnumPoints (if full unit
%              vectors are given) set of direction cosines in 2D to turn
%              into polar angles.
%   systemType An optional parameter specifying the axis from which the
%              angles are measured. Possible values are
%              0 (The default if omitted or an empty matrix is passed) The
%                azimuth angle is counterclockwise from the x axis.
%              1 The azimuth angle is measured clockwise from the y axis.
%
%OUTPUTS: azimuth The 1XN set of polar angles corresponding to u.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
    systemType=0;
end

u=uVals(1,:);
if(size(uVals,1)==1)%If it does not come with v
    switch(systemType)
        case 0
            azimuth=acos(u);
        case 1
            azimuth=asin(u);
        otherwise
            error('Invalid system type specified.')
    end
else
    v=uVals(2,:);
    switch(systemType)
        case 0
            azimuth=atan2(v,u);
        case 1
            azimuth=atan2(u,v);
        otherwise
            error('Invalid system type specified.')
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
