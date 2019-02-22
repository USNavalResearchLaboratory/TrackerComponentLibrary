function u=polAng2U2D(azimuth,systemType,includeV)
%%POLANG2U2D Convert 2D polar angles into equivalent direction cosine
%            values. The "direction cosine" u is just the x coordinate of a
%            unit vector from the receiver to the target in the coordinate
%            system at the receiver. Optionally, a full [u;v] unit vector
%            in 2D can be returned if includeV is true.
%
%INPUTS: azimuth A 1XN or NX1 set of azimuthal values.
%     systemType An optional parameter specifying the axis from which the
%                angles are measured. Possible values are
%                0 (The default if omitted or an empty matrix is passed)
%                  The azimuth angle is counterclockwise from the x axis.
%                1 The azimuth angle is measured clockwise from the y axis.
%       includeV An optional boolean value indicating whether a second
%                direction cosine component should be included. The u
%                direction cosine is one parts of a 2D unit vector.
%                Generally, one might assume that the target is in front of
%                the sensor, so the second component would be positive and
%                is not needed. However, the second component can be
%                included if ambiguity exists. The default if this
%                parameter is omitted or an empty matrix is passed is
%                false.
%
%OUTPUTS: u A 1XN (without v) or 2XN( with v) set of direction cosines in
%           2D corresponding to the specified azimuthal values.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(includeV))
   includeV=false; 
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

%Make a row vector.
azimuth=azimuth(:).';

switch(systemType)
    case 0
        if(includeV)
            u=[cos(azimuth);sin(azimuth)];
        else
            u=cos(azimuth);
        end
    case 1
        if(includeV)
            u=[sin(azimuth);cos(azimuth)];
        else
            u=sin(azimuth);
        end
    otherwise
        error('Invalid system type specified.')
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
