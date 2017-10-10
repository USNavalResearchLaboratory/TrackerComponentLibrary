function [zone,invalidLat]=findUTMZone(latLon,hemisphereFlip,includeExceptions)
%%FINDUTMZONE Universal transverse Mercator (UTM) coordinates are broken
%             into 60 zones that (when exceptions are not considered) are
%             all 6 degrees in width in longitude. In each zone,
%             coordinates are given according to a transverse mercator
%             projection that is centered at the center of that zone. This
%             function takes a point in latitude and longitude and returns
%             the number of the zone from 1 to 60 (also including -1 to
%             -60 if hemisphereFlip=true and the point is in the southern
%             hemisphere).
%
%INPUTS: latLon A 2X1 point in ellipsoidal [latitude;longtude] coordinates
%               given in radians. Latitude must be from -pi/2 to pi/2.
%               Longitude can be any value. The zone returned will be
%               relevant as applied to the reference ellipsoid on which the
%               latitude and longitude are based (e.g. the WGS-84
%               ellipsoid). If a single scalar value is passed, it will be
%               assumed to be just a longitude and a latitude of 0 will be
%               assigned.
% hemisphereFlip If this optional parameter is given, then for points with
%               a latitude below zero (in the southern hemisphere), the
%               zone number will be negative. Otherwise, all zone numbers
%               will be positive. Negative numbers are used in [1].
%               However, many others do not use them. The default if this
%               parameter is omitted or an empty matrix is passed is true.
% includeExceptions UTM coordinates are typically implemented such that all
%               zones are 6 degrees in width with exceptions near Norway
%               and near Svalbard, a Norwegian archipelago in the arctic.
%               For latitudes up there, the zones are modifes so as not to
%               split landmasses in inconvenient areas. The default if this
%               parameter is omitted or an empty matrix is passed is true.
%
%OUTPUTS: zone This is the zone number of the given point in UTM
%              coordinates.
%   invalidLat UTM zones are only defined for latitudes from -80 up to but
%              not including +84 degrees. If the latitude of the given
%              point falls outside this range, then this is true.
%              Otherwise, this is false.
%
%The rules governing how the zones are divided are given in Sections 7.4
%and 7.5 of [1] and are implemented here.
%
%REFERENCES:
%[1] Office of Geomatics, "National geospatial-intelligence agency
%    standardization document: Implementation practice: The universal
%    grids and the transverse mercator and polar stereographic
%    map projections," 25 Mar. 2014. [Online]. Available:
%    http://earth-info.nga.mil/GandG/publications/NGA_SIG_0012_2_0_0_UTMUPS/NGA.SIG.0012_2.0.0_UTMUPS.pdf
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(hemisphereFlip))
   hemisphereFlip=true; 
end

if(nargin<3||isempty(includeExceptions))
    includeExceptions=true; 
end

if(isscalar(latLon))%If only longitude is provided.
   latLon=[0;latLon]; 
end

phi=latLon(1);%Latitude
phiDeg=phi*(180/pi);

%If the latitude value is outside of the valid range.
if(phiDeg>=84||phiDeg<-80)
    invalidLat=true;
else
    invalidLat=false;
end

lambda=latLon(2);%Longitude
lambdaDeg=wrapRange(lambda*(180/pi),-180,180);

%This will never be 61, since the wrapRange function above does not include
%the 180 at the top of the range.
zone=floor((lambdaDeg+180)/6)+1;
if(phi<0&&hemisphereFlip)
    zone=-zone;
end

if(includeExceptions&&phi>0)
    if(zone==31&&56<=phiDeg&&phiDeg<64&&lambdaDeg>=3)
        zone=32;
    elseif(zone==32&&phiDeg>=72)
        if(lambdaDeg<9)
            zone=31; 
        else
            zone=33;
        end
    elseif(zone==34&&phiDeg>=72)
        if(lambdaDeg<21)
            zone=33;
        else
            zone=35;
        end
    elseif(zone==36&&phiDeg>=72)
        if(lambdaDeg<33)
            zone=35; 
        else
            zone=37;
        end
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
