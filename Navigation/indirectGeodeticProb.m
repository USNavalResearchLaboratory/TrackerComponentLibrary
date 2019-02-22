function [azStart,dist,azEnd]=indirectGeodeticProb(latLonStart,latLonEnd,a,f)
%%INDIRECTGEODETICPROB Solve the indirect geodetic problem. That is, given
%                      two points on an ellipsoidal Earth, find the
%                      initial bearing and distance one must travel to
%                      take the shortest (geodesic) path between the
%                      points.
%
%INPUTS:latLonStart A 2XN matrix of the N initial points given in geodetic
%                  latitude and longitude in radians of the format
%                  [latitude;longitude] for each column (point). The
%                  latitude must be between -pi/2 and pi/2 radians and the
%                  longitude between -pi and pi radians.
%        latLonEnd The 2XN matrix of final points for the geodesic path
%                  given in geodetic latitude  and longitude in radians.
%                  latLonEnd has the same format at latLonStart.
%                a The semi-major axis of the reference ellipsoid (in
%                  meters). If this argument is omitted, the value in
%                  Constants.WGS84SemiMajorAxis is used.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted, the value in
%                  Constants.WGS84Flattening is used.
%
%OUTPUTS: azStart The NX1 forward azimuth at the starting points in
%                 radians East of true North on the reference ellipsoid.
%                 This is the initial heading one would travel to go
%                 between latLonStart and latLonEnd for each of the point
%                 pairs
%            dist The NX1 vector of geodetic distances between the
%                 starting and stopping points in meters.
%           azEnd The forward azimuth at the ending point in radians East
%                 of true North on the reference ellipsoid.
%
%The function is essentially a Matlab interface for the implementation in
%GeographicLib, which is documented in [1], [2], and [3]. GeographicLib
%can be downloaded from
%http://geographiclib.sourceforge.net
%Though a native Matlab version of the relevant function in GepgraphicLib
%exists, it is rather slow. hence the need for this interface to the
%compiled version.
%
%The algorithm can be compiled for use in Matlab using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[azStart,dist,azEnd]=indirectGeodeticProb(latLonStart,latLonEnd);
%or if something other than the WGS84 reference ellipsoid is used
%[azStart,dist,azEnd]=indirectGeodeticProb(latLonStart,latLonEnd,a,f);
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-45, Jan. 2013. [Online].
%    Available: http://arxiv.org/pdf/1109.4448.pdf
%[2] C. F. F. Karney. (2013, 31 Aug.) Addenda and errata for papers on
%    geodesics. [Online].
%    Available: http://geographiclib.sourceforge.net/geod-addenda.html
%[3] C. F. F. Karney. (2011, 7 Feb.) Geodesics on an ellipsoid of
%    revolution.
%    [Online]. Available: http://arxiv.org/pdf/1102.1215.pdf
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

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
