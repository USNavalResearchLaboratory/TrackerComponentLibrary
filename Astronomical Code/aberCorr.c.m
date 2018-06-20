function vAberr=aberCorr(vOrig,obsVel,sunDist)
%%ABERCORR Rotate vectors for stellar aberration, approximately including
%          the effects of the gravitational potential of the Sun if the
%          distances between the observers and the Sun are provided.
%
%INPUTS: vOrig A 3XN set of N vectors pointing from the observers to the
%              objects being observed without corruption due to
%              aberration. The units of the vector do not matter.
%              Normally, one would expect vOrig to be composed of unit
%              vectors as aberration affects the pointing direction of the
%              vectors.
%       obsVel The 3XN matrix of the velocity of each observer with
%              respect to the stellar coordinate system origin in meters
%              per second. For example, if using the barycentric celestial
%              reference system (BCRS), they are the set of velocity
%              vectors of the observers with respect to the barycenter.
%      sunDist An optional 3XN vector of the scalar distances in meters
%              from the Sun to the observer. If this parameter is
%              included, an approximate correction due to the
%              gravitational potential of the Sun is included.
%
%OUTPUTS: vAberr The 3XN set of N vectors rotated to deal with aberration
%                effects.
%
%If all parameters are provided, the function is essentially a wrapper for
%the function iauAb in the International Astronomical Union's Standards of
%Fundamental Astronomy library with an adjustment so that the vector in
%question need not have unit magnitude. If the distance to the Sun is
%omitted, then the standard special relativistic aberration correction
%without any gravitational effects is applied. The standard special
%relativistic aberration correction is described in Chapter 7 of [1]
%
%Note that if the vectors provided are meant to be apparent distances,
%rather than just unit vectors representing directions, the transformation
%does not adjust the magnitudes of the vectors to account for special
%relativistic contraction of space due to the motion of the observer with
%respect to the origin. The magnitudes of the output vectors equal the
%magnitudes of the input vectors.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%vAberr=aberCorr(vOrig,obsVel,sunDist);
%or
%vAberr=aberCorr(vOrig,obsVel);
%
%REFERENCES:
%[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
%    Books, 2013.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
