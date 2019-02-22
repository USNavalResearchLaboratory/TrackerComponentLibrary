function vDefl=lightDeflectCorr(vObsSource,posObs,MSolar,xBody,deflecLimit)
%%LIGHTDEFLECTCORR Rotate vectors to account for light-deflection due to
%             one or more solar system bodies using an approximate formula
%             to deal with the light-time problem without resorting to ray
%             tracing.
%
%INPUTS: vObsSource A 3XN set of vectors pointing from the observer to the 
%              objects being observed without corruption due to light 
%              deflection by the Sun or a planet These are direction 
%              vectors in the International Celestial Reference System 
%              (ICRS). Normally, one would expect vObsSource to be unit
%              vectors.
%      posObs The 3X1 Cartesian position vector of the observer in the
%             Barycentric Celestial Reference System (BCRS) at the time of
%             observation. The units are meters.
%      MSolar The numBodX1 vector of the masses of the gravitating bodies
%             that will deflect the light. The masses are solar masses.
%             That is, the ratio of the mass of the body to the mass of
%             the Sun. Thus, the input is unitless.
%       xBody A 6XnumBodies set of position and velocity state vectors in
%             the BCRS for the defecting bodies at the time of
%             observation. The units are meters and meters per second.
% deflecLimit An optional numBodX1 set of deflection limiter parameters in
%             radians that are phi^2/2 for each body, where phi is the
%             angular distance between the source and the thing that is to
%             deflect the light of the source. As phi goes below the
%             threshold, the deflection applied goes to zero. If this is
%             omitted, the default value of 1e-10 for all of the inputs
%             is used.
%
%OUTPUTS: uAberr The 3XN set of N vectors rotated to deal with the
%            deflection due to the specified astronomical bodies.
%             
%This function is primarily a wrapper for the approximate light deflection
%formula used in the iauLdn in the International Astronomical Union's
%(IAU's) Standards of Fundamental Astronomy (SOFA) library.
%
%Example solar mass values and deflection limit values taken from the
%IAU's documentation for the iauLdn function are
%Body       Solar Mass     Deflection Limit
%Sun        1.0            6e-6
%Jupiter    0.00095435     3e-9
%Saturn     0.00028574     3e-10
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%vDefl=lightDeflectCorr(vObsSource,posObs,MSolar,xBody);
%or
%vDefl=lightDeflectCorr(vObsSource,posObs,MSolar,xBody,deflecLimit);
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
