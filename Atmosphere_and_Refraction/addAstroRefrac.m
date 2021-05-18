function [z0,deltaZ]=addAstroRefrac(algorithm,plhObs,zTrue,Rh,P,T,wl)
%%ADDASTROREFRAC Add the effects of refraction to the true zenith angle of
%                an object outside the atmosphere using low-precision
%                atmospheric models for an observer near the surface of
%                the Earth viewing an object using a narrowband sensor.
%
%INPUTS: algorithm This specified the algorithm used. The possible values
%                  are:
%                  0: Use a numerical integration method from Chapter 7.2
%                     of [1] and from [2]. This algorithm uses the inputs
%                     plhObs,z0, Rh, P, T, and wl and makes use of the
%                     Sinclair atmospheric model. This algorithm is the
%                     most precise of all of the methods. If a point is too
%                     far below the horizon, z0 and deltaZ will be empty
%                     matrices.
%                  1: Use the simple formula of [3] for the refraction
%                     experienced by an observer at sea level viewing light
%                     with a wavelength of 0.574 micrometers (yellow). This
%                     algorithm uses the inputs z0, Rh, P, T. This
%                     algorithm is only valid for refraction-corrupted
%                     zenith distances below 70 degrees. Values that are
%                     too high will lead to z0 and deltaZ being empty
%                     matrices.
%                  2: Use the algorithm from the International Astronomical
%                     Union's (IAU) standards of fundamental astronomy
%                     library. This algorithm uses the inputs z0,Rh,P,T,wl
%                     for an observer at sea level. This algorithm will
%                     produce results for all positive zenith distances,
%                     though the results might not be very good for large
%                     values.
%           plhObs The WGS-84 ellipsoidal latitude, longitude and height
%                  of the observer in meters. Only algorithm 0 uses this
%                  parameter, and only the latitude and height above the
%                  reference ellipsoid are used (the latitude does not
%                  matter). The height is treated as an approximate height
%                  above mean sea level. Given the precision of the
%                  algorithm, the fact that geoid undulations are ignored
%                  probably does not matter.For algorithms 1 and 2, an
%                  empty matrix can be passed.
%            zTrue A vector or matrix of true (refraction-free) positive
%                  zenith "distances" (in radians) of the celestial object
%                  being observed. A zenith distance is the observed angle
%                  of an object down from the gravitational vertical (the
%                  zenith) at the observer's location. This value is
%                  inaccurate as one approaches pi/2 radians (the horizon).
%                  zTrue must be greater than zero. Given the precision of
%                  the algorithms that are available, an angle with respect
%                  to the vertical defined by the WGS-84 reference
%                  ellipsoid could probably be substituted for an angle
%                  with respect to the true gravitational vertical.
%               Rh The relative humidity at the observer (between 0 and 1).
%                  If this parameter is omitted or an empty matrix is
%                  passed, then Constants.standardRelHumid is used.
%                P The atmospheric pressure at the observer in Pascals
%                  (N/m^2). If this parameter is omitted or an empty matrix
%                  is passed, then Constants.standardAtmosphericPressure is
%                  used.
%                T The air temperature at the observer in degrees Kelvin.
%                  If this parameter is omitted or an empty matrix is
%                  passed, then Constants.standardTemp is used.
%               wl The wavelength at which the observation is made in units
%                  of meters. If this parameter is omitted or an empty
%                  matrix is passed, then a wavelength of 0.574 micrometers
%                  is used, which is in the visible spectrum (a rather
%                  yellow color).
%
%OUTPUTS: z0  An NX1 vector of the zTrue values in radians with atmospheric
%             refraction added.
%      deltaZ The refraction value that was applied. z0=zTrue-deltaZ.
%
%The function removeAstroRefrac computes deltaZ for the problem where z0 is
%known and zTrue is unknown. This function iterates the solution a few
%times to try to solve the inverse problem. A fixed 20 iterations are used,
%which is generally sufficient to ensure convergence to within working
%precision limits for all of the algorithms.
%
%REFERENCES:
%[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
%    Books, 2013.
%[2] C. Y. Hohenkerk and A. T. Sinclair, "The computation of angular
%    atmospheric refraction at large zenith angles," United Kingdom
%    Hydrographic Office, HM Nautical Almanac, Tech. Rep. 63, Apr. 1985.
%    http://astro.ukho.gov.uk/data/tn/naotn63.pdf
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<7||isempty(wl))
       wl=0.574e-6; 
    end

    if(nargin<6||isempty(T))
       T=Constants.standardTemp; 
    end
    
    if(nargin<5||isempty(P))
        P=Constants.standardAtmosphericPressure;
    end
    
    if(nargin<4||isempty(Rh))
        Rh=Constants.standardRelHumid;
    end
    
    %The initial estimate of deltaZ is given by solving the inverse problem
    %at zTrue
    [~,deltaZ]=removeAstroRefrac(algorithm,plhObs,zTrue,Rh,P,T,wl);
    
    if(~isempty(deltaZ))
        numIter=20;
        for curIter=1:numIter
            [~,deltaZ]=removeAstroRefrac(algorithm,plhObs,zTrue-deltaZ,Rh,P,T,wl);
            if(isempty(deltaZ))
                %This can arise if the observation ends up too far
                %underground.
                break
            end
        end
    end
    if(isempty(deltaZ))
        z0=[];
    else
        z0=zTrue-deltaZ;
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
