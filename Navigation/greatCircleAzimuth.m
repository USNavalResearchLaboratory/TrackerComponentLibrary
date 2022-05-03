function [azStart,azEnd]=greatCircleAzimuth(latLonStart,latLonEnd,algorithm)
%%GREATCIRCLEAZIMUTH On a spherical Earth, find the initial bearing to
%       travel from latLonStart to latLonEnd. Also find the bearing with
%       which one will arrive at the endpoint.
%
%INPUTS: latLonStart The 2X1 initial point given in latitude and longitude
%                    in radians in the format [latitude;longitude]
%                    (on a reference sphere, latitude is spherical
%                    elevation; longitude is spherical azimuth).
%          latLonEnd The 2X1 final point given in latitude and longitude in
%                    radians in the format [latitude;longitude].
%          algorithm Optionally, select the formula used to solve the
%                    problem. This can affect the numeric stability.
%                    Possible values are:
%                    0, 1 Simply call indirectGreatCircleProb with the
%                       algorithm of the same number.
%                    2 (The default if omitted or an empty matrix is
%                       passed). Use formula 5-4b in [1]. This tends to be
%                       simpler and numerically stabler than the others.
%
%OUTPUTS: azStart The scalar forward azimuth at the starting point in
%                 radians East of true North on the reference sphere.
%           azEnd The forward azimuth at the ending point in radians East
%                 of true North on the reference sphere.
%
%EXAMPLE:
%This example just shows how to call the function between two
%widely-separated points.
% latLon1=deg2rad([19.7241;-155.0868]);
% latLon2=deg2rad([48.8566;2.3522]);
% [azStart,azEnd]=greatCircleAzimuth(latLon1,latLon2)
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=2;
end

if(algorithm==0||algorithm==1)
    [azStart,~,azEnd]=indirectGreatCircleProb(latLonStart,latLonEnd,[],[],algorithm);
elseif(algorithm==2)%Use a simple formula.
    phi0=latLonStart(1,:);
    lambda0=latLonStart(2,:);
    phi=latLonEnd(1,:);
    lambda=latLonEnd(2,:);

    deltaLambda=lambda-lambda0;
    sinDeltaLambda=sin(deltaLambda);
    cosDeltaLambda=cos(deltaLambda);
    sinPhi=sin(phi);
    cosPhi=cos(phi);
    sinPhi0=sin(phi0);
    cosPhi0=cos(phi0);
    
    azStart=atan2(cosPhi.*sinDeltaLambda,cosPhi0.*sinPhi-sinPhi0.*cosPhi.*cosDeltaLambda);
    
    if(nargout>1)
        azEnd=atan2(-cosPhi0.*sinDeltaLambda,cosPhi.*sinPhi0-sinPhi.*cosPhi0.*cosDeltaLambda);
        azEnd=wrapRange(azEnd+pi,-pi,pi);
    end
else
    error('Unknown algorithm specified.')
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
