function latLonPts=azEquidistantProj2Ellipse(xyPts,latLonRef,a,f,algorithm)
%%AZEQUIDISTANTPROJ2ELLIPSE Given points as [x;y] values in a azimuthal
%       equidistant projection, convert the points to latitudes and
%       longitudes on a reference ellipsoid (or sphere). A third height
%       coordinate can also be provided, which just gets appended to the
%       output.
%
%INPUTS: xyPts A 2XN set of the azimuthal equidistant projection points
%              about latLonRef to convert. Alternatively, if heights are
%              given, this can be a 3XN set of points.
%    latLonRef A 2X1 [latitude;longitude] reference point in radians about
%              which the projection is taken.
%            a The semi-major axis of the reference ellipsoid (in meters).
%              If this argument is omitted or an empty matrix is passed,
%              the value in Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%    algorithm If f=0, then this optinal input selects the algorithm to be
%              used. Possibel values are:
%              0 Use the function directGreatCircleProb.
%              1 (The default if omitted or an empty matrix is passed) Use
%                the formulae in Chapter 25 of [1]. This might have numeric
%                issues if the reference point is very close to but not
%                exactly at a pole.
%
%OUTPUTS: latLonPts The 2XN set of converted [latitude;longitude] points in
%                   radians. If heights were given in xyPts, then this is a
%                   3XN set of converted [latitude;longitude;height] points
%                   with the third row the same as in xyPts.
%
%The conversion is described in Chapter 25 of [1], where expressions for
%a spherical Earth are given. However, we only use those formulae for
%algorithm 1 with f=0. The norm of each xy point corresponds to a distance
%across the surface of the curved Earth and using the inverse tangent
%function, one can obtain a launch azimuth in radians East of North. From
%there, the latitude and longitude of the point can be obtained using the
%directGreatCircleProb function for a spherical Earth or the
%directGeodeticProb for an ellipsoidal Earth.
%
%EXAMPLE:
%Here, a value is converted into azimuthal equidistant coordinates on a
%spherical Earth and is then converted back. The relative error of the
%reverse conversion under algorithms 0 and 1 are shown and one can see the
%results are within finite precision limitiations.
% latLonPt=deg2rad([24.266997;-152.244802]);
% latLonRef=deg2rad([21.295516;-158.002386]);
% rE=osculatingSpher4LatLon(latLonRef);
% xyPt=ellips2AzEquidistantProj(latLonPt,latLonRef,rE,0);
% latLonPt0=azEquidistantProj2Ellipse(xyPt,latLonRef,rE,0,0);
% latLonPt1=azEquidistantProj2Ellipse(xyPt,latLonRef,rE,0,1);
% RelErr0=max(abs((latLonPt0-latLonPt)./latLonPt))
% RelErr1=max(abs((latLonPt1-latLonPt)./latLonPt))
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=1;
end

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

hasHeight=(size(xyPts,1)==3);
numPts=size(xyPts,2);
if(hasHeight)
    latLonPts=zeros(3,numPts);
else
    latLonPts=zeros(2,numPts);
end
if(f==0)
    %Under a spherical Earth approximation.
    if(algorithm==0)
        %Just call directGreatCircleProb.
        for k=1:numPts
            distVal=norm(xyPts(1:2,k));
            az=atan2(xyPts(1,k),xyPts(2,k));
            
            latLonPts(1:2,k)=directGreatCircleProb(latLonRef,az,distVal,a);
        end
    else
        %Use the formulae in Chapter 25 of [1].
        phi1=latLonRef(1);
        lambda0=latLonRef(2);
        cosPhi1=cos(phi1);
        sinPhi1=sin(phi1);

        for k=1:numPts
            x=xyPts(1,k);
            y=xyPts(2,k);
            rho=sqrt(x^2+y^2);%Equation 2-18.
            c=rho/a;%Equation 25-15.
            sinC=sin(c);
            cosC=cos(c);
            %Equation 20-14.
            phi=asin(cosC*sinPhi1+(y*sinC*cosPhi1/rho));
            if(phi1==pi/2)
                %Equation 20-16.
                lambda=lambda0+atan2(x,-y);
            elseif(phi1==-pi/2)
                %Equation 20-17.
                lambda=lambda0+atan2(x,y);
            else
                %Equation 20-15.
                lambda=lambda0+atan2(x*sinC,(rho*cosPhi1*cosC-y*sinPhi1*sinC));
            end
            latLonPts(1:2,k)=[phi;lambda];
        end

    end
else
    %Under an ellipsoidal Earth approximation.
    for k=1:numPts
        distVal=norm(xyPts(1:2,k));
        az=atan2(xyPts(1,k),xyPts(2,k));
        
        latLonPts(1:2,k)=directGeodeticProb(latLonRef,az,distVal,a,f);
    end
end

if(hasHeight)
    latLonPts(3,:)=xyPts(3,:);
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
