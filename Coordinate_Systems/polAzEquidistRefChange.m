function pAzPts=polAzEquidistRefChange(pAzPts,latLonRefOld,latLonRefNew,a,f,algorithm,desAngSpan)
%%POLAZEQUIDISTREFCHANGE Points given in polar azimuthal equidistant
%   coordinates are defined with respect to a particular reference location
%   on the surface of the reference ellipsoid or sphere. Here, given one
%   set of polar azimuthal equidistant coordinates, obtain the coordinates
%   when the reference points has been changed. This can be useful for
%   expressing rays traced from one sensor in the coordinate system of
%   another sensor.
%
%INPUTS: pAzPts A 2XN set of the [ground distance; heading] points, with the
%              heading given in radians East of North, to convert.
%              Alternatively, if heights are also given, this can be a 3XN
%              set of points with the height being the third dimension.
% latLonRefOld A 2X1 [latitude;longitude] reference point in radians about
%              which pAzPts was computed.
% latLonRefNew A 2X1 [latitude;longitude] reference point in radians with
%              respect to which pAzPts should be converted.
%            a The semi-major axis of the reference ellipsoid (in meters).
%              If this argument is omitted or an empty matrix is passed,
%              the value in Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%    algorithm If f is not zero, then this parameter is not used. If f=0,
%              then this parameter selects the algorithm to use. Possible
%              values are:
%              0 (The default if omitted or an empty matrix is used) Use a
%                formula specific for the transformation on a sphere.
%              1 Call the polarAzEquidistProj2Ellipse function and then the
%                ellips2PolarAzEquidistProj function.
%   desAngSpan A 2X1 range of angular values [min Val;max Val] in which it
%              is desired that the converted azimuthal values should
%              attempt to be kept (by adding or subtracing either 0 or 2*pi
%              [but not a multiple thereof] to the converted value to get
%              it as close as possible to the span). If omitted or an empty
%              matrix is passed, no attempt is made to keep the azimuth
%              within any particular bounds.
%
%OUTPUTS: pAzPts A 2XN (or 3XN) set of [ground distance;heading;(height)]
%                points in polar azimuthal equidistant coordinates taken
%                with respect to latLonRefNew.
%
%If f is 0 and algorithm=0, then the coordinate system change is occurring
%on a sphere and a direct conversion is used. Otherwise, this function just
%calls polarAzEquidistProj2Ellipse to convert the points into latitude and
%longitude coordinates and then uses ellips2PolarAzEquidistProj to convert
%them into polar azimuthal equidistant coordinates with a different
%reference point.
%
%EXAMPLE 1:
%This example shows that the relative error of this conversion is on the
%order of finite precision limitatios when compared to simply using the
%desired alternate reference point from the start.
% latLonPt=[deg2rad([20.631756;-155.551818;10]),deg2rad([21.295516;-158.002386;20])];
% latLonRef1=deg2rad([20.906029;-157.078918]);
% latLonRef2=deg2rad([35;-125]);
% pAzPts1=ellips2PolarAzEquidistProj(latLonPt,latLonRef1);
% pAzPts2=ellips2PolarAzEquidistProj(latLonPt,latLonRef2);
% pAzPtsConv=polAzEquidistRefChange(pAzPts1,latLonRef1,latLonRef2);
% RelErr=max(max(abs((pAzPtsConv-pAzPts2)./pAzPts2)))
%
%EXAMPLE 2:
%This example generates a number of random values and then computes the
%error in the conversion compared to directly converting it. The maximum
%relative error magnitude is displayed and is reasonable given the
%scale of the distances across the surface of the Earth and finite
%precision limits.
% rE=Constants.WGS84MeanRadius;
% numMCRuns=10e3;
% maxDiff=0;
% for curRun=1:numMCRuns
%     latLonRef1=[UniformD.rand(1,[-pi/2;pi/2]);
%                 UniformD.rand(1,[-pi;pi])];
%     latLonRef2=[UniformD.rand(1,[-pi/2;pi/2]);
%                 UniformD.rand(1,[-pi;pi])];
%     latLonPt=[UniformD.rand(1,[-pi/2;pi/2]);
%               UniformD.rand(1,[-pi;pi])];    
% 
%     pAzPts1=ellips2PolarAzEquidistProj(latLonPt,latLonRef1,rE,0);
%     pAzPts2=ellips2PolarAzEquidistProj(latLonPt,latLonRef2,rE,0);
%     pAzConv=polAzEquidistRefChange(pAzPts1,latLonRef1,latLonRef2,rE,0);
% 
%     diffVal=(pAzPts2-pAzConv);
%     diffVal(2)=wrapRange(diffVal(2),-pi,pi);
%     diffVal=diffVal./pAzPts2;%Make a relative error.
%     maxDiff=max(maxDiff,norm(diffVal));
% end
% maxDiff
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7)
    desAngSpan=[];
end

if(nargin<6||isempty(algorithm))
    algorithm=0;
end

if(nargin<5||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<4||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(algorithm~=0&&algorithm~=1)
    error('Unknown algorithm specified.')
end

if(f==0&&algorithm==0)
    [azStart,azEnd]=greatCircleAzimuth(latLonRefOld,latLonRefNew);
    azOffset1=pi/2-azStart;
    azOffset2=pi/2-azEnd;
    lambdaDelta=greatCircleDistance(latLonRefOld,latLonRefNew,1);
    
    %Transform into the system with both sensors on the equator.
    pAzPts(2,:)=pAzPts(2,:)+azOffset1;
    
    p=pAzPts(1,:);
    Az=pAzPts(2,:);
    rE=a;
    pAzPts(1:2,:)=[rE.*acos(cos(p./rE).*cos(lambdaDelta)+sin(p./rE).*sin(Az).*sin(lambdaDelta));
                   atan2(sin(p./rE).*cos(lambdaDelta).*sin(Az)-cos(p./rE).*sin(lambdaDelta),sin(p./rE).*cos(Az))];
    
    %Convert back to the global coordinate system.
    pAzPts(2,:)=pAzPts(2,:)-azOffset2;    
else
    convLatLon=polarAzEquidistProj2Ellipse(pAzPts,latLonRefOld,a,f);
    pAzPts=ellips2PolarAzEquidistProj(convLatLon,latLonRefNew,a,f);
end

if(~isempty(desAngSpan))
    N=size(pAzPts,2);
    
    minVal=desAngSpan(1);
    maxVal=desAngSpan(2);

    for k=1:N
        curAz=pAzPts(2,k);

        if(minVal>curAz)
            %See if adding 2*pi places it closer to the span.
            azTest=curAz+2*pi;
            offsetOld=minVal-curAz;
        elseif(maxVal<curAz)
            %See if subtracting 2*pi places it closer to the span.
            azTest=curAz-2*pi;
            offsetOld=curAz-maxVal;
        else
            %It is in the span.
            continue;
        end

        %If it is closer to the span than before
        offsetNew=min(abs(minVal-azTest),abs(maxVal-azTest));

        if(offsetNew<offsetOld)
            pAzPts(2,k)=azTest;
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
