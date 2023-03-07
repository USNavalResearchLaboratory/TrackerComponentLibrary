function J=polAzEquiDistRefChangeCrossGrad(pAzPt,latLonRefOld,latLonRefNew,rE)
%%POLAZEQUIDISTREFCHANGECROSSGRAD Given a point in polar azimuthal
%      coordinates with respect to a reference position on a spherical
%      Earth, find the gradient of the point in a polar azimuthal
%      coordinate system with a difference reference point taken with
%      respect to the point in the original system. This is the gradient of
%      the polAzEquidistRefChange function.
%
%INPUTS: pAzPts A 2XN set of the [ground distance; heading] points, with
%              the heading given in radians East of North, to convert.
%              Alternatively, if heights are also given, this can be a 3XN
%              set of points with the height being the third dimension.
% latLonRefOld A 2X1 [latitude;longitude] reference point in radians about
%              which pAzPts was computed.
% latLonRefNew A 2X1 [latitude;longitude] reference point in radians with
%              respect to which the function to be differentied is defined.
%           rE The radius of the reference sphere. If this argument is
%              omitted or an empty matrix is passed, the value in
%              Constants.WGS84MeanRadius is used.
%
%OUTPUTS: J A 2X2XN matrix holding the gradient of the transformation
%           evaluated at each point in pAzPts (or a 3X3XN matrix if height
%           is also provided). If one says that a point in pAzPts=[p;az;h]
%           and a converted point is [p1;az1;h1], the gradients in 3D are
%           ordered:
%           [dp1dp,   dp1dAz,  dp1dh;
%            dAz1dp, dAz1dAz, dAz1dh;
%            dh1dp,   dh1dAz,  dh1dh];          
% 
%The gradient is the analytic gradient of the spherical Earth conversion
%given in the function polAzEquidistRefChange. Due to a singularity, the
%gradient is not valid if the target is exactly on the opposite side of the
%Earth from the first sensor.
%
%EXAMPLE:
%Given three random points, the accuracy of this function is compared to
%numeric differentiation. The relative error is on the order of what one
%might expect due to finite precision limitations.
% rE=Constants.WGS84MeanRadius;
% latLonRef1Old=[UniformD.rand(1,[-pi/2;pi/2]);
%                UniformD.rand(1,[-pi;pi])];
% latLonRefNew=[UniformD.rand(1,[-pi/2;pi/2]);
%               UniformD.rand(1,[-pi;pi])];
% latLonPt=[UniformD.rand(1,[-pi/2;pi/2]);
%           UniformD.rand(1,[-pi;pi])];   
% pAzPt=ellips2PolarAzEquidistProj(latLonPt,latLonRef1Old,rE,0);
% J=polAzEquiDistRefChangeCrossGrad(pAzPt,latLonRef1Old,latLonRefNew,rE);
% f=@(x)polAzEquidistRefChange(x,latLonRef1Old,latLonRefNew,rE,0);
% JNumDiff=numDiff(pAzPt,f,2);
% RelErr=max(max(abs((J-JNumDiff)./JNumDiff)))
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(rE))
    rE=Constants.WGS84MeanRadius; 
end

numPoints=size(pAzPt,2);
if(all(latLonRefOld==latLonRefNew))
    %For the special case of no change, force the result to be exact (avoid
    %finite precision issues).
    numDim=size(pAzPt,1);%2 or 3.
    J=repmat(eye(numDim,numDim),[1,1,numPoints]);
    return;
end

azStart=greatCircleAzimuth(latLonRefOld,latLonRefNew);
azOffset1=pi/2-azStart;
lambda0=greatCircleDistance(latLonRefOld,latLonRefNew,1);

if(size(pAzPt,1)==3)
    J=zeros(3,3,numPoints);
    J(3,3,:)=1;
else
    J=zeros(2,2,numPoints);
end

%Transform into the system with both sensors on the equator.
pAzPt(2,:)=pAzPt(2,:)+azOffset1;

p=pAzPt(1,:);
az=pAzPt(2,:);

cosPre=cos(p/rE);
sinPre=sin(p/rE);

cosLambda0=cos(lambda0);
sinLambda0=sin(lambda0);
sinAz=sin(az);
cosAz=cos(az);

denom2=1-(cosPre.*cosLambda0+sinPre.*sinAz.*sinLambda0).^2;
denom1=sqrt(denom2);

dp0dp=(cosLambda0.*sinPre-cosPre.*sinAz.*sinLambda0)./denom1;
dp0dAz=-rE.*sinPre.*cosAz.*sinLambda0./denom1;

dAz0dp=cosAz.*sinLambda0./(rE.*denom2);
dAz0dAz=sinPre.*(cosLambda0.*sinPre-cosPre.*sinAz.*sinLambda0)./denom2;

J(1,1,:)=dp0dp;
J(2,1,:)=dAz0dp;
J(1,2,:)=dp0dAz;
J(2,2,:)=dAz0dAz;

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
