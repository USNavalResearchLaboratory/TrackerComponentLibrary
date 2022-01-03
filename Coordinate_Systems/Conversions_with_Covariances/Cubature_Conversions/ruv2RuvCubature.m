function [zConv,RConv]=ruv2RuvCubature(z,SR,useHalfRange,zTx1,zRx1,M1,zTx2,zRx2,M2,includeW,xi,w)
%%RUV2RUVCUBATURE Use cubature integration to approximate the moments of a
%         measurement converted from bistatic r-u-v (or r-u-v-w) 
%         coordinates with respect to one bistatic radar pair into those 
%         for another. This function can be useful for converting
%         bistatically received measurements into the coordinate system of
%         the transmitter so that they can be gated with the transmit beam
%         of a radar. If r-u-v coordinates are used, the target is assumed
%         to be in front of the receiver (local z coordinate is positive).
%         r-u-v coordinates consist of a bistatic range and direction
%         cosines at the receiver. The "direction cosines" u and v are just
%         the x and y coordinates of a unit vector from the receiver to the
%         target in the coordinate system at the receiver. This basically
%         assumes that the boresight direction of the receiver is the z 
%         axis. Assuming the target is in front of the receiver, the third
%         unit vector coordinate is not needed. However, r-u-v-w
%         coordinates include the third component.
%
%INPUTS: z A 3XN matrix of vectors with elements [r;u;v], where r is the
%          bistatic range from the transmitter to the target to the
%          receiver, and u and v are direction cosines. Each u,v pair
%          should have a magnitude less than or equal to one. If the
%          magnitude is greater than one, then the pair is normalized
%          before conversion to avoid imaginary results. It is assumed
%          that all directions are in front of the receiver (the third
%          component of the direction unit vector is positive).
%       SR The 3X3XN lower-triangular square root of the measurement
%          covariance matrices for the measurements. If all of the
%          matrices are the same, then this can just be a single matrix.
% useHalfRange A scalar boolean value or a 2X1 or 1X2 vector of boolean
%          values specifying whether the bistatic range value should be
%          divided by two. This normally comes up when operating in
%          monostatic mode, so that the range reported is a one-way range.
%          The default if an empty matrix is provided is false.
%          useHalfRange(1) applies to the first (source) bistatic channel
%          and useHalfRange(2) to the second (destination) bistatic
%          channel. If a scalar is passed, then both values are taken to
%          be the same.
%     zTx1 The 3X1 [x;y;z] location vector of the transmitter originating
%          the measurements z in global Cartesian coordinates. zTx1 can
%          have more than 3 rows; additional rows are ignored.
%     zRx1 The 3X1 [x;y;z] location vector of the receiver in Cartesian
%          coordinates that produced the measurements z. zRx1 can
%          have more than 3 rows; additional rows are ignored.
%       M1 A 3X3 rotation matrix to go from the alignment of the global
%          coordinate system to the local alignment of the original
%          sensor. The z vector of the local coordinate system of the
%          sensor is the pointing direction of the sensor.
% zTx2,zRx2,M2 These are the same as zTx1,zRx1, and M1, but for the
%          coordinate system into which the bistatic r-u-v(-w)
%          measurements should be converted.
% includeW A value indicating whether a third direction cosine component
%          should be included in the output and how to perform the
%          averaging. This helps for determining whether the converted
%          measurement is located behind the receiving radar. The u and v
%          direction cosines are two parts of a 3D unit vector. Generally,
%          one might assume that the target is in front of the sensor, so
%          the third component would be positive and is not needed.
%          However, the third component can be included if ambiguity
%          exists. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Do not
%             include the third direction component.
%          1 Include the third direction component and perform the
%            averaging over the directions as a mean direction. A linear
%            covariance matrix is computed with respect to the mean
%            direction values.
%          2 Include the third component and just perform weighted linear
%            averaging. The three directional components will no longer
%            have unit magnitude.
%       xi A 3 X numCubaturePoints (or 4X numCubaturePoints for r-u-v-w)
%          matrix of cubature points for the numeric integration. If this
%          and the final parameter are omitted or empty matrices are
%          passed, then fifthOrderCubPoints is used to generate cubature
%          points.
%       w  A numCubaturePoints X 1 vector of the weights associated with
%          the cubature points.
%
%OUTPUTS: zConv The approximate mean of the PDF of the measurement
%               in bistatic [r;u;v] (or [r;u;v;w]) coordinates for each
%               measurement in the coordinate system of the second bistatic
%               path. This is a 3XN or 4XN matrix.
%         RConv The approximate 3X3XN (or 4X4XN) covariance
%               matrices of the PDF of the bistatic [r;u;v] converted
%               measurements. If the w component is included, so that the
%               output is 4X4XN, then the matrices should be assumed to be
%               singular.
%
%Given just r-u-v coordinates, the conversion is performed using linear
%sums as in [1]. If a very large directional uncertainty exists, then the
%use of a linear covariance matrix for directional data is inappropriate.
%
%Note that when given only r-u-v and not r-u-v-w coordinates, an implicit
%assumption is that the measurement is not extremely close to being outside
%of the receiver's viewing area. For many practical radars having around a
%+/-60 degree field of view and angular measurement standard deviations
%less than one degress (once u-v deviations are converted to angles), then
%this is not a problem at all. However, given a measurement very close to
%+/- 90 degrees from the boresight, when cubature points are added for the
%integration, this can push u^2+v^2 above 1. At such extremes, however, the
%Gaussian approximation tends to be invalid anyway. Of course, utilizing
%r-u-v-w coordinates (with an often singular covariance matrix) can cause
%this issue to be avoided.
%
%EXAMPLE:
%Here, we verify that the covariance matrices of the converted measurements
%are consistent.
% pointCart=1e4*[1.173024396049598;
%               -4.844345843776918;
%                3.969150579064054];
% 
% RRuv=[50^2,   0,      0;
%        0,   1e-3^2,   0;
%        0,   0,      1e-3^2];
% SRuv1=chol(RRuv,'lower');
% 
% xTx1=[0;0;0];
% xRx1=1e4*[1;2;0];
% xTx2=xTx1;
% xRx2=xTx1;
% M1=eye(3,3);
% M2=M1;
% 
% useHalfRange=false;
% 
% ruv1True=Cart2Ruv(pointCart,useHalfRange,xTx1,xRx1,M1);
% ruv2True=Cart2Ruv(pointCart,useHalfRange,xTx2,xRx2,M2);
% 
% %Now, we convert a noisy ruv1True into the coordinate system for the
% %second bistatic path and check its consistency with ruv2True.
% numRuns=1000;
% includeW=false; 
%
% zRuv1=bsxfun(@plus,ruv1True,SRuv1*randn(3,numRuns));
% [zConv,RConv]=ruv2RuvCubature(zRuv1,SRuv1,useHalfRange,xTx1,xRx1,M1,xTx2,xRx2,M2,includeW);
% calcNEES(ruv2True,zConv,RConv)
% %One will see that the NEES is near 1 indicating that the conversion is
% %consistent.
% %If we include the w component in the output, then we have to discard the
% %redundancy to assess the NEES, because the covariance matrix is
% %otherwise singular.
% includeW=true;
% [zConv,RConv]=ruv2RuvCubature(zRuv1,SRuv1,useHalfRange,xTx1,xRx1,M1,xTx2,xRx2,M2,includeW);
% calcNEES(ruv2True,zConv(1:3,:),RConv(1:3,1:3,:))
% %Again, one will get an NEES near 1.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%May 2015 David Karnick, with major changes by David F. Crouse February
%2017, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(z,2);

if(includeW)
    numDim=4;
else
    numDim=3;
end

if(isscalar(useHalfRange))
    useHalfRange=[useHalfRange;useHalfRange];
end

if(nargin<11||isempty(xi))
    [xi,w]=fifthOrderCubPoints(3);
end

if(nargin<10||isempty(includeW))
    includeW=false;
end

if(size(SR,3)==1)
    SR=repmat(SR,[1,1,numMeas]);
end

zConv=zeros(numDim,numMeas);
RConv=zeros(numDim,numDim,numMeas);
if(includeW==0||includeW==2)
        h=@(x)ruv2Ruv(x,useHalfRange,zTx1,zRx1,M1,zTx2,zRx2,M2,includeW);
        for curMeas=1:numMeas
            [zConv(:,curMeas), RConv(:,:,curMeas)]=calcCubPointMoments(z(:,curMeas),SR(:,:,curMeas),h,xi,w);
        end
elseif(includeW==1)
        numCub=length(w);
        for curMeas=1:numMeas
            %Transform the cubature points to match the given Gaussian.
            cubPoints=transformCubPoints(xi,z(:,curMeas),SR(:,:,curMeas));

            %Convert all of the points.
            cubPointsConv=ruv2Ruv(cubPoints,useHalfRange,zTx1,zRx1,M1,zTx2,zRx2,M2,includeW);

            %Given r-u-v-w coordinates, the mean direction is just the
            %linear average projected onto the unit sphere.

            %Deal with any cubature points that might have had zero
            %magnitude.
            sel=any(~isfinite(cubPointsConv(2:end,:)),1);
            cubPointsConv(2,sel)=1;
            cubPointsConv(3,sel)=0;
            cubPointsConv(4,sel)=0;

            zMean=sum(bsxfun(@times,cubPointsConv,w'),2);
            %Normalize the result
            zMean(2:end)=zMean(2:end)/norm(zMean(2:end));
            %Deal with the case where the points are so spread out that the
            %linear mean is zero.
            if(any(~isfinite(zMean(2:end))))
                zMean(2:end)=[1;0;0];
            end

            zConv(:,curMeas)=zMean;
            %Now, just compute a standard 3D linear covariance matrix.
            RCov=zeros(4,4);
            for curCub=1:numCub
                diff=cubPointsConv(:,curCub)-zMean(:);

                RCov=RCov+w(curCub)*(diff*diff');
            end
            RConv(:,:,curMeas)=RCov;
        end
else
    error('Unknown value of includeW provided.')
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
