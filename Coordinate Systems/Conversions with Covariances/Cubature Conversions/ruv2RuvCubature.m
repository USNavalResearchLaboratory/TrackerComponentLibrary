function [zConv,RConv]=ruv2RuvCubature(z,SR,useHalfRange,zTx1,zRx1,M1,zTx2,zRx2,M2,xi,w)
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
%INPUTS: z  A 3XN matrix of vectors with elements [r;u;v], where r is the
%           bistatic range from the transmitter to the target to the
%           receiver, and u and v are direction cosines. Each u,v pair
%           should have a magnitude less than or equal to one. If the
%           magnitude is greater than one, then the pair is normalized
%           before conversion to avoid imaginary results. Alternatively,
%           one can pass a 4XN matrix of [r;u;v;w] vectors where [u;v;w]
%           form a full unit vector in the receiver's local 3D Cartesian
%           coordinates.
%       SR  The 3X3XN (or 4X4XN or r-u-v-w coordinates) lower-triangular
%           square root of the measurement covariance matrices for each
%           measurement. When dealing with r-u-v-w coordinates, this is
%           generally going to be singular (positive semi-definite). The
%           cholSemiDef function can be used to find the lower-triangular
%           square root of such a matrix. If all of the matrices are the
%           same, then this can just be a single matrix.
%useHalfRange A boolean value specifying whether the bistatic range value
%           should be divided by two. This normally comes up when operating
%           in monostatic mode, so that the range reported is a one-way
%           range. The default if an empty matrix is provided is false.
%      zTx1 The 3X1 [x;y;z] location vector of the transmitter originating
%           the measurements z in global Cartesian coordinates. zTx1 can
%           have more than 3 rows; additional rows are ignored.
%      zRx1 The 3X1 [x;y;z] location vector of the receiver in Cartesian
%           coordinates that produced the measurements z. zRx1 can
%           have more than 3 rows; additional rows are ignored.
%        M1 A 3X3 rotation matrix to go from the alignment of the global
%           coordinate system to the local alignment of the original
%           sensor. The z vector of the local coordinate system of the
%           sensor is the pointing direction of the sensor.
% zTx2,zRx2,M2 These are the same as aTx1,zRx1, and M1, but for the
%           coordinate system into which the bistatic r-u-v(-w)
%           measurements should be converted.
%        xi A 3 X numCubaturePoints (or 4X numCubaturePoints for r-u-v-w)
%           matrix of cubature points for the numeric integration. If this
%           and the final parameter are omitted or empty matrices are
%           passed, then fifthOrderCubPoints is used to generate cubature
%           points.
%        w  A numCubaturePoints X 1 vector of the weights associated with
%           the cubature points.
%
%OUTPUTS:zConv The approximate mean of the PDF of the the measurement
%              in bistatic [r;u;v] (or [r;u;v;w]) coordinates for each
%              measurement in the coordinate system of the second bistatic
%              path. This is a 3XnumMeas or 4XnumMeas matrix.
%        RConv The approximate 3X3XnumMeas (or 4X4XnumMeas) covariance
%              matrices of the PDF of the bistatic [r;u;v] converted
%              measurements.
%
%Given just r-u-v coordinates, the conversion is performed using linear
%sums as in [1]. On the other hand, when given r-u-v-w coordinates, the
%directional components of the transformed cubature points are projected 
%onto the unit sphere prior to computing the mean and covariance matrix.
%The directional components of the mean (which is computed in a linear
%manner) is projected onto the unit sphere to provide a valid direction
%vector and that is the vector used for computing a linear covariance
%matrix to return. These ad-hoc fixed for the directional nature of
%r-u-v-w coordinates are conceptually similar to how angular coordinates
%are handled in [2]. If a very large directional uncertainty exists, then
%the use of a linear covariance matrix for directional data is
%inappropriate.
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
%r-u-v-w coordiantes (with an often singular covariance matrix) can cause
%this issue to be avoided.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%May 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(z,1);
numMeas=size(z,2);

if(nargin<10||isempty(xi))
    [xi,w]=fifthOrderCubPoints(numDim);
end

if(size(SR,3)==1)
    SR=repmat(SR,[1,1,numMeas]);
end

numCub=length(w);
zConv=zeros(numDim,numMeas);
RConv=zeros(numDim,numDim,numMeas);
for curMeas=1:numMeas
    %Transform the cubature points to match the given Gaussian.
    cubPoints=transformCubPoints(xi,z(:,curMeas),SR(:,:,curMeas));

    %Convert all of the points.
    cubPointsConv=ruv2Ruv(cubPoints,useHalfRange,zTx1,zRx1,M1,zTx2,zRx2,M2);
    
    if(numDim==3)
        %Just extract the first two moments of the measurement using linear
        %means if it is not an r-u-v-w measurement.
        [zConv(:,curMeas),RConv(:,:,curMeas)]=calcMixtureMoments(cubPointsConv,w);
    else
        %Given r-u-v-w coordinates, the mean direction is just the linear
        %average, but it must be normalized. We will first project to
        %cubature points onto the unit sphere.
        cubPointsConv(2:end,:)=bsxfun(@rdivide,cubPointsConv(2:end,:),sqrt(sum(cubPointsConv(2:end,:).^2,1)));
        
        %Deal with any cubature points that might have had zero magnitude.
        sel=any(~isfinite(cubPointsConv(2:end,:)),1);
        cubPointsConv(2,sel)=1;
        cubPointsConv(3,sel)=0;
        cubPointsConv(4,sel)=0;
        
        zMean=sum(bsxfun(@times,cubPointsConv,w'),2);
        zMean(2:end)=zMean(2:end)/sum(zMean(2:end));
        %Deal with the case where the points are so spread out that the
        %linear mean is zero.
        if(any(~isfinite(zMean(2:end))))
            zMean(2:end)=[1;0;0];
        end
        
        zConv=zMean;
        %Now, just compute a standard linear covariance matrix.
        RCov=zeros(4,4);
        for curCub=1:numCub
            diff=cubPointsConv(:,curCub)-zMean;
            
            RCov=RCov+w(curCub)*(diff*diff');
        end
        RConv(:,:,curMeas)=RCov;
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
