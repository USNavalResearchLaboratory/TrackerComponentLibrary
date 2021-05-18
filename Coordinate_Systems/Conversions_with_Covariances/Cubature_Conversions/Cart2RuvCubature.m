function [zRuv, RRuv]=Cart2RuvCubature(z,SR,useHalfRange,zTx,zRx,M,includeW,xi,w)
%CART2RUVCUBATURE Use cubature integration to approximate the moments of
%                 measurements converted from Cartesian coordinates into
%                 bistatic r-u-v coordinates. For a two-way monostatic
%                 conversion, set zTx=[0;0;0]; to make the transmitter and
%                 receiver collocated.
%
%INPUTS: z A 3XnumMeas matrix of Cartesian points in global [x;y;z]
%          Cartesian coordinates that are to be converted.
%       SR The 3X3XnumMeas lower-triangular square root of the measurement
%          covariance matrices for each measurement. If all of the matrices
%          are the same, then a single 3X3 matrix can be passed.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided is false.
%      zTx The 3X1 [x;y;z] location vector of the transmitter in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%      zRx The 3X1 [x;y;z] location vector of the receiver in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%        M A 3X3 rotation matrix to go from the alignment of the global
%          coordinate system to the local alignment of the receiver. The z
%          vector of the local coordinate system of the receiver is the
%          pointing direction of the receiver. If this matrix is omitted,
%          then the identity matrix is used.
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
%            include the third direction component.
%          1 Include the third direction component and perform the
%            averaging over the directions as a mean direction. A linear
%            covariance matrix is computed with respect to the mean
%            direction values.
%          2 Include the third component and just perform weighted linear
%            averaging. The three directional components will no longer
%            have unit magnitude.
%       xi A 3 X numCubaturePoints matrix of cubature points for the
%          numeric integration. If this and the final parameter are omitted
%          or empty matrices are passed, then fifthOrderCubPoints is used
%          to generate cubature points.
%        w A numCubaturePoints X 1 vector of the weights associated
%          with the cubature points.
%
%OUTPUTS: zRuv The approximate means of the PDF of the measurements
%              in bistatic [r;u;v] coordinates or [r;u;v;w] coordinates if
%              includeW is true. This is a 3XnumMeas  or a 4XnumMeas
%              matrix.
%         RRuv The approximate 3X3XnumMeas or 4X4XnumMeas covariance
%              matrices of the PDF of the bistatic [r;u;v] or [r;u;v;w]
%              converted measurements. This is a 3X3XnumMeas or a
%              4X4XnumMeas hypermatrix. If includeW is true, then these
%              matrices should be assumed to be singular.
%
%Details of the conversion are given in [1].
%
%EXAMPLE:
%Here, we verify that the covariance matrices of the converted measurements
%are consistent.
% pointCart=1e4*[1.173024396049598;
%               -4.844345843776918;
%                3.969150579064054];
% 
% RCart=[50^2,     0,  -25;
%           0,  50^2,    0;
%         -25,     0, 50^2];
% SCart=chol(RCart,'lower');
% 
% zTx=[0;0;0];
% zRx=1e4*[1;2;0];
% M=eye(3,3);
% 
% useHalfRange=false;
% 
% ruvTrue=Cart2Ruv(pointCart,useHalfRange,zTx,zRx,M);
% 
% %Now, we convert the noisy Cartesian point into r-u-v coordinates and
% %make sure that the NEES is consistent.
% numRuns=1000;
% includeW=false; 
% 
% zCart=bsxfun(@plus,pointCart,SCart*randn(3,numRuns));
% [zConv,RConv]=Cart2RuvCubature(zCart,SCart,useHalfRange,zTx,zRx,M,includeW);
% calcNEES(ruvTrue,zConv,RConv)
% %One will see that the NEES is near 1 indicating that the conversion is
% %consistent.
% %If we include the w component in the output, then we have to discard the
% %redundancy
% includeW=true;
% [zConv,RConv]=Cart2RuvCubature(zCart,SCart,useHalfRange,zTx,zRx,M,includeW);
% calcNEES(ruvTrue,zConv(1:3,:),RConv(1:3,1:3,:))
% %One find the NEES to be near 1, which indicates covariance consistency.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%May 2015 David Karnick, with major changes by David F. Crouse February
%2017, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8||isempty(xi))
    [xi,w]=fifthOrderCubPoints(3);
end

if(nargin<7||isempty(includeW))
   includeW=false; 
end

if(nargin<6||isempty(M))
    M=eye(3);
end

if(nargin<5||isempty(zRx))
    zRx=zeros(3,1);
end

if(nargin<4||isempty(zTx))
    zTx=zeros(3,1);
end

if(nargin<3||isempty(useHalfRange))
    useHalfRange=false;
end

numMeas=size(z,2);
if(size(SR,3)==1)
    SR=repmat(SR,[1,1,numMeas]);
end

hasW=(includeW~=0);

if(hasW)
    numDim=4;
else
    numDim=3;
end

zRuv=zeros(numDim,numMeas);
RRuv=zeros(numDim,numDim,numMeas);
if(includeW==0||includeW==2)
    h=@(x)Cart2Ruv(x,useHalfRange,zTx,zRx,M,hasW);

    for curMeas=1:numMeas
        [zRuv(:,curMeas), RRuv(:,:,curMeas)]=calcCubPointMoments(z(:,curMeas),SR(:,:,curMeas),h,xi,w);
    end
else
    numCub=size(xi,2);
    for curMeas=1:numMeas
        %Transform the cubature points to match the given Gaussian.
        cubPoints=transformCubPoints(xi,z(:,curMeas),SR(:,:,curMeas));
        %Convert all of the points.
        cubPointsConv=Cart2Ruv(cubPoints,useHalfRange,zTx,zRx,M,hasW);

        %Given r-u-v-w coordinates, the mean direction is just the linear
        %average projected onto the unit sphere.

        %Deal with any cubature points that might have had zero magnitude.
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

        zRuv(:,curMeas)=zMean;
        %Now, just compute a standard 3D linear covariance matrix.
        RCov=zeros(4,4);
        for curCub=1:numCub
            diff=cubPointsConv(:,curCub)-zMean(:);

            RCov=RCov+w(curCub)*(diff*diff');
        end
        RRuv(:,:,curMeas)=RCov;
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
