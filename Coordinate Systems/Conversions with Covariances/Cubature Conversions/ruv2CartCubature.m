function [zCart, RCart]=ruv2CartCubature(z,SR,useHalfRange,zTx,zRx,M,xi,w)
%RUV2CARTCUBATURE Use cubature integration to approximate the moments of
%                 measurements converted from bistatic r-u-v coordinates
%                 into Cartesian coordinates. For a two-way monostatic
%                 conversion, set zTx=[0;0;0]; to make the transmitter and
%                 receiver collocated.
%
%INPUTS:      z A 3XnumMeas matrix of numMeas vectors to convert. Each
%               has elements [r;u;v], where r is the bistatic range from
%               the transmitter to the target to the receiver, and u and v
%               are direction cosines.
%            SR The 3X3 lower-triangular square roots of the
%               measurement covariance matrices for the measurements.
% useHalfRange  A boolean value specifying whether the bistatic range value
%               has been divided by two. This normally comes up when
%               operating in monostatic mode, so that the range reported is
%               a one-way range. The default if this parameter is not
%               provided is false.
%           zTx The 3X1 [x;y;z] location vector of the transmitter in
%               global Cartesian coordinates.  If this parameter is omitted
%               or an empty matrix is passed, then the receiver is placed
%               at the origin.
%           zRx The 3X1 [x;y;z] location vector of the receiver in global
%               Cartesian coordinates. If this parameter is omitted or an
%               empty matrix is passed, then the receiver is placed at the
%               origin.
%             M A 3X3 rotation matrix to go from the alignment of the
%               global coordinate system to the local alignment of the
%               receiver. The z vector of the local coordinate system of
%               the receiver is the pointing direction of the receiver. If
%               this matrix is omitted, then the identity matrix is used.
%            xi A 3XnumCubaturePoints matrix of cubature points for the
%               numeric integration. If this and the final parameter are
%               omitted or empty matrices are passed, then
%               fifthOrderCubPoints is used to generate cubature points.
%             w A numCubaturePointsX1 vector of the weights associated
%               with the cubature points.
%
%OUTPUTS:   zCart The approximate means of the PDF of the Cartesian
%                 converted measurements in [x;y;z] Cartesian coordinates
%                 for each measurement. This is a 3XnumMeas matrix.
%           RCart The approximate 3X3XnumMeas set of covariance matrices of
%                 the PDFs of the Cartesian converted measurements.
%
%Details of the conversion are given in [1].
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
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(xi))
    [xi,w]=fifthOrderCubPoints(3);
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

%Perform the conversion
numCubaturePoints=size(xi,2);

numMeas=size(z,2);

if(size(SR,3)==1)
    SR=repmat(SR,[1,1,numMeas]);
end

zCart=zeros(3,numMeas);
RCart=zeros(3,3,numMeas);
for curMeas=1:numMeas
    %Transform the cubature points to match the given Gaussian.
    cubPoints=transformCubPoints(xi,z(:,curMeas),SR(:,:,curMeas));

    %Convert all of the points into Cartesian space
    for curPoint=1:numCubaturePoints
        cubPoints(:,curPoint)=ruv2Cart(cubPoints(:,curPoint),useHalfRange,zTx,zRx,M);
    end

    %Extract the first two moments of the transformed points.
    [zCart(:,curMeas),RCart(:,:,curMeas)]=calcMixtureMoments(cubPoints,w);
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
