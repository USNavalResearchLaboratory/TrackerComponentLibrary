function [zRuv, RRuv] = Cart2RuvCubature(z,SR,useHalfRange,zTx,zRx,M,xi,w)
%CART2RUVCUBATURE Use cubature integration to approximate the moments of
%                 measurements converted from Cartesian coordinates into
%                 bistatic r-u-v coordinates. For a two-way monostatic
%                 conversion, set zTx=[0;0;0]; to make the transmitter and
%                 receiver collocated.
%
%INPUTS:    z   A 3XnumMeas matrix of Cartesian points in global [x;y;z]
%               Cartesian coordinates that are to be converted.
%           SR  The 3X3XnumMeas lower-triangular square root of the
%               measurement covariance matrices for each measurement.
% useHalfRange  A boolean value specifying whether the bistatic range value
%               should be divided by two. This normally comes up when
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
%           xi  A 3 X numCubaturePoints matrix of cubature points for the
%               numeric integration. If this and the final parameter are
%               omitted or empty matrices are passed, then
%               fifthOrderCubPoints is used to generate cubature points.
%           w   A numCubaturePoints X 1 vector of the weights associated
%               with the cubature points.
%
%OUTPUTS: zRuv The approximate means of the PDF of the the measurements
%              in bistatic [r;u;v] coordinates. This is a 3XnumMeas matrix.
%         RRuv The approximate 3X3XnumMeas covariance matrices of the
%              PDF of the bistatic [r;u;v] converted measurements. This is
%              a 3X3XnumMeas hypermatrix.
%
%Details of the conversion are given in [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%May 2015 David Karnick, Naval Research Laboratory, Washington D.C.
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

numMeas=size(z,2);
if(size(SR,3)==1)
    SR=repmat(SR,[1,1,numMeas]);
end

zRuv=zeros(3,numMeas);
RRuv=zeros(3,3,numMeas);
for curMeas=1:numMeas
    %Transform the cubature points to match the given Gaussian.
    cubPoints=transformCubPoints(xi,z(:,curMeas),SR(:,:,curMeas));

    %Convert all of the points into RUV space
    cubPointsRUV=Cart2Ruv(cubPoints,useHalfRange,zTx,zRx,M);
    
    %Extract the first two moments of the transformed points.
    [zRuv(:,curMeas),RRuv(:,:,curMeas)]=calcMixtureMoments(cubPointsRUV,w);
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
