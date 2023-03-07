function [zCart,RCart]=ruv2CartStdRefracCubature(zRUVBiased,SR,useHalfRange,zTx,zRx,M,Ns,xi,w,ce,rE,spherCent,xMax)
%%RUV@CARTSTDREFRACCUBATURE Use cubature integration to approximate the 
%              moments of measurements converted from refraction-corrupted
%              bistatic r-u-v coordinates into Cartesian coordinates using
%              a standard exponential atmospheric model. For a two-way
%              monostatic conversion, set zTx=[0;0;0]; to make the
%              transmitter and receiver collocated.
%
%INPUTS: z A 3XnumMeas matrix of numMeas vectors to convert. Each has
%          elements [r;u;v], where r is the bistatic range from the
%          transmitter to the target to the receiver, and u and v are
%          direction cosines.
%       SR The 3X3XnumMeas lower-triangular square roots of the measurement
%          covariance matrices for the measurements. If all of the matrices
%          are the same, then this can just be a single 3X3 matrix.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided (or an
%          empty matrix is provided) is false.
%      zTx The 3X1 [x;y;z] location vector of the transmitter in global
%          Cartesian coordinates.  If this parameter is omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%      zRx The 3X1 [x;y;z] location vector of the receiver in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%        M A 3X3 rotation matrix to go from the alignment of the global
%          coordinate system to the local alignment fo the receiver. The z
%          vector of the local coordinate system of the receiver is the
%          pointing direction of the receiver. If this matrix is omitted,
%          then the identity matrix is used.
%       Ns The atmospheric refractivity reduced to the reference sphere.
%          Note that the refractivity is (n-1)*1e6, where n is the index
%          of refraction. The function reduceStdRefrac2Spher can be used
%          to reduce a refractivity to the surface of a reference
%          ellipsoid. This function does not allow different
%          refractivities to be used as the transmitter and receiver. If
%          this parameter is omitted or an empty matrix is passed, a
%          default value of 313 is used.
%       xi A 3 X numCubaturePoints matrix of cubature points for the
%          numeric integration. If this and the next parameter are omitted
%          or empty matrices are passed, then fifthOrderCubPoints is used
%          to generate cubature points.
%        w A numCubaturePoints X 1 vector of the weights associated with
%          the cubature points.
%       ce The optional decay constant of the exponential model. The
%          refractivity N at height h is N=Ns*exp(-ce*(h-h0)) where h0 is
%          the reference height (in this function, the height of the
%          reference ellipsoid surface is used). ce is related to the
%          change in refractivity at an elevation of 1km based on the
%          refractivity at sea level as
%          ce=log(Ns/(Ns+DeltaN))/1000;%Units of inverse meters.
%          where the change in refractivity for a change in elevation of
%          1km is DeltaN=-multConst*exp(expConst*Ns); In [1], standard
%          values for the two constants are expConst=0.005577; and
%          multConst=7.32; If ce is omitted or an empty matrix is passed,
%          the value based on the standard model is used.
% rE,spherCent The radius of the Earth to use for the spherical Earth
%           approximation used in the model and also the offset between the
%           global model and the local spherical model. It is assumed that
%           zC,zTx,and zRx are all given in the global model and will need
%           to be transformed to the local model to the used. If rE is
%           omitted or an empty matrix is passed, then the default of
%           [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx)) is
%           used. The defaults here mean that a WGS-84 reference ellipsoid
%           is approximated by the local osculating sphere.
%     xMax This function traces the ray to a maximum displacement in the
%          local tangent plane prior (or vertically for nearly vertical
%          directions) to attempting to find a particular range. This is an
%          optional parameter specifying the maximum distance to trace the
%          ray. The default if this parameter is omitted or an empty matrix
%          is passed is 1000e3 (1000km).
%
%OUTPUTS: zCart The approximate means of the PDF of the measurements
%               in global Cartesian [x;y;z] coordinates. This is a
%               3XnumMeas matrix.
%         RCart The approximate 3X3XnumMeas covariance matrices of the
%               PDF of the Cartesian converted measurements. This is a
%               3X3XnumMeas hypermatrix.
%
%The basic cubature conversion approach is detailed in [1] and [2]. The
%standard exponential measurement model is from [3] and the model is
%discussed in more detail in the comments to ruv2CartStdRefrac.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using 3D monostatic and bistatic
%    measurements in refractive environments," IEEE Aerospace and
%    Electronic Systems Magazine, vol. 29, no. 8, Part II, pp. 54-75, Aug.
%    2014.
%[2] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[3] B. R. Bean and G. D. Thayer, CRPL Exponential Reference Atmosphere.
%    Washington, D.C.: U. S. Department of Commerce, National Bureau of
%    Standards, Oct. 1959. [Online]. Available:
%    http://digicoll.manoa.hawaii.edu/techreports/PDF/NBS4.pdf
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(zRUVBiased,2);

if(nargin<3||isempty(useHalfRange))
    useHalfRange=false;
end

if(nargin<4||isempty(zTx))
   zTx=zeros(3,1); 
end

if(nargin<5||isempty(zRx))
   zRx=zeros(3,1); 
end

if(nargin<6||isempty(M))
   M=eye(3,3); 
end

if(nargin<7||isempty(Ns))
   Ns=313;
end

if(nargin<8||isempty(xi))
    [xi,w]=fifthOrderCubPoints(3);
end

if(nargin<10||isempty(ce))
    expConst=0.005577;
    multConst=7.32;

    %The change in refractivity at an elevation of 1km based on the
    %refractivity on the surface of the Earth.
    DeltaN=-multConst*exp(expConst*Ns);
    ce=log(Ns/(Ns+DeltaN))/1000;%Units of inverse meters.
end

if(nargin<11||isempty(rE))
    %Use the radius of the Earth that is the radius of the osculating
    %sphere at the location of the observed. This will be the radius used
    %in the local spherical Earth approximation for computing atmospheric
    %refraction. This uses the WGS-84 reference ellipsoid.
    [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx));
end

if(nargin<13||isempty(xMax))
    xMax=1000e3;%1000 kilometer assumed maximum x displacement.
end

if(size(SR,3)==1)
    SR=repmat(SR,[1,1,numMeas]);
end

zCart=zeros(3,numMeas);
RCart=zeros(3,3,numMeas);
for curMeas=1:numMeas
    %Transform the cubature points to match the given Gaussian.
    cubPoints=transformCubPoints(xi,zRUVBiased(:,curMeas),SR(:,:,curMeas));

    %Convert all of the points into Cartesian space
    cubPoints=ruv2CartStdRefrac(cubPoints,useHalfRange,zTx,zRx,M,Ns,ce,rE,spherCent,xMax);

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
