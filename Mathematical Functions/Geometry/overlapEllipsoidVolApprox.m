function vol=overlapEllipsoidVolApprox(A,z0,gammaVal,numSamp)
%%OVERLAPELLIPSOIDVOLAPPROX Determine the approximate volume of overlapping
%                   ellipsoids using a Monte Carlo sampling method. The ith
%                   ellipsoid is defined as the shape
%                   (z0(:,i)-z)'*A(:,:,i)*(z0(:,i)-z)=gammaVal
%
%INPUTS:  A A numDimXnumDimXN set of N positive definite matrices that
%           specify the size and shape of the ellipses or ellipsoids, where
%           a point zp is on the ith ellipse/ ellipsoid if
%           (zp-z(:,i))'*A(:,:,i)*(zp-z(:,i))=gammaVal.
%        z0 The numDimXN set of centers of the ellipsoids.
%  gammaVal A parameter specifying the size of the ellipse/ ellipsoid.
%           gammaVal must be positive. To specify a probability region of
%           probReg as is commonly used in tracking where A is a Gaussian
%           covariance matrix, one can get gammaVal from
%           ChiSquareD.invCDF(probReg,size(A,1)) The default if this
%           parameter is omitted or an empty matrix is passed is 1.
%   numSamp Optionally, the number of samples to use. If this parameter is
%           omitted or an empty matrix is passed, the default of 2e4 is
%           used.
%
%OUTPUTS: vol The approximate volume of the overlapping ellipsoids.
%
%The algorithm simply determines the smallest hyperrectangle that contains
%all of the ellipsoids. It then uniformly samples within the rectangle. The
%ratio of the points that fall in the ellipsoids to the total number
%samples times the volume of the rectangle is the approximate volume of the
%ellipsoids, which might be overlapping.
%
%EXAMPLE:
%Here, we have two circles in 2D.
% A=zeros(2,2,2);
% A(:,:,1)=inv(10*eye(2));
% A(:,:,2)=inv(10*eye(2));
% z0=[[5;0],[0;0]];
% gammaVal=16.2235;
% %Visualize how the circles overlap.
% figure(1)
% clf
% drawEllipse(z0,A,gammaVal)
% %Determine the volume
% vol=overlapEllipsoidVolApprox(A,z0,gammaVal)
%One should get a volume that is greater than the volume of any one circle,
%but less than their combined volumes.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(numSamp))
   numSamp=2e4; 
end

if(nargin>2&&~isempty(gammaVal))
   A=A/gammaVal;
end

numEllipses=size(z0,2);
numDim=size(z0,1);

invA=zeros(numDim,numDim,numEllipses);
for curEl=1:numEllipses
    invA(:,:,curEl)=inv(A(:,:,curEl));
end

ellipsVols=zeros(numEllipses,1);
for curEl=1:numEllipses
     ellipsVols(curEl)=ellipsoidVolume(A(:,:,curEl));
end

%Sorting the ellipsoids in descending order of volume can speed up the
%algorithms when generating a large number of samples.
[~,idx]=sort(ellipsVols,'descend');
z0=z0(:,idx);
A=A(:,:,idx);
invA=invA(:,:,idx);

%First, determine the minimum and maximum values in each dimension
%so as to form a bounding box around all of the ellipses.
minVals=Inf(numDim,1);
maxVals=-Inf(numDim,1);
for curEl=1:numEllipses
    curMin=z0(:,curEl)-diag(invA(:,:,curEl));
    curMax=z0(:,curEl)+diag(invA(:,:,curEl));

    minVals=min(minVals,curMin);
    maxVals=max(maxVals,curMax);
end

%Next, we draw uniform samples from the box and determine the
%fraction that gate to find the volume of the ellipses.
span=maxVals-minVals;

VB=prod(span);%The volume of the box.

gateCount=0;
for curSamp=1:numSamp
    z=minVals+span.*rand(numDim,1);

    for curEl=1:numEllipses
        diff=z-z0(:,curEl);

        if(diff'*A(:,:,curEl)*diff<=1)
            gateCount=gateCount+1;
            break;
        end
    end
end

%The approximate volume.
vol=VB*(gateCount/numSamp);

%The volume is not less than the volume of the largest ellipsoid.
vol=max(vol,max(ellipsVols));
%The volume is not more than the sum of the volumes of all ellipsoids.
vol=min(vol,sum(ellipsVols));
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
