function [x,P]=twoPointDiffInit(T,z,R,q)
%TWOPOINTDIFFINIT State vector and covariance matrix initialization using a
%                 two-point differencing method. This assumes that a linear
%                 dynamic model (change in position equals velocity times
%                 time plus noise) is used and the target state consists
%                 of only position and velocity components in 2D or 3D.
%                 Optionally, process noise can be taken into account.
%
%INPUTS: T The time delta between the measurements, in seconds. The
%          transition matrix is assumed to have the form of a matrix with
%          position and velocity components as in FPolyKal.
%        z The zDimX2XN matrix of Cartesian position measurements. zDim can
%          be 2 or 3. The second measurement is assumed to arrive at time T
%          after the first. The N dimensions allows for multiple
%          simultaneous conversions.
%        R The zDimXzDimX2XN hypermatrix of measurement covariance
%          matrices. If just a zDimXzDim matrix is passed, then it is
%          assumed that both measurements have the same measurement matrix
%          for all N conversions. If just a zDimXzDimX2 matrix is passed,
%          then it is assumed that the same set of R matrices is used in
%          each of the N conversions.
%        q The optional process noise parameter. If this parameter is
%          omitted or an empty matrix is passed, then it is assumed there
%          is no process noise. A discretized dynamic model is used, so q
%          has units of length^2/time^3 and is the power spectral density
%          parameters as in the function QPolyKal.
%
%OUTPUTS: x The (2*zDim)XN initialized state vectors consisting of position
%           and velocity with all position components coming before any of
%           the velocity components.
%         P The (2*zDim)X(2*zDim)XN covariance matrixces, one for each of
%           the N initiations.
%
%Two point differencing assuming a linear dynamic model without process
%noise is given in Equations 39 and 40 in [1]. The inclusion of process
%noise assuming a discretized dynamic model is given in Equation 56 of [1].
%
%REFERENCES:
%[1] Mallick, M.,La Scala, B., "Comparison of single-point and two-point
%    difference track initiation algorithms using position measurements". 
%    Acta Automatica Sinica, 2008.
%
%March 2016 AJ Rivera, Naval Rearch Laboratory, Washington D.C.
%Process noise added May 2016 by David F. Crouse, Naval Rearch Laboratory,
%Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(q))
    q=0;
end

zDim=size(z,1);
if(size(R,3)==1)
   R=repmat(R,[1,1,2]);
end

N=size(z,3);
if(size(R,4)==1)
   R=repmat(R,[1,1,1,N]);
end

xDim=2*zDim;
x=zeros(xDim,N);
P=zeros(xDim,xDim,N);

posIdx=1:zDim;
velIdx=(zDim+1):(2*zDim);

for idx=1:N
    %Equation 39 in [1].
    %The position components
    x(posIdx,idx)=z(:,2,idx);
    %The velocity components
    x(velIdx,idx)=(z(:,2,idx)-z(:,1,idx))/T;

    %Equation 40 in [1].
    P(posIdx,posIdx,idx)=R(:,:,2,idx);
    P(posIdx,velIdx,idx)=R(:,:,2,idx)/T;
    P(velIdx,posIdx,idx)=R(:,:,2,idx)/T;
    P(velIdx,velIdx,idx)=sum(R(:,:,:,idx),3)/T^2;

    %The bias due to process noise as given in Equation 56 of [1].
    P(velIdx,velIdx,idx)=P(velIdx,velIdx,idx)+(1/3)*q*T*eye(zDim,zDim);
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
