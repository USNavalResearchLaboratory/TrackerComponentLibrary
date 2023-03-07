function vEst=RROnlyStaticVelEst(rr,xTx,xRx,zTar)
%%RRONLYSTATICVELEST Perform unweighted least-squares estimation of a
%           target velocity vector given only range rate measurements from
%           different bistatic channels (or from multiple receivers if the
%           target is the transmitter). This function can work in 2D and
%           3D space and it will produce a least-squares estimate when more
%           than the minimum number of range rate measurements needed is
%           given. This uses a non-relativistic model for the range rates
%           and ignores possible atmospheric effects. 
%
%INPUTS: rr A numMeasX1 or 1XnumMeas vector of range rates. numMeas>=3 for
%           the velocity vector to be observable in 3D and numMeas>=2 in
%           2D.
%       xTx A 6XnumMeas matrix (in 3D or a 4XnumMeas matrix in 2D) of
%           stacked transmitter position and velocity vectors. If the
%           target is the transmitter, then an empty matrix should be
%           passed. xTx(:,n)=[x;y;z;xDot;yDot;zDot] is the state (3
%           position and 3 velocity components) corresponding to the nth
%           range rate measurement. If all of the transmitters are in the
%           same spot, then a single 6X1 or 4X1 vector can be passed.
%       xRx A 6XnumMeas matrix (in 3D or a 4XnumMeas matrix in 2D) of
%           stacked transmitter position and velocity vectors.
%           states1(:,n)=[x;y;z;xDot;yDot;zDot] is the state (3 position
%           and 3 velocity components) corresponding to the nth range rate
%           measurement. If all of the receivers are in the same spot, then
%           a single 6X1 or 4X1 vector can be passed.
%      zTar The 3X1 (in 3D) or 2X1 (in 2D) Cartesian position of the
%           target.
%
%OUTPUTS: vEst The 3X1 (in 3D) or 2X1 (in 2D) least squares Cartesian
%              velocity estimate of the target.
%
%This implements the least squares velocity vector estimation procedure of
%Equation 41 in Section IV E of [1]. The case of the transmitter being the
%target is handled specially (to get rid fo the singulaity).
%
%EXAMPLE 1:
%This is just a simple example showing that for three range rate values,
%one can get back the orignal velocity value in 3D.
% zTar=[0;40e3;40e3];
% vTar=[400;-200;100];
% xTx1=[100;10e3;3e3;50;50;-50];
% xTx2=[0;0;0;0;0;-20];
% xTx3=[10e3;10e3;3e3;100;-100;100];
% 
% xRx1=[-10e3;0;3e3;100;100;100];
% xRx2=[0;10e3;30;-80;-200;-20];
% xRx3=xRx2;
% 
% states1=[xTx1,xTx2,xTx3];
% states2=[xRx1,xRx2,xRx3];
% 
% xTar=[zTar;vTar];
% 
% useHalfRange=false;
% rr=zeros(3,1);
% rr(1)=getRangeRate(xTar,useHalfRange,xTx1,xRx1);
% rr(2)=getRangeRate(xTar,useHalfRange,xTx2,xRx2);
% rr(3)=getRangeRate(xTar,useHalfRange,xTx3,xRx3);
% 
% vEst=RROnlyStaticVelEst(rr,states1,states2,zTar)
%One will see that vEst=vzTar.
%
%EXAMPLE 2:
%This is the same as example 1, except the demonstration is now occurring
%in 2D.
% zTar=[0;40e3];
% vTar=[400;-200];
% xTx1=[100;10e3;50;50];
% xTx2=[0;0;0;0];
% xTx3=[10e3;10e3;100;-100];
% 
% xRx1=[-10e3;0;100;100];
% xRx2=[0;10e3;-80;-200];
% xRx3=xRx2;
% 
% states1=[xTx1,xTx2,xTx3];
% states2=[xRx1,xRx2,xRx3];
% 
% xTar=[zTar;vTar];
% 
% useHalfRange=false;
% rr=zeros(3,1);
% rr(1)=getRangeRate(xTar,useHalfRange,xTx1,xRx1);
% rr(2)=getRangeRate(xTar,useHalfRange,xTx2,xRx2);
% rr(3)=getRangeRate(xTar,useHalfRange,xTx3,xRx3);
% 
% vEst=RROnlyStaticVelEst(rr,states1,states2,zTar)
%One will see that vEst=vzTar.
%
%EXAMPLE 3:
%In this example, the range rate is one-way from the target (e.g. the
%target is an emitter). The relative error of the velocity estimate with
%the truth (noiseless scenario) is within finite precision limits.
% zTar=randn(3,1);
% vTar=randn(3,1);
% numMeas=5;
% xRx=randn(6,numMeas);
% xTar=[zTar;vTar];
% useHalfRange=false;
% rr=zeros(numMeas,1);
% for k=1:numMeas
%     rr(k)=getRangeRate(xTar,useHalfRange,xTar,xRx(:,k));
% end
% vEst=RROnlyStaticVelEst(rr,[],xRx,zTar);
% RelErr=max(abs((vTar-vEst)./vTar))
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=length(rr);
rr=rr(:);

if(size(xTx,2)==1)
    xTx=repmat(xTx,1,numMeas);
end

if(size(xRx,2)==1)
    xRx=repmat(xRx,1,numMeas);
end

posDim=length(zTar);

zRx=xRx(1:posDim,:);
vRx=xRx((posDim+1):(2*posDim),:);

h=bsxfun(@minus,zTar,zRx);
hNorm=sqrt(sum(h.*h,1));
h=bsxfun(@rdivide,h,hNorm);

if(~isempty(xTx))
    %The target is not the transmitter.
    zTx=xTx(1:posDim,:);
    vTx=xTx((posDim+1):(2*posDim),:);
    hi=bsxfun(@minus,zTar,zTx);
    hiNorm=sqrt(sum(hi.*hi,1));
    hi=bsxfun(@rdivide,hi,hiNorm);

    rDotB=rr+sum(h.*vRx,1)'+sum(hi.*vTx,1)';
    Hv=bsxfun(@plus,h',hi');
else
    %The target is the transmitter.
    rDotB=rr+sum(h.*vRx,1)';
    Hv=h';
end

vEst=lsqminnorm(Hv,rDotB);

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
