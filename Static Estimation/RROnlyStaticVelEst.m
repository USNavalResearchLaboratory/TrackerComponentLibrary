function [vEst]=RROnlyStaticVelEst(rr,states1,state2,zTar)
%%RRONLYSTATICVELEST  Perform least-squares estimation of a target velocity
%                     vector given only range rate measurements from one
%                     receiver and multiple transmitters or from one
%                     transmitter and multiple receivers. This uses a
%                     non-relativistic model for the range rates and
%                     ignores possible atmospheric effects. 
%
%INPUTS:    rr  A numMeas X 1 vector of range rates. numMeas >=3 for the
%               velocity vector to be observable.
%           states1 A 6X numMeas matrix of stacked transmitter position and
%                   velocity vectors if multiple transmitter and one
%                   receiver are used or a matrix of stacked receiver
%                   positions and velocities if one transmitter and
%                   multiple receivers are used.states(:,n) is the state
%                   (3 position and 3 velocity components) corresponding to
%                   the Nth range rate measurement.
%           state2  A 6X1 vector of the receiver state (Cartesian position
%                   and velocity componets), if multiple transmitters and
%                   one receiver are used or a vector of the transmitter
%                   state if one transmitter and multiple receivers are
%                   used.
%           zTar    The 3X1 Cartesian position of the target.
%
%OUTPUTS:   vEst   The 3X1 least squares Cartesian velocity estimate of the
%                  target.
%
%This implements the least squares velocity vector estimation procedure of
%Section IV E of [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zTx=states1(1:3,:);
zRx=state2(1:3);
vTx=states1(4:6,:);
vRx=state2(4:6);

h=zTar-zRx;
h=h/norm(h);

hi=bsxfun(@minus,zTar,zTx);
hiNorm=sqrt(sum(hi.*hi,1));
hi=bsxfun(@rdivide,hi,hiNorm);

rDotB=rr+h'*vRx+sum(hi.*vTx,1)';
Hv=bsxfun(@plus,h',hi');

vEst=pinv(Hv)*rDotB;

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
