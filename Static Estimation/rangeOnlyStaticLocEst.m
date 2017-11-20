function [xEst,PTaylor,PCRLB]=rangeOnlyStaticLocEst(rBi,zLoc1,zLoc2,RCov)
%%RANGEONLYSTATICLOCEST  Given multiple bistatic range-only measurements
%                        consisting of one receiver and multiple
%                        transmitters (or vice versa), find the location of
%                        a target in 3D. Covariance matrices for the
%                        estimate are also available. A minimum of 3
%                        receivers (or transmitters) is necessary.
%                        
%INPUTS: rBi A numMeas X 1 vector of bistatic range measurements where
%            numMeas>3 to ensure observability.
%      zLoc1 A 3 X numMeas vector of transmitter location vectors, if one
%            receiver and multiple transmitters are used. If one
%            transmitter and multiple receivers are used, then zLoc1
%            should be a 3 X numMeas vector of receiver locations.
%      zLoc2 A 3X1 vector of the receiver locations, if multiple
%            transmitters are used, or a 3X1 vector of the transmitter
%            location if multiple receivers are used.
%       RCov If a covariance matrix for the estimate is desired on the
%            output, then RCov, a measurement covariance matrix must be
%            provided. Otherwise, RCov can be omitted as the algorithm
%            implicitly assumes that all of the measurements have the same
%            accuracy.
%
%OUTPUTS: xEst The Cartesian target location estimate. If onlt the minimum
%              number of receivers/ transmitters is provided (3), then
%              there are two columns, each corresponding to one of the two
%              solutions.
%      PTaylor A covariance matrix for the estimate based on a Taylor
%              series expansion. then this is a 3D matrix, with the extra
%              dimension selecting the solution for which the covariance
%              matrix is valid.
%        PCRLB A covariance matrix (or matrices for two solutions) for the
%              estimates) based on the Cramér-Rao lower bound.
%
%The target localization algorithm is based on the spherical intersection
%method of [1], where the Taylor series expansion covariance matrix for the
%estimate is also derived. The PCRLB covariance estimate is standard
%assuming additive Gaussian noise.
%
%REFERENCES:
%[1] Malanowski, M, and Kulpa, K., "Two methods for target localization in
%    multistatic passive radar," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 48, no. 1, pp. 572-580, Jan. 2012.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Move the receiver to the origin.
zTx=bsxfun(@plus,zLoc1,zLoc2);

measDim=length(rBi);

S=zTx';
SStar=pinv(S);

z=0.5*(bsxfun(@minus,sum(S.*S,2),rBi.^2));
a=SStar*z;
b=SStar*rBi;

%The paper suggests just taking the real part to add robustness to noise in
%the measurements.
rootTerm=real(sqrt(4*(a'*b)^2-4*((b'*b)-1)*(a'*a)));
denom=2*(b'*b-1);
rMono1=(-2*a'*b-rootTerm)/denom;
rMono2=(-2*a'*b+rootTerm)/denom;

%Compute the error norm with both solutions and choose the one with the
%smallest error if more than the minimum number of receivers is given.
xEst1=a+b*rMono1;
diff=bsxfun(@minus,xEst1,zTx);
d1=norm(rBi-norm(xEst1)-sqrt(sum(diff.*diff,1)'));

xEst2=a+b*rMono2;
diff=bsxfun(@minus,xEst2,zTx);
d2=norm(rBi-norm(xEst2)-sqrt(sum(diff.*diff,1)'));

if(measDim==3)
    xEst(:,1)=xEst1;
    xEst(:,2)=xEst2;
    Rt1=rMono1;
    %xEst2 is returned as the second solution.
    Rt2=rMono2;
else
    if(d1<d2)
        xEst=xEst1;
        Rt1=rMono1;
    else
        xEst=xEst2;
        Rt1=rMono2;
    end
    xEst2=[];%There is a unique solution.
end

%If a measurement covariance matrix is desired
if(nargout>1)
    numSol=size(xEst,2);
    %Allocate space
    PTaylor=zeros(3,3,numSol);
    PCRLB=zeros(3,3,numSol);
    
    Delta=S-rBi*xEst1'/norm(xEst1);
    Gamma=diag(rBi);

    dxdr=lsqminnorm(Delta,(eye(measDim)*Rt1-Gamma));
    PTaylor(:,:,1)=dxdr*RCov*dxdr';
    PCRLB(:,:,1)=pinv(dxdr*pinv(RCov)*dxdr');
    
    %If a second solution is available, then save it.
    if(numSol>1)
        Delta=S-rBi*xEst2'/norm(xEst2);
        Gamma=diag(rBi);

        dxdr=lsqminnorm(Delta,(eye(measDim)*Rt2-Gamma));
        PTaylor(:,:,2)=dxdr*RCov*dxdr';
        PCRLB(:,:,2)=pinv(dxdr*pinv(RCov)*dxdr');
    end
end

%Adjust for the receiver not being at the origin.
xEst=bsxfun(@minus,xEst,zLoc2);

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
