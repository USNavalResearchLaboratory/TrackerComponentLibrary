function rr=getRangeRate(xTar,useHalfRange,xTx,xRx,numDim)
%GETRANGERATE Obtain the bistatic range rate of targets in the absence
%             of refraction under non-relativistic mechanics, ignoring
%             atmospheric effects, when the transmitter, target and
%             receiver are all moving.
%
%INPUTS:    xTar    The Cartesian state of the targets in 2D or 3D
%                   Cartesian space. If 3D, then xTar is a 6XN matrix, (or 
%                   where the first 3 components are position and the
%                   second 3 are velocity. Otherwise, xTar is a 4XN matrix,
%                   where the first two components are position and the
%                   second two velocity. In 2D and 3D, it is allowed for
%                   the state to have extra rows. In such an instance, the 
%                   extra rows are ignored. This allowed states with
%                   extra components, such as acceleration, to be passed.
%      useHalfRange A boolean value specifying whether the bistatic range
%                   value should be divided by two, which means that the
%                   range rate is divided by two. This normally comes up
%                   when operating in monostatic mode, so that the range
%                   reported is a one-way range. The default if this
%                   parameter is not provided is false.
%           xTx     An xTxDimXN matrix of the states of the transmitters
%                   consisting of stacked 3D position and velocity
%                   components. Other components will be ignored. If this
%                   parameter is omitted, the transmitters are assumed to
%                   be stationary at the origin. If only a single vector is
%                   passed, then the transmitter state is assumed the same
%                   for all of the target states being converted.
%           xRx     An xRxDimXN matrix of the states of the receivers
%                   consisting of stacked 3D position and velocity
%                   components. Other components will be ignored. If this
%                   parameter is omitted, the receivers are assumed to be
%                   stationary at the origin. If only a single vector is
%                   passed, then the receiver state is assumed the same for
%                   all of the target states being converted.
%          numDim   An optional value specifying the number of Cartesian
%                   dimensions of the states. This can be 2 or 3. If this
%                   parameter is omitted, then the default value of 3 is
%                   used.
%
%OUTPUTS:   rr      The 1 X N bistatic range rates of the targets. If
%                   useHalfRange=true, then the range rate is halved to
%                   reflect a halved range.
%
%This assumes that the target state that is provided has the same
%dimensionality as the states of the transmitter and receiver and that all
%of the states are in Cartesian coordinates.
%
%A derivation of this non-relativistic approximation is given in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%July 2012 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5)
        numDim=3;
    end

    N=size(xTar,2);
    
    if(nargin<4)
        xRx=zeros(2*numDim,N);
    elseif(size(xRx,2)==1)
        xRx=repmat(xRx,[1,N]);
    end

    if(nargin<3)
        xTx=zeros(2*numDim,N);
    elseif(size(xTx,2)==1)
        xTx=repmat(xTx,[1,N]);
    end
    
    if(nargin<2)
        useHalfRange=false;
    end
    
    selPos=1:numDim;
    selVel=(numDim+1):(2*numDim);
    
    zt=xTar(selPos,:);
    zr=xRx(selPos,:);
    li=xTx(selPos,:);

    vt=xTar(selVel,:);
    vr=xRx(selVel,:);
    vi=xTx(selVel,:);

    dtr=zt-zr;
    dtl=zt-li;

    dtrRat=bsxfun(@rdivide,dtr,sqrt(sum(dtr.*dtr,1)));
    dtlRat=bsxfun(@rdivide,dtl,sqrt(sum(dtl.*dtl,1)));
    
    rr=sum((dtrRat+dtlRat).*vt,1)-sum(dtrRat.*vr,1)-sum(dtlRat.*vi,1);
    
    if(useHalfRange)
        rr=rr/2;
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
