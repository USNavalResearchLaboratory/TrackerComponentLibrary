function rr=getRangeRate(xTar,useHalfRange,xTx,xRx)
%GETRANGERATE Obtain the bistatic range rates of targets in the absence
%             of refraction under non-relativistic mechanics, ignoring
%             atmospheric effects, when the transmitter, target and
%             receiver are all moving. If the states of the target and
%             transmitter are the same, then it is assumed that the
%             target is the emitter.
%
%INPUTS: xTar The Cartesian states of the targets in 1D, 2D, or 3D space.
%             If 3D, then xTar is a 6XN matrix, (or where the first 3
%             components are position and the second 3 are velocity. If
%             xTar is a 4XN matrix, where the first two components are
%             position and the second two velocity. Components are ordered
%             position (x,y,z) and velocity (xDot,yDot,zDot). The same
%             pattern applies in 1D.
% useHalfRange A boolean value specifying whether the bistatic range value
%             should be divided by two, which means that the range rate is
%             divided by two. This normally comes up when operating in
%             monostatic mode, so that the range reported is a one-way
%             range. The default if an empty matrix is passed is false.
%         xTx An xTxDimXN matrix of the states of the transmitters
%             consisting of stacked 1D, 2D, or 3D position and velocity
%             components. If this parameter is omitted, the transmitters
%             are assumed to be stationary at the origin. If only a single
%             vector is passed, then the transmitter state is assumed the
%             same for all of the target states being converted.
%         xRx An xRxDimXN matrix of the states of the receivers consisting
%             of stacked 1D, 2D, or 3D position and velocity components. If
%             this parameter is omitted, the receivers are assumed to be
%             stationary at the origin. If only a single vector is passed,
%             then the receiver state is assumed the same for all of the
%             target states being converted.
%
%OUTPUTS: rr The 1XN bistatic range rates of the targets. If
%            useHalfRange=true, then the range rate is halved to reflect a
%            halved range.
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

    numDim=size(xTar,1)/2;

    N=size(xTar,2);
    
    if(nargin<4||isempty(xRx))
        xRx=zeros(2*numDim,N);
    elseif(size(xRx,2)==1&&N~=1)
        xRx=repmat(xRx,[1,N]);
    end

    if(nargin<3||isempty(xTx))
        xTx=zeros(2*numDim,N);
    elseif(size(xTx,2)==1&&N~=1)
        xTx=repmat(xTx,[1,N]);
    end
    
    if(nargin<2||isempty(useHalfRange))
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

    targetIsTrans=all(xTar==xTx);%If the target is the transmitter.

    %One could just do the following three lines
    % dtrRat=dtr./sqrt(sum(dtr.*dtr,1));
    % dtlRat=dtl./sqrt(sum(dtl.*dtl,1));
    % rr=sum((dtrRat+dtlRat).*vt,1)-sum(dtrRat.*vr,1)-sum(dtlRat.*vi,1);
    %However, this appears to be slower than explicitly writing out the
    %multiplications, which we do for the 2D and 3D cases.
    if(numDim==2)
        if(targetIsTrans)
            dtlRat=zeros(2,N);
        else 
            dtlRat=bsxfun(@rdivide,dtl,sqrt(dtl(1,:).*dtl(1,:)+dtl(2,:).*dtl(2,:)));
        end

        dtrRat=bsxfun(@rdivide,dtr,sqrt(dtr(1,:).*dtr(1,:)+dtr(2,:).*dtr(2,:)));
        rr=(dtrRat(1,:)+dtlRat(1,:)).*vt(1,:)+(dtrRat(2,:)+dtlRat(2,:)).*vt(2,:)...
            -dtrRat(1,:).*vr(1,:)-dtrRat(2,:).*vr(2,:)...
            -dtlRat(1,:).*vi(1,:)-dtlRat(2,:).*vi(2,:);
    elseif(numDim==3)
        if(targetIsTrans)
            dtlRat=zeros(3,N);
        else
            dtlRat=bsxfun(@rdivide,dtl,sqrt(dtl(1,:).*dtl(1,:)+dtl(2,:).*dtl(2,:)+dtl(3,:).*dtl(3,:)));
        end
        
        dtrRat=bsxfun(@rdivide,dtr,sqrt(dtr(1,:).*dtr(1,:)+dtr(2,:).*dtr(2,:)+dtr(3,:).*dtr(3,:)));
        rr=(dtrRat(1,:)+dtlRat(1,:)).*vt(1,:)+(dtrRat(2,:)+dtlRat(2,:)).*vt(2,:)+(dtrRat(3,:)+dtlRat(3,:)).*vt(3,:)...
            -dtrRat(1,:).*vr(1,:)-dtrRat(2,:).*vr(2,:)-dtrRat(3,:).*vr(3,:)...
            -dtlRat(1,:).*vi(1,:)-dtlRat(2,:).*vi(2,:)-dtlRat(3,:).*vi(3,:);
    else%1D
        if(targetIsTrans)
            dtlRat=zeros(numDim,N);
        else
            dtlRat=dtl./sqrt(sum(dtl.*dtl,1));
        end
        dtrRat=dtr./sqrt(sum(dtr.*dtr,1));
        rr=sum((dtrRat+dtlRat).*vt,1)-sum(dtrRat.*vr,1)-sum(dtlRat.*vi,1);
    end

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
