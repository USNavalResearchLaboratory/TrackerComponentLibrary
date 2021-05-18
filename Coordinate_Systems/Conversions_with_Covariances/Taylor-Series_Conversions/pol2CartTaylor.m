function [zCart,RCart]=pol2CartTaylor(polPoint,R,zRx,M,algorithm)
%%POL2CARTTAYLOR Approximate the first and second Cartesian moments of an
%                additive Gaussian noise-corrupted measurement in polar
%                coordinates. The polar angle (azimuth) is measured from
%                the x-axis counterclockwise. This function approaches the
%                conversion using traditional methods that use Taylor
%                series expansions. The function pol2CartCubature can have
%                better performance when the measurement noise is high and
%                it can make use of cross terms in the covariance matrix.
%
%INPUTS: polPoint A 2XN set of N polar points of the form [range;azimuth]
%               that are to be converted into Cartesian coordinates. The
%               range is one-way monostatic; the azimuth is in radians.
%             R The 2X2XN covariance matrices associated with polPoint. If
%               all of the matrices are the same, then this can just be a
%               single 2X2 matrix. Note that all of the algorithms ignore
%               cross terms in the covariance matrices.
%           zRx The optional 2X1 [x;y] location vector of the receiver in
%               global Cartesian coordinates. If omitted or an empty matrix
%               is passed, the receiver is placed at the origin.
%             M A 2X2 rotation matrix to go from the alignment of the
%               global coordinate system to that at the receiver. If
%               omitted or an empty matrix is passed, then it is assumed
%               that the local coordinate system is aligned with the global
%               and M=eye(2) --the identity matrix is used.
%     algorithm An optional parameter specifying the algorithm to use.
%               Possible values are
%               0 The "standard conversion" from Ch. 10.4.3 of [1].
%               1 The additively debiased conversion from [2].
%               2 (The default if omitted or an empty matrix is passed) The
%                 multiplicative unbiased conversion from [3]
%               3 The modified multiplicative unbiased conversion from [4].
%
%OUTPUTS: zCart The 2XN approximate means of the PDF of the polar
%               measurements converted to [x;y] Cartesian coordinates.
%         RCart The 2X2XN set of approximate covariance matrices for the N
%               estimates.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%[2] D. Lerro and Y. Bar-Shalom, "Tracking with debiased consistent
%    converted measurements versus EKF," IEEE Transactions on Aerospace
%    and Electronic Systems, vol. 29, no. 3, pp. 1015-1022, Jul. 1993.
%[3] M. Longbin, S. Xiaoquan, Z. Yiyu, S. Z. Kang, and Y. Bar-Shalom,
%    "Unbiased converted measurements for tracking," IEEE Transactions
%    on Aerospace and Electronic Systems, vol. 34, no. 3, pp. 1023-1027,
%    Jul. 1998.
%[4] Z. Duan, C. Han, and X. R. Li, "Comments on 'unbiased converted
%    measurements for tracking'," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 40, no. 4, pp. 1374-1377, Oct. 2004.
%[5] S. V. Bordonaro, P. Willett, and Y. Bar-Shalom, "Tracking with
%    converted position and Doppler measurements," in Proceedings of SPIE:
%    Signal and Data Processing of Small Targets, vol. 8137, San Diego, CA,
%    21 Aug. 2011.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPoints=size(polPoint,2);

if(size(R,3)==1)
    R=repmat(R,[1,1,numPoints]);
end

if(nargin<3||isempty(zRx))
    zRx=zeros(2,1);
end

if(nargin<4||isempty(M))
	M=eye(2,2); 
end

if(nargin<5||isempty(algorithm))
    algorithm=2;
end

if(size(zRx,2)==1)
    zRx=repmat(zRx,[1,numPoints]);
end

%For the inverse rotation (local to global).
M=M';

zCart=zeros(2,numPoints);
RCart=zeros(2,2,numPoints);
switch(algorithm)
    case 0%The "standard conversion" from Ch. 10.4.3 of [1].
        for curPoint=1:numPoints
            r=polPoint(1,curPoint);
            theta=polPoint(2,curPoint);
            sigmaR2=R(1,1,curPoint);
            sigmaTheta2=R(2,2,curPoint);
        
            sT=sin(theta);
            cT=cos(theta);

            %Equation 10.4.3-3
            x=r*cT;
            y=r*sT;
            zCart(:,curPoint)=M*[x;y]+zRx(:,curPoint);
            
            %Equation 10.4.3-4
            RCart(1,1,curPoint)=r^2*sigmaTheta2*sT^2+sigmaR2*cT^2;
            %Equation 10.4.3-6
            RCart(1,2,curPoint)=(sigmaR2-r^2*sigmaTheta2)*sT*cT;
            RCart(2,1,curPoint)=RCart(1,2,curPoint);
            %Equation 10.4.3-5
            RCart(2,2,curPoint)=r^2*sigmaTheta2*cT^2+sigmaR2*sT^2;
            
            RCart(:,:,curPoint)=M*RCart(:,:,curPoint)*M';
        end
    case 1%The additively debiased conversion from [2].
        for curPoint=1:numPoints
            r=polPoint(1,curPoint);
            theta=polPoint(2,curPoint);
            sigmaR2=R(1,1,curPoint);
            sigmaTheta2=R(2,2,curPoint);
            
            sT=sin(theta);
            cT=cos(theta);

            expTerm=exp(-sigmaTheta2);
            expTermd2=exp(-sigmaTheta2/2);
            expTerm2=exp(-2*sigmaTheta2);

            %Equation 12 and 14
            x=r*cT-r*cT*(expTerm-expTermd2);
            y=r*sT-r*sT*(expTerm-expTermd2);
            zCart(:,curPoint)=M*[x;y]+zRx(:,curPoint);

            deltaCosh=cosh(2*sigmaTheta2)-cosh(sigmaTheta2);
            deltaCosh2=2*cosh(2*sigmaTheta2)-cosh(sigmaTheta2);
            deltaSinh=sinh(2*sigmaTheta2)-sinh(sigmaTheta2);
            deltaSinh2=2*sinh(2*sigmaTheta2)-sinh(sigmaTheta2);

            %Equation 13a
            RCart(1,1,curPoint)=r^2*expTerm2*(cT^2*deltaCosh+sT^2*deltaSinh)+sigmaR2*expTerm2*(cT^2*deltaCosh2+sT^2*deltaSinh2);
            %Equation 13c
            RCart(1,2,curPoint)=sT*cT*expTerm2^2*(sigmaR2+(r^2+sigmaR2)*(1-1/expTerm));
            RCart(2,1,curPoint)=RCart(1,2,curPoint);
            %Equation 13b
            RCart(2,2,curPoint)=r^2*expTerm2*(sT^2*deltaCosh+cT^2*deltaSinh)+sigmaR2*expTerm2*(sT^2*deltaCosh2+cT^2*deltaSinh2);
            
            RCart(:,:,curPoint)=M*RCart(:,:,curPoint)*M';
        end
    case 2%The multiplicative unbiased conversion from [3]; it is also in
          %[1]. 
        for curPoint=1:numPoints
            r=polPoint(1,curPoint);
            theta=polPoint(2,curPoint);
            sigmaR2=R(1,1,curPoint);
            sigmaTheta2=R(2,2,curPoint);

            sT=sin(theta);
            cT=cos(theta);
            c2T=cos(2*theta);
            s2T=sin(2*theta);

            %Equation 12 in [2].
            b1Inv=exp(sigmaTheta2/2);
            %Equation 12 in [2].
            b2=b1Inv^(-4);

            %Equation 5 in [2].
            x=b1Inv*r*cT;
            y=b1Inv*r*sT;
            zCart(:,curPoint)=M*[x;y]+zRx(:,curPoint);

            %Equation 7a in [2].
            RCart(1,1,curPoint)=(b1Inv^2-2)*r^2*cT^2+(1/2)*(r^2+sigmaR2)*(1+b2*c2T);
            %Equation 7c in [2].
            RCart(1,2,curPoint)=b1Inv^2*r^2*cT*sT-2*r^2*cT*sT+(1/2)*(r^2+sigmaR2)*b2*s2T;
            RCart(2,1,curPoint)=RCart(1,2,curPoint);
            %Equation 7b in [2].
            RCart(2,2,curPoint)=(b1Inv^2-2)*r^2*sT^2+(1/2)*(r^2+sigmaR2)*(1-b2*c2T);
            
            RCart(:,:,curPoint)=M*RCart(:,:,curPoint)*M';
        end
    case 3%The modified multiplicative unbiased conversion from [4]. We use
        %the equations given in [5], which are a bit clearer as they
        %actually give the converted values rather than requiring the user
        %to subtract. Equations numbers are from [5].
        
        for curPoint=1:numPoints
            r=polPoint(1,curPoint);
            theta=polPoint(2,curPoint);
            sigmaR2=R(1,1,curPoint);
            sigmaTheta2=R(2,2,curPoint);
            
            sT=sin(theta);
            s2T=sin(2*theta);
            cT=cos(theta);
            c2T=cos(2*theta);

            expTerm=exp(-sigmaTheta2);
            expTermd2=exp(-sigmaTheta2/2);
            expTerm2=exp(-2*sigmaTheta2);

            %Equation 14
            x=expTermd2*r*cT;
            %Equation 15
            y=expTermd2*r*sT;
            zCart(:,curPoint)=M*[x;y]+zRx(:,curPoint);

            %Equation 16
            RCart(1,1,curPoint)=(1/2)*(r^2+sigmaR2)*(1+c2T*expTerm2)-expTerm*r^2*cT^2;
            %Equation 18
            RCart(1,2,curPoint)=(1/2)*(r^2+sigmaR2)*s2T*expTerm2-expTerm*r^2*cT*sT;
            RCart(2,1,curPoint)=RCart(1,2,curPoint);
            %Equation 17
            RCart(2,2,curPoint)=(1/2)*(r^2+sigmaR2)*(1-c2T*expTerm2)-expTerm*r^2*sT^2;
            
            RCart(:,:,curPoint)=M*RCart(:,:,curPoint)*M';
        end
    otherwise
        error('Unknown algorithm specified')     
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
