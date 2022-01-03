function [zCart,RCart]=spher2CartTaylor(zMeas,R,zRx,M,algorithm)
%%SPHER2CARTTAYLOR Approximate the Cartesian moments of a Gaussian
%                 noise-corrupted measurement in monostatic spherical
%                 coordinates. This function approaches the conversion
%                 using traditional methods that use Taylor series
%                 expansions. The function spher2CartCubature can have
%                 better performance when the measurement noise is high and
%                 it can make use of cross terms in the covariance matrix;
%                 it also supports bistatic measurements and an alternative
%                 definition of the spherical coordinate system. In this
%                 function, azimuth is measured counterclockwise from the
%                 x-axis in the x-y plane. Elevation is measured up from
%                 the x-y plane (towards the z-axis).
%
%INPUTS:zMeas One or more points given in terms of range, azimuth and
%             elevation, with the angles in radians To convert N points,
%             zMeas is a 3XN matrix with each column having the format
%             [range;azimuth; elevation].
%           R The 3X3XN covariance matrices associated with polPoint. If
%             all of the matrices are the same, then this can just be a
%             single 3X3 matrix.
%         zRx The 3XN [x;y;z] location vectors of the receivers in
%             Cartesian coordinates. If this parameter is omitted, then the
%             receivers are assumed to be at the origin. If only a single
%             vector is passed, then the receiver location is assumed the
%             same for all of the target states being converted. zRx can
%             have more than 3 rows; additional rows are ignored. If
%             monostatic or no range values are provided, an empty matrix
%             can be passed.
%           M A 3X3XN hypermatrix of the rotation matrices to go from the
%             alignment of the global coordinate system to that at the
%             receiver. The z-axis of the local coordinate system of the
%             receiver is the pointing direction of the receiver. If
%             omitted, then it is assumed that the local coordinate system
%             is aligned with the global and M=eye(3) --the identity matrix
%             is used. If only a single 3X3 matrix is passed, then it is
%             assumed to be the same for all of the N conversions.
%   algorithm An optional parameter specifying the algorithm to use.
%             Possible values are
%             0 The multiplicative unbiased conversion from [1].
%             1 (The default if omitted or an empty matrix is passed) The
%               modified multiplicative unbiased conversion from [2].
%
%OUTPUTS: zCart The approximate means of the PDF of the Cartesian converted
%               measurements in [x;y;z] Cartesian coordinates for each
%               measurement. This is a 3XN matrix.
%         RCart The approximate 3X3XN set of covariance matrices of the
%               PDFs of the Cartesian converted measurements.
%
%REFERENCES:
%[1] M. Longbin, S. Xiaoquan, Z. Yiyu, S. Z. Kang, and Y. Bar-Shalom,
%    "Unbiased converted measurements for tracking," IEEE Transactions
%    on Aerospace and Electronic Systems, vol. 34, no. 3, pp. 1023-1027,
%    Jul. 1998.
%[2] Z. Duan, C. Han, and X. R. Li, "Comments on 'unbiased converted
%    measurements for tracking'," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 40, no. 4, pp. 1374-1377, Oct. 2004.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=1;
end

if(nargin<3||isempty(zRx))
    zRx=zeros(3,1);
end

if(nargin<4||isempty(M))
    M=eye(3,3);
end

numMeas=size(zMeas,2);

if(size(R,3)==1)
    R=repmat(R,[1,1,numMeas]);
end

if(size(M,3)==1)
    M=repmat(M,[1,1,numMeas]);
end

if(size(zRx,2)==1)
    zRx=repmat(zRx,[1,numMeas]);
end

zCart=zeros(3,numMeas);
RCart=zeros(3,3,numMeas);

switch(algorithm)
    case 0%The multiplicative unbiased conversion from [1].
        for curMeas=1:numMeas
            r=zMeas(1,curMeas);
            beta=zMeas(2,curMeas);
            eta=zMeas(3,curMeas);

            %This ignores any cross terms.
            sigmaR2=R(1,1,curMeas);
            sigmaBeta2=R(2,2,curMeas);
            sigmaEta2=R(3,3,curMeas);
            
            %Transpose equals inverse of a rotation matrix.
            MInvCur=M(:,:,curMeas)';
            
            cosBeta=cos(beta);
            cos2Beta=cos(2*beta);
            sinBeta=sin(beta);
            sin2Beta=sin(2*beta);
            cosEta=cos(eta);
            cos2Eta=cos(2*eta);
            sinEta=sin(eta);
            sin2Eta=sin(2*eta);

            %From Equation 12
            lambdaBeta=exp(-sigmaBeta2/2);
            lambdaBetaInv=1/lambdaBeta;
            lambdaBetaPrime=lambdaBeta^4;
            lambdaEta=exp(-sigmaEta2/2);
            lambdaEtaInv=1/lambdaEta;
            lambdaEtaPrime=lambdaEta^4;

            %Equation 9a
            x=lambdaBetaInv*lambdaEtaInv*r*cosBeta*cosEta;
            %Equation 9b
            y=lambdaBetaInv*lambdaEtaInv*r*sinBeta*cosEta;
            %Equation 9c
            z=lambdaEtaInv*r*sinEta;
            zCart(:,curMeas)=MInvCur*[x;y;z]+zRx(:,curMeas);

            %Equation 11a
            RCart(1,1,curMeas)=(1/4)*lambdaBetaInv^2*lambdaEtaInv^2*(r^2+2*sigmaR2)*(1+lambdaBetaPrime^2*cos2Beta)*(1+lambdaEtaPrime^2*cos2Eta)-(1/4)*(r^2+sigmaR2)*(1+lambdaBetaPrime*cos2Beta)*(1+lambdaEtaPrime*cos2Eta);
            %Equation 11b
            RCart(2,2,curMeas)=(1/4)*lambdaBetaInv^2*lambdaEtaInv^2*(r^2+2*sigmaR2)*(1-lambdaBetaPrime^2*cos2Beta)*(1+lambdaEtaPrime^2*cos2Eta)-(1/4)*(r^2+sigmaR2)*(1-lambdaBetaPrime*cos2Beta)*(1+lambdaEtaPrime*cos2Eta);
            %Equation 11c
            RCart(3,3,curMeas)=(1/2)*lambdaEtaInv^2*(r^2+2*sigmaR2)*(1-lambdaEtaPrime^2*cos2Eta)-(1/2)*(r^2+sigmaR2)*(1-lambdaEtaPrime*cos2Eta);
            %Equation 11d
            RCart(1,2,curMeas)=(1/4)*lambdaBetaInv^2*lambdaEtaInv^2*lambdaBetaPrime^2*(r^2+2*sigmaR2)*sin2Beta*(1+lambdaEtaPrime^2*cos2Eta)-(1/4)*lambdaBetaPrime*(r^2+sigmaR2)*sin2Beta*(1+lambdaEtaPrime*cos2Eta);
            RCart(2,1,curMeas)=RCart(1,2,curMeas);
            %Equation 11e
            RCart(1,3,curMeas)=(1/2)*lambdaBeta*lambdaEtaInv^2*lambdaEtaPrime^2*(r^2+2*sigmaR2)*cosBeta*sin2Eta-(1/2)*lambdaBeta*lambdaEtaPrime*(r^2+sigmaR2)*cosBeta*sin2Eta;
            RCart(3,1,curMeas)=RCart(1,3,curMeas);
            %Equation 11f
            RCart(2,3,curMeas)=(1/2)*lambdaBeta*lambdaEtaInv^2*lambdaEtaPrime^2*(r^2+2*sigmaR2)*sinBeta*sin2Eta-(1/2)*lambdaBeta*lambdaEtaPrime*(r^2+sigmaR2)*sinBeta*sin2Eta;
            RCart(3,2,curMeas)=RCart(2,3,curMeas);
            
            RCart(:,:,curMeas)=MInvCur*RCart(:,:,curMeas)*MInvCur';
        end
    case 1%The modified multiplicative unbiased conversion from [2].
        for curMeas=1:numMeas
            r=zMeas(1,curMeas);
            beta=zMeas(2,curMeas);
            eta=zMeas(3,curMeas);

            %This ignores any cross terms.
            sigmaR2=R(1,1,curMeas);
            sigmaBeta2=R(2,2,curMeas);
            sigmaEta2=R(3,3,curMeas);

            %Transpose equals inverse of a rotation matrix.
            MInvCur=M(:,:,curMeas)';

            cosBeta=cos(beta);
            cos2Beta=cos(2*beta);
            sinBeta=sin(beta);
            sin2Beta=sin(2*beta);
            cosEta=cos(eta);
            cos2Eta=cos(2*eta);
            sinEta=sin(eta);
            sin2Eta=sin(2*eta);
        
            %From Equation 12
            lambdaBeta=exp(-sigmaBeta2/2);
            lambdaBetaPrime=lambdaBeta^4;
            lambdaEta=exp(-sigmaEta2/2);
            lambdaEtaPrime=lambdaEta^4;

            x=lambdaBeta*lambdaEta*r*cosBeta*cosEta;
            y=lambdaBeta*lambdaEta*r*sinBeta*cosEta;
            z=lambdaEta*r*sinEta;
            zCart(:,curMeas)=MInvCur*[x;y;z]+zRx(:,curMeas);

            RCart(1,1,curMeas)=-lambdaBeta^2*lambdaEta^2*r^2*cosBeta^2*cosEta^2+(1/4)*(r^2+sigmaR2)*(1+lambdaBetaPrime*cos2Beta)*(1+lambdaEtaPrime*cos2Eta);
            RCart(2,2,curMeas)=-lambdaBeta^2*lambdaEta^2*r^2*sinBeta^2*cosEta^2+(1/4)*(r^2+sigmaR2)*(1-lambdaBetaPrime*cos2Beta)*(1+lambdaEtaPrime*cos2Eta);
            RCart(3,3,curMeas)=-lambdaEta^2*r^2*sinEta^2+(1/2)*(r^2+sigmaR2)*(1-lambdaEtaPrime*cos2Eta);
            RCart(1,2,curMeas)=-lambdaBeta^2*lambdaEta^2*r^2*sinBeta*cosBeta*cosEta^2+(1/4)*(r^2+sigmaR2)*lambdaBetaPrime*sin2Beta*(1+lambdaEtaPrime*cos2Eta);
            RCart(2,1,curMeas)=RCart(1,2,curMeas);
            RCart(1,3,curMeas)=-lambdaBeta*lambdaEta^2*r^2*cosBeta*sinEta*cosEta+(1/2)*(r^2+sigmaR2)*lambdaBeta*lambdaEtaPrime*cosBeta*sin2Eta;
            RCart(3,1,curMeas)=RCart(1,3,curMeas);
            RCart(2,3,curMeas)=-lambdaBeta*lambdaEta^2*r^2*sinBeta*sinEta*cosEta+(1/2)*(r^2+sigmaR2)*lambdaBeta*lambdaEtaPrime*sinBeta*sin2Eta;
            RCart(3,2,curMeas)=RCart(2,3,curMeas);

            RCart(:,:,curMeas)=MInvCur*RCart(:,:,curMeas)*MInvCur';
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
