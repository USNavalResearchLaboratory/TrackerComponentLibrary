function [zCart,RCart]=monostatRuv2CartTaylor(zMeas,R,zRx,M,algorithm)
%MONOSTATRUV2CARTTAYLOR Approximate the Cartesian moments of a Gaussian
%                 noise-corrupted measurement in monostatic r-u-v
%                 coordinates. This function approaches the conversion
%                 using traditional methods that use Taylor series
%                 expansions. The function ruv2CartCubature can have better
%                 performance when the measurement noise is high and it can
%                 make use of cross terms in the covariance matrix.
%
%INPUTS:      z A 3XnumMeas matrix of numMeas vectors to convert. Each
%               has elements [r;u;v], where r is the one-way monostatic
%               range from the target to the receiver, and u and v are
%               direction cosines.
%             R The 3X3XnumMeas measurement covariance matrices for the
%               measurements. If all of the matrices are the same, then
%               this can just be a single 3X3 matrix.
%           zRx The 3XN [x;y;z] location vector of the receivers in
%               Cartesian coordinates.  If this parameter is omitted, then
%               the receivers are assumed to be at the origin. If only a
%               single vector is passed, then the receiver location is
%               assumed the same for all of the target states being
%               converted.
%            M  A 3X3XN hypermatrix of the rotation matrices to go from the
%               alignment of the global coordinate system to that at the
%               receiver. The z-axis of the local coordinate system of the
%               receiver is the pointing direction of the receiver. If 
%               omitted, then it is assumed that the local coordinate 
%               system is aligned with the global and M=eye(3) --the
%               identity matrix is used. If only a single 3X3 matrix is
%               passed, then is is assumed to be the same for all of the N
%               conversions.
%     algorithm An optional parameter specifying the algorithm to use.
%               Possible values are
%               0 The CM1 conversion from [1] (uses a first-order Taylor
%                 series approximation).
%               1 (The default if omitted or an empty matrix is passed) The
%                 CM2 conversion from [1] (uses a second-order Taylor
%                 series approximation).
%
%OUTPUTS:   zCart   The approximate means of the PDF of the Cartesian
%                   converted measurements in [x;y;z] Cartesian coordinates
%                   for each measurement. This is a 3XnumMeas matrix.
%           RCart   The approximate 3X3XnumMeas set of covariance
%                   matrices of the PDFs of the Cartesian converted
%                   measurements.
%
%REFERENCES:
%[1] X. Tian and Y. Bar-Shalom, "Coordinate conversion and tracking
%    for very long range radars," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 45, no. 3, pp. 1073-1088, Jul. 2009.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(zMeas,2);

if(nargin<5||isempty(algorithm))
    algorithm=1;
end

if(nargin<4||isempty(M))
    M=eye(3);
end

if(nargin<3||isempty(zRx))
    zRx=zeros(3,1);
end

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
    case 0%The CM1 algorithm from [1]. This uses a first-order Taylor
          %series expansion.
        for curMeas=1:numMeas
            r=zMeas(1,curMeas);
            u=zMeas(2,curMeas);
            v=zMeas(3,curMeas);

            %Cross terms are neglected
            sigmaR2=R(1,1,curMeas);
            sigmaU2=R(2,2,curMeas);
            sigmaV2=R(3,3,curMeas);
            
            %Transpose equals inverse of a rotation matrix.
            MInvCur=M(:,:,curMeas)';

            %First derivatives with respect to z.
            temp1=abs(sqrt(1-u^2-v^2));%Abs added to deal with bad inputs.
            dfzdr=temp1;
            dfzdu=-r*u/temp1;
            dfzdv=-r*v/temp1;

            %Equation 15
            x=r*u;
            %Equation 16
            y=r*v;
            %Equation 17
            z=r*temp1;

            %The bias removal algorithm breaks down at extreme angles. This
            %should not occur in practice. It is assumed that the target is
            %is in front of the radar, so the z component should be
            %positive.
            z=max(z,0);
            zCart(:,curMeas)=MInvCur*[x;y;z]+zRx(:,curMeas);

            %Equation 18
            R11=u^2*sigmaR2+r^2*sigmaU2;
            %Equation 19
            R22=v^2*sigmaR2+r^2*sigmaV2;
            %Equation 20
            R33=dfzdr^2*sigmaR2+dfzdu^2*sigmaU2+dfzdv^2*sigmaV2;
            %Equation 21
            R12=u*v*sigmaR2;
            %Equation 22
            R13=u*dfzdr*sigmaR2+r*dfzdu*sigmaU2;
            %Equation 23
            R23=v*dfzdr*sigmaR2+r*dfzdv*sigmaV2;

            RCart(:,:,curMeas)=MInvCur*[R11, R12, R13;
                                        R12, R22, R23;
                                        R13, R23, R33]*MInvCur';
        end 
    case 1%The CM2 algorithm from [1]. This uses a second-order Taylor
          %series expansion.
        for curMeas=1:numMeas
            r=zMeas(1,curMeas);
            u=zMeas(2,curMeas);
            v=zMeas(3,curMeas);

            %Cross terms are neglected
            sigmaR2=R(1,1,curMeas);
            sigmaU2=R(2,2,curMeas);
            sigmaV2=R(3,3,curMeas);
            
            %Transpose equals inverse of a rotation matrix.
            MInvCur=M(:,:,curMeas)';
            
            %First derivatives with respect to z.
            temp1=abs(sqrt(1-u^2-v^2));%Abs  added to deal with bad inputs.
            dfzdr=temp1;
            dfzdu=-r*u/temp1;
            dfzdv=-r*v/temp1;

            %Second derivatives with respect to z.
            temp2=(1-u^2-v^2)^(3/2);
            d2fzdu=r*(v^2-1)/temp2;
            d2fzdv=r*(u^2-1)/temp2;
            d2fzdrdu=-u/temp1;
            d2fzdrdv=-v/temp1;
            d2fzdudv=-r*u*v/temp2;

            cz=(1/2)*d2fzdu*sigmaU2+(1/2)*d2fzdv*sigmaV2;
            
            %Calculate the debiased converted measurements.
            %Equation 28, taking the mean (eliminates w terms)
            x=r*u;
            %Equation 29, taking the mean (eliminates w terms)
            y=r*v;
            %Equation 27
            z=r*temp1-cz;

            %The bias removal algorithm breaks down at extreme angles. This
            %should not occur in practice. It is assumed that the target is in
            %front of the radar, so the z component should be positive.
            z=max(z,0);
            zCart(:,curMeas)=MInvCur*[x;y;z]+zRx(:,curMeas);

            %Now for the derivatives.
            %Equation 34
            R11=r^2*sigmaU2+u^2*sigmaR2+sigmaR2*sigmaU2;
            %Equation 35
            R22=r^2*sigmaV2+v^2*sigmaR2+sigmaR2*sigmaV2;
            %Equation 33
            R33=cz^2+dfzdr^2*sigmaR2+dfzdu^2*sigmaU2+dfzdv^2*sigmaV2+cz*(d2fzdu*sigmaU2+d2fzdv*sigmaV2)+(3/4)*d2fzdu*sigmaU2^2+(3/4)*d2fzdv*sigmaV2^2+d2fzdrdu^2*sigmaR2*sigmaU2+d2fzdrdv^2*sigmaR2*sigmaV2+d2fzdudv^2*sigmaU2*sigmaV2;
            %Equation 36
            R12=u*v*sigmaR2;
            %Equation 37
            R13=dfzdu*r*sigmaU2+dfzdr*u*sigmaR2;
            %Equation 38
            R23=dfzdv*r*sigmaV2+dfzdr*v*sigmaR2;

            RCart(:,:,curMeas)=MInvCur*[R11, R12, R13;
                                        R12, R22, R23;
                                        R13, R23, R33]*MInvCur';
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
