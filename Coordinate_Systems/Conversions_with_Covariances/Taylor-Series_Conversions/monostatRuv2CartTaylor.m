function [zCart,RCart]=monostatRuv2CartTaylor(zMeas,R,useHalfRange,zRx,M,algorithm)
%MONOSTATRUV2CARTTAYLOR Approximate the Cartesian mean and covariance
%                 matrix of a Gaussian noise-corrupted measurement in
%                 monostatic r-u-v coordinates. This function approximates
%                 the conversion using traditional methods that use Taylor
%                 series expansions. The function ruv2CartCubature can have
%                 better performance when the measurement noise is high and
%                 it can make use of cross terms in the covariance matrix.
%
%INPUTS: z A 3XnumMeas matrix of numMeas vectors to convert. Each has
%          elements [r;u;v], where r is the one-way monostatic range from
%          the target to the receiver, and u and v are direction cosines.
%        R The 3X3XnumMeas measurement covariance matrices for the
%          measurements. If all of the matrices are the same, then this can
%          just be a single 3X3 matrix.
% useHalfRange A boolean value specifying whether the bistatic range value
%          has been divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided is false.
%      zRx The 3XN [x;y;z] location vector of the receivers in Cartesian
%          coordinates. If this parameter is omitted, then the receivers
%          are assumed to be at the origin. If only a single vector is
%          passed, then the receiver location is assumed the same for all
%          of the target states being converted.
%        M A 3X3XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. The z-axis of the local coordinate system of the
%          receiver is the pointing direction of the receiver. If omitted,
%          then it is assumed that the local coordinate system is aligned
%          with the global and M=eye(3) --the identity matrix is used. If
%          only a single 3X3 matrix is passed, then is=t is assumed to be
%          the same for all of the N conversions.
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are
%          0 The CM1 conversion from [1] (uses a first-order Taylor series
%            approximation).
%          1 (The default if omitted or an empty matrix is passed) The CM2
%            conversion from [1] (uses a second-order Taylor series
%            approximation) as given in [2], which has corrections.
%          2 Use the CM3 algorithm of [2] (uses a third-order Taylor series
%            approximation).
%
%OUTPUTS: zCart The approximate means of the PDF of the Cartesian converted
%               measurements in [x;y;z] Cartesian coordinates for each
%               measurement. This is a 3XnumMeas matrix.
%         RCart The approximate 3X3XnumMeas set of covariance matrices of
%               the PDFs of the Cartesian converted measurements.
%
%In algorithms 1 and 2, the sign of cz has been corrected compared to [2].
%
%EXAMPLE:
%Here, we compare the NEES of the CM1, CM2, and CM3 algorithms with
%cubature integration. One will see that the CM1 algorithm performs poorly
%whereas the CM2 algorithm is a bit better and the CM3 and cubature
%integration algorithms have identical NEES performance.
% numMCRuns=5000;
% numUVals=10;
% 
% r=100e3;%True one-way range.
% useHalfRange=true;%Indicated one-way range.
% v=0;%True V value.
% sigmaR=1;
% sigmaU=1e-3;
% sigmaV=5e-3;
% SR=diag([sigmaR;sigmaU;sigmaV]);
% R=SR*SR';%Measurement covariance matrix.
% uVals=linspace(-0.9,0.9,numUVals);
% 
% zCart1=zeros(3,numUVals,numMCRuns);
% RCart1=zeros(3,3,numUVals,numMCRuns);
% zCart2=zeros(3,numUVals,numMCRuns);
% RCart2=zeros(3,3,numUVals,numMCRuns);
% zCart3=zeros(3,numUVals,numMCRuns);
% RCart3=zeros(3,3,numUVals,numMCRuns);
% zCart4=zeros(3,numUVals,numMCRuns);
% RCart4=zeros(3,3,numUVals,numMCRuns);
% 
% zCartTrue=zeros(3,numUVals);
% for curU=1:numUVals
%     u=uVals(curU);
%     zTrue=[r;u;v];
%     zCartTrue(:,curU)=ruv2Cart(zTrue,useHalfRange);
%     for curRun=1:numMCRuns
%         z=zTrue+SR*randn(3,1);
%         
%         [zCart1(:,curU,curRun),RCart1(:,:,curU,curRun)]=monostatRuv2CartTaylor(z,R,useHalfRange,[],[],0);
%         [zCart2(:,curU,curRun),RCart2(:,:,curU,curRun)]=monostatRuv2CartTaylor(z,R,useHalfRange,[],[],1);
%         [zCart3(:,curU,curRun),RCart3(:,:,curU,curRun)]=monostatRuv2CartTaylor(z,R,useHalfRange,[],[],2);
%         [zCart4(:,curU,curRun),RCart4(:,:,curU,curRun)]=ruv2CartCubature(z,SR,useHalfRange);
%     end
% end
% NEES1=calcNEES(zCartTrue,zCart1,RCart1,true);
% NEES2=calcNEES(zCartTrue,zCart2,RCart2,true);
% NEES3=calcNEES(zCartTrue,zCart3,RCart3,true);
% NEES4=calcNEES(zCartTrue,zCart4,RCart4,true);
% 
% figure(1)
% clf
% hold on
% plot(uVals,NEES1,'-b','linewidth',2);
% plot(uVals,NEES2,'-g','linewidth',4);
% plot(uVals,NEES3,'--r','linewidth',4);
% plot(uVals,NEES4,'-.c','linewidth',2);
% h1=xlabel('u');
% h2=ylabel('NEES');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% %Draw the 99% confidence interval bounds.
% bounds=getNEESConfBounds(0.99,3,numMCRuns);
% plot([uVals(1),uVals(end)],[bounds(1),bounds(1)],'-k','linewidth',2)
% plot([uVals(1),uVals(end)],[bounds(2),bounds(2)],'-k','linewidth',2)
% legend('CM1','CM2','CM3','Cubature','99% Bounds')
%
%REFERENCES:
%[1] X. Tian and Y. Bar-Shalom, "Coordinate conversion and tracking
%    for very long range radars," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 45, no. 3, pp. 1073-1088, Jul. 2009.
%[2] J. R. Cookson, Z. T. Chance, and L. F. Urbano, "Consistent Estimation
%    for Very Long Range Radars," Lincoln Laboratory, Lexington, MA., Tech
%    Rep. 1184, 1 December 2017.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(zMeas,2);

if(nargin<6||isempty(algorithm))
    algorithm=1;
end

if(nargin<5||isempty(M))
    M=eye(3);
end

if(nargin<4||isempty(zRx))
    zRx=zeros(3,1);
end

if(nargin<3||isempty(useHalfRange))
    useHalfRange=false;
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

if(useHalfRange==false)
    zMeas(1)=zMeas(1)/2;
    D=diag([1/2;1;1]);
    for curMeas=1:numMeas
        R(:,:,curMeas)=D*R(:,:,curMeas)*D';
    end
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
            w=sqrt(abs(1-u^2-v^2));%Abs added to deal with bad inputs.
            dfzdr=w;
            dfzdu=-r*u/w;
            dfzdv=-r*v/w;

            %Equation 15
            x=r*u;
            %Equation 16
            y=r*v;
            %Equation 17
            z=r*w;
            
            r2=r*r;
            u2=u*u;
            v2=v*v;

            %The bias removal algorithm breaks down at extreme angles. This
            %should not occur in practice. It is assumed that the target is
            %is in front of the radar, so the z component should be
            %positive.
            z=max(z,0);
            zCart(:,curMeas)=MInvCur*[x;y;z]+zRx(:,curMeas);

            %Equation 18
            R11=u2*sigmaR2+r2*sigmaU2;
            %Equation 19
            R22=v2*sigmaR2+r2*sigmaV2;
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
    case 1%The CM2 algorithm from [1], as corrected in [2]. This uses a
          %second-order Taylor series expansion.
        for curMeas=1:numMeas
            r=zMeas(1,curMeas);
            u=zMeas(2,curMeas);
            v=zMeas(3,curMeas);

            %Cross terms are neglected
            sigmaR2=R(1,1,curMeas);
            sigmaU2=R(2,2,curMeas);
            sigmaU4=sigmaU2*sigmaU2;
            sigmaV2=R(3,3,curMeas);
            sigmaV4=sigmaV2*sigmaV2;
            
            r2=r*r;
            u2=u*u;
            u3=u2*u;
            u4=u3*u;
            v2=v*v;
            v3=v2*v;
            v4=v3*v;
            
            %Transpose equals inverse of a rotation matrix.
            MInvCur=M(:,:,curMeas)';
        
            %Equation 22 in [2].
            w=sqrt(abs(1-u^2-v^2));%Abs added to deal with bad inputs.
            w2=w*w;
            w3=w2*w;
            w4=w3*w;
            w5=w4*w;
            w6=w5*w;

            %Equation 23 in [2].
            cz=-(1/2)*(r/w+r*u2/w3)*sigmaU2-(1/2)*(r/w+r*v2/w3)*sigmaV2;
            
            %Equation 21 in [2].
            x=r*u;
            y=r*v;
            z=r*w+cz;%Sign corrected versus [2].
            
            %The bias removal algorithm breaks down at extreme angles. This
            %should not occur in practice. It is assumed that the target is
            %in front of the radar, so the z component should be positive.
            z=max(z,0);
            zCart(:,curMeas)=MInvCur*[x;y;z]+zRx(:,curMeas);
            
            %Equation 24 in [2].
            sigmaXX2=u2*sigmaR2+r2*sigmaU2+sigmaR2*sigmaU2;
            %Equation 25 in [2].
            sigmaYY2=v2*sigmaR2+r2*sigmaV2+sigmaR2*sigmaV2;
            %Equation 26 in [2].
            sigmaZZ2=w2*sigmaR2+(r2/w2)*(u2*sigmaU2+v2*sigmaV2)+(u2/w2)*sigmaR2*sigmaU2+(v2/w2)*sigmaR2*sigmaV2...
                +(r2/2)*(u2*v2/w6-(u2+v2)/w4-1/w2)*sigmaU2*sigmaV2...
                +(r2/2)*(u4/w6+2*u2/w4+1/w2)*sigmaU4+(r2/2)*(v4/w6+2*v2/w4+1/w2)*sigmaV4;
            %Equation 27 in [2].
            sigmaXY2=u*v*sigmaR2;
            %Equation 28 in [2].
            sigmaXZ2=u*w*sigmaR2-(r2*u/w)*sigmaU2-(u/w)*sigmaR2*sigmaU2;
            %Equation 29 in [2].
            sigmaYZ2=v*w*sigmaR2-(r2*v/w)*sigmaV2-(v/w)*sigmaR2*sigmaV2;
            
            RCart(:,:,curMeas)=MInvCur*[sigmaXX2, sigmaXY2, sigmaXZ2;
                                        sigmaXY2, sigmaYY2, sigmaYZ2;
                                        sigmaXZ2, sigmaYZ2, sigmaZZ2]*MInvCur';   
        end
    case 2%The CM3 algorithm of [2].
        for curMeas=1:numMeas
            r=zMeas(1,curMeas);
            u=zMeas(2,curMeas);
            v=zMeas(3,curMeas);

            %Cross terms are neglected
            sigmaR2=R(1,1,curMeas);
            sigmaU2=R(2,2,curMeas);
            sigmaU4=sigmaU2*sigmaU2;
            sigmaU6=sigmaU4*sigmaU2;
            sigmaV2=R(3,3,curMeas);
            sigmaV4=sigmaV2*sigmaV2;
            sigmaV6=sigmaV4*sigmaV2;
            
            %Transpose equals inverse of a rotation matrix.
            MInvCur=M(:,:,curMeas)';
        
            %Equation 31 in [2].
            w=sqrt(abs(1-u^2-v^2));%Abs added to deal with bad inputs.
            w2=w*w;
            w3=w2*w;
            w4=w3*w;
            w5=w4*w;
            w6=w5*w;
            w7=w6*w;
            w8=w7*w;
            w9=w8*w;
            w10=w9*w;
            
            r2=r*r;
            u2=u*u;
            u3=u2*u;
            u4=u3*u;
            u6=u4*u2;
            v2=v*v;
            v3=v2*v;
            v4=v3*v;
            v6=v4*v2;
            
            %Equation 32 in [2].
            cz=-(1/2)*(r/w+r*u2/w3)*sigmaU2-(1/2)*(r/w+r*v2/w3)*sigmaV2;
            
            %Equation 30 in [2].
            x=r*u;
            %Equation 30 in [2].
            y=r*v;
            %Equation 30 in [2] (the sign of cz was wrong and has been
            %corrected).
            z=r*w+cz;
            
            %The bias removal algorithm breaks down at extreme angles. This
            %should not occur in practice. It is assumed that the target is
            %in front of the radar, so the z component should be positive.
            z=max(z,0);
            zCart(:,curMeas)=MInvCur*[x;y;z]+zRx(:,curMeas);
            
            %Equation 33 in [2].
            sigmaXX2=u2*sigmaR2+r2*sigmaU2+sigmaR2*sigmaU2;
            %Equation 34 in [2].
            sigmaYY2=v2*sigmaR2+r2*sigmaV2+sigmaR2*sigmaV2;
            %Equation 35 in [2].
            sigmaZZ2=w2*sigmaR2+(r2/w2)*(u2*sigmaU2+v2*sigmaV2)-sigmaR2*sigmaU2-sigmaR2*sigmaV2+(7*r2*u2*v2/w6+r2*(u2+v2)/w4)*sigmaU2*sigmaV2...
                     +(1/2)*(7*r2*u4/w6+8*r2*u2/w4+r2/w2)*sigmaU4...
                     +(1/2)*(7*r2*v4/w6+8*r2*v2/w4+r2/w2)*sigmaV4...
                     +(3/4)*(u4/w6+2*u2/w4+1/w2)*sigmaR2*sigmaU4...
                     +(3/4)*(v4/w6+2*v2/w4+1/w2)*sigmaR2*sigmaV4...
                     +(1/4)*(45*r2*u4*v2/w10+(36*r2*u2*v2+6*r2*u4)/w8+(6*r2*u2+3*r2*v2)/w6)*sigmaU4*sigmaV2...
                     +(1/4)*(45*r2*u2*v4/w10+(36*r2*u2*v2+6*r2*v4)/w8+(6*r2*v2+3*r2*u2)/w6)*sigmaU2*sigmaV4...
                     +(15/4)*(r2*u6/w10+2*r2*u4/w8+r2*u2/w6)*sigmaU6...
                     +(15/4)*(r2*v6/w10+2*r2*v4/w8+r2*v2/w6)*sigmaV6...
                     +(1/2)*(3*u2*v2/w6+(u2+v2)/w4+1/w2)*sigmaR2*sigmaU2*sigmaV2;
            %Equation 36 in [2].
            sigmaXY2=u*v*sigmaR2;
            %Equation 37 in [2].
            sigmaXZ2=u*w*sigmaR2-(r2*u/w)*sigmaU2-(1/2)*(u3/w3+3*u/w)*sigmaR2*sigmaU2-(1/2)*(u*v2/w3+u/w)*sigmaR2*sigmaV2...
                     -(1/2)*(3*r2*u*v2/w5+r2*u/w3)*sigmaU2*sigmaV2-(3/2)*(r2*u3/w5+r2*u/w3)*sigmaU4;
            %Equation 38 in [2].
            sigmaYZ2=v*w*sigmaR2-(r2*v/w)*sigmaV2-(1/2)*(u2*v/w3+v/w)*sigmaR2*sigmaU2-(1/2)*(v3/w3+3*v/w)*sigmaR2*sigmaV2...
                     -(1/2)*(3*r2*u2*v/w5+r2*v/w3)*sigmaU2*sigmaV2-(3/2)*(r2*v3/w5+r2*v/w3)*sigmaV4;
         
            RCart(:,:,curMeas)=MInvCur*[sigmaXX2, sigmaXY2, sigmaXZ2;
                                        sigmaXY2, sigmaYY2, sigmaYZ2;
                                        sigmaXZ2, sigmaYZ2, sigmaZZ2]*MInvCur';    
        end
  case 3%The CM2 algorithm from [1]. This uses a second-order Taylor
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
            temp1=sqrt(abs(1-u^2-v^2));%Abs added to deal with bad inputs.
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
            %should not occur in practice. It is assumed that the target is
            %in front of the radar, so the z component should be positive.
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
