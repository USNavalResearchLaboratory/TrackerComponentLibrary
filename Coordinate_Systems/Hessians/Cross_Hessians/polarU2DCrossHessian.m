function H=polarU2DCrossHessian(uList,systemType)
%%POLARU2DCROSSHESSIAN Given the direction cosine value u in 2D, obtain the
%              second derivative of the polar azimuth angle with respect to
%              u. This is the derivative of the output of u2PolAng2D with
%              respect to its input.
%
%INPUTS: uList A 1XnumPoints (for only u) or a 2XnumPoints (if full unit
%              vectors are given) set of direction cosines in 2D.
%   systemType An optional parameter specifying the axis from which the
%              angles are measured. Possible values are
%              0 (The default if omitted or an empty matrix is passed) The
%                azimuth angle is counterclockwise from the x axis.
%              1 The azimuth angle is measured clockwise from the y axis.
%
%OUTPUTS: H A 1X1XnumPoints or 2X2XnumPoints (if both u and v are provided)
%           set of second derivatives of the azimuth angle. If this is a
%           1X1XnumPoints matrix, then the second derivatives are with
%           respect to u. If this is a 2X2XnumPoints matrix, then the
%           seocnd derivatives are
%           [d2u/dudu,d2ududv,
%            d2u/dvdu,d2udvdv];
%
%Note that u and v are not independent. Thus, in the comparison to finite
%differencing below (Example 2), when we offset u, we change v and vice
%versa. This dependence means that unlike normal, each H matrix is NOT
%symmetric and it is singular. It also means that the off-diagonal terms
%are generally not very useful.
%
%EXAMPLE 1:
%Here, we verify that the derivatives returned by this function are about
%equal to those returned via numeric differentiation (forward
%differencing).
% points=[0.1,0.2,-0.2,0,-0.9];%u points
% systemType=0;
% epsVal=1e-9;
% 
% J=polarU2DCrossGrad(points,systemType);
% J1=polarU2DCrossGrad(points+epsVal,systemType);
% HNumDiff=(J1-J)/epsVal;
% H=polarU2DCrossHessian(points,systemType);
% max(abs(HNumDiff(:)-H(:)))
%One will see that the difference is on the order of 5e-7, which is a good
%agreement.
%
%EXAMPLE 2:
%This is the same as the first example, except both a u and a v component
%are provided, so derivatives with respect to each term are provided. The
%relative errors are O(1e-6) or O(1-e7), which is good agreement.
% points=[0.1,0.2,-0.2,-0.1,-0.9];%u
% points(2,:)=sqrt(1-points(1,:).^2);%v
% %Do the same thing, but switch the sign of v.
% points=[points,[points(1,:);-points(2,:)]];
% numPts=size(points,2);
% systemType=0;
% epsVal=1e-9;
% 
% J=polarU2DCrossGrad(points,systemType);
% zDiff(1,:)=points(1,:)+epsVal;
% zDiff(2,:)=sign(points(2,:)).*sqrt(1-zDiff(1,:).^2);
% J1=polarU2DCrossGrad(zDiff,systemType);
% HNumDiffDeltaU=reshape((J1-J)/epsVal,[2,1,numPts]);
% zDiff(2,:)=points(2,:)+epsVal;
% zDiff(1,:)=sign(points(1,:)).*sqrt(1-zDiff(2,:).^2);
% J1=polarU2DCrossGrad(zDiff,systemType);
% HNumDiffDeltaV=reshape((J1-J)/epsVal,[2,1,numPts]);
% 
% HNumDiff=[HNumDiffDeltaU,HNumDiffDeltaV];
% H=polarU2DCrossHessian(points,systemType);
% RelErr=max(abs(HNumDiff(:)-H(:))./HNumDiff(:))
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
    systemType=0; 
end

hasV=size(uList,1)>1;

N=size(uList,2);

numDim=1+hasV;
H=zeros(numDim,numDim,N);
for curPoint=1:N
    u=uList(1,curPoint);
    
    if(hasV)
        v=uList(2,curPoint);
    else
        v=sqrt(1-u^2);
    end

    switch(systemType)
        case 0
            H(1,1,curPoint)=-u./v.^3;
            if(hasV)
                H(1,2,curPoint)=1./v.^2;
                H(2,1,curPoint)=-1/u^2;
                H(2,2,curPoint)=v./u.^3;
            end
        case 1
            H(1,1,curPoint)=u./v.^3;
            if(hasV)
                H(1,2,curPoint)=-1./v.^2;
                H(2,1,curPoint)=1/u^2;
                H(2,2,curPoint)=-v./u.^3;
            end
        otherwise
            error('Invalid system type specified.')
    end
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
