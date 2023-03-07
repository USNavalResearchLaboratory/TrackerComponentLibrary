function H=rangeHessian(x,useHalfRange,lTx,lRx)
%%RANGEHESSIAN Determine the Hessian matrix (a matrix of second-order
%          partial derivatives) of a bistatic range rate measurement with
%          respect to Cartesian position. Atmospheric and other propagation
%          effects are not taken into account. The transmitter can be
%          collocated with the target  (in which case it is assumed that
%          they move together).
%
%INPUTS: x A numPosDimXN set of N target position vectors of the form
%          [x], [x;y] or [x;y;z].
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided is false.
%      lTx The numPosDimX1 position vector of the transmitter. If this
%          parameter is omitted or an empty matrix is passed, then the
%          transmitter is assumed to be at the origin.
%      lRx The numPosDimX1 position vector of the receiver. If this
%          parameter is omitted or an empty matrix is passed, then the
%          receiver is assumed to be at the origin.
%
%OUTPUTS: H A numPosDimXnumPosDimXN set of N Hessian matrices, one for each
%           point in x. The ordering of the derivatives in H(:,:,i) is
%           [d^2r/(dxdx), d^2r/(dxdy), d^2r/(dxdz);
%            d^2r/(dydx), d^2r/(dydy), d^2r/(dydz);
%            d^2r/(dzdx), d^2r/(dzdy), d^2r/(dzdz)];
%           note that each matrix is symmetric (i.e.
%           d^2r/(dydx)=d^2r/(dxdy) ).
%
%Expressions for bistatic range are given in [1]. The Hessian comes from
%simple differentiation of the expression.
%
%EXAMPLE 1:
%Here, we verify that a numerically differentiated Hessian is consistent
%with the analytic one produced by this function.
% x=[100;-1000;500];
% lRx=[500;20;-400];
% lTx=[-1000;10;2];
% useHalfRange=true;
% epsVal=1e-5;
% 
% gradVal=rangeGradient(x,useHalfRange,lTx,lRx);
% gradValdX=rangeGradient(x+[epsVal;0;0],useHalfRange,lTx,lRx);
% gradValdY=rangeGradient(x+[0;epsVal;0],useHalfRange,lTx,lRx);
% gradValdZ=rangeGradient(x+[0;0;epsVal],useHalfRange,lTx,lRx);
% HNumDiff=zeros(3,3);
% HNumDiff(1,1)=(gradValdX(1)-gradVal(1))/epsVal;
% HNumDiff(1,2)=(gradValdX(2)-gradVal(2))/epsVal;
% HNumDiff(2,1)=HNumDiff(1,2);
% HNumDiff(1,3)=(gradValdX(3)-gradVal(3))/epsVal;
% HNumDiff(3,1)=HNumDiff(1,3);
% HNumDiff(2,2)=(gradValdY(2)-gradVal(2))/epsVal;
% HNumDiff(2,3)=(gradValdY(3)-gradVal(3))/epsVal;
% HNumDiff(3,2)=HNumDiff(2,3);
% HNumDiff(3,3)=(gradValdZ(3)-gradVal(3))/epsVal;
% 
% H=rangeHessian(x,useHalfRange,lTx,lRx);
% 
% max(abs((H(:)-HNumDiff(:))./H(:)))
%The relative error will be on the order of 9e-7, indicating good agreement
%between the numerical Hessian matrix and the actual Hessian matrix.
%
%EXAMPLE 2:
%This example is similar to the first one, except the transmitter is on the
%target.
% x=[100;-1000;500];
% lRx=[500;20;-400];
% useHalfRange=false;
% epsVal=1e-4;
% 
% gradVal=rangeGradient(x,useHalfRange,x,lRx);
% gradValdX=rangeGradient(x+[epsVal;0;0],useHalfRange,x+[epsVal;0;0],lRx);
% gradValdY=rangeGradient(x+[0;epsVal;0],useHalfRange,x+[0;epsVal;0],lRx);
% gradValdZ=rangeGradient(x+[0;0;epsVal],useHalfRange,x+[0;0;epsVal],lRx);
% HNumDiff=zeros(3,3);
% HNumDiff(1,1)=(gradValdX(1)-gradVal(1))/epsVal;
% HNumDiff(1,2)=(gradValdX(2)-gradVal(2))/epsVal;
% HNumDiff(2,1)=HNumDiff(1,2);
% HNumDiff(1,3)=(gradValdX(3)-gradVal(3))/epsVal;
% HNumDiff(3,1)=HNumDiff(1,3);
% HNumDiff(2,2)=(gradValdY(2)-gradVal(2))/epsVal;
% HNumDiff(2,3)=(gradValdY(3)-gradVal(3))/epsVal;
% HNumDiff(3,2)=HNumDiff(2,3);
% HNumDiff(3,3)=(gradValdZ(3)-gradVal(3))/epsVal;
% 
% H=rangeHessian(x,useHalfRange,x,lRx);
% max(abs((H(:)-HNumDiff(:))./H(:)))
%The relative error will be on the order of 9.5e-8, indicating good
%agreement between the numerical Hessian matrix and the actual Hessian
%matrix.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%June 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);

if(nargin<2||isempty(useHalfRange))
   useHalfRange=false; 
end

if(nargin<4||isempty(lRx))
    lRx=zeros(numDim,1); 
end

if(nargin<3||isempty(lTx))
    lTx=zeros(numDim,1); 
end

N=size(x,2);
H=zeros(numDim,numDim,N);
for curX=1:N
    deltaRx=x(:,curX)-lRx;
    normDeltRx=norm(deltaRx);
    deltaTx=x(:,curX)-lTx;
    normDeltTx=norm(deltaTx);
    
    H(1,1,curX)=-deltaRx(1)^2/normDeltRx^3+1/normDeltRx;

    if(normDeltTx~=0)
        %If the target is not collocated with the transmitter.
        H(1,1,curX)=H(1,1,curX)-deltaTx(1)^2/normDeltTx^3+1/normDeltTx;
    end
    
    if(numDim>1)
        H(2,2,curX)=-deltaRx(2)^2/normDeltRx^3+1/normDeltRx;
        H(1,2,curX)=-deltaRx(1)*deltaRx(2)/normDeltRx^3;

        if(normDeltTx~=0)
            H(2,2,curX)=H(2,2,curX)-deltaTx(2)^2/normDeltTx^3+1/normDeltTx;
            H(1,2,curX)=H(1,2,curX)-deltaTx(1)*deltaTx(2)/normDeltTx^3;
        end

        H(2,1,curX)=H(1,2,curX);

        if(numDim>2)
            H(3,3,curX)=-deltaRx(3)^2/normDeltRx^3+1/normDeltRx;
            H(1,3,curX)=-deltaRx(1)*deltaRx(3)/normDeltRx^3;
            H(2,3,curX)=-deltaRx(2)*deltaRx(3)/normDeltRx^3;

            if(normDeltTx~=0)
                H(3,3,curX)=H(3,3,curX)-deltaTx(3)^2/normDeltTx^3+1/normDeltTx;
                H(1,3,curX)=H(1,3,curX)-deltaTx(1)*deltaTx(3)/normDeltTx^3;
                H(2,3,curX)=H(2,3,curX)-deltaTx(2)*deltaTx(3)/normDeltTx^3;
            end

            H(3,1,curX)=H(1,3,curX);
            H(3,2,curX)=H(2,3,curX);
        end
    end
end

if(useHalfRange)
    H=H/2; 
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
