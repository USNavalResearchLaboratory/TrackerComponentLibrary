function [dpinvAdx,d2pinvAdxdy]=pseudoInverseDerivative(A,dAdx,d2Adxdy,pinvA)
%%PSEUDOINVERSEDERIAVTIVE Given an mXn matrix A and its derivatives with
%       respect to ones or more parameters x, this function finds the
%       first derivatives of the pseudoinverse of A, pinv(A). Additionally,
%       if one provides second derivatives of A with respect to the
%       parameters in x, this function can return the second deriavtives of
%       pinv(A).
%
%INPUTS: A An mXn matrix.
%     dAdx The mXnXN set of partial derivatives of A with respect to the
%          elements of the length N vector x.
%  d2Adxdy If the second output of this function is not desired, then this
%          input may be omitted. Otherwise, this is an mXnXNXN set of
%          second partial derivatives of the matrix A with respect to the
%          elements of a length n parameter vector x. d2Adxdy 
%    pinvA Optionally, one can pass the nXm psuedoinverse of A. If this is
%          omitted or an empty matrix is passed, then it is just computed
%          here as pinv(A).
%
%OUTPUTS: dpinvAdx The nXmXN set of derivatives of pinv(A) with respect to
%                  the elements of a vector x. dpinvAdx(:,:,i) holds the
%                  derivative of pinv(A) with respect to x(i).
%      d2pinvAdxdy The nXmXNXN set of second derivative of pinv(A) with
%                  respect to the elements of x. d2pinvAdxdy(:,:,,i,j) holds
%                  the second derivative of pinv(A) with respect to x(i)
%                  and x(j).
%
%This implements Equation 4.12 of [1]. Note that there are discontinuities
%in pseudoinverse function. This solution is based on an assumed local
%continuity at the point where the derivative is taken.
%
%EXAMPLE:
%This example considers the relative error of the first and second partial
%derivatives returnes by this function and those obtained using finite
%differencing. One can see that the results agree to more than 6 decimal
%places, which is reasonably good.
% C0=[-14, -9;
%       1,  0;
%     -14, -8];
% C1=[100, -179;
%      94, -170;
%      95, -171];
% C2=[74,   52;
%    168,   118;
%     50,   35];
% m=size(C0,1);
% n=size(C0,2);
% 
% AFun=@(x,y)C0*y^3+C1*x*y^2+C2*x^2;
% AFunDerivx=@(x,y)(C1*y^2+2*C2*x);
% AFunDerivy=@(x,y)(3*C0*y^2+2*C1*x*y);
% AFunDerivxx=@(x,y)reshape(2*C2,[1,1,m,n]);
% AFunDerivyy=@(x,y)reshape(6*C0*y+2*C1*x,[1,1,m,n]);
% AFunDerivxy=@(x,y)reshape(2*C1*y,[1,1,m,n]);
% AMat1stDerivs=@(x,y)cat(3,AFunDerivx(x,y),AFunDerivy(x,y));
% AMat2ndDerivs=@(x,y)permute([AFunDerivxx(x,y),AFunDerivxy(x,y);
%                             AFunDerivxy(x,y),AFunDerivyy(x,y)],[3,4,1,2]);
% 
% x=3;
% y=2;
% A=AFun(x,y);
% pinvA=pinv(A);
% dA=AMat1stDerivs(x,y);
% %First, verify that the first derivatives returned by this function are
% %consistent with finite differencing.
% J=pseudoInverseDerivative(A,dA,[],pinvA);
% epsVal=1e-8;
% f=@(x,y)pinv(AFun(x,y));
% pinvAdX=f(x+epsVal,y);
% pinvAdY=f(x,y+epsVal);
% JNumDiff=cat(3,(pinvAdX-pinvA)/epsVal,(pinvAdY-pinvA)/epsVal);
% RelError=max(abs((J(:)-JNumDiff(:))./JNumDiff(:)))
% 
% %Next, verify that the second derivatives returned by this function are
% %consistent with finite differencing.
% d2A=AMat2ndDerivs(x,y);
% [~,H]=pseudoInverseDerivative(A,dA,d2A,pinvA);
% AdX=AFun(x+epsVal,y);
% pinvAdX=pinv(A);
% dAdX=AMat1stDerivs(x+epsVal,y);
% AdY=AFun(x,y+epsVal);
% pinvAdY=pinv(A);
% dAdY=AMat1stDerivs(x,y+epsVal);
% JdX=pseudoInverseDerivative(AdX,dAdX,[],pinvAdX);
% JdY=pseudoInverseDerivative(AdY,dAdY,[],pinvAdY);
% 
% HNumDiff=zeros(n,m,2,2);
% HNumDiff(:,:,1,1)=(JdX(:,:,1)-J(:,:,1))/epsVal;
% HNumDiff(:,:,1,2)=(JdX(:,:,2)-J(:,:,2))/epsVal;
% HNumDiff(:,:,2,1)=HNumDiff(:,:,1,2);
% HNumDiff(:,:,2,2)=(JdY(:,:,2)-J(:,:,2))/epsVal;
% RelErr=max(abs((H(:)-HNumDiff(:))./HNumDiff(:)))
%
%REFERENCES:
%[1] G. H. Golub and V. Pereyra, "The differentiation of pseudo-inverses
%    and nonlinear least squares problems whose variables separate," SIAM
%    Journal on Numerical Analysis, vol. 10, no. 2, pp. 413-432, Apr. 1973.
%
%April 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(pinvA))
    pinvA=pinv(A);
end

numVar=size(dAdx,3);

m=size(A,1);
n=size(A,2);
I1=eye(m,m);
I2=eye(n,n);

diff1=(I1-A*pinvA);
diff2=(I2-pinvA*A);

dpinvAdx=zeros(n,m,numVar);
for k=1:numVar
    dpinvAdx(:,:,k)=-pinvA*dAdx(:,:,k)*pinvA+(pinvA*pinvA')*(dAdx(:,:,k)')*diff1+diff2*(dAdx(:,:,k)')*(pinvA'*pinvA);
end

if(nargout>1)
    %Compute the other first derivative.
    d2pinvAdxdy=zeros(n,m,numVar,numVar);

    for k1=1:numVar
        for k2=k1:numVar
            d2pinvAdxdy(:,:,k1,k2)=-dpinvAdx(:,:,k2)*dAdx(:,:,k1)*pinvA-pinvA*d2Adxdy(:,:,k1,k2)*pinvA-pinvA*dAdx(:,:,k1)*dpinvAdx(:,:,k2) ...
                    +(dpinvAdx(:,:,k2)*pinvA'+pinvA*dpinvAdx(:,:,k2)')*(dAdx(:,:,k1)')*diff1+pinvA*pinvA'*(d2Adxdy(:,:,k1,k2)')*diff1-pinvA*pinvA'*(dAdx(:,:,k1)')*(dAdx(:,:,k2)*pinvA+A*dpinvAdx(:,:,k2))...
                    -(dpinvAdx(:,:,k2)*A+pinvA*dAdx(:,:,k2))*(dAdx(:,:,k1)')*(pinvA')*pinvA+diff2*(d2Adxdy(:,:,k1,k2)')*(pinvA')*pinvA+diff2*(dAdx(:,:,k1)')*(dpinvAdx(:,:,k2)'*pinvA+pinvA'*dpinvAdx(:,:,k2));
            d2pinvAdxdy(:,:,k2,k1)=d2pinvAdxdy(:,:,k1,k2);
        end
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
