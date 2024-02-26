function [dpinvAdx,dpinvAdy,d2pinvAdxdy]=pseudoInverseDerivative(A,dAdx,dAdy,d2Adxdy,pinvA)
%%PSEUDOINVERSEDERIAVTIVE Given an mXn matrix A and its derivative with
%       respect to a parameter x, this function finds the derivative of the
%       pseudoinverse of A, pinv(A). Additionally, if one provides the
%       derivative of A with respect to y (which could be the same as x or
%       different) as well as the second derivative of A with respect to x
%       and y this function can return the second derivative of pinv(A)
%       with respect to x and y.
%
%INPUTS: A An mXn matrix.
%     dAdx The mXn derivative of A with respect to a scalar x.
%     dAdy The mXn derivative of A with respect to a scalar y (only needed
%          if the second and third outputs are desired).
%  d2Adxdy The second derivative of A with respect to x and y (only needed
%          if the second and third outputs are desired).
%    pinvA Optionally, one can pass the nXm psuedoinverse of A. If this is
%          omitted or an empty matrix is passed, then it is just computed
%          here as pinv(A).
%
%OUTPUTS: dpinvAdx The nXm derivative of pinv(A) with respect to x.
%         dpinvAdy The nXm derivative of pinv(A) with respect to y.
%      d2pinvAdxdy The nXm second derivative of pinv(A) with respect to x
%                  and y.
%
%This implements Equation 4.12 of [1]. Note that there are discontinuities
%in pseudoinverse function. This solution is bsed on an assumed local
%continuity at the point where the derivative is taken.
%
%EXAMPLE 1:
%This example compares the first deriavtive results to finite differencing
%with a function that is quadratic in x. One can see that the results agree
%to more than 9 decimal places, indicating good agreement.
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
% AFun=@(x)(C0+C1*x+C2*x^2);
% AFunDeriv=@(x)(C1+2*C2*x);
% 
% x=3;
% A=AFun(x);
% dAdx=AFunDeriv(x);
% dpinvAdx=pseudoInverseDerivative(A,dAdx);
% 
% f=@(x)vec(pinv(AFun(x)));
% dpinvAdxNumDiff=reshape(numDiff(x,f,m*n),[n,m]);
% RelErr=max(max(abs((dpinvAdxNumDiff-dpinvAdx)./dpinvAdxNumDiff)))
%
%EXAMPLE 2:
%In this example, the second derivative with respect to x and y is
%considered. It is shown to agree with finite differencing to more than 8
%decimal places, which is good agreement.
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
% AFun=@(x,y)(C0*y^3+C1*x*y^2+C2*x^2);
% AFunDerivx=@(x,y)(C1*y^2+2*C2*x);
% AFunDerivy=@(x,y)(3*C0*y^2+2*C1*x*y);
% AFunDerivxy=@(x,y)2*C1*y;
% 
% x=3;
% y=2;
% A=AFun(x,y);
% dAdx=AFunDerivx(x,y);
% dAdy=AFunDerivy(x,y);
% d2Adxdy=AFunDerivxy(x,y);
% 
% [~,~,d2pinvAdxdy]=pseudoInverseDerivative(A,dAdx,dAdy,d2Adxdy);
% f=@(y)vec(pseudoInverseDerivative(AFun(x,y),AFunDerivx(x,y)));
% d2pinvAdxdyNumDiff=reshape(numDiff(y,f,m*n),[n,m]);
% RelErr=max(max(abs((d2pinvAdxdyNumDiff-d2pinvAdxdy)./d2pinvAdxdyNumDiff)))
%
%EXAMPLE 3:
%In this example, the second derivative with respect to just x is
%considered. It is shown to agree with finite differencing to more than 8
%decimal places, which is good agreement.
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
% AFun=@(x,y)(C0+C1*x+C2*x^3);
% AFunDerivx=@(x,y)(C1+3*C2*x^2);
% AFunDerivx2=@(x,y)6*C2*x;
% 
% x=3;
% A=AFun(x);
% dAdx=AFunDerivx(x);
% d2Adx2=AFunDerivx2(x);
% [~,~,d2pinvAdxdy]=pseudoInverseDerivative(A,dAdx,dAdx,d2Adx2);
% f=@(x)vec(pseudoInverseDerivative(AFun(x),AFunDerivx(x)));
% d2pinvAdx2NumDiff=reshape(numDiff(x,f,m*n),[n,m]);
% RelErr=max(max(abs((d2pinvAdx2NumDiff-d2pinvAdxdy)./d2pinvAdx2NumDiff)))
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

m=size(A,1);
n=size(A,2);
I1=eye(m,m);
I2=eye(n,n);

diff1=(I1-A*pinvA);
diff2=(I2-pinvA*A);

dpinvAdx=-pinvA*dAdx*pinvA+(pinvA*pinvA')*(dAdx')*diff1+diff2*(dAdx')*(pinvA'*pinvA);

if(nargout>1)
    %Compute the sother first derivative.
    dpinvAdy=-pinvA*dAdy*pinvA+(pinvA*pinvA')*(dAdy')*diff1+diff2*(dAdy')*(pinvA'*pinvA);

    if(nargout>2)
        %Compute the second derivative.
        d2pinvAdxdy=-dpinvAdy*dAdx*pinvA-pinvA*d2Adxdy*pinvA-pinvA*dAdx*dpinvAdy ...
                    +(dpinvAdy*pinvA'+pinvA*dpinvAdy')*(dAdx')*diff1+pinvA*pinvA'*(d2Adxdy')*diff1-pinvA*pinvA'*(dAdx')*(dAdy*pinvA+A*dpinvAdy)...
                    -(dpinvAdy*A+pinvA*dAdy)*(dAdx')*(pinvA')*pinvA+diff2*(d2Adxdy')*(pinvA')*pinvA+diff2*(dAdx')*(dpinvAdy'*pinvA+pinvA'*dpinvAdy);
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
