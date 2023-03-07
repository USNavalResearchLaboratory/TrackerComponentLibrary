function HTrans=HessianChainRule(Hf,Hg,Jf,Jg)
%%HESSIANCHAINRULE Suppose that one has a (scalar or vector) function f(x),
%    where x=[x1;...;xn] is a vector, whose gradient matrix and Hessian
%    (matrix or 3D matrix) can be evaluated. However, one wants the
%    second derivatives not of f(x) but of f(g(x)), where g is a function
%    that takes x and returns n outputs. This function converts the Hessian
%    of f given the gradient and Hessian of g.
%
%INPUTS: Hf An xDimXxDimXfDim matrix of second derivatives of the function
%           f(x) with respect to x evaluated at the point g(x). Hf(:,:,i)
%           is a matrix of second derivatives of the ith component of f.
%           For example, if f is a function of [x;y;z] Cartesian position,
%           then Hf(:,:,i) holds second derivatives in the order:
%             [d^2/(dxdx), d^2/(dxdy), d^2/(dxdz);
%              d^2/(dydx), d^2/(dydy), d^2/(dydz);
%              d^2/(dzdx), d^2/(dzdy), d^2/(dzdz)];
%        Hg An xDimXxDimXxDim matrix of second derivatives of the function
%           g evaluated at x. Hg(:,:,i) is the second derivative matrix for
%           the ith element of g. the ordering of the secon deriavtives is
%           the same as in Hf.
%        Jf An fDimXxDim matrix of first derivatives of the f function
%           evaluated at g(x). The rows are the elements of f and the
%           columns the element of x with which the derivative has been
%           taken. 
%        Jg An xDimXxDim matrix of derivatives of with respect to the
%           elements of x, evaluated at x. The rows select the dimension of
%           the output of g and the columns select with which element of x
%           the derivative is taken.
%
%OUTPUTS: HTrans The matrix of second derivatives of the function
%           f(g(x)) evaluated at x. The ordering of the elements if the
%           same is an Hf.
%
%This function is implemented using the identities in [1].
%
%EXAMPLE 1:
%In this example, we consider a rotated coordinate system compared to the
%desired coordinate system. The function calcSpherHessian has been coded to
%allow one to specify the rotation directly. Here, we consider compare the
%explicit Hessian with the rotation to computing the Hessian without the
%rotation and then transforming the Hessian. The relative error indicates
%agreement on the order of finite precision limitations.
% M=randRotMat(3);
% x=[12;-8;16];
% xTrans=M*x;
% HTrue=calcSpherHessian(x,0,true,[],[],M);
% %The Hessian not accounting for the rotation.
% Hf=calcSpherHessian(xTrans);
% %The Jacobian not accounting for the rotation.
% Jf=calcSpherJacob(xTrans);
% Hg=zeros(3,3,3);
% Jg=M;
% HTrans=HessianChainRule(Hf,Hg,Jf,Jg);
% RelErr=max(abs((HTrans(:)-HTrue(:))./HTrue(:)))
%
%EXAMPLE 2:
%Here, we use the example of a transofrmation of Rosenbrack;'s function
%from [1] and show that we get the correct result. It is a bivariate
%polynomial function, so an explicit solution for compairson is easily
%obtained. The asbolute difference in this case between the exact and the
%transformed solution is zero.
%f(x,y)=(1-x)^2+100*(x^2-y)^2;
%We use the transformation.
%g1(x,y)=x;
%g2(x,y)=x^2-y;
%
% Jf=@(xy)[-2*(1-xy(1))+400*xy(1)*(xy(1)^2-xy(2)), -200*(xy(1)^2-xy(2))];
% Hf=@(xy)[2+800*xy(1)^2+400*(xy(1)^2-xy(2)),-400*xy(1);
%         -400*xy(1),                 200];
% Jg=@(xy)[1,       0;
%          2*xy(1),-1];
% Hg=@(xy)cat(3, zeros(2,2),[2,0;0,0]);
% HfTransExact=[2,0;
%               0,200];
% xy=[15;64];
% HfTrans=HessianChainRule(Hf(xy),Hg(xy),Jf(xy),Jg(xy));
% AbsDiff=HfTrans-HfTransExact
%
%REFERENCES:
%[1] M. Skorski, "Chain rules for Hessian and higher derivatives made
%    easy by tensor calculus,"" arXiv, 29 Nov. 2019. [Online]. Available:
%    https://arxiv.org/abs/1911.13292
%
%May 2022 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDimOut=size(Hf,3);
numDimIn=size(Hf,1);

HTrans=zeros(numDimIn,numDimIn,numDimOut);
for k=1:numDimOut
    HTrans(:,:,k)=Jg'*Hf(:,:,k)*Jg;
    
    for k2=1:numDimIn
        HTrans(:,:,k)=HTrans(:,:,k)+Jf(k,k2)*Hg(:,:,k2);
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
