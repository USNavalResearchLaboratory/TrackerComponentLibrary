function H=HessianOfAffineTransFun(H,M)
%HESSIANOFAFFINETRANSFUN Let H be the Hessian hyper matrix of a vector
%           function h(x). That is a 3D matrix of second derivatives of the
%           function h such with respect to the elements of x. This takes H
%           and returns H transformed such that it is the Hessian matrix of
%           M*h(x) instead of h(x). That is, obtain the Hessian of an
%           arbitrary affine transformation of the original function (e.g.
%           M*h(x)+c, where c can be any constant vector). Note that if one
%           wants the Hessian of h(M*x+c) as opposed to of M*h(x)+c, then
%           one should use the function HessianChainRule.
%
%INPUTS: H An xDimXxDimXhDim hypermatrix such that H(:,:,1) is the Hessian
%          matrix for the first component of h(x), H(:,:,2) is the Hessian
%          matrix for the second component of h(x), etc. The ordering of
%          the derivatives in each matrix is:
%          [d^2/(dx(1)dx(1)), d^2/(dx(1)dx(2)), d^2/(dx(1)dx(3)),...
%           d^2/(dx(2)dx(1)), d^2/(dx(2)dx(2)), d^2/(dx(2)dx(3)),...
%           d^2/(dx(3)dx(1)), d^2/(dx(3)dx(2)), d^2/(dx(3)dx(3)),...
%          note that each matrix in the hypermatrix is symmetric.
%        M An hDimXhDim matrix.
%
%OUTPUTS: H The Hessian of M*h(x).
%
%EXAMPLE:
%Here, we consider the measurement function h(x) being Cart2Sphere. The
%Jacobian of M*h(x) for any h(x) is just M*J, where J is the Jacobian
%matrix of j(x). Here, we compute the Hessian matrix of M*h(x) via numeric
%differentiation of the Jacobian and also by transforming the analysic
%Hessian using HessianOfAffineTransFun. The Relative error between them is
%about 3e-10, which indicates agreement within reasonable finite precision
%limits.
% M=magic(3);
% x=[12;-8;16];
% s=[];
% s.type="()";
% HNumDiff=zeros(3,3,3);
% for i1=1:3
%     s.subs={':',i1};
%     J=@(x)subsref(M*calcSpherJacob(x),s);
%     HNumDiff(:,:,i1)=numDiff(x,J,3);
% end
% HNumDiff=permute(HNumDiff,[3,2,1]);
% H=calcSpherHessian(x);
% H=HessianOfAffineTransFun(H,M);
% RelErr=max(abs((HNumDiff(:)-H(:))./HNumDiff(:)))
%
%January 2021 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    H=permute(H,[3,2,1]);
    numK=size(H,3);
    for k=1:numK
        H(:,:,k)=M*H(:,:,k);
    end
    H=permute(H,[3,2,1]); 
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
