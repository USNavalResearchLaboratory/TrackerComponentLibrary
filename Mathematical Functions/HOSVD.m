function [U,S]=HOSVD(X)
%%HOSVD Evaluate the high-order singular value decomposition (HOSVD) of a
%       tensor, also known as a Tucker decomposition. The Tucker
%       decomposition can have multiple solutions, and the HOSVD is just
%       one possible solution, though not necessarily the best one. it is,
%       however, simple to compute. The HOSVD is akin to a generalization
%       of the concept of singular value decomposition.  Though only the
%       left-singular vectors are used, the full matrix can still be
%       recreated. Thus, when run with a matrix input, this produces a
%       slightly different type of decomposition than the standard svd
%       command.
%
%INPUTS: X  The hypermatrix (tensor) whose Tucker decomposition is desired.
%           This should have 2 or more modes. For example, a 4-mode
%           hypermatrix would be addressable as X(i1,i2,i3,i4). The matrix
%           can be real or complex.
%
%OUTPUTS: U A cell array of orthogonal matrices from the decomposition.
%         S The core tensor of the decomposition. See the example below to
%         see how the outputs are used to reconstruct X.
%
%The algorithm is written to be consistent with the default n-way
%matricization operation of tensor2Mat. An example of how the decomposition
%can be used to decompose and then reconstruct a tensor is as follows:
% X=zeros(3,2,3);
% X(1,1,1)=1;
% X(1,1,2)=1;
% X(2,1,1)=1;
% X(2,1,2)=-1;
% X(2,1,3)=2;
% X(3,1,1)=2;
% X(3,1,3)=2;
% X(1,2,1)=2;
% X(1,2,2)=2;
% X(2,2,1)=2;
% X(2,2,2)=-2;
% X(2,2,3)=4;
% X(3,2,1)=4;
% X(3,2,3)=4;
% 
% [U,S]=HOSVD(X);
% 
% XRec=S;
% for n=1:length(U)
%     XRec=nModeProd(XRec,U{n},n);
% end
%
%One will observe that XRec is equal (within finite precision bounds) 
%to the original tensor X. 
%
%The algorithm is as in [1], though the definition of the n-way
%matricization is consistent with the more commonly used version, as in
%[2].
%
%REFERENCES:
%[1] L. de Lathauwer, B. de Moore, and J. Vandewalle, "A multilinear
%    singular value decomposition," SIAM Journal on Matrix Analysis and
%    Applications, vol. 21, no. 4, pp. 1253-1278, 2000.
%[2] R. G. Kolda and B. W. Bader, "Tensor decompositions and applications,"
%    SIAM Review, vol. 51, no. 3, pp. 455?500, 2009.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of dimensions of S
N=length(size(X));

U=cell(N,1);

S=X;
for curDim=1:N
    AMode=tensor2Mat(X,curDim);

    [UCur,~,~]=svd(AMode,'econ');
    U{curDim}=UCur;
    S=nModeProd(S,UCur',curDim);
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
