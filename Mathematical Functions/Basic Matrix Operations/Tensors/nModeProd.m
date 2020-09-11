function Y=nModeProd(X,U,n)
%%NMODEPROD Given a real or complex N-mode tensor A, whose dimensions are
%           I_1XI_2...XI_N, and a JXI_n matrix U, compute the n-mode
%           matrix product of the tensor and the matrix.
%           
%INPUTS: X The hypermatrix (tensor) that is to be involved in the n-mode
%          product. This should have 2 or more modes. For example, a
%          4-mode hyper matrix would be addressable as X(i1,i2,i3,i4).
%          The hypermatrix can be real or complex.
%        U The matrix by which the tensor is to be multiplied. The
%          dimensionality of the matrix is JXI_n, where J can be any
%          positive integer, but I_n has to match the  dimensionality of
%          the nth mode of X.
%        n The integer mode over which the multiplication is to be
%          performed. This ranges from 1 to N, the number of modes in X.
%
%OUTPUTS: Y The n-mode product of X and U. This has size
%           I_1 X...XI_{n-1} X J X I_{n+1} X... XI_N. That is, it has as
%           many modes as X, but the dimensionality of the nth mode has
%           changed to J.
%
%The n-mode product of tensors is defined in [1] and [2]. The
%implementation is using the tensor to matrix conversions as described in
%[3].
%
%EXAMPLE 1:
% X(:,:,1)=[1,4,7,10;
%           2,5,8,11;
%           3,6,9,12];
%       
% X(:,:,2)=[13,16,19,22;
%           14,17,20,23;
%           15,18,21,24];
% U=[1,3,5;
%    2,4,6];
% n=1;
% Y=nModeProd(X,U,n)
%Would produce the results
% Y(:,:,1)=[22    49    76   103;
%           28    64   100   136];
% Y(:,:,2)=[130   157   184   211;
%           172   208   244   280]
%
%EXAMPLE 2:
%This is the same as Example 3.5 in [2]. We have an orthonormal matrix
%(C'*C=eye(2,2)), and we multiply and then rerecover the original vector
%demonstrating  that for an orthonormal matrix C, that
%X=nModeProd(Y,C,n) can be reversed as Y=nModeProd(X,C',n)
% Y(:,:,1)=[1,4,7,10;
%           2,5,8,11;
%           3,6,9,12];
% Y(:,:,2)=[13,16,19,22;
%           14,17,20,23;
%           15,18,21,24];
% C=[0.58,0;
%    0.58,-0.71;
%    0.58,0.71];
% C(:,1)=C(:,1)/norm(C(:,1));
% C(:,2)=C(:,2)/norm(C(:,2));
% X=nModeProd(Y,C,3);
% Z=nModeProd(X,C',3);
% max(abs(Z(:)-Y(:)))
%One will see that Z and Y are equal within reasonable finite precision
%bounds.
%
%REFERENCES:
%[1] R. G. Kolda and B. W. Bader, "Tensor decompositions and applications,"
%    SIAM Review, vol. 51, no. 3, pp. 455-500, 2009.
%[2] T. G. Kolda, "Multilinear operators for higher-order decompositions,"
%    Sandia National Laboratories, Tech. Rep. SAND2006-2081, Apr. 2006.
%    [Online]. Available: http://www.sandia.gov/~tgkolda/pubs/pubfiles/SAND2006-2081.pdf
%[3] J. Salmi, A. Richter, and V. Koivunen, "Sequential unfolding SVD for
%    tensors with applications in array signal processing," IEEE
%    Transactions on Signal Processing, vol. 57, no. 12, pp. 4719-4733,
%    Dec. 2009.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

MVals=size(X);
N=length(MVals);%The number of ways of the matrix.

Rn=size(U,1);
AU=reshape(U*tensor2Mat(X,n,[(n+1):1:N,1:1:(n-1)]),[Rn,MVals([(n+1):1:N,1:1:(n-1)])]);
Y=ipermute(AU,[n:1:N,1:1:(n-1)]);

%Note that the above three lines are equivalent to the following commented
%out section, which might be easier to understand:
% J=size(U,1);
% YDims=MVals;
% YDims(n)=J;
% Y=zeros(YDims);
% 
% sumDims=MVals;
% sumDims(n)=1;
% 
% xIdx=repmat({':'},[1,N]);
% yIdx=repmat({':'},[1,N]);
% for j=1:J
%     yIdx{n}=j;
%     sumVal=zeros(sumDims);
%     for i=1:MVals(n)
%         xIdx{n}=i;
%         sumVal=sumVal+X(xIdx{:})*U(j,i);
%     end
%     
%     Y(yIdx{:})=sumVal;
% end

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