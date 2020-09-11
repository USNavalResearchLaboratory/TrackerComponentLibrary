function Y=tensor2Mat(X,param2,C)
%%TENSOR2MAT Unfold a (real or complex) tensor into a matrix. This can
%            either be done as an n-way matricization using a standard
%            ordering of the modes, or one can explicitly specify the
%            ordering of the modes. The tensor is just a hypermatrix. The
%            unfolding produces a 2D matrix as the output. The
%            transformation of tensors into matrices plays a role in
%            numerous multilinear algorithms, such as the multilinear
%            singular value decomposition.
%
%INPUTS: X The hypermatrix (tensor) that is to be unfolded. This should
%          have 2 or more modes. For example, a 4-mode hyper matrix would
%          be addressable as X(i1,i2,i3,i4). The matrix can be real or
%          complex.
% param2, C If the algorithm is called with only two parameters, then
%          param2=n is the number n of the mode of the unfolding and an
%          n-mode matricization is performed. The mode number ranges
%          from 1 through N --the number of modes of the matrix X (the
%          number of indices needed to address an element in X). The
%          permutation of the modes to the columns is as in [1], which is
%          [1:1:(n-1),(n+1):1:N] and the number of rows of the matrix Y
%          returned is equal to the number of elements in the nth way. On
%          the other hand, some authors, such as in [2], use a different
%          permutation. In such an instance, the explicit ordering can be
%          specified using param2=R and also providing C as in [3]. R
%          contains the indices of the modes that are to be mapped to rows
%          and C contains the indices of modes that are to be mapped to
%          columns. For example, for a 4-mode matrix one might use R=[2;3],
%          C=[4;1], which makes the dimensionality of the rows of the
%          output equal to the sum of the dimensionalities of
%          modes 2 and 3. 
%           
%OUTPUT: Y A 2D matrix containing the elements of the hypermatrix X
%          arranged according to the inputs.
%
%The algorithm can be called as
%Y=tensor2Mat(X,n);
%for a standard n-mode matricization or as
%Y=tensor2Mat(X,R,C);
%for general ordering to be used. The standard n-mode matricization is
%consistent with the orderings used in the nModeProd function.
%
%As an example, consider the 3X2X3 hypermatrix
% A=zeros(3,2,3);
% A(1,1,1)=1;
% A(1,1,2)=1;
% A(2,1,1)=1;
% A(2,1,2)=-1;
% A(2,1,3)=2;
% A(3,1,1)=2;
% A(3,1,3)=2;
% A(1,2,1)=2;
% A(1,2,2)=2;
% A(2,2,1)=2;
% A(2,2,2)=-2;
% A(2,2,3)=4;
% A(3,2,1)=4;
% A(3,2,3)=4;
% %A standard 1-mode matrix unfolding is 
% Y=tensor2Mat(A,1)
%Providing the result
% Y=[1     2     1     2     0     0;
%    1     2    -1    -2     2     4;
%    2     4     0     0     2     4];
%On the other hand, the example in [2] gives a different answer, because
%they use a different ordering. That is, they rearrange, the columns in the
%unfolding differently. Here, one can get the same answer using
% Y=tensor2Mat(A,1,[3,2])
%which returns
% Y=[1     1     0     2     2     0;
%    1    -1     2     2    -2     4;
%    2     0     2     4     0     4];
%
%Introductions to tensor operations are in [3] and [4].
%
%REFERENCES:
%[1] J. Salmi, A. Richter, and V. Koivunen, "Sequential unfolding SVD for
%    tensors with applications in array signal processing," IEEE
%    Transactions on Signal Processing, vol. 57, no. 12, pp. 4719-4733,
%    Dec. 2009.
%[2] L. de Lathauwer, B. de Moore, and J. Vandewalle, "A multilinear
%    singular value decomposition," SIAM Journal on Matrix Analysis and
%    Applications, vol. 21, no. 4, pp. 1253-1278, 2000.
%[3] T. G. Kolda, "Multilinear operators for higher-order decompositions,"
%    Sandia National Laboratories, Tech. Rep. SAND2006-2081, Apr. 2006.
%    [Online]. Available: http://www.sandia.gov/~tgkolda/pubs/pubfiles/SAND2006-2081.pdf
%[4] R. G. Kolda and B. W. Bader, "Tensor decompositions and applications,"
%    SIAM Review, vol. 51, no. 3, pp. 455-500, 2009.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

I=size(X);%The dimensionality of each way of X.

%If a standard n-mode matricization is being performed
if(nargin<3)
    d=length(I);
    C=[1:1:(param2-1),(param2+1):1:d];
end

R=param2;
J=prod(I(R));
K=prod(I(C));
Y=reshape(permute(X,[R(:);C(:)]),J,K);%Convert X to the matrix Y.

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
