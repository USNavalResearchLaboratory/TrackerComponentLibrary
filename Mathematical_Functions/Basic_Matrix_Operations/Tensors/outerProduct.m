function C=outerProduct(varargin)
%%OUTPERPRODUCT Take the outer product of two or more vectors, matrices or
%        hypermatrices (tensors). For two tensors, the output of this
%        function is the product such that
%        C(ia(1),ia(2),...,ia(ndA),ib(1),ib(2),...,ib(ndB))=
%                                 A(ia(1),...,ia(ndA))*B(ib(1),...,ib(ndB))
%        where ia is a vector of indices for (hyper)matrix A and ib is a
%        vector of indices for (hyper)matrix B  This defines a tensor outer
%        product. Trailing singleton dimensions are not used when computing
%        the indexation of C. This is sometimes just called a tensor
%        product. The function can be given multiple matrices. When given,
%        for example 4, A1,A2,A3,and A4, the order of operations is
%        A1*(A2*(A3*A4)) where in that case, the * represents the tensor
%        outper product.
%
%INPUTS: varargin This is either a sequence of matrices or hypermatrices,
%                 as in A, B, C, etc. or a cell array containing all of the
%                 matrices.
%
%OUTPUTS: C The outer product of the matrices given on the input.
%
%EXAMPLE 1:
%Two vectors become a matrix and the outer product with a third will turns
%it into a hyper matrix.
% A=[1;2;3;4];
% B=[5;6;7;8;3];
% C=[-1;1];
% D=outerProduct(A,B)
% %D will be a 4X5 matrix such that D(ia,ib)=A(ia)*B(ib)
% E=outerProduct(D,C)
% %E will be a 4X5X2 matrix such that E(ia,ib,ic)=A(ia)*B(ib)*C(ic)
%
%EXAMPLE 2:
%Non-trailing singleton dimensions will remain singleton in the output
%(They can be suppressed using the squeeze function on the output).
% A=[1;2;3;4];
% B=[5,6,7,8,3];
% C=outerProduct(A,B)
%C will be a 4X1X5 matrix, because B is 1X5.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(iscell(varargin{1}))
    %Assume that a cell array was passed instead of multiple inputs.
    varargin=varargin{1};
end

numMat=length(varargin);

if(numMat==2)
    A=varargin{1};
    B=varargin{2};
    C=twoMatOuterProduct(A,B);
else
    C=varargin{numMat};
    
    for curMat=(numMat-1):-1:1
        C=twoMatOuterProduct(varargin{curMat},C);
    end
end
end

function C=twoMatOuterProduct(A,B)

%Get the dimensions of the inputs and remove trailing ones that are
%appended automatically by Matlab for vectors.
dims1=size(A);
dims2=size(B);

if(numel(dims1)==2&&dims1(2)==1)
    dims1=dims1(1); 
end

if(numel(dims2)==2&&dims2(2)==1)
    dims2=dims2(1); 
end

C=reshape(reshape(A,[],1)*reshape(B,1,[]),[dims1,dims2]);

end
