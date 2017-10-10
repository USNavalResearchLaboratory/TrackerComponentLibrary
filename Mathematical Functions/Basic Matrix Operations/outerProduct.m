function C=outerProduct(A,B)
%%OUTPERPRODUCT Take the outer product of two vectors, matrices or
%        hypermatrices (tensors). The output of this function is the
%        product such that
%        C(ia(1),ia(2),...,ia(ndA),ib(1),ib(2),...,ib(ndB))=
%                                 A(ia(1),...,ia(ndA))*B(ib(1),...,ib(ndB))
%        where ia is a vector of indices for (hyper)matrix A and ib is a
%        vector of indices for (hyper)matrix B  This defines a tensor outer
%        product. Trailing singleton dimensions are not used when computing
%        the indexation of C.
%
%INPUTS: A A scalar, vector, matrix or hypermatrix.
%        B A scalar, vector, matrix or hypermatrix.
%
%OUTPUTS: C The outer product of A and B.
%
%EXAMPLE 1:
%In the simplest case, two vectors become a matrix and the outer product
%with a third will turn it into a hyper matrix.
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