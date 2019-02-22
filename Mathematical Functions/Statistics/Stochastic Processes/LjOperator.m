function Lj=LjOperator(b,pfpx)
%%LJOPERATOR This function implements the Lj operation that is used in
%            Chapter 10 of [1] as applied to a single vector or matrix
%            argument and is defined in Chapter 5.3 of [1]. This is a
%            differential operator under Itô calculus. It is:
%            L^j=sum_{k=1}^d b(k,j) partial/(partial x_k)
%            where partial is the nabla (partial derivative) operator. This
%            does not implement multiple Lj operators, that is L^{j1}L^{j2}
%            applied to function.
%
%INPUTS: b The dXm diffusion  matrix of the stochastic differential
%          equation of the form dx=a*dt+b*dW at the time under
%          consideration, where x is the value being estimated and a and b
%          typically depend on x.
%     pfpx The d1Xd or d1Xd2Xd matrix of partial derivatives of the
%          function to which the operator L_j is being applied with respect
%          to the elements of x.
%
%OUTPUTS: Lj A d1Xm or d1Xd2Xm (depending on pfpx) matrix that results from
%            the application of the Lj operator.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

fDims=size(pfpx);

d=size(b,1);%State dimensionality. a is dX1.
m=size(b,2);%Noise dimensionality. b is dXm.

numFDims=ndims(pfpx)-1;%This should be 1 or 2.

if(numFDims>1)
    dF=prod(fDims(1:numFDims));%Total dimensionality of f.
    
    pfpx=reshape(pfpx,[dF,d]);
else
    dF=fDims(1);
end
    
Lj=zeros(dF,m);%Lj(k,j)
for j=1:m
    for k=1:d
        Lj(:,j)=Lj(:,j)+b(k,j)*pfpx(:,k);%Derivative with respect to x_k
    end
end

if(numFDims>1)
    Lj=reshape(Lj,[fDims(1:numFDims),m]);
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
