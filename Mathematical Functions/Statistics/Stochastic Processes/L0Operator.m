function L0=L0Operator(a,b,pfpt,pfpx,p2fpxpx)
%%L0OPERATOR This function implements the L0 operation that is used in
%         Chapter 10 of [1] as applied to a single vector or matrix
%         argument and is defined in Chapter 5.3 of [1]. This is a
%         differential operator under Itô calculus. It is:
%         L^0=partial/(partial t)+sum_{k=1}^d a(k)*partial/(partialx(k))+
%             (1/2)*sum_{k=1}^d sum_{l=1}^d sum_{j=1}^m
%                       b(k,j)*b(l,k)*partial^2/(partial x(k)*partial x(l)
%         where partial is the nabla (partial derivative) operator.
%
%INPUTS: a The dX1 drift matrix of the stochastic differential equation of
%          the form dx=a*dt+b*dW at the time under consideration, where x
%          is the value being estimated and a and b typically depend on x.
%        b The dXm diffusion matrix.
%     pfpt The d1 or d1Xd2 partial derivatives of the elements of the
%          function f (onto which the operator is being applied) with
%          respect to t.
%     pfpx The d1Xd or d1Xd2Xd partial derivatives of the elements of the
%          function f with respect to the elements of the dX1 state x.
%  p2fpxpx The d1XdXd or d1Xd2XdXd matrix of second partial derivatives of
%          the elements of f with respect to the elements of x.
%
%OUTPUTS: L0 The d1X1 or d1Xd2 (depending on dfdx) matrix that results from
%            the application of the L0 operator. 
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

fDims=size(pfpx);

d=size(a,1);%State dimensionality. a is dX1.
m=size(b,2);%Noise dimensionality. b is dXm.

numFDims=ndims(pfpx)-1;

if(numFDims>1)
    dF=prod(fDims(1:numFDims));%Total dimensionality of f.
    
    pfpx=reshape(pfpx,[dF,d]);
    pfpt=pfpt(:);
    
    p2fpxpx=reshape(p2fpxpx,[dF,d,d]);
else
    dF=fDims(1);
end

L0=pfpt;%First term in the sum.

%Second term in the sum.
for k=1:d
    L0=L0+a(k)*pfpx(:,k);
end

term3=zeros(dF,1);
for k=1:d
    for l=1:d
        for j=1:m
            term3=term3+b(k,j)*b(l,j)*p2fpxpx(:,k,l);
        end
    end
end
L0=L0+term3/2;

if(numFDims>1)
    L0=reshape(L0,fDims(1:numFDims));
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
