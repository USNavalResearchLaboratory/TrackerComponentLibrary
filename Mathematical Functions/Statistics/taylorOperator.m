function L=taylorOperator(a,b,order,dfdx,dfdt,d2fdx)
%%TAYLOROPERATOR Apply the L operator to the drift or diffusion functions
%                of a nonlinear continuous-time random process. This
%                operator can be used to simplify some stochastic Taylor
%                approximations.
%
%INPUTS: a A xDimX1 vector of the current drift vector.
%        b A xDimxdColDim matrix of the current diffusion matrix.
%    order The order of the operator to apply. order==0 for the L0
%          operator, and anything else for the Lj operator. If order is an
%          empty matrix, all Lj outputs will be returned.
%     dfdx A xDimXxDimXdim3 matrix of the current derivative of the
%          function with respect to the state, where dim3 differentiates
%          between the drift and diffusion functions. If not supplied or an
%          empty matrix is passed, this is assumed to be zero for all
%          values.
%     dfdt A xDimXdim3 matrix of the current derivative of the function
%          with respect to the time, where dim3 differentiates between the
%          drift and diffusion functions. If not supplied or an empty
%          matrix is passed, this is assumed to be zero for all values.
%    d2fdx A xDimXxDimXxDimXdim3 hyper matrix of the current second
%          derivative of the function with respect to the state, where dim3
%          differentiates between the drift and diffusion functions. The
%          value at point (m,k,l,j) represents d2f(m,j)/dx(k)dx(l). If not
%          supplied or an empty matrix is passed, this is assumed to be
%          zero for all values.
%
%OUTPUTS: L The result of applying the L operator to either the drift or
%           diffusion function. If an order is provided, the size of L is
%           size xDimXdim3. If order is an empty matrix, L is size
%           xDimXdColDimXdim3, where dim3 differentiates between the drift
%           and diffusion functions
%
%This operator is used in Chapter 10.1 of [1] and is defined in Chapter
%5.3. The L0 operator is defined in Equation 3.1 of Chapter 5.3 and the Lj
%operator is defined in Equation 3.2 of Chapter 5.3. These are under Itô
%calculus, not Stratonovich.
%
%The operators perform partial derivatives. Using dp to represent the
%partial derivaitve operator, we have 
% L0=dp/(dp t)+sum_{k=1}^xDim a(k) dp/(dp x(k))+(1/2)sum_{k,l=1}^xDim sum_{j=1}^dColDim B(k,j)*b(l,j) dp^2/(dp x(k) dp x(l))
% Lj=sum_{k=1}^xDim b(k,j) dp/(dp x(k))
%where a and b are the drift and diffusion functions evaluated at the
%current point.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%April 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(a,1);
dColDim=size(b,2);

if(nargin<4||isempty(dfdx))
    dfdx=zeros(xDim,xDim);
end
dim3=size(dfdx,3);

if(order==0)
    if(nargin<5||isempty(dfdt))
        dfdt=zeros([xDim dim3]);
    end
    if(nargin<6||isempty(d2fdx))
        d2fdx=zeros([xDim xDim xDim dim3]);
    end
    
    term1=zeros(xDim,dim3);
    for j=1:dim3
        term1(:,j)=dfdx(:,:,j)*a;
    end
    term2=zeros(xDim,dim3);
    for k=1:xDim
        for l=1:xDim
            for j=1:dColDim
                term2(k,:)=term2(k,:)+(1/2)*b(k,j)*b(l,j)*sum(reshape(d2fdx(:,k,l,:),[xDim dim3]),1);
            end
        end
    end
    L=dfdt+term1+term2; 
else   
    L=zeros(xDim,dim3,dColDim);
    for j=1:dim3
        L(:,j,:)=dfdx(:,:,j)*b;
    end
    
    if(~isempty(order))
        L=L(:,:,order);
    end
    if(dim3==1)
        L=squeeze(L);
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
