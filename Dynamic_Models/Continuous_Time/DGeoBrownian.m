function [DVal,DJacob,DHess,pDpt]=DGeoBrownian(x,DMat)
%%DGEOBROWNIAN This provides the diffusion matrix of the multivariate
%              geometric Brownian motion process, also known as the Black-
%              Scholes models. The stochastic differential equation for
%              the ith component of the model is:
%              dx(i)=a(i)*x(i) dt+sum_{j=1}^m DMat(i,j)x(i) dw(j) 
%              where dw is the differential of a Wiener process and a is a
%              vector of constants. This diffusion function goes with
%              the drift function aGeoBrownian.
%
%INPUTS: x The xDimX1 target state vector.
%     DMat The xDimXm diffusion coefficient matrix.
%
%OUTPUTS: DVal The xDimXm diffusion matrix.
%       DJacob The xDimXmXxDim hypermatrix such that DJacob(:,:,k) is the
%              partial derivative of DVal with respect to the kth element
%              of x.
%        DHess The xDimXmXxDim hypermatrix such that DHess(:,:,k1,k2) is
%              the second partial derivative of DVal with respect to
%              elements k1 and k2 of x (all zeros in this case).
%         pDpt The xDimXm matrix of the partial derivative of DVal with
%              respect to t (all zeros in this case).
%
%The model is given in Chapter 2.4 of [1].
%
%REFERENCES:
%[1] E. Platen and N. Bruti-Liberati, Numerical Solution of Stochastic
%    Differential Equations with Jumps in Finance. Berlin: Springer-Verlag,
%    2010.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

DVal=bsxfun(@times,x,DMat);

if(nargout>1)
    d=size(DMat,1);
    m=size(DMat,2);
    DJacob=zeros(d,m,d);
    for k=1:d
        DJacob(k,:,k)=DMat(k,:);
    end

    if(nargout>2)
        DHess=zeros(d,m,d,d);
        if(nargout>3)
           pDpt=zeros(d,m); 
        end
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
